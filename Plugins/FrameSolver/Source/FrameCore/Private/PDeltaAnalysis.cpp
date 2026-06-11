#include "FrameCore/PDeltaAnalysis.h"
#include "PreparedSystemImpl.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace frame {

namespace {

// Free-free (nf x nf) sparse submatrix of a full N x N sparse matrix via the free-DOF map.
// (Same scatter the buckling path uses; kept local so each analysis TU is self-contained.)
SpMat reduceFF(const SpMat& Afull, const std::vector<int>& fmap, int nf) {
    std::vector<Triplet> t;
    t.reserve(static_cast<size_t>(Afull.nonZeros()));
    for (int c = 0; c < Afull.outerSize(); ++c)
        for (SpMat::InnerIterator it(Afull, c); it; ++it) {
            const int r = it.row();
            if (fmap[r] >= 0 && fmap[c] >= 0) t.emplace_back(fmap[r], fmap[c], it.value());
        }
    SpMat R(nf, nf);
    R.setFromTriplets(t.begin(), t.end());
    R.makeCompressed();
    return R;
}

// Positive-definiteness test of an LDLT-factored matrix from its diagonal D. SimplicialLDLT
// reports info()==Success even for an indefinite matrix (it does not pivot for stability), so a
// non-positive / near-zero entry in D is how a singular-or-indefinite K_T is detected — mirrors
// the mechanism-detection test in assembleAndFactor. tol is relative to max|D| (scale-invariant).
bool ldltPositiveDefinite(const LDLTSolver& s, real relTol) {
    const VecX D = s.vectorD();
    real maxAbs = 0;
    for (int i = 0; i < D.size(); ++i) maxAbs = std::max(maxAbs, std::abs(D(i)));
    const real tol = relTol * std::max<real>(1, maxAbs);
    for (int i = 0; i < D.size(); ++i)
        if (!(D(i) > tol)) return false;
    return true;
}

// Recover the FULL second-order state from a reduced free-DOF displacement uf. Scatters in the
// prescribed (support) values, rebuilds the load vector to report reactions, and recovers member
// end forces AT the second-order displacement (= the second-order internal forces). The reactions
// use the LINEAR elastic stiffness K_e (K*u - F): an honest "nodal force residual at the deflected
// shape", not a K_T equilibrium (the geometric term is a fictitious stress-stiffening operator, not
// an applied force). finalState mirrors solveLoad's tail so consumers get an ordinary SolveResult.
SolveResult recoverState(const PreparedSystem::Impl& S, const FrameModel& model, const VecX& uf) {
    SolveResult R;
    const int N = S.N;
    R.u.assign(static_cast<size_t>(std::max(0, N)), 0.0);
    R.reactions.assign(static_cast<size_t>(std::max(0, N)), 0.0);

    std::vector<real> presc(static_cast<size_t>(N), 0.0);
    for (size_t k = 0; k < model.nodes.size(); ++k)
        for (int d = 0; d < 6; ++d)
            if (model.nodes[k].fixed[d]) presc[static_cast<size_t>(gdof((int)k, d))] = model.nodes[k].prescribed[d];

    VecX u = VecX::Zero(N);
    for (int g = 0; g < N; ++g) u(g) = (S.fmap[g] >= 0) ? uf(S.fmap[g]) : presc[static_cast<size_t>(g)];
    for (int g = 0; g < N; ++g) R.u[static_cast<size_t>(g)] = u(g);

    VecX F = VecX::Zero(N);
    for (const auto& nl : model.nodalLoads) {
        const int ni = model.nodeIndex(nl.node);
        if (ni < 0) continue;
        for (int d = 0; d < 6; ++d) F(gdof(ni, d)) += nl.comp[d];
    }
    for (const auto& el : S.elems) el->addEquivalentNodalLoads(F);
    const VecX Rv = S.K * u - F;
    for (int g = 0; g < N; ++g) R.reactions[static_cast<size_t>(g)] = Rv(g);

    R.memberForces.resize(model.members.size());
    R.shellForces.resize(model.shells.size());
    for (size_t e = 0; e < model.members.size(); ++e) R.memberForces[e].member = model.members[e].id;
    for (size_t s = 0; s < model.shells.size(); ++s)  R.shellForces[s].shell   = model.shells[s].id;
    for (const auto& el : S.elems) el->recover(u, R);
    R.pivotMargin = S.pivotMargin;
    return R;
}

}  // namespace

PDeltaResult runPDelta(const FrameModel& model, const PDeltaOptions& opts) {
    PDeltaResult R;

    // 1) factorize once + first-order linear solve (axial forces freeze Kg).
    PreparedSystem ps = assembleAndFactor(model, opts.solve);
    const PreparedSystem::Impl& S = *ps.impl;
    if (S.singular) { R.finalState.singular = true; R.finalState.diagnostic = S.diagnostic; return R; }

    const SolveResult lin = solveLoad(ps, model);
    if (lin.singular) { R.finalState = lin; return R; }
    const int N = S.N, nf = S.nf;

    std::vector<real> axial(model.members.size(), 0.0);
    for (size_t e = 0; e < lin.memberForces.size() && e < axial.size(); ++e)
        axial[e] = lin.memberForces[e].endI.N;   // compression-positive

    // 2) geometric stiffness Kg from the frozen axial forces; reduce to the free block.
    std::vector<Triplet> gtrips;
    for (const auto& el : S.elems) el->assembleGeometric(gtrips, axial);
    SpMat Kg(N, N);
    Kg.setFromTriplets(gtrips.begin(), gtrips.end());
    Kg.makeCompressed();
    const SpMat Kgff = reduceFF(Kg, S.fmap, nf);

    // P = 0 (no member in compression -> empty Kg): the second-order state IS the linear state.
    // Return it bit-for-bit (first step is an immediate convergence with zero increment).
    if (Kgff.nonZeros() == 0) {
        R.converged = true; R.iterations = 0; R.lastIncrement = 0; R.finalState = lin;
        return R;
    }

    const SpMat Kff = reduceFF(S.K, S.fmap, nf);

    // Fixed external force on the free DOFs, recovered from the linear solve:  K_ff u_lin_ff = F_ff
    // (since u_lin_ff = ldlt0.solve(F_ff)). One sparse mat-vec; automatically carries the nodal
    // loads, the distributed equivalent loads AND the prescribed -K_fc u_c term. Both paths use
    // THIS same F_ff, so frozen-vs-reference agreement is independent of any rebuild round-off.
    VecX ulinff = VecX::Zero(nf);
    for (int g = 0; g < N; ++g) if (S.fmap[g] >= 0) ulinff(S.fmap[g]) = lin.u[static_cast<size_t>(g)];
    const VecX Fff = Kff * ulinff;

    // ---- REFERENCE PATH: one fresh LDLT of K_T = K_e + Kg (Wilson & Habibullah 1987) ----
    if (opts.refactorPath) {
        const SpMat KT = Kff + Kgff;
        LDLTSolver ldltT; ldltT.compute(KT);
        if (ldltT.info() != Eigen::Success || !ldltPositiveDefinite(ldltT, opts.solve.pivotTol)) {
            R.diverged = true;
            R.finalState.singular = true;
            R.finalState.diagnostic = "P-Delta reference K_T not positive-definite (P exceeds critical)";
            return R;
        }
        const VecX uff = ldltT.solve(Fff);
        if (ldltT.info() != Eigen::Success || !uff.allFinite()) {
            R.diverged = true;
            R.finalState.singular = true;
            R.finalState.diagnostic = "P-Delta reference solve produced non-finite displacements";
            return R;
        }
        R.converged = true; R.iterations = 1;
        R.lastIncrement = 0;
        R.finalState = recoverState(S, model, uff);
        return R;
    }

    // ---- FROZEN PATH: pseudo-load iteration reusing the existing K_e factorization ----
    //   x_{k+1} = K_e^-1 ( F_ff - Kg_ff x_k ),   x_0 = u_lin_ff
    // Converges (rate ~ P/P_cr) for P < P_cr. Two [NEW CODE] safeguards harden the bare iteration
    // the research prototype used:
    //   (a) PROTECTED extrapolation: Aitken geometric extrapolation gated by a stability window
    //       (0<rho<0.95 AND |rho_k - rho_{k-1}|<0.2) with UNDO — if the step after an extrapolation
    //       does not reduce ||du||, that extrapolation is reverted and acceleration is disabled for
    //       the rest of this solve (so a mis-estimated rho can never do worse than bare iteration).
    //   (b) DIVERGENCE detector: a 20-step sliding window catches the slow growth (~x1.05/step) that
    //       a single-step ratio test misses past P_cr; maxIter is the final backstop.
    VecX x = ulinff;
    VecX prevD;                 // previous increment (for the contraction-ratio estimate)
    real prevDn  = -1.0;        // ||previous increment||
    real prevRho = -2.0;        // previous contraction ratio (sentinel: none yet)
    bool accel   = opts.accelerate;

    bool pendingUndo = false;   // did the LAST accepted step come from an extrapolation?
    VecX undoX;                 // the un-extrapolated x of that step, to revert to
    real preExtrapDn = 0.0;     // its increment norm, to test "did the next step reduce ||du||?"

    std::vector<real> dHist;
    dHist.reserve(static_cast<size_t>(std::max(1, opts.maxIter)));

    real lastRel = 0.0;
    for (int m = 1; m <= opts.maxIter; ++m) {
        const VecX xNext = S.ldlt.solve(Fff - Kgff * x);
        VecX d   = xNext - x;
        const real dn = d.norm();
        const real xn = std::max<real>(1e-300, xNext.norm());
        lastRel = dn / xn;

        if (lastRel < opts.tolU) {
            R.converged = true; R.iterations = m; R.lastIncrement = lastRel;
            R.finalState = recoverState(S, model, xNext);
            return R;
        }

        // (a) UNDO: the previous step was an extrapolation; if this step did not reduce ||du||,
        //     revert to the un-extrapolated point and fall back to bare iteration.
        if (pendingUndo) {
            pendingUndo = false;
            if (!(dn < preExtrapDn)) {
                x = undoX; accel = false;
                prevD = VecX(); prevDn = -1.0; prevRho = -2.0;
                continue;
            }
        }

        // (b) DIVERGENCE: within a 20-step window the increment rose above the window start AND the
        //     very first increment. For P<P_cr the increments contract geometrically (never trips);
        //     for P>P_cr they grow (trips). maxIter below is the backstop if growth is non-monotone.
        dHist.push_back(dn);
        if (dHist.size() >= 20) {
            const real winFront = dHist[dHist.size() - 20];
            if (dn > winFront && dn > dHist.front()) {
                R.diverged = true; R.iterations = m; R.lastIncrement = lastRel;
                R.finalState.singular = true;
                R.finalState.diagnostic = "P-Delta frozen iteration diverged (P exceeds critical)";
                return R;
            }
        }

        // protective geometric extrapolation
        real rho = -2.0;
        VecX xAdv = xNext;
        if (accel && prevD.size() == d.size() && prevDn > 0) {
            const real denom = prevD.dot(prevD);
            if (denom > 0) rho = d.dot(prevD) / denom;
            if (rho > 0 && rho < 0.95 && prevRho > -2.0 && std::abs(rho - prevRho) < 0.2) {
                undoX = xNext; preExtrapDn = dn; pendingUndo = true;
                xAdv = xNext + (rho / (1.0 - rho)) * d;
                d = xAdv - x;     // keep the increment bookkeeping consistent with the advanced point
            }
        }
        if (rho > -2.0) prevRho = rho;
        prevD = d; prevDn = d.norm();
        x = xAdv;
    }

    // maxIter reached: no convergence, no divergence verdict. Report the best-effort state.
    R.converged = false; R.iterations = opts.maxIter; R.lastIncrement = lastRel;
    R.finalState = recoverState(S, model, x);
    return R;
}

}  // namespace frame
