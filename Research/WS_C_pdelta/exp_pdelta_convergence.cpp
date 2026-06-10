// exp_pdelta_convergence.cpp — WS-C / N3 frozen-factor P-Delta experiment (scratch, NOT engine).
//
// Cantilever column (base fixed, tip axial compression P + lateral H), meshed into nElem
// beams. Three solutions compared:
//   exact   : beam-column stability-function tip deflection  delta = H (tan(kL) - kL) / (P k),
//             k = sqrt(P/EI)   (limit H L^3 / 3EI as P -> 0)
//   ref     : discretized second-order solve  (K_ff + Kg_ff(N)) u = F   (fresh LDLT on K_T;
//             N frozen at the first-order axial = P, standard Th.II practice)
//   frozen  : N3 pseudo-load iteration  u <- K_ff^-1 (F - Kg_ff u)  reusing the ORIGINAL
//             LDLT only (no refactor), plain and Aitken-accelerated.
//
// Sweep P/Pcr in {0, .1, .3, .5, .7, .9, .95}; divergence demo at 1.05; mesh 8 vs 16.

#include "research_common.h"
#include <cstring>

using namespace research;

namespace {

struct FrozenResult { VecX u; int iters = 0; bool converged = false; bool diverged = false; };

// plain (aitken=false) or Aitken-Delta^2 accelerated (aitken=true) pseudo-load iteration
FrozenResult frozenIterate(const LDLTSolver& ldlt0, const SpMat& Kgff, const VecX& Fff,
                           int maxIter, real tol, bool aitken) {
    FrozenResult R;
    VecX x = ldlt0.solve(Fff);          // first-order start
    VecX dPrev;                          // previous increment (for Aitken)
    real normPrev = -1;
    for (int m = 1; m <= maxIter; ++m) {
        VecX xNext = ldlt0.solve(Fff - Kgff * x);
        VecX d = xNext - x;
        const real dn = d.norm(), xn = std::max<real>(1e-300, xNext.norm());
        if (dn / xn < tol) { R.u = xNext; R.iters = m; R.converged = true; return R; }
        if (normPrev > 0 && dn > 4.0 * normPrev && dn / xn > 1.0) {   // geometric growth
            R.u = xNext; R.iters = m; R.diverged = true; return R;
        }
        if (aitken && dPrev.size() == d.size() && dPrev.size() > 0) {
            const real denom = dPrev.dot(dPrev);
            if (denom > 0) {
                const real rho = d.dot(dPrev) / denom;     // estimated contraction ratio
                if (rho > 0 && rho < 0.999) {
                    xNext = xNext + (rho / (1.0 - rho)) * d;   // geometric-series extrapolation
                    d = xNext - x;                              // keep bookkeeping consistent
                }
            }
        }
        dPrev = d; normPrev = dn;
        x = xNext;
    }
    R.u = x; R.iters = maxIter;
    return R;
}

void runCase(int nElem, real frac) {
    const real L = 6000.0, side = 200.0, H = 1000.0, E = 200000.0;
    const real I = side * side * side * side / 12.0;
    const real PI = 3.14159265358979323846;
    const real Pcr = PI * PI * E * I / (4.0 * L * L);
    const real P = frac * Pcr;

    FrameModel m = makeColumn(nElem, L, side, P, H);
    assertNodalOnly(m, "exp_pdelta");
    PreparedSystem ps = assembleAndFactor(m);
    const PreparedSystem::Impl& S = *ps.impl;
    if (S.singular) { std::printf("[pdelta] frac=%.2f FATAL singular\n", frac); return; }

    SolveResult lin = solveLoad(ps, m);
    std::vector<real> axial(m.members.size(), 0.0);
    for (size_t e = 0; e < lin.memberForces.size(); ++e) axial[e] = lin.memberForces[e].endI.N;

    std::vector<Triplet> gtrips;
    for (const auto& el : S.elems) el->assembleGeometric(gtrips, axial);
    SpMat Kg(S.N, S.N);
    Kg.setFromTriplets(gtrips.begin(), gtrips.end());
    const SpMat Kgff = reduceFF(Kg, S.fmap, S.nf);
    const SpMat Kff  = reduceFF(S.K, S.fmap, S.nf);
    const VecX  Fff  = reduceVec(nodalLoadVector(m), S.fmap, S.nf);

    const int tipDof = S.fmap[static_cast<size_t>(gdof(nElem, Ux))];

    // exact (P > 0) or Euler-Bernoulli closed form (P == 0)
    real tipExact;
    if (P > 0) {
        const real k = std::sqrt(P / (E * I));
        tipExact = H * (std::tan(k * L) - k * L) / (P * k);
    } else {
        tipExact = H * L * L * L / (3.0 * E * I);
    }

    // ref: K_T = Kff + Kgff, fresh LDLT
    real tipRef = 0; bool refOk = false;
    {
        SpMat KT = Kff + Kgff;
        LDLTSolver ldltT; ldltT.compute(KT);
        if (ldltT.info() == Eigen::Success) {
            const VecX uref = ldltT.solve(Fff);
            tipRef = uref(tipDof); refOk = true;
        }
    }

    // frozen plain + aitken
    const FrozenResult fp = frozenIterate(S.ldlt, Kgff, Fff, 5000, 1e-14, false);
    const FrozenResult fa = frozenIterate(S.ldlt, Kgff, Fff, 5000, 1e-14, true);

    const real tipLin = lin.u[static_cast<size_t>(gdof(nElem, Ux))];
    const real B1 = (frac < 1.0) ? 1.0 / (1.0 - frac) : -1.0;

    std::printf("[pdelta] nElem=%d frac=%.2f itersPlain=%d conv=%d div=%d itersAitken=%d convA=%d | tipLin=%.6g tipFrozen=%.6g tipRef=%.6g tipExact=%.6g | relRefVsExact=%.3e relFrozenVsRef=%.3e ampActual=%.4f B1=%.4f\n",
                nElem, frac, fp.iters, fp.converged ? 1 : 0, fp.diverged ? 1 : 0,
                fa.iters, fa.converged ? 1 : 0,
                static_cast<double>(tipLin), static_cast<double>(fp.u(tipDof)),
                static_cast<double>(tipRef), static_cast<double>(tipExact),
                refOk ? std::abs(tipRef - tipExact) / std::abs(tipExact) : -1.0,
                (refOk && fp.converged) ? std::abs(fp.u(tipDof) - tipRef) / std::abs(tipRef) : -1.0,
                static_cast<double>(fp.u(tipDof) / tipLin), B1);
}

} // namespace

int main() {
    for (real f : { 0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95 }) runCase(8, f);
    runCase(16, 0.3);     // mesh refinement: discretization error should drop ~4x
    runCase(8, 1.05);     // beyond Pcr: frozen iteration must report divergence
    std::printf("[done]\n");
    return 0;
}
