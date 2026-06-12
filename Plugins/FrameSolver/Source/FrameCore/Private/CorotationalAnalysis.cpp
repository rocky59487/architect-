#include "FrameCore/CorotationalAnalysis.h"
#include "FrameEigen.h"
#include "FrameCore/FrameModel.h"
#include "FrameCore/SolveResult.h"

#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

// S9 -- PLANAR co-rotational beam (geometric nonlinearity), XY plane, bending about global Z.
//
// Each member carries a 3-DOF natural deformation [u_bar, theta_bar_i, theta_bar_j] in its CURRENT
// co-rotated chord frame; the local stiffness is the small-strain Euler-Bernoulli beam (EA, EIz). The
// 3x6 strain matrix B maps the planar global DOFs [ux,uy,rz]x2 to the natural deformation, the internal
// force is B^T q_l (virtual-work consistent), and the tangent is the material part B^T k_l B plus the
// co-rotational geometric stiffness
//     Kg = (N/Ln) q q^T + ((Mi+Mj)/Ln^2) (p q^T + q p^T),
//   p = [-c,-s,0, c,s,0]^T  (= d Ln / d d),   q = [s,-c,0,-s,c,0]^T  (= Ln * d beta / d d).
// Both p,q derivations and Kg are the standard 2D corotational beam (Crisfield Vol.1 Ch.7); Kg is
// symmetric. Planar rotations about the single Z axis commute, so the nodal angle is read directly from
// the Rz displacement DOF -- no SO(3) update needed (3D CR is S9b). A rigid in-plane rotation gives
// u_bar=0 and theta_bar=theta_z-(beta-beta0)=0, hence zero internal force (the defining CR property;
// exercised by the rotational-invariance oracle).

namespace frame {
namespace {

constexpr real kPi = 3.14159265358979323846;

// Per-member planar co-rotational element (initial invariants + unwrap state).
struct CrBeam {
    int      e      = -1;
    MemberId id     = 0;
    int      ni     = -1, nj = -1;     // node indices
    int      gmap[6] = { 0,0,0,0,0,0 };// global DOF of [uxi,uyi,rzi,uxj,uyj,rzj]
    real     L0     = 0;               // initial length
    real     beta0  = 0;               // initial chord angle
    real     EA     = 0, EIz = 0;
    real     betaPrev = 0;             // last chord angle (atan2 unwrap continuity)
    // recovered natural state (for member-force output)
    real     N = 0, Mi = 0, Mj = 0, Ln = 0;
};

// Build the planar internal force fe(6) and tangent Ke(6x6) at the current displacement u. Also returns
// the (unwrapped) chord angle so the caller can advance betaPrev, and stores N/Mi/Mj on the element.
void crCompute(CrBeam& b, const FrameModel& M, const std::vector<real>& u,
               Eigen::Matrix<real, 6, 1>& fe, Eigen::Matrix<real, 6, 6>& Ke, real& betaOut) {
    const Vec3 Xi = M.nodes[(size_t)b.ni].pos, Xj = M.nodes[(size_t)b.nj].pos;
    const real uxi = u[(size_t)b.gmap[0]], uyi = u[(size_t)b.gmap[1]], tzi = u[(size_t)b.gmap[2]];
    const real uxj = u[(size_t)b.gmap[3]], uyj = u[(size_t)b.gmap[4]], tzj = u[(size_t)b.gmap[5]];

    const real x1 = Xi.x + uxi, y1 = Xi.y + uyi;
    const real x2 = Xj.x + uxj, y2 = Xj.y + uyj;
    const real dx = x2 - x1, dy = y2 - y1;
    const real Ln = std::max<real>(1e-300, std::sqrt(dx * dx + dy * dy));   // guard a fully collapsed chord

    real beta = std::atan2(dy, dx);
    while (beta - b.betaPrev >  kPi) beta -= 2 * kPi;       // unwrap relative to previous chord angle
    while (beta - b.betaPrev < -kPi) beta += 2 * kPi;
    betaOut = beta;

    const real c = std::cos(beta), s = std::sin(beta);
    const real ubar  = Ln - b.L0;
    const real alpha = beta - b.beta0;                      // chord rigid rotation
    const real tbi   = tzi - alpha;                         // local deformational rotations
    const real tbj   = tzj - alpha;

    const real ea = b.EA  / b.L0;
    const real ei = b.EIz / b.L0;
    const real N  = ea * ubar;                              // tension-positive natural axial
    const real Mi = ei * (4.0 * tbi + 2.0 * tbj);
    const real Mj = ei * (2.0 * tbi + 4.0 * tbj);
    b.N = N; b.Mi = Mi; b.Mj = Mj; b.Ln = Ln;

    // B (3x6) rows: u_bar, theta_bar_i, theta_bar_j over [uxi,uyi,rzi,uxj,uyj,rzj].
    Eigen::Matrix<real, 3, 6> B; B.setZero();
    B(0, 0) = -c;      B(0, 1) = -s;      B(0, 3) = c;       B(0, 4) = s;
    B(1, 0) = -s / Ln; B(1, 1) = c / Ln;  B(1, 2) = 1.0;     B(1, 3) = s / Ln; B(1, 4) = -c / Ln;
    B(2, 0) = -s / Ln; B(2, 1) = c / Ln;                     B(2, 3) = s / Ln; B(2, 4) = -c / Ln; B(2, 5) = 1.0;

    Eigen::Matrix<real, 3, 1> ql; ql << N, Mi, Mj;
    fe = B.transpose() * ql;

    Eigen::Matrix<real, 3, 3> kl; kl.setZero();
    kl(0, 0) = ea;
    kl(1, 1) = 4.0 * ei; kl(1, 2) = 2.0 * ei;
    kl(2, 1) = 2.0 * ei; kl(2, 2) = 4.0 * ei;
    Ke = B.transpose() * kl * B;                            // material tangent

    Eigen::Matrix<real, 6, 1> p, q;
    p << -c, -s, 0, c, s, 0;
    q <<  s, -c, 0, -s, c, 0;
    Ke += (N / Ln) * (q * q.transpose());                  // co-rotational geometric stiffness
    Ke += ((Mi + Mj) / (Ln * Ln)) * (p * q.transpose() + q * p.transpose());
}

// Positive-definite / non-singular test of an LDLT factor from its diagonal D (mirrors the
// assembleAndFactor mechanism test and PDeltaAnalysis::ldltPositiveDefinite). tol is relative to max|D|.
bool ldltPosDef(const LDLTSolver& solver, real relTol) {
    const VecX D = solver.vectorD();
    real maxAbs = 0;
    for (int i = 0; i < D.size(); ++i) maxAbs = std::max(maxAbs, std::abs(D(i)));
    const real tol = relTol * std::max<real>(1, maxAbs);
    for (int i = 0; i < D.size(); ++i)
        if (!(D(i) > tol)) return false;
    return true;
}

}  // namespace

CorotationalResult runCorotational(const FrameModel& model, const CorotationalOptions& opts) {
    CorotationalResult R;

    // A rejected / failed run still returns a WELL-FORMED SolveResult: u and reactions are zero-filled
    // to the model DOF count and memberForces are id-tagged, so a consumer that reads finalState.u
    // (e.g. the CLI printState) never indexes an empty vector. singular=true + a diagnostic mark it.
    auto reject = [&](const char* diag) -> CorotationalResult {
        CorotationalResult RR;
        const int n = std::max(0, model.dofCount());
        RR.finalState.singular = true;
        RR.finalState.diagnostic = diag;
        RR.finalState.u.assign((size_t)n, 0.0);
        RR.finalState.reactions.assign((size_t)n, 0.0);
        RR.finalState.memberForces.resize(model.members.size());
        for (size_t e = 0; e < model.members.size(); ++e) RR.finalState.memberForces[e].member = model.members[e].id;
        return RR;
    };
    // Fill a partial (e.g. diverged) finalState from a best-effort displacement snapshot.
    auto fillPartial = [&](CorotationalResult& RR, const std::vector<real>& uu) {
        RR.finalState.u = uu;
        RR.finalState.reactions.assign(uu.size(), 0.0);
        RR.finalState.memberForces.resize(model.members.size());
        for (size_t e = 0; e < model.members.size(); ++e) RR.finalState.memberForces[e].member = model.members[e].id;
    };

    // --- scope guards (honest: REJECT rather than silently mis-handle; all fill finalState.u) ---
    for (const auto& sh : model.shells)
        if (sh.active) return reject("co-rotational large-displacement is beam-column only; model contains shells");
    std::string why;
    if (!model.validate(why)) return reject(why.c_str());
    // nodal force loads only -- a member UDL would otherwise be silently dropped
    for (const auto& udl : model.memberUDLs)
        for (const auto& mem : model.members)
            if (mem.id == udl.member && mem.active)
                return reject("co-rotational v1 supports nodal force loads only; member UDLs not yet supported");
    // no prescribed support displacement
    for (const auto& nd : model.nodes)
        for (int d = 0; d < 6; ++d)
            if (nd.fixed[(size_t)d] && nd.prescribed[(size_t)d] != 0.0)
                return reject("co-rotational v1 does not support prescribed support displacements");
    // formed plastic hinges are not honoured (would be computed at FULL stiffness -> silently wrong),
    // consistent with the UDL rejection policy (S10 N-M hinges land after S9b, R4)
    for (const auto& h : model.hinges)
        for (const auto& mem : model.members)
            if (mem.id == h.member && mem.active)
                return reject("co-rotational v1 does not honour formed plastic hinges");
    // PLANAR geometry: every active member must lie in the global XY plane. The formulation uses only
    // x,y, so a member with out-of-plane z-extent would be SILENTLY computed as its xy-projection.
    for (const auto& mem : model.members) {
        if (!mem.active) continue;
        const int ni = model.nodeIndex(mem.i), nj = model.nodeIndex(mem.j);
        if (ni < 0 || nj < 0) continue;
        const Vec3 pi = model.nodes[(size_t)ni].pos, pj = model.nodes[(size_t)nj].pos;
        const real dx = pj.x - pi.x, dy = pj.y - pi.y, dz = pj.z - pi.z;
        const real L3 = std::sqrt(dx * dx + dy * dy + dz * dz);
        if (std::fabs(dz) > 1e-9 * std::max<real>(1.0, L3))
            return reject("co-rotational v1 requires all members in the global XY plane (out-of-plane z-extent); 3D CR is S9b");
    }

    const int N = model.dofCount();
    if (N <= 0) return reject("empty model");

    // --- free-DOF map ---
    std::vector<int> fmap((size_t)N, -1);
    int nf = 0;
    for (size_t k = 0; k < model.nodes.size(); ++k)
        for (int d = 0; d < 6; ++d)
            if (!model.nodes[k].fixed[(size_t)d]) fmap[(size_t)gdof((int)k, d)] = nf++;
    if (nf == 0) return reject("fully constrained model");

    // --- build planar CR elements from active members ---
    std::vector<CrBeam> elems;
    elems.reserve(model.members.size());
    for (size_t e = 0; e < model.members.size(); ++e) {
        const Member& m = model.members[e];
        if (!m.active) continue;
        const int ni = model.nodeIndex(m.i), nj = model.nodeIndex(m.j);
        if (ni < 0 || nj < 0) continue;
        CrBeam b;
        b.e = (int)e; b.id = m.id; b.ni = ni; b.nj = nj;
        b.gmap[0] = gdof(ni, Ux); b.gmap[1] = gdof(ni, Uy); b.gmap[2] = gdof(ni, Rz);
        b.gmap[3] = gdof(nj, Ux); b.gmap[4] = gdof(nj, Uy); b.gmap[5] = gdof(nj, Rz);
        const Vec3 Xi = model.nodes[(size_t)ni].pos, Xj = model.nodes[(size_t)nj].pos;
        const real dx = Xj.x - Xi.x, dy = Xj.y - Xi.y;
        b.L0 = std::sqrt(dx * dx + dy * dy);
        b.beta0 = std::atan2(dy, dx);
        b.betaPrev = b.beta0;
        const Material& mat = model.materials[(size_t)m.matIdx];
        const Section&  sec = model.sections[(size_t)m.secIdx];
        b.EA = mat.E * sec.A; b.EIz = mat.E * sec.Iz;
        elems.push_back(b);
    }

    // --- external nodal force vector (full). UDL / prescribed reserved for a later CR revision. ---
    VecX Fext = VecX::Zero(N);
    for (const auto& nl : model.nodalLoads) {
        const int ni = model.nodeIndex(nl.node);
        if (ni < 0) continue;
        for (int d = 0; d < 6; ++d) Fext((Eigen::Index)gdof(ni, d)) += nl.comp[(size_t)d];
    }
    VecX Fext_f = VecX::Zero(nf);
    for (int g = 0; g < N; ++g) if (fmap[(size_t)g] >= 0) Fext_f(fmap[(size_t)g]) = Fext((Eigen::Index)g);

    // --- Newton-Raphson with load stepping ---
    std::vector<real> u((size_t)N, 0.0);
    const int steps = std::max(1, opts.loadSteps);
    int totalIters = 0;
    real lastRel = 0;

    for (int s = 1; s <= steps; ++s) {
        const real lambda = real(s) / real(steps);
        bool stepConverged = false;

        for (int it = 1; it <= std::max(1, opts.maxIter); ++it) {
            // assemble internal force + tangent at current u
            VecX fint = VecX::Zero(N);
            std::vector<Triplet> trips;
            trips.reserve(elems.size() * 36);
            std::vector<real> betaNow(elems.size());
            for (size_t i = 0; i < elems.size(); ++i) {
                Eigen::Matrix<real, 6, 1> fe; Eigen::Matrix<real, 6, 6> Ke; real bo;
                crCompute(elems[i], model, u, fe, Ke, bo);
                betaNow[i] = bo;
                for (int a = 0; a < 6; ++a) {
                    fint((Eigen::Index)elems[i].gmap[a]) += fe(a);
                    for (int c = 0; c < 6; ++c)
                        trips.emplace_back(elems[i].gmap[a], elems[i].gmap[c], Ke(a, c));
                }
            }
            for (size_t i = 0; i < elems.size(); ++i) elems[i].betaPrev = betaNow[i];   // advance unwrap

            // residual on free DOFs: r = lambda*Fext - fint
            VecX rf = VecX::Zero(nf);
            for (int g = 0; g < N; ++g)
                if (fmap[(size_t)g] >= 0) rf(fmap[(size_t)g]) = lambda * Fext((Eigen::Index)g) - fint((Eigen::Index)g);

            // reduce tangent to free-free
            std::vector<Triplet> tf;
            tf.reserve(trips.size());
            for (const auto& t : trips) {
                const int fr = fmap[(size_t)t.row()], fc = fmap[(size_t)t.col()];
                if (fr >= 0 && fc >= 0) tf.emplace_back(fr, fc, t.value());
            }
            SpMat Kff((Eigen::Index)nf, (Eigen::Index)nf);
            Kff.setFromTriplets(tf.begin(), tf.end());
            Kff.makeCompressed();

            LDLTSolver ldlt; ldlt.compute(Kff);
            if (ldlt.info() != Eigen::Success || !ldltPosDef(ldlt, opts.solve.pivotTol)) {
                R.diverged = true;
                R.finalState.singular = true;
                R.finalState.diagnostic = "co-rotational tangent not positive-definite (limit point; snap-through needs arc-length)";
                R.loadStepsCompleted = s - 1;
                R.totalIterations = totalIters;
                R.lastResidual = lastRel;
                fillPartial(R, u);
                return R;
            }

            const VecX duf = ldlt.solve(rf);
            if (ldlt.info() != Eigen::Success || !duf.allFinite()) {
                R.diverged = true;
                R.finalState.singular = true;
                R.finalState.diagnostic = "co-rotational solve produced non-finite increment";
                R.loadStepsCompleted = s - 1;
                R.totalIterations = totalIters;
                fillPartial(R, u);
                return R;
            }

            // update u
            real un = 0;
            for (int g = 0; g < N; ++g)
                if (fmap[(size_t)g] >= 0) { u[(size_t)g] += duf(fmap[(size_t)g]); un += u[(size_t)g] * u[(size_t)g]; }
            ++totalIters;

            const real dn = duf.norm();
            const real Frn = (lambda * Fext_f).norm();
            const real relR = rf.norm() / std::max<real>(1e-300, Frn);
            const real relU = dn / std::max<real>(1.0, std::sqrt(un));
            lastRel = relR;
            if (relR < opts.tolR || relU < opts.tolU) { stepConverged = true; break; }
        }

        R.loadStepsCompleted = stepConverged ? s : (s - 1);
        if (!stepConverged) {
            // hit maxIter without a verdict: report best-effort state, not converged.
            R.converged = false;
            R.totalIterations = totalIters;
            R.lastResidual = lastRel;
            // still recover so the caller sees the partial state
            break;
        }
        if (s == steps) R.converged = true;
    }
    R.totalIterations = totalIters;
    R.lastResidual = lastRel;

    // --- recover finalState (SolveResult) ---
    SolveResult& SR = R.finalState;
    SR.u.assign((size_t)N, 0.0);
    SR.reactions.assign((size_t)N, 0.0);
    for (int g = 0; g < N; ++g) SR.u[(size_t)g] = u[(size_t)g];

    // recompute internal force at the final state for reactions + member forces
    VecX fint = VecX::Zero(N);
    for (size_t i = 0; i < elems.size(); ++i) {
        Eigen::Matrix<real, 6, 1> fe; Eigen::Matrix<real, 6, 6> Ke; real bo;
        crCompute(elems[i], model, u, fe, Ke, bo);
        for (int a = 0; a < 6; ++a) fint((Eigen::Index)elems[i].gmap[a]) += fe(a);
    }
    for (int g = 0; g < N; ++g) SR.reactions[(size_t)g] = fint((Eigen::Index)g) - Fext((Eigen::Index)g);

    SR.memberForces.resize(model.members.size());
    SR.shellForces.clear();
    for (size_t e = 0; e < model.members.size(); ++e) SR.memberForces[e].member = model.members[e].id;
    for (const CrBeam& b : elems) {
        const real V = (b.Mi + b.Mj) / std::max<real>(1e-300, b.Ln);   // transverse shear (current length)
        MemberForcePair& mp = SR.memberForces[(size_t)b.e];
        mp.member = b.id;
        // local end forces (planar): N compression-positive, in-plane shear Vy, bending Mz.
        mp.endI = MemberEndForces{ -b.N,  V, 0, 0, 0, b.Mi };
        mp.endJ = MemberEndForces{ -b.N, -V, 0, 0, 0, b.Mj };
    }
    SR.singular = false;
    return R;
}

}  // namespace frame
