#include "FrameCore/CorotationalAnalysis.h"
#include "FrameEigen.h"
#include "ElementStiffness.h"          // localAxes (initial element frame from refVec)
#include "FrameCore/FrameModel.h"
#include "FrameCore/SolveResult.h"

#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

// S9b -- 3D GENERAL co-rotational beam (geometric nonlinearity): arbitrary spatial orientation, torsion +
// biaxial bending + finite SO(3) rotation. Generalises the S9 PLANAR driver. Each member co-rotates with
// its CURRENT chord AND current section triad; the local strain stays small (small-strain large-rotation
// CR). The local stiffness is the same Euler-Bernoulli beam (EA/EIy/EIz/GJ) used by localStiffness12.
//
// KEY DIFFERENCE FROM S9: 3D finite rotations DO NOT COMMUTE, so the per-node rotation is carried as a
// matrix R_node in SO(3) (initial I) and updated by a SPATIAL increment R_node <- exp(skew(dtheta))*R_node
// after each NR step (avoids the total-rotation-vector 2.pi singularity, Battini 2002). The planar single-
// axis "accumulate Rz directly" of S9 is replaced by this. The element frame E=[e1|e2|e3] is rebuilt every
// iteration from the deformed geometry (chord) + the mean nodal triad (Gram-Schmidt), the local deformation
// is the nodal triad's rotation relative to E (logSO3), the internal force is f_int = T^T*pb (virtual-work
// consistent), and the tangent is the material part T^T*Kl*T plus the axial geometric stiffness Ksigma1.
//
// A rigid body motion R_g gives q_I=q_J=R_g*E0 and E=R_g*E0, so E^T*q = E0^T*E0 = const -> zero local
// deformation -> pb=0 -> f_int=0 (the defining CR property; exercised by the arbitrary-axis rotational-
// invariance oracle). In the planar limit (members in XY, rotation about Z only) every quantity reduces to
// the S9 planar formulation bit-for-bit modulo float path (the planar elastica / P-Delta oracles still pass).
//
// Tangent (WS_F2 section 6 + the S9-verified durable lesson that the NR converged solution depends ONLY on
// f_int, the tangent only on convergence speed): the tangent is T^T*Kl*T + Ksigma1 (the strict axial
// geometric term, reducing to S9's (N/Ln)q q^T, symmetric so the same SimplicialLDLT path as S9/PDelta is
// reused). The analytical spin/moment geometric terms (OpenSees Ksigma2/3) are NOT added -- the main term
// already converges the elastica alpha=1..10 (F50/F51), and S9c's opts.consistentTangent provides a numerical
// FD consistent tangent for quadratic NR. The analytical Ksigma2/3 stays a future convergence-speed optimisation.

namespace frame {
namespace {

// kPi (FrameEigen.h) and ldltPositiveDefinite (FrameEigen.h) are shared — no file-local copies.

// ----------------------------- SO(3) helpers (Mat3, hand-written, no Eigen/Geometry) -----------------------------
// skew(w): the cross-product matrix S with S*x = w x x.
inline Mat3 skew(const Vec3e& w) {
    Mat3 S; S << 0, -w(2), w(1),
                 w(2), 0, -w(0),
                 -w(1), w(0), 0;
    return S;
}
// vee of (R - R^T): the axial vector of the skew part, = 2 sin(theta) * axis.
inline Vec3e veeAsym(const Mat3& R) {
    return Vec3e(R(2, 1) - R(1, 2), R(0, 2) - R(2, 0), R(1, 0) - R(0, 1));
}
// Rodrigues exp: rotation matrix for rotation vector w (axis*angle). Series near 0 to avoid 0/0.
inline Mat3 expSO3(const Vec3e& w) {
    const real th = w.norm();
    const Mat3 S = skew(w);
    if (th < 1e-7) return Mat3::Identity() + S + 0.5 * S * S;           // sin/th->1, (1-cos)/th^2->1/2
    const real a = std::sin(th) / th;
    const real b = (1.0 - std::cos(th)) / (th * th);
    return Mat3::Identity() + a * S + b * (S * S);
}
// log of a rotation matrix -> rotation vector (axis*angle). Handles theta->0 (series) and theta->pi.
inline Vec3e logSO3(const Mat3& R) {
    real cth = 0.5 * (R.trace() - 1.0);
    cth = std::max<real>(-1.0, std::min<real>(1.0, cth));
    const real th = std::acos(cth);
    if (th < 1e-7) {
        // log ~ (1/2)(R-R^T) axial, with the (theta/(2 sin theta))->1/2 prefactor (series exact to O(th^2))
        return 0.5 * (1.0 + th * th / 6.0) * veeAsym(R);
    }
    if (th < kPi - 1e-5) {
        return (th / (2.0 * std::sin(th))) * veeAsym(R);
    }
    // theta near pi: sin(theta)~0, recover the axis from the symmetric part B = (R + I)/2 = axis axis^T.
    // The 1e-5 switch widens the safety band: the narrow theta in (pi-2e-6, pi-1e-6) seam where both
    // branches dip to ~1e-4 is harmless -- the rotation vector here is OUTPUT-only (every R used in the
    // computation is built by expSO3, exact everywhere), and no oracle or engineering rotation reaches it.
    Mat3 B = 0.5 * (R + Mat3::Identity());
    int k = 0; real best = B(0, 0);
    if (B(1, 1) > best) { best = B(1, 1); k = 1; }
    if (B(2, 2) > best) { best = B(2, 2); k = 2; }
    Vec3e axis;
    const real d = std::sqrt(std::max<real>(0.0, B(k, k)));
    axis(k) = d;
    if (d > 1e-12) { for (int i = 0; i < 3; ++i) if (i != k) axis(i) = B(k, i) / d; }
    axis.normalize();
    // fix the sign from the (vanishing but directional) skew part
    const Vec3e s = veeAsym(R);
    if (s.dot(axis) < 0) axis = -axis;
    return th * axis;
}

// Per-member 3D co-rotational element (initial invariants + recovered state).
struct CrBeam3D {
    int      e  = -1;
    MemberId id = 0;
    int      ni = -1, nj = -1;
    int      gmap[12] = { 0,0,0,0,0,0,0,0,0,0,0,0 };  // global DOF of [uI(3),thI(3),uJ(3),thJ(3)]
    real     L0 = 0;
    Mat3     E0col = Mat3::Identity();                // initial element frame (COLUMNS = local x,y,z in global)
    real     EA = 0, EIy = 0, EIz = 0, GJ = 0;
    // recovered natural state (member-force output)
    real     Nax = 0, MzI = 0, MzJ = 0, MyI = 0, MyJ = 0, Tx = 0, Ln = 0;
};

// Translation of a node from the working displacement vector.
inline Vec3e nodeTrans(const std::vector<real>& u, int n) {
    return Vec3e(u[(size_t)gdof(n, Ux)], u[(size_t)gdof(n, Uy)], u[(size_t)gdof(n, Uz)]);
}

// Build the 3D internal force fe(12) and tangent Ke(12x12) at the current translations u + nodal triads
// Rnode. Stores the recovered natural forces on the element. The tangent is the material part T^T*Kl*T plus
// the axial geometric stiffness Ksigma1 (the strict axial term that reduces to S9's (N/Ln)q q^T, verified
// by the F50c P-Delta degeneration oracle). The full spin/moment corrections (Ksigma2/3, OpenSees
// getKs2Matrix) are NOT added: the NR converged solution depends only on f_int (which is complete), the
// tangent only on convergence speed, and the main term already converges the elastica alpha=1..10 (S9c).
void crCompute3D(CrBeam3D& b, const FrameModel& M, const std::vector<real>& u,
                 const std::vector<Mat3>& Rnode,
                 Eigen::Matrix<real, 12, 1>& fe, Eigen::Matrix<real, 12, 12>& Ke) {
    const Vec3e Xi = toE(M.nodes[(size_t)b.ni].pos) + nodeTrans(u, b.ni);
    const Vec3e Xj = toE(M.nodes[(size_t)b.nj].pos) + nodeTrans(u, b.nj);
    const Vec3e d  = Xj - Xi;
    const real  Ln = std::max<real>(1e-300, d.norm());
    const Vec3e e1 = d / Ln;

    // current nodal triads (columns = current section axes in global)
    const Mat3 qI = Rnode[(size_t)b.ni] * b.E0col;
    const Mat3 qJ = Rnode[(size_t)b.nj] * b.E0col;

    // mean triad (geodesic midpoint of qI, qJ), then Gram-Schmidt its axes onto the chord e1
    const Mat3 qm = expSO3(0.5 * logSO3(qJ * qI.transpose())) * qI;
    const Vec3e r1 = qm.col(0), r2 = qm.col(1), r3 = qm.col(2);
    const real denom = 1.0 + r1.dot(e1);                       // ~2 normally; ~0 only at a 180 flip
    const Vec3e e1r1 = e1 + r1;
    Vec3e e2 = r2 - (r2.dot(e1) / denom) * e1r1; e2.normalize();
    Vec3e e3 = r3 - (r3.dot(e1) / denom) * e1r1; e3.normalize();
    Mat3 E; E.col(0) = e1; E.col(1) = e2; E.col(2) = e3;

    // local deformational rotations: nodal triad rotation relative to the element frame
    const Vec3e thI = logSO3(E.transpose() * qI);              // [twist_x, bend_y, bend_z] at I
    const Vec3e thJ = logSO3(E.transpose() * qJ);
    const real ubar = Ln - b.L0;

    // natural forces pb = Kl * v, v = [ubar, tz_I, tz_J, ty_I, ty_J, tx]
    const real eaL = b.EA / b.L0, eiyL = b.EIy / b.L0, eizL = b.EIz / b.L0, gjL = b.GJ / b.L0;
    const real tzI = thI(2), tzJ = thJ(2), tyI = thI(1), tyJ = thJ(1);
    const real tx  = thJ(0) - thI(0);                          // relative twist about the beam axis
    const real N   = eaL * ubar;                               // axial (tension positive)
    const real MzI = eizL * (4.0 * tzI + 2.0 * tzJ);
    const real MzJ = eizL * (2.0 * tzI + 4.0 * tzJ);
    const real MyI = eiyL * (4.0 * tyI + 2.0 * tyJ);
    const real MyJ = eiyL * (2.0 * tyI + 4.0 * tyJ);
    const real Tx  = gjL * tx;
    b.Nax = N; b.MzI = MzI; b.MzJ = MzJ; b.MyI = MyI; b.MyJ = MyJ; b.Tx = Tx; b.Ln = Ln;

    Eigen::Matrix<real, 6, 1> pb; pb << N, MzI, MzJ, MyI, MyJ, Tx;

    // transformation T (6x12): rows = basic [ubar, tz_I, tz_J, ty_I, ty_J, tx], cols = [uI,thI,uJ,thJ].
    // Derived from the first-order frame variation; reduces EXACTLY to the S9 planar B in the XY/Z limit.
    Eigen::Matrix<real, 6, 12> T; T.setZero();
    auto setBlk = [&](int row, int c0, const Vec3e& v) { for (int i = 0; i < 3; ++i) T(row, c0 + i) = v(i); };
    const Vec3e e2L = e2 / Ln, e3L = e3 / Ln;
    setBlk(0, 0, -e1);   setBlk(0, 6,  e1);                                  // ubar
    setBlk(1, 0,  e2L);  setBlk(1, 3,  e3);  setBlk(1, 6, -e2L);            // tz_I
    setBlk(2, 0,  e2L);  setBlk(2, 9,  e3);  setBlk(2, 6, -e2L);            // tz_J
    setBlk(3, 0, -e3L);  setBlk(3, 3,  e2);  setBlk(3, 6,  e3L);            // ty_I
    setBlk(4, 0, -e3L);  setBlk(4, 9,  e2);  setBlk(4, 6,  e3L);            // ty_J
    setBlk(5, 3, -e1);   setBlk(5, 9,  e1);                                  // tx (relative twist)

    fe = T.transpose() * pb;

    // material tangent
    Eigen::Matrix<real, 6, 6> Kl; Kl.setZero();
    Kl(0, 0) = eaL;
    Kl(1, 1) = 4.0 * eizL; Kl(1, 2) = 2.0 * eizL; Kl(2, 1) = 2.0 * eizL; Kl(2, 2) = 4.0 * eizL;
    Kl(3, 3) = 4.0 * eiyL; Kl(3, 4) = 2.0 * eiyL; Kl(4, 3) = 2.0 * eiyL; Kl(4, 4) = 4.0 * eiyL;
    Kl(5, 5) = gjL;
    Ke = T.transpose() * Kl * T;

    // Ksigma1: axial geometric stiffness (chord stress stiffening), A = (I - e1 e1^T)/Ln. MUST add.
    const Mat3 A = (Mat3::Identity() - e1 * e1.transpose()) / Ln;
    const Mat3 NA = N * A;
    Ke.block(0, 0, 3, 3) += NA; Ke.block(0, 6, 3, 3) -= NA;
    Ke.block(6, 0, 3, 3) -= NA; Ke.block(6, 6, 3, 3) += NA;
}

}  // namespace

CorotationalResult runCorotational(const FrameModel& model, const CorotationalOptions& opts) {
    CorotationalResult R;
    const int N = model.dofCount();

    // A rejected / failed run still returns a WELL-FORMED SolveResult (zero-filled to the model DOF count,
    // member forces id-tagged) so a consumer that reads finalState.u never indexes an empty vector.
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
    // S9c: member UDLs (-> equivalent nodal loads) and prescribed support displacements (-> lambda-ramped
    // Dirichlet BC) are now SUPPORTED; the guards below still reject what CR genuinely cannot do.
    for (const auto& h : model.hinges)
        for (const auto& mem : model.members)
            if (mem.id == h.member && mem.active)
                return reject("co-rotational does not honour formed plastic hinges (S10, after S9b)");
    // tension-only / member-end releases are not honoured in CR (would be computed at full stiffness)
    for (const auto& mem : model.members) {
        if (!mem.active) continue;
        if (mem.tensionOnly) return reject("co-rotational does not honour tension-only members");
        for (int k = 0; k < 12; ++k) if (mem.release[(size_t)k])
            return reject("co-rotational does not honour member-end releases");
    }
    // NOTE: S9b removes the S9 out-of-plane-z guard -- members may now have any spatial orientation.

    if (N <= 0) return reject("empty model");

    // --- free-DOF map ---
    std::vector<int> fmap((size_t)N, -1);
    int nf = 0;
    for (size_t k = 0; k < model.nodes.size(); ++k)
        for (int d = 0; d < 6; ++d)
            if (!model.nodes[k].fixed[(size_t)d]) fmap[(size_t)gdof((int)k, d)] = nf++;
    if (nf == 0) return reject("fully constrained model");

    // --- build 3D CR elements from active members ---
    std::vector<CrBeam3D> elems;
    elems.reserve(model.members.size());
    for (size_t e = 0; e < model.members.size(); ++e) {
        const Member& m = model.members[e];
        if (!m.active) continue;
        const int ni = model.nodeIndex(m.i), nj = model.nodeIndex(m.j);
        if (ni < 0 || nj < 0) continue;
        CrBeam3D b;
        b.e = (int)e; b.id = m.id; b.ni = ni; b.nj = nj;
        for (int d = 0; d < 6; ++d) { b.gmap[d] = gdof(ni, d); b.gmap[6 + d] = gdof(nj, d); }
        const Vec3 pi = model.nodes[(size_t)ni].pos, pj = model.nodes[(size_t)nj].pos;
        b.L0 = norm(pj - pi);
        b.E0col = localAxes(pi, pj, m.refVec).transpose();     // rows->cols: columns are local axes in global
        const Material& mat = model.materials[(size_t)m.matIdx];
        const Section&  sec = model.sections[(size_t)m.secIdx];
        b.EA = mat.E * sec.A; b.EIy = mat.E * sec.Iy; b.EIz = mat.E * sec.Iz; b.GJ = mat.G * sec.J;
        elems.push_back(b);
    }

    // --- external nodal force vector ---
    VecX Fext = VecX::Zero(N);
    for (const auto& nl : model.nodalLoads) {
        const int ni = model.nodeIndex(nl.node);
        if (ni < 0) continue;
        for (int d = 0; d < 6; ++d) Fext((Eigen::Index)gdof(ni, d)) += nl.comp[(size_t)d];
    }
    // S9c: member UDL -> INITIAL-configuration equivalent nodal loads (same consistent Qf as BeamColumnElement),
    // added to F_ext and scaled by lambda. Small-strain equivalent; deformation-tracking (follower) UDL -> later.
    for (const auto& udl : model.memberUDLs) {
        for (size_t e = 0; e < model.members.size(); ++e) {
            const Member& m = model.members[e];
            if (!m.active || m.id != udl.member) continue;
            const int ni = model.nodeIndex(m.i), nj = model.nodeIndex(m.j);
            if (ni < 0 || nj < 0) continue;
            const Vec3 pi = model.nodes[(size_t)ni].pos, pj = model.nodes[(size_t)nj].pos;
            const real L = norm(pj - pi);
            const real wx = -udl.w_local.x, wy = -udl.w_local.y, wz = -udl.w_local.z;
            Vec12 Qf = Vec12::Zero();
            Qf(0) += wx * L / 2;          Qf(6)  += wx * L / 2;
            Qf(1) += wy * L / 2;          Qf(7)  += wy * L / 2;
            Qf(5) += wy * L * L / 12.0;   Qf(11) += -wy * L * L / 12.0;
            Qf(2) += wz * L / 2;          Qf(8)  += wz * L / 2;
            Qf(4) += -wz * L * L / 12.0;  Qf(10) += wz * L * L / 12.0;
            const Vec12 peq = -transform12(localAxes(pi, pj, m.refVec)).transpose() * Qf;   // global equivalent load
            for (int d = 0; d < 6; ++d) { Fext((Eigen::Index)gdof(ni, d)) += peq(d); Fext((Eigen::Index)gdof(nj, d)) += peq(6 + d); }
        }
    }
    VecX Fext_f = VecX::Zero(nf);
    for (int g = 0; g < N; ++g) if (fmap[(size_t)g] >= 0) Fext_f(fmap[(size_t)g]) = Fext((Eigen::Index)g);

    // --- shared assembly / increment-application / recovery (load-stepping AND arc-length use these) ---
    // Internal force (always analytical) + tangent Kff. useFD -> numerical CONSISTENT tangent
    // (Kff[:,j] = d f_int/d u_j by forward difference -> quadratic NR), else analytical T^T Kl T + Ksigma1.
    auto assemble = [&](const std::vector<real>& U, const std::vector<Mat3>& Rn, VecX& fint, SpMat& Kff, bool useFD) {
        fint = VecX::Zero(N);
        std::vector<Triplet> tf; tf.reserve(elems.size() * 144);
        for (auto& el : elems) {
            Eigen::Matrix<real, 12, 1> fe; Eigen::Matrix<real, 12, 12> Ke;
            crCompute3D(el, model, U, Rn, fe, Ke);
            for (int a = 0; a < 12; ++a) {
                fint((Eigen::Index)el.gmap[a]) += fe(a);
                if (useFD) continue;
                const int fr = fmap[(size_t)el.gmap[a]];
                if (fr < 0) continue;
                for (int c = 0; c < 12; ++c) { const int fc = fmap[(size_t)el.gmap[c]]; if (fc >= 0) tf.emplace_back(fr, fc, Ke(a, c)); }
            }
        }
        if (useFD) {                                   // numerical consistent tangent (oracle / small models)
            VecX f0 = VecX::Zero(nf);
            for (int g = 0; g < N; ++g) if (fmap[(size_t)g] >= 0) f0(fmap[(size_t)g]) = fint((Eigen::Index)g);
            const real eps = 1e-7;
            for (size_t k = 0; k < model.nodes.size(); ++k) {
                for (int d = 0; d < 6; ++d) {
                    const int g = gdof((int)k, d), col = fmap[(size_t)g];
                    if (col < 0) continue;
                    std::vector<real> Uc = U; std::vector<Mat3> Rc = Rn;
                    if (d < 3) Uc[(size_t)g] += eps;
                    else { Vec3e ax(0, 0, 0); ax(d - 3) = eps; Rc[(size_t)k] = expSO3(ax) * Rc[(size_t)k]; }
                    VecX fp = VecX::Zero(N);
                    for (auto& el : elems) {
                        Eigen::Matrix<real, 12, 1> fe; Eigen::Matrix<real, 12, 12> Ke;
                        crCompute3D(el, model, Uc, Rc, fe, Ke);
                        for (int a = 0; a < 12; ++a) fp((Eigen::Index)el.gmap[a]) += fe(a);
                    }
                    for (int r = 0; r < N; ++r) { const int rr = fmap[(size_t)r]; if (rr >= 0) { const real v = (fp((Eigen::Index)r) - f0(rr)) / eps; if (v != 0.0) tf.emplace_back(rr, col, v); } }
                }
            }
        }
        Kff = SpMat((Eigen::Index)nf, (Eigen::Index)nf);
        Kff.setFromTriplets(tf.begin(), tf.end());
        Kff.makeCompressed();
        if (useFD) { SpMat Kt = SpMat(Kff.transpose()); Kff = (Kff + Kt) * 0.5; Kff.makeCompressed(); }   // symmetrise for LDLT
    };
    auto applyInc = [&](const VecX& duf, std::vector<real>& U, std::vector<Mat3>& Rn) {
        for (size_t k = 0; k < model.nodes.size(); ++k) {
            for (int d = 0; d < 3; ++d) { const int g = gdof((int)k, d); if (fmap[(size_t)g] >= 0) U[(size_t)g] += duf(fmap[(size_t)g]); }
            Vec3e dth; for (int d = 0; d < 3; ++d) { const int g = gdof((int)k, Rx + d); dth(d) = (fmap[(size_t)g] >= 0) ? duf(fmap[(size_t)g]) : 0.0; }
            if (dth.squaredNorm() > 0) Rn[(size_t)k] = expSO3(dth) * Rn[(size_t)k];
        }
    };
    // S9c: prescribed support displacements as a lambda-ramped Dirichlet BC (prescribed DOFs are FIXED, not
    // in the free set; NR never touches them). Translation set directly; rotation via R_node = expSO3(lambda*presc).
    bool anyPrescribed = false;
    for (const auto& nd : model.nodes) for (int d = 0; d < 6; ++d) if (nd.fixed[(size_t)d] && nd.prescribed[(size_t)d] != 0.0) anyPrescribed = true;
    auto applyPrescribed = [&](real lam, std::vector<real>& U, std::vector<Mat3>& Rn) {
        for (size_t k = 0; k < model.nodes.size(); ++k) {
            const Node& nd = model.nodes[k];
            for (int d = 0; d < 3; ++d) if (nd.fixed[(size_t)d] && nd.prescribed[(size_t)d] != 0.0) U[(size_t)gdof((int)k, d)] = lam * nd.prescribed[(size_t)d];
            Vec3e pth(0, 0, 0); bool anyRot = false;
            for (int d = 0; d < 3; ++d) if (nd.fixed[(size_t)(Rx + d)] && nd.prescribed[(size_t)(Rx + d)] != 0.0) { pth(d) = lam * nd.prescribed[(size_t)(Rx + d)]; anyRot = true; }
            if (anyRot) Rn[(size_t)k] = expSO3(pth);
        }
    };
    auto recoverState = [&](const std::vector<real>& U, const std::vector<Mat3>& Rn, real lam) {
        SolveResult& SR = R.finalState;
        SR.u.assign((size_t)N, 0.0);
        SR.reactions.assign((size_t)N, 0.0);
        for (size_t k = 0; k < model.nodes.size(); ++k) {
            SR.u[(size_t)gdof((int)k, Ux)] = U[(size_t)gdof((int)k, Ux)];
            SR.u[(size_t)gdof((int)k, Uy)] = U[(size_t)gdof((int)k, Uy)];
            SR.u[(size_t)gdof((int)k, Uz)] = U[(size_t)gdof((int)k, Uz)];
            const Vec3e w = logSO3(Rn[(size_t)k]);
            SR.u[(size_t)gdof((int)k, Rx)] = w(0);
            SR.u[(size_t)gdof((int)k, Ry)] = w(1);
            SR.u[(size_t)gdof((int)k, Rz)] = w(2);
        }
        VecX fint = VecX::Zero(N);
        for (auto& el : elems) {
            Eigen::Matrix<real, 12, 1> fe; Eigen::Matrix<real, 12, 12> Ke;
            crCompute3D(el, model, U, Rn, fe, Ke);
            for (int a = 0; a < 12; ++a) fint((Eigen::Index)el.gmap[a]) += fe(a);
        }
        for (int g = 0; g < N; ++g) SR.reactions[(size_t)g] = fint((Eigen::Index)g) - lam * Fext((Eigen::Index)g);
        SR.memberForces.resize(model.members.size());
        SR.shellForces.clear();
        for (size_t e = 0; e < model.members.size(); ++e) SR.memberForces[e].member = model.members[e].id;
        for (const CrBeam3D& b : elems) {
            const real Ln = std::max<real>(1e-300, b.Ln);
            const real Vy = (b.MzI + b.MzJ) / Ln;
            const real Vz = -(b.MyI + b.MyJ) / Ln;
            MemberForcePair& mp = SR.memberForces[(size_t)b.e];
            mp.member = b.id;
            mp.endI = MemberEndForces{ -b.Nax,  Vy,  Vz, -b.Tx, b.MyI, b.MzI };
            mp.endJ = MemberEndForces{ -b.Nax, -Vy, -Vz,  b.Tx, b.MyJ, b.MzJ };
        }
        SR.singular = false;
    };

    // --- S9c: Crisfield cylindrical arc-length (snap-through; load factor lambda is an unknown) ---
    if (opts.useArcLength) {
        int mdof = opts.monitorDof;
        if (mdof < 0) {                          // auto: last node's first free translational DOF
            const int ln = (int)model.nodes.size() - 1;
            for (int d : { Uy, Uz, Ux }) { const int g = gdof(ln, d); if (fmap[(size_t)g] >= 0) { mdof = g; break; } }
            if (mdof < 0) mdof = 0;
        }
        std::vector<real> uu((size_t)N, 0.0);
        std::vector<Mat3> RR(model.nodes.size(), Mat3::Identity());
        real lambda = 0.0, Dl = opts.arcLength;
        VecX dutPrev = VecX::Zero(nf);
        bool firstStep = true;
        int totalIters = 0;
        for (int step = 1; step <= std::max(1, opts.arcSteps); ++step) {
            VecX fint; SpMat Kff;
            assemble(uu, RR, fint, Kff, opts.consistentTangent);
            LDLTSolver ldlt; ldlt.compute(Kff);
            VecX dut = ldlt.solve(Fext_f);       // reference tangent
            if (ldlt.info() != Eigen::Success || !dut.allFinite()) {
                R.diverged = true; R.finalState.singular = true;
                R.finalState.diagnostic = "arc-length predictor tangent singular (sharp limit point)";
                R.loadStepsCompleted = step - 1;
                break;
            }
            if (Dl <= 0) Dl = dut.norm() / std::max(1, opts.loadSteps);   // auto arc-length from first tangent
            const real dutn = std::sqrt(std::max<real>(1e-300, dut.squaredNorm()));
            const real sgn = firstStep ? 1.0 : (dutPrev.dot(dut) >= 0 ? 1.0 : -1.0);   // GSP sign (avoid back-track)
            const real Dlam0 = sgn * Dl / dutn;
            VecX Du = Dlam0 * dut;
            applyInc(Du, uu, RR);
            real lam = lambda + Dlam0, lastRel = 0;
            if (anyPrescribed) applyPrescribed(lam, uu, RR);
            bool conv = false;
            for (int it = 1; it <= std::max(1, opts.maxIter); ++it) {   // corrector NR on the constraint sphere
                assemble(uu, RR, fint, Kff, opts.consistentTangent);
                VecX rf = VecX::Zero(nf);
                for (int g = 0; g < N; ++g) if (fmap[(size_t)g] >= 0) rf(fmap[(size_t)g]) = lam * Fext((Eigen::Index)g) - fint((Eigen::Index)g);
                ldlt.compute(Kff);
                if (ldlt.info() != Eigen::Success) break;
                const VecX dub = ldlt.solve(rf);
                const VecX dt  = ldlt.solve(Fext_f);
                if (!dub.allFinite() || !dt.allFinite()) break;
                const real a = dt.squaredNorm();                 // cylindrical constraint: a dl^2 + b dl + c = 0
                const VecX Dub = Du + dub;
                const real b = 2.0 * Dub.dot(dt);
                const real c = Dub.squaredNorm() - Dl * Dl;
                real disc = b * b - 4.0 * a * c; if (disc < 0) disc = 0;
                const real sq = std::sqrt(disc);
                const real dl1 = (-b + sq) / (2.0 * std::max<real>(1e-300, a));
                const real dl2 = (-b - sq) / (2.0 * std::max<real>(1e-300, a));
                const real cos1 = Du.dot(Dub + dl1 * dt), cos2 = Du.dot(Dub + dl2 * dt);
                const real dlam = (cos1 >= cos2) ? dl1 : dl2;     // root closest to the incoming direction
                const VecX du = dub + dlam * dt;
                applyInc(du, uu, RR);
                Du += du; lam += dlam;
                if (anyPrescribed) applyPrescribed(lam, uu, RR);
                ++totalIters;
                const real Frn = std::max<real>(1e-300, (lam * Fext_f).norm());
                lastRel = rf.norm() / Frn;
                if (lastRel < opts.tolR || du.norm() < opts.tolU * std::max<real>(1.0, Du.norm())) { conv = true; break; }
            }
            lambda = lam;
            dutPrev = dut; firstStep = false;   // GSP sign uses the dot of CONSECUTIVE raw tangent solutions
                                                // dut^(i).dut^(i-1): it flips negative exactly when Kt^-1*F
                                                // reverses at the limit point, so the path descends past it.
                                                // (Storing sgn*dut instead makes this case climb monotonically
                                                // and miss the snap-through -- verified: it breaks F52.)
            R.totalIterations = totalIters;
            R.lastResidual = lastRel;
            if (!conv) { R.converged = false; R.loadStepsCompleted = step - 1; recoverState(uu, RR, lambda); return R; }
            R.pathLambda.push_back(lambda);
            R.pathDisp.push_back(uu[(size_t)mdof]);
            R.loadStepsCompleted = step;
        }
        if (!R.diverged) R.converged = true;
        recoverState(uu, RR, lambda);
        return R;
    }

    // --- driver state: translations in u (accumulated), rotations in Rnode (SO(3), spatial update) ---
    std::vector<real> u((size_t)N, 0.0);
    std::vector<Mat3> Rnode(model.nodes.size(), Mat3::Identity());
    const int steps = std::max(1, opts.loadSteps);
    int totalIters = 0;
    real lastRel = 0;

    for (int s = 1; s <= steps; ++s) {
        const real lambda = real(s) / real(steps);
        if (anyPrescribed) applyPrescribed(lambda, u, Rnode);
        bool stepConverged = false;

        for (int it = 1; it <= std::max(1, opts.maxIter); ++it) {
            VecX fint; SpMat Kff;
            assemble(u, Rnode, fint, Kff, opts.consistentTangent);

            VecX rf = VecX::Zero(nf);
            for (int g = 0; g < N; ++g)
                if (fmap[(size_t)g] >= 0) rf(fmap[(size_t)g]) = lambda * Fext((Eigen::Index)g) - fint((Eigen::Index)g);

            LDLTSolver ldlt; ldlt.compute(Kff);
            if (ldlt.info() != Eigen::Success || !ldltPositiveDefinite(ldlt, opts.solve.pivotTol)) {
                R.diverged = true;
                R.finalState.singular = true;
                R.finalState.diagnostic = "co-rotational tangent not positive-definite (limit point; snap-through needs arc-length)";
                R.loadStepsCompleted = s - 1; R.totalIterations = totalIters; R.lastResidual = lastRel;
                fillPartial(R, u);
                return R;
            }
            const VecX duf = ldlt.solve(rf);
            if (ldlt.info() != Eigen::Success || !duf.allFinite()) {
                R.diverged = true;
                R.finalState.singular = true;
                R.finalState.diagnostic = "co-rotational solve produced non-finite increment";
                R.loadStepsCompleted = s - 1; R.totalIterations = totalIters;
                fillPartial(R, u);
                return R;
            }

            // scatter increment; accumulate translations, spatial-update nodal rotations
            std::vector<real> du((size_t)N, 0.0);
            for (int g = 0; g < N; ++g) if (fmap[(size_t)g] >= 0) du[(size_t)g] = duf(fmap[(size_t)g]);
            real un = 0;
            for (size_t k = 0; k < model.nodes.size(); ++k) {
                u[(size_t)gdof((int)k, Ux)] += du[(size_t)gdof((int)k, Ux)];
                u[(size_t)gdof((int)k, Uy)] += du[(size_t)gdof((int)k, Uy)];
                u[(size_t)gdof((int)k, Uz)] += du[(size_t)gdof((int)k, Uz)];
                const Vec3e dth(du[(size_t)gdof((int)k, Rx)], du[(size_t)gdof((int)k, Ry)], du[(size_t)gdof((int)k, Rz)]);
                if (dth.squaredNorm() > 0) Rnode[k] = expSO3(dth) * Rnode[k];   // spatial left-multiply
                for (int d = 0; d < 3; ++d) { const real t = u[(size_t)gdof((int)k, d)]; un += t * t; }
            }
            ++totalIters;

            const real dn = duf.norm();
            const real Frn = (lambda * Fext_f).norm();
            const real relR = rf.norm() / std::max<real>(1e-300, Frn);
            const real relU = dn / std::max<real>(1.0, std::sqrt(un));
            lastRel = relR;
            if (relR < opts.tolR || relU < opts.tolU) { stepConverged = true; break; }
        }

        R.loadStepsCompleted = stepConverged ? s : (s - 1);
        if (!stepConverged) { R.converged = false; R.totalIterations = totalIters; R.lastResidual = lastRel; break; }
        if (s == steps) R.converged = true;
    }
    R.totalIterations = totalIters;
    R.lastResidual = lastRel;

    // --- recover finalState at lambda = 1 (shared helper; load-stepping and arc-length use the same) ---
    recoverState(u, Rnode, 1.0);
    return R;
}

}  // namespace frame
