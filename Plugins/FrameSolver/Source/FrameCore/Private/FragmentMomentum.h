#pragma once
//
// Fragment momentum handoff (S2 N4): build a detached fragment's consistent mass on its OWN node
// set and extract linear / angular momentum from the global velocity field at the detach instant,
// then fill FragmentCluster::vel / angVel for the physics (UE5 Chaos) handoff. Header-inline (the
// SparseEigsolver.h pattern) so BOTH the dynamic-collapse driver and the deep audit exercise the
// identical code. Eigen-carrying -> Private only.
//
// Honest scope: p (linear) and the TRANSVERSE angular momentum match the FragmentCluster closed
// form exactly under rigid motion; the OWN-AXIS angular momentum uses the cluster's slender-rod
// inertia (zero about a member's own axis), so the FE section polar term is dropped (negligible
// for slender members) -- the audit reports that gap rather than asserting it is zero.
//
#include "FrameCore/Connectivity.h"
#include "FrameCore/FrameModel.h"
#include "FrameCore/SolveOptions.h"
#include "FrameEigen.h"
#include "IElement.h"
#include "BeamColumnElement.h"
#include "MITC4ShellElement.h"
#include <algorithm>
#include <vector>
#include <string>

namespace frame {

// Fragment-local consistent mass (6*nn x 6*nn, nodes remapped to fc.nodes order). A node shared by
// two fragment elements is summed once (node-keyed remap), so multi-element fragments are NOT
// double-counted (the C1 trap).
inline MatX fragmentConsistentMass(const FragmentCluster& fc, const FrameModel& work) {
    const int nn = (int)fc.nodes.size();
    auto localOf = [&](NodeId id) -> int {
        const auto it = std::lower_bound(fc.nodes.begin(), fc.nodes.end(), id);   // fc.nodes ascending
        return (it != fc.nodes.end() && *it == id) ? (int)(it - fc.nodes.begin()) : -1;
    };
    MatX M = MatX::Zero(6 * nn, 6 * nn);
    auto stamp = [&](IElement& el) {
        std::vector<Triplet> gt; el.assembleMass(gt);
        for (const auto& t : gt) {
            const int gr = (int)t.row(), gc = (int)t.col();
            const int lr = localOf(work.nodes[(size_t)(gr / 6)].id);
            const int lc = localOf(work.nodes[(size_t)(gc / 6)].id);
            if (lr >= 0 && lc >= 0) M(6 * lr + gr % 6, 6 * lc + gc % 6) += t.value();
        }
    };
    std::string why;
    for (MemberId mid : fc.members) {
        const int e = work.memberIndex(mid);
        if (e >= 0) { BeamColumnElement el(e); if (el.prepare(work, SolveOptions{}, why)) stamp(el); }
    }
    for (int sid : fc.shells) {
        const int s = work.shellIndex(sid);
        if (s >= 0) { MITC4ShellElement el(s); if (el.prepare(work, SolveOptions{}, why)) stamp(el); }
    }
    return M;
}

// Linear (p) and angular (L, about fc.com, GLOBAL axes) momentum under the global velocity v_N.
inline void fragmentMomentum(const FragmentCluster& fc, const FrameModel& work, const VecX& v_N,
                             Vec3& p, Vec3& L) {
    const int nn = (int)fc.nodes.size();
    const MatX M = fragmentConsistentMass(fc, work);
    VecX vf = VecX::Zero(6 * nn);
    for (int a = 0; a < nn; ++a) {
        const int gni = work.nodeIndex(fc.nodes[(size_t)a]);
        for (int d = 0; d < 6; ++d) vf(6 * a + d) = v_N(gdof(gni, d));
    }
    const VecX Mvf = M * vf;
    auto transMode = [&](int axis) { VecX T = VecX::Zero(6 * nn); for (int a = 0; a < nn; ++a) T(6 * a + axis) = 1; return T; };
    auto rotMode = [&](const Vec3& ek) {
        VecX T = VecX::Zero(6 * nn);
        for (int a = 0; a < nn; ++a) {
            const Vec3 rel = work.nodes[(size_t)work.nodeIndex(fc.nodes[(size_t)a])].pos - fc.com;
            const Vec3 tr  = cross(ek, rel);
            T(6 * a + 0) = tr.x; T(6 * a + 1) = tr.y; T(6 * a + 2) = tr.z;
            T(6 * a + 3) = ek.x; T(6 * a + 4) = ek.y; T(6 * a + 5) = ek.z;
        }
        return T;
    };
    p = Vec3{ transMode(0).dot(Mvf), transMode(1).dot(Mvf), transMode(2).dot(Mvf) };
    L = Vec3{ rotMode({ 1, 0, 0 }).dot(Mvf), rotMode({ 0, 1, 0 }).dot(Mvf), rotMode({ 0, 0, 1 }).dot(Mvf) };
}

// Fill fc.vel = p/mass and fc.angVel = I^{-1} L (symmetric pseudo-inverse: a collinear fragment's
// zero own-axis inertia, whose own-axis L is ~0 too, resolves to zero angular velocity there).
// fc.mass == 0 -> leave vel/angVel at their zero defaults.
inline void fillFragmentVelocity(FragmentCluster& fc, const FrameModel& work, const VecX& v_N) {
    if (!(fc.mass > 0) || fc.nodes.empty()) return;
    Vec3 p, L; fragmentMomentum(fc, work, v_N, p, L);
    fc.vel = p * (real(1) / fc.mass);
    Mat3 I; I << fc.inertia[0], fc.inertia[3], fc.inertia[4],
                 fc.inertia[3], fc.inertia[1], fc.inertia[5],
                 fc.inertia[4], fc.inertia[5], fc.inertia[2];
    Eigen::SelfAdjointEigenSolver<Mat3> es(I);
    const Vec3e ev = es.eigenvalues();
    const Vec3e Ll = es.eigenvectors().transpose() * Vec3e(L.x, L.y, L.z);
    const real emax = ev.cwiseAbs().maxCoeff();
    Vec3e wl;
    for (int k = 0; k < 3; ++k) wl(k) = (ev(k) > real(1e-12) * emax) ? Ll(k) / ev(k) : real(0);
    const Vec3e w = es.eigenvectors() * wl;
    fc.angVel = Vec3(w(0), w(1), w(2));
}

// Fragment rigid-body kinetic energy 1/2 m|vel|^2 + 1/2 angVel^T I angVel.
inline real fragmentKE(const FragmentCluster& fc) {
    const real lin = real(0.5) * fc.mass * dot(fc.vel, fc.vel);
    const Vec3 w = fc.angVel;
    const real Iww = fc.inertia[0] * w.x * w.x + fc.inertia[1] * w.y * w.y + fc.inertia[2] * w.z * w.z
                   + 2 * (fc.inertia[3] * w.x * w.y + fc.inertia[4] * w.x * w.z + fc.inertia[5] * w.y * w.z);
    return lin + real(0.5) * Iww;
}

} // namespace frame
