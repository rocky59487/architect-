#include "BeamColumnElement.h"
#include "ElementStiffness.h"
#include "FrameCore/FrameModel.h"
#include "FrameCore/SolveOptions.h"
#include "FrameCore/SolveResult.h"

#include <string>

namespace frame {

bool BeamColumnElement::prepare(const FrameModel& model, const SolveOptions& opts, std::string& why) {
    const Member& mem = model.members[static_cast<size_t>(e_)];
    id_ = mem.id;

    const int ni = model.nodeIndex(mem.i), nj = model.nodeIndex(mem.j);
    const Vec3 pi = model.nodes[ni].pos, pj = model.nodes[nj].pos;
    L_ = norm(pj - pi);

    const Mat3 Rm = localAxes(pi, pj, mem.refVec);
    // Timoshenko (shear-flexible) only when explicitly enabled AND the section
    // carries shear areas; otherwise the Euler-Bernoulli element, bit-for-bit.
    if (opts.useTimoshenko && mem.sec->Asy > 0 && mem.sec->Asz > 0) {
        kl_ = localStiffness12T(mem.mat->E, mem.mat->G, mem.sec->A,
                                mem.sec->Iy, mem.sec->Iz, mem.sec->J, L_,
                                mem.sec->Asy, mem.sec->Asz);
    } else {
        kl_ = localStiffness12(mem.mat->E, mem.mat->G, mem.sec->A,
                               mem.sec->Iy, mem.sec->Iz, mem.sec->J, L_);
    }
    T_  = transform12(Rm);
    // consistent mass (rho kg/m^3 -> tonne/mm^3 via 1e-12). Not condensed by releases
    // (releases are a static device; modal analysis uses the unreleased element).
    ml_ = localMass12(mem.mat->rho * 1.0e-12, mem.sec->A, mem.sec->Iy, mem.sec->Iz, L_);
    for (int d = 0; d < 6; ++d) { dofs_[d] = gdof(ni, d); dofs_[6 + d] = gdof(nj, d); }

    // ---- fixed-end forces Qf (local) from this member's UDLs. Same per-element
    // accumulation order as the original monolithic loop, so the result is identical. ----
    Qf_ = Vec12::Zero();
    for (const auto& udl : model.memberUDLs) {
        if (udl.member != mem.id) continue;
        // Qf = -(consistent load vector): peq = -T^T*Qf applies w_local in its TRUE
        // direction, and Q = k*T*u + Qf is the textbook member end-force recovery.
        const real wx = -udl.w_local.x, wy = -udl.w_local.y, wz = -udl.w_local.z;
        Vec12 q = Vec12::Zero();
        // axial uniform (split half/half)
        q(0) += wx * L_ / 2; q(6) += wx * L_ / 2;
        // transverse in local y (couples v=1,7 ; rz=5,11)
        q(1) += wy * L_ / 2; q(7) += wy * L_ / 2;
        q(5) += wy * L_ * L_ / 12.0; q(11) += -wy * L_ * L_ / 12.0;
        // transverse in local z (couples w=2,8 ; ry=4,10)
        q(2) += wz * L_ / 2; q(8) += wz * L_ / 2;
        q(4) += -wz * L_ * L_ / 12.0; q(10) += wz * L_ * L_ / 12.0;
        Qf_ += q;
    }

    // ---- member-end release condensation (on kl AND Qf together) ----
    if (opts.enableReleases) {
        bool any = false;
        for (int k = 0; k < 12; ++k) if (mem.release[k]) { any = true; break; }
        if (any && !condenseReleases(kl_, Qf_, mem.release)) {
            why = "member-end release leaves a singular released sub-block "
                  "(free mechanism, e.g. both torsional ends released) at member "
                  + std::to_string(static_cast<long long>(mem.id));
            return false;
        }
    }
    return true;
}

void BeamColumnElement::assemble(std::vector<Triplet>& trips) const {
    const Mat12 kg = T_.transpose() * kl_ * T_;
    for (int a = 0; a < 12; ++a)
        for (int b = 0; b < 12; ++b)
            trips.emplace_back(dofs_[a], dofs_[b], kg(a, b));
}

void BeamColumnElement::assembleMass(std::vector<Triplet>& trips) const {
    const Mat12 mg = T_.transpose() * ml_ * T_;
    for (int a = 0; a < 12; ++a)
        for (int b = 0; b < 12; ++b)
            trips.emplace_back(dofs_[a], dofs_[b], mg(a, b));
}

void BeamColumnElement::addEquivalentNodalLoads(VecX& F) const {
    if (Qf_.isZero(0)) return;
    const Vec12 peq = -(T_.transpose() * Qf_);   // P_equiv = -T^T Qf
    for (int a = 0; a < 12; ++a) F(dofs_[a]) += peq(a);
}

void BeamColumnElement::recover(const VecX& u, SolveResult& R) const {
    Vec12 ue;
    for (int a = 0; a < 12; ++a) ue(a) = u(dofs_[a]);
    const Vec12 Q = kl_ * (T_ * ue) + Qf_;   // Q = k_local (T u_e) + Qf

    MemberForcePair& mp = R.memberForces[static_cast<size_t>(e_)];
    mp.member = id_;
    // end i : N compression-positive = +Q[0]
    mp.endI = { Q(0), Q(1), Q(2), Q(3), Q(4), Q(5) };
    // end j : N compression-positive = -Q[6]
    mp.endJ = { -Q(6), Q(7), Q(8), Q(9), Q(10), Q(11) };
}

} // namespace frame
