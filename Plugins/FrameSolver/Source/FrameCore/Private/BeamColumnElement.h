#pragma once
//
// The 3D Euler-Bernoulli (optionally Timoshenko) beam-column as an IElement.
// Wraps the existing free functions (localStiffness12 / transform12 /
// condenseReleases / end-force recovery) so the refactor that routed the solver
// through IElement is bit-for-bit identical to the old hard-coded loop.
//
#include "IElement.h"
#include "FrameEigen.h"

namespace frame {

class BeamColumnElement final : public IElement {
public:
    explicit BeamColumnElement(int memberIndex) : e_(memberIndex) {}

    int  localDof() const override { return 12; }
    bool prepare(const FrameModel& model, const SolveOptions& opts, std::string& why) override;
    void assemble(std::vector<Triplet>& trips) const override;
    void assembleMass(std::vector<Triplet>& trips) const override;
    void addEquivalentNodalLoads(VecX& F) const override;
    void recover(const VecX& u, SolveResult& R) const override;

private:
    int      e_  = -1;                 // member index in model.members
    MemberId id_ = 0;                  // member id (for the result mapping)
    Mat12    kl_ = Mat12::Zero();      // local stiffness (possibly release-condensed)
    Mat12    ml_ = Mat12::Zero();      // local consistent mass (for modal analysis)
    Mat12    T_  = Mat12::Zero();      // transform: u_local = T u_global
    real     L_  = 0;
    int      dofs_[12] = { 0 };        // global DOF indices [node_i 6][node_j 6]
    Vec12    Qf_ = Vec12::Zero();      // fixed-end forces (possibly release-condensed)
};

} // namespace frame
