#pragma once
//
// MITC4 Reissner-Mindlin flat-shell facet as an IElement (engine seam #1, the same
// hook BeamColumnElement uses). Each facet is membrane (plane stress + drilling Rz)
// + plate bending (w, Rx, Ry) with MITC4 assumed natural transverse shear strain
// (Bathe-Dvorkin) to defeat shear locking, assembled into a 24-DOF (4 nodes x 6)
// local stiffness and rotated into 3D. The solver treats it exactly like a beam:
// build -> prepare -> assemble -> addEquivalentNodalLoads -> recover.
//
#include "IElement.h"
#include "FrameEigen.h"

namespace frame {

class MITC4ShellElement final : public IElement {
public:
    explicit MITC4ShellElement(int shellIndex) : s_(shellIndex) {}

    int  localDof() const override { return 24; }
    bool prepare(const FrameModel& model, const SolveOptions& opts, std::string& why) override;
    void assemble(std::vector<Triplet>& trips) const override;
    void assembleMass(std::vector<Triplet>& trips) const override;
    void addEquivalentNodalLoads(VecX& F) const override;
    void recover(const VecX& u, SolveResult& R) const override;

private:
    using Mat24 = Eigen::Matrix<real, 24, 24>;
    using Vec24 = Eigen::Matrix<real, 24, 1>;

    int   s_   = -1;                 // shell index in model.shells
    int   id_  = 0;                  // shell id (result mapping)
    Mat24 kl_  = Mat24::Zero();      // local 24x24 stiffness (membrane + bending + drilling)
    Mat24 ml_  = Mat24::Zero();      // local 24x24 consistent mass (for modal analysis)
    Mat24 T_   = Mat24::Zero();      // transform: u_local = T u_global  (blockdiag R x8)
    int   dofs_[24] = { 0 };         // global DOF indices, 6 per corner node
    Vec24 Qf_  = Vec24::Zero();      // equivalent nodal loads from transverse pressure (local)

    // Facet geometry / material cache (filled in prepare, reused by recover).
    Mat3  R_   = Mat3::Identity();   // rows = facet local axes in global coords (R v_g = v_l)
    real  xl_[4] = { 0, 0, 0, 0 };   // corner x in the facet local 2D frame (mm)
    real  yl_[4] = { 0, 0, 0, 0 };   // corner y in the facet local 2D frame (mm)
    real  t_   = 0;                  // thickness
    real  E_   = 0, nu_ = 0, G_ = 0; // material
};

} // namespace frame
