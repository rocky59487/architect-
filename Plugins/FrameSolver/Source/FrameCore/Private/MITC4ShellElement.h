#pragma once
//
// MITC4 Reissner-Mindlin flat-shell facet as an IElement (engine seam #1, the same
// hook BeamColumnElement uses). Each facet is membrane (plane stress + drilling Rz)
// + plate bending (w, Rx, Ry) with MITC4 assumed natural transverse shear strain
// (Bathe-Dvorkin) to defeat shear locking, assembled into a 24-DOF (4 nodes x 6)
// local stiffness and rotated into 3D. The solver treats it exactly like a beam:
// build -> prepare -> assemble -> addEquivalentNodalLoads -> recover.
//
// KNOWN ELEMENT-LEVEL TRAIT (documented, not a bug): the local 24x24 stiffness has the
// 6 true rigid-body zero modes PLUS one inherent low-energy (near-zero, non-rigid) plate-
// bending mode (a w-dominated mode, ~5e-10 relative to the largest eigenvalue, present even
// on a regular square). This is a standard MITC4 plate-bending trait, NOT distortion-induced
// and NOT the drilling DOF. Adjacent elements constrain it on assembly, so the drilling/patch/
// benchmark gates are all non-singular; recover() and all results are unaffected. Disclosed
// for honesty; eliminating it would need a different stabilization (high-risk, no present need).
//
// recover() samples stress resultants at the element CENTRE (single Gauss point 0,0): the
// returned {Mxx,Myy,Mxy,Qx,Qy,Nxx,Nyy,Nxy} are centre values (the element average for a
// linearly-varying field), NOT nodal/peak values. No nodal extrapolation is performed.
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

    // Audit-only: the assembled local 24x24 stiffness (valid after prepare()). Lets an
    // out-of-engine check verify the element-level eigen-spectrum (the 6 rigid-body zero
    // modes + the inherent low-energy plate-bending mode disclosed above). This is a
    // Private header, so it does NOT touch the POD public API.
    const Eigen::Matrix<real, 24, 24>& localKForAudit() const { return kl_; }

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

    // S8 opt-in membrane/plate formulation selectors (read from SolveOptions in prepare()).
    // Both false -> original MITC4 (bilinear Q4 membrane + assumed-shear plate), bit-for-bit.
    bool  useQM6_ = false;           // 8a: QM6 incompatible-mode membrane
    bool  useDKQ_ = false;           // 8b: DKQ discrete-Kirchhoff thin-plate bending
    // Note: QM6 needs NO recover cache. B_inc(0,0)=0 (dP1/dxi=-2xi, dP2/deta=-2eta vanish at
    // the centre), and ShellElementForces.N is a CENTRE value, so the centre membrane strain
    // is Bm(0,0)*d with zero bubble contribution -- the existing recover path is already
    // consistent and simply benefits from the QM6-improved displacement solution. (A future
    // per-corner membrane output WOULD need the condensed bubble DOF; documented if added.)
    // DKQ recover branches on useDKQ_ (B_DKQ plate moments, Qx=Qy=0).
};

} // namespace frame
