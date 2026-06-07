#pragma once
//
// Internal element seam (#1 in the upgrade roadmap). This interface traffics in
// Eigen stiffness/triplets, so it lives in Private — FrameCore's PUBLIC API stays
// a pure POD boundary (FrameModel / SolveResult). The solver never hard-codes the
// 12x12 beam: it builds a transient std::vector<std::unique_ptr<IElement>> and
// iterates. BeamColumnElement is the first implementation; MITC shells (a later
// milestone) slot in behind the same interface with a different localDof().
//
#include "FrameEigen.h"
#include <vector>
#include <string>

namespace frame {

struct FrameModel;     // POD model (Public)
struct SolveOptions;   // POD options (Public)
struct SolveResult;    // POD result (Public)

struct IElement {
    virtual ~IElement() = default;

    // Number of LOCAL DOFs the element contributes (12 for a 3D beam-column).
    virtual int localDof() const = 0;

    // Build local stiffness + fixed-end forces and condense any released DOFs.
    // Must be called before assemble/addEquivalentNodalLoads/recover. Returns
    // false and fills `why` when the element is ill-posed (e.g. a released
    // sub-block is singular -> a free mechanism), so the solver can flag the
    // model singular with a clear diagnostic.
    virtual bool prepare(const FrameModel& model, const SolveOptions& opts, std::string& why) = 0;

    // Append this element's global-DOF stiffness contribution as triplets
    // (K = T^T kl T scattered to the element's global DOFs).
    virtual void assemble(std::vector<Triplet>& trips) const = 0;

    // Append this element's global-DOF CONSISTENT MASS contribution as triplets
    // (M = T^T ml T). Default empty so elements without mass (or analyses that don't
    // need it) are unaffected; BeamColumnElement and MITC4ShellElement override it.
    virtual void assembleMass(std::vector<Triplet>& trips) const { (void)trips; }

    // Add this element's equivalent nodal loads (P_equiv = -T^T Qf) into global F.
    virtual void addEquivalentNodalLoads(VecX& F) const = 0;

    // Recover element end forces from the full global displacement vector u and
    // store them in R (beam: into R.memberForces[memberIndex]).
    virtual void recover(const VecX& u, SolveResult& R) const = 0;
};

} // namespace frame
