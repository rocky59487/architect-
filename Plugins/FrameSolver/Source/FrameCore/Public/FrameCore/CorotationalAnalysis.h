#pragma once
#include "FrameCore/FrameModel.h"
#include "FrameCore/SolveResult.h"
#include "FrameCore/SolveOptions.h"

namespace frame {

// POD options for runCorotational (S9). Load-controlled Newton-Raphson over `loadSteps` equal
// increments of a single proportional load factor lambda: 0 -> 1 applied to all nodal force loads.
// POD/std only -- no Eigen, no UE -- so this stays on the pure public boundary (same as PDeltaOptions).
struct CorotationalOptions {
    int  loadSteps = 10;     // equal lambda increments 0->1 (more steps reach larger displacement
                             //   before an NR step over-rotates; elastica alpha=10 wants ~10).
    int  maxIter   = 50;     // NR iterations per load step (cap; non-convergence at a NON-limit point
                             //   means raise loadSteps, not a bug).
    real tolR      = 1e-9;   // convergence: ||residual_free|| / ||lambda*F_ext_free||   (force residual)
    real tolU      = 1e-12;  // convergence: ||du_free|| / max(1,||u_free||)             (displacement)
    SolveOptions solve;      // pivotTol passthrough (useTimoshenko reserved; planar v1 is Euler-Bernoulli).
};

// POD result. converged / diverged are mutually exclusive on a healthy run; both false means a step
// hit maxIter without a limit-point verdict (raise loadSteps/maxIter). diverged = a limit point was
// reached (tangent K_T not positive-definite under load control -> snap-through needs arc-length).
struct CorotationalResult {
    bool        converged = false;
    bool        diverged  = false;   // limit point: K_T not positive-definite (snap-through; S9 stops)
    int         loadStepsCompleted = 0;   // lambda increments fully equilibrated (== loadSteps on success)
    int         totalIterations    = 0;   // summed NR iterations across all completed steps
    real        lastResidual       = 0;   // last ||residual_free|| / ||lambda*F_ext_free||
    SolveResult finalState;          // large-displacement state at lambda=1 (u, reactions, member forces);
                                     // singular flag forwarded on a failed / limit-point / invalid run.
};

// Co-rotational large-displacement analysis (geometric nonlinearity). Newton-Raphson, load-controlled.
//
// SCOPE (v1, honest): PLANAR co-rotational beam in the global XY plane (bending about global Z, using
// the section Iz, Euler-Bernoulli). Each member rotates with its CURRENT chord; the local strain stays
// small (small-strain large-rotation CR). Because all rotations are about the single global Z axis they
// commute, so the nodal rotation is accumulated directly in the Rz displacement DOF -- the spatial
// incremental-rotation-vector machinery (R_node, exp/log SO(3)) needed for genuine 3D rotation is NOT
// required here and is reserved for S9b (3D CR; OpenSees CorotCrdTransf3d formulas archived in
// docs/research/WS_F2_corot3d_opensees.md). The caller MUST restrain the out-of-plane DOFs (Uz, Rx, Ry)
// at every node (the planar fixtures do); an unrestrained out-of-plane DOF has no CR stiffness and the
// solve is reported singular.
//
// NO snap-through: a limit point is reported diverged (not tracked); arc-length (Riks/Crisfield) reserved.
// The main oracle (a transverse end-loaded cantilever's elastica) is monotone, so NR reaches alpha=1..10.
//
// The caller's model is NEVER mutated (internal working copy; safe to call concurrently, same contract as
// runPDelta / runProgressiveCollapse). A model containing shells is rejected (beam-column only).
FRAMECORE_API CorotationalResult runCorotational(const FrameModel& model,
                                                 const CorotationalOptions& opts = {});

}  // namespace frame
