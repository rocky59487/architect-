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
// SCOPE (S9b, honest): 3D GENERAL co-rotational beam of ARBITRARY spatial orientation -- torsion + biaxial
// bending (section Iy/Iz/J, Euler-Bernoulli) + finite SO(3) rotation. Each member co-rotates with its
// CURRENT chord and section triad; the local strain stays small (small-strain large-rotation CR, NOT a
// geometrically-exact Reissner beam; the two agree under small strain). Because 3D finite rotations do not
// commute, each node carries a rotation matrix R_node in SO(3) (initial I) updated by a SPATIAL increment
// R_node <- exp(skew(dtheta))*R_node after each NR step (avoids the total-rotation-vector 2.pi singularity).
// In the planar limit (members in XY, rotation about Z) it reduces to the S9 planar formulation -- the
// planar elastica / rigid-rotation / P-Delta degeneration oracles still pass. The tangent is T^T Kl T +
// Ksigma1 (the strict axial geometric term); the full spin/moment corrections (Ksigma2/3) are not added
// (they only accelerate convergence, which is already reached for the elastica alpha=1..10; -> S9c).
//
// NO snap-through: a limit point is reported diverged (not tracked); arc-length (Riks/Crisfield) reserved.
// Nodal force loads ONLY: member UDLs / prescribed support displacements / formed plastic hinges /
// tension-only members / member-end releases are REJECTED (singular + diagnostic, never silently wrong).
//
// The caller's model is NEVER mutated (internal working copy; safe to call concurrently, same contract as
// runPDelta / runProgressiveCollapse). A model containing shells is rejected (beam-column only).
FRAMECORE_API CorotationalResult runCorotational(const FrameModel& model,
                                                 const CorotationalOptions& opts = {});

}  // namespace frame
