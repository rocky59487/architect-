#pragma once
#include "FrameCore/FrameModel.h"
#include "FrameCore/SolveResult.h"
#include "FrameCore/SolveOptions.h"

namespace frame {

// POD options for runCorotational (S9-S9c). Default load-control uses `loadSteps` equal increments
// of a single proportional load factor lambda: 0 -> 1 applied to nodal and equivalent member loads.
// POD/std only -- no Eigen, no UE -- so this stays on the pure public boundary (same as PDeltaOptions).
struct CorotationalOptions {
    int  loadSteps = 10;     // equal lambda increments 0->1 (more steps reach larger displacement
                             //   before an NR step over-rotates; elastica alpha=10 wants ~10).
    int  maxIter   = 50;     // NR iterations per load step (cap; non-convergence at a NON-limit point
                             //   means raise loadSteps, not a bug).
    real tolR      = 1e-9;   // convergence: ||residual_free|| / ||lambda*F_ext_free||   (force residual)
    real tolU      = 1e-12;  // convergence: ||du_free|| / max(1,||u_free||)             (displacement)
    SolveOptions solve;      // pivotTol passthrough (useTimoshenko reserved; CR uses Euler-Bernoulli beams).
    // --- S9c (all default off -> S9b behaviour bit-for-bit; none enter modelFingerprint) ---
    bool useArcLength      = false;  // Crisfield cylindrical arc-length (snap-through); ignores loadSteps
    real arcLength         = 0;      // arc-length increment Dl (0 -> auto from the first tangent / loadSteps)
    int  arcSteps          = 50;     // max arc-length increments
    int  monitorDof        = -1;     // global DOF recorded in result.pathDisp (-1 -> tip translation auto)
    bool consistentTangent = false;  // numerical (finite-difference) consistent tangent -> quadratic NR
};

// POD result. converged / diverged are mutually exclusive on a healthy run; both false means a step
// hit maxIter without a limit-point verdict (raise loadSteps/maxIter). diverged = a limit point was
// reached (tangent K_T not positive-definite under load control -> snap-through needs arc-length).
struct CorotationalResult {
    bool        converged = false;
    bool        diverged  = false;   // load-control: limit point reached (snap-through). Arc-length tracks past it.
    int         loadStepsCompleted = 0;   // lambda increments fully equilibrated (== loadSteps on success)
    int         totalIterations    = 0;   // summed NR iterations across all completed steps
    real        lastResidual       = 0;   // last ||residual_free|| / ||lambda*F_ext_free||
    SolveResult finalState;          // load-control: state at lambda=1; arc-length: final tracked lambda;
                                     // singular flag forwarded on a failed / limit-point / invalid run.
    // --- S9c arc-length path (empty unless useArcLength): load factor + monitored displacement per step ---
    std::vector<real> pathLambda;    // load factor lambda at each completed arc-length increment
    std::vector<real> pathDisp;      // displacement of the monitor DOF at each increment (snap-through curve)
};

// Co-rotational large-displacement analysis (geometric nonlinearity). Newton-Raphson, load-controlled.
//
// SCOPE (S9c, honest): 3D GENERAL co-rotational beam of ARBITRARY spatial orientation -- torsion + biaxial
// bending (section Iy/Iz/J, Euler-Bernoulli) + finite SO(3) rotation. Each member co-rotates with its
// CURRENT chord and section triad; the local strain stays small (small-strain large-rotation CR, NOT a
// geometrically-exact Reissner beam; the two agree under small strain). Because 3D finite rotations do not
// commute, each node carries a rotation matrix R_node in SO(3) (initial I) updated by a SPATIAL increment
// R_node <- exp(skew(dtheta))*R_node after each NR step (avoids the total-rotation-vector 2.pi singularity).
// In the planar limit (members in XY, rotation about Z) it reduces to the S9 planar formulation -- the
// planar elastica / rigid-rotation / P-Delta degeneration oracles still pass.
//
// TANGENT: default = T^T Kl T + Ksigma1 (strict axial geometric term; already converges elastica 1..10).
// opts.consistentTangent -> numerical FD consistent tangent K_t = d f_int/d u (quadratic NR). The analytical
// spin/moment corrections (OpenSees Ksigma2/3) are not added (FD covers consistency; analytical -> future).
//
// LOADS (S9c): nodal forces, member UDLs (-> equivalent nodal loads), and prescribed support displacements
// (lambda-ramped Dirichlet BC) are all SUPPORTED. Formed plastic hinges / tension-only members / member-end
// releases are REJECTED (singular + diagnostic, never silently wrong; S10 covers N-M hinges).
//
// SNAP-THROUGH (S9c): opts.useArcLength -> Crisfield cylindrical arc-length tracks the load-displacement
// path PAST a limit point (result.pathLambda/pathDisp). Under load control a limit point is reported
// diverged (not tracked). Snap-back (load+displacement double reversal) needs spherical arc-length -> future.
//
// The caller's model is NEVER mutated (internal working copy; safe to call concurrently, same contract as
// runPDelta / runProgressiveCollapse). A model containing shells is rejected (beam-column only).
FRAMECORE_API CorotationalResult runCorotational(const FrameModel& model,
                                                 const CorotationalOptions& opts = {});

}  // namespace frame
