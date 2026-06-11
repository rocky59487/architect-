#pragma once
#include "FrameCore/FrameModel.h"
#include "FrameCore/SolveResult.h"
#include "FrameCore/SolveOptions.h"
#include <vector>

namespace frame {

// POD options for runTensionOnly (S4). POD/std only — no Eigen, no UE. Members flagged
// Member::tensionOnly are iteratively deactivated while they carry compression (a slack
// cable / out-of-plane-buckling X-brace) and — unless allowReactivation is false —
// re-activated when their ends pull apart again.
struct TensionOnlyOptions {
    int  maxIter           = 32;     // active-set iteration cap before reporting non-convergence
    bool allowReactivation = true;   // false = MONOTONE (deactivate-only): conservative but guarantees
                                     // finite termination in <= n_TO steps (no flip-flop possible)
    real axialTol          = 0.0;    // deactivate when N > axialTol (N is compression-positive);
                                     // 0 = machine-level (any compression slackens the member)
    SolveOptions solve;              // threaded into the baseline assembleAndFactor / re-solves
};

// POD result. `cycled` means the active-set iteration revisited a state (a flip-flop): the driver
// then auto-restarts with allowReactivation=false, so `finalState` is the MONOTONE solution (a
// valid, order-dependent-conservative equilibrium that satisfies the tension-only sign constraint).
struct TensionOnlyResult {
    bool                  converged  = false;
    bool                  cycled     = false;  // cycle guard tripped -> result is the monotone fallback
    int                   iterations = 0;
    SolveResult           finalState;          // converged second-state result (singular flag forwarded)
    std::vector<MemberId> slack;               // tension-only members left DEACTIVATED at convergence
};

// Tension-only (cable / brace eliminator) analysis. Active-set fixed-point iteration: solve, then
// for each tension-only member deactivate it if it reads compression / re-activate it if its ends
// elongate, until no member flips. The inner re-solves REUSE the baseline factorization through a
// ReSolveSession (each flip is an exact rank-6 Woodbury update — the natural client of the S1
// ladder). The converged displacements equal a model with the slack members OMITTED bit-for-bit.
//
// Honest scope: an active-set fixed point has no general convergence guarantee (the LCP view), so
// the driver pairs a transition-hash cycle guard with a monotone (deactivate-only) fallback to
// guarantee FINITE TERMINATION. The criterion is the axial-force SIGN (no pre-tension / no slack
// length / no cable sag). A cycled solution is order-dependent and reported via `cycled`.
FRAMECORE_API TensionOnlyResult runTensionOnly(const FrameModel& model,
                                               const TensionOnlyOptions& opts = {});

}  // namespace frame
