#pragma once
#include "FrameCore/FrameSolver.h"
#include "FrameCore/ModalResult.h"
#include <vector>
#include <string>

namespace frame {

struct ModalDynamicsOptions {
    real dt     = 0.0;    // time step (s)
    int  nSteps = 0;      // number of steps
    real zeta   = 0.05;   // modal (proportional) damping ratio
};

struct ModalTimeHistory {
    bool                           singular = false;
    std::string                    diagnostic;
    real                           dt = 0;
    std::vector<std::vector<real>> u;   // u[step] = 6N displacement at that time (step 0 = t=0)
};

// Modal-superposition TRANSIENT response to a STEP nodal load (the model's current nodalLoads
// applied suddenly at t=0 to a system at rest). Each modal coordinate is integrated by
// Newmark-beta (average acceleration: beta=1/4, gamma=1/2 — unconditionally stable); the
// physical response is u(t) = sum_n phi_n q_n(t). Cost is O(nModes) per step, so a handful of
// precomputed modes gives a real-time (UE5) sway / vibration response without re-factorizing.
FRAMECORE_API ModalTimeHistory solveModalStepResponse(
    const PreparedSystem& prepared, const FrameModel& model, const ModalResult& modes,
    const ModalDynamicsOptions& opts);

}  // namespace frame
