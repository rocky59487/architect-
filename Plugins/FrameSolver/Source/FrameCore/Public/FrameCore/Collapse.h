#pragma once
#include "FrameCore/FrameModel.h"
#include "FrameCore/SolveOptions.h"
#include "FrameCore/Connectivity.h"
#include "FrameCore/ISectionStrength.h"
#include <string>
#include <vector>

namespace frame {

// C2 progressive-collapse driver options.
struct CollapseOptions {
    // Dynamic load factor (GSA LSP sudden-removal amplification), applied ONCE to all FORCE
    // loads in the working copy (nodalLoads, memberUDLs, shellPressures). Prescribed
    // displacement VALUES are NOT scaled (they are kinematic, not forces). Set 1.0 to run the
    // baked loads as-is.
    real dlf = 2.0;

    // The governing element is removed only while its D/C exceeds this; at or below it the
    // run terminates Stable. This is the screening D/C of ElasticAllowable (allowable-stress
    // ratio), NOT a code acceptance check.
    real removeThreshold = 1.0;

    // Total recorded steps INCLUDING the step-0 baseline (>= 1). Exhausting it while events
    // are still occurring terminates MaxSteps.
    int maxSteps = 256;

    // Damage scenario applied at step 0 (GSA-style "suddenly remove these members"). May be
    // empty: step 0 is then just the baseline solve. Unknown ids -> Invalid.
    std::vector<MemberId> initialRemovals;

    SolveOptions solve;   // pivotTol / enableReleases / useTimoshenko passthrough
};

enum class CollapseOutcome {
    Stable,     // every screened active element has D/C <= removeThreshold
    Collapsed,  // the grounded remainder went singular (mechanism) after fragment cleanup,
                // or no active element remains grounded. The engine does NOT distinguish
                // local from global mechanism (see diagnostic).
    MaxSteps,   // step budget exhausted while events were still occurring
    Invalid     // invalid model / options (see diagnostic)
};

// One driver step. The event lists are what was APPLIED at the START of this step (step 0
// carries initialRemovals and may be empty; every later step removes exactly one element).
// The state fields describe the model AFTER removal -> fragment cleanup -> re-solve.
struct CollapseStep {
    int step = 0;

    std::vector<MemberId> removedMembers;   // applied at the start of this step
    FailMode mode = FailMode::None;         // governing mode that selected this step's event
                                            // (None for the scenario-imposed step 0)
    real triggerRatio = 0;                  // the D/C that condemned it (0 at step 0)

    // Fragments cut loose by this step's removal: each is wholly deactivated, its nodes are
    // temporarily pinned at zero (they read u = 0 from this step on) and their nodal loads
    // leave the model with the debris. Hand these to the physics layer.
    std::vector<FragmentCluster> detached;

    bool solved = false;       // false on the terminal step of a Collapsed run
    real maxDC = 0;            // worst screening D/C over active members after the re-solve
    real safetyFactor = 0;     // 1/maxDC (+inf when maxDC==0; 0 when nothing screenable),
                               // mirroring DemandSummary semantics
    real pivotMargin = 0;      // criticality margin of this step's factorization (C4)
    std::vector<real> u;       // 6N displacement snapshot for UE replay (~6N*8 bytes/step)
};

struct CollapseHistory {
    CollapseOutcome outcome = CollapseOutcome::Invalid;
    std::string diagnostic;
    std::vector<CollapseStep> steps;
};

// C2 progressive-collapse driver: a GSA-style Linear Static Procedure run as sequential
// linear analysis. Per step: apply the pending removal -> connectivity cleanup (detached
// debris is deactivated wholesale, its nodes pinned, its loads dropped; bare free nodes are
// pinned) -> fresh assembleAndFactor + solveLoad -> screen D/C -> remove the governing
// member while it exceeds removeThreshold. Terminates Stable / Collapsed / MaxSteps.
//
// The caller's model is NEVER mutated (internal working copy; same contract as
// reactionInfluenceLine -- safe to call concurrently). Loads are whatever is baked in the
// model: call addSelfWeight()/compose the collapse combination first.
//
// Honest boundaries (LSP-grade): linear elastic between events; no inertia or dynamics
// beyond the scalar dlf; no plastic redistribution, membrane or catenary action (those make
// real structures both shed and pick up load nonlinearly -- expect conservative collapse
// extents, literature places LSP at roughly +/-30%); members are screened and removed,
// shells participate in stiffness/mass/debris but have no failure criterion yet (stage 3d).
// Removal order is deterministic: worst D/C, ties broken by smallest member id.
FRAMECORE_API CollapseHistory runProgressiveCollapse(const FrameModel& model,
                                                     const CollapseOptions& opts = {});

} // namespace frame
