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
    std::vector<int>      initialShellRemovals;   // same, for shell facet ids (stage 3d)

    // Event-to-event plastic hinges (stage 4b). When true, a HINGE-CAPABLE member (fy > 0 and
    // Zy/Zz > 0) fails in bending DUCTILELY: its bending no longer triggers removal; instead a
    // hinge forms when |M| >= Mp = fy*Z on an end/axis (the driver inserts the PlasticHinge AND
    // its node-side moment, see Hinge.h), and the run goes on until a hinge mechanism (LDLT
    // singular -> Collapsed) or quiescence. Its BRITTLE modes shrink to the SEPARABLE ratios --
    // pure axial N/A, shear, torsion (the combined-fibre sM +/- sN allowable screen stays
    // reporting-only, since its bending content is ductile here); members without fy/Z fall
    // back to full brittle removal, shells stay brittle (von Mises). Honest boundary:
    // sequential linear analysis, NOT true elastoplasticity -- no hinge unloading/reversal,
    // uniaxial Mp, no N-M interaction.
    bool plasticHinges = false;

    // N-M interaction for the plastic hinges above (stage 4b / S10). Only meaningful when
    // plasticHinges == true. When set, each end's hinge threshold and the residual moment it
    // carries use the AXIALLY-REDUCED plastic moment Mp_eff(N) = fy*Z * max(0,1-(N/Ny)^2),
    // Ny = fy*A (reducedPlasticMoment, NMInteraction.h) -- a hinge forms earlier under axial
    // load. Default false -> fixed Mp = fy*Z, bit-for-bit identical to the stage-4b behaviour
    // (N is fed as 0 -> no reduction). Honest boundary: RECTANGULAR is exact (the textbook
    // plastic N-M envelope, EC3 EN1993-1-1 sec 6.2.9; AISC H1.1 is a more conservative bilinear
    // design check, NOT this envelope), CIRCULAR is conservative; Mp_eff is FROZEN at formation
    // (the released end then recovers M = 0, so later steps never re-reduce it); UNIAXIAL N
    // reduction only, NOT My-Mz biaxial coupling; sequential linear analysis, NOT true
    // elastoplasticity (no unloading, no N-M tangent coupling).
    bool nmInteraction = false;

    SolveOptions solve;   // pivotTol / enableReleases / useTimoshenko passthrough
};

enum class CollapseOutcome {
    Stable,     // no event triggers: every brittle ratio is at/below removeThreshold and (in
                // hinge mode) every |M|/Mp < 1. NOTE in hinge mode the REPORTED maxDC may sit
                // above the threshold while Stable: bending of a hinge-capable member is
                // ductile and waits for Mp (the allowable screen sits BELOW fy by design).
    Collapsed,  // the grounded remainder went singular (mechanism -- including a hinge
                // mechanism in hinge mode) after fragment cleanup, or no active element
                // remains grounded. The engine does NOT distinguish local from global
                // mechanism (see diagnostic).
    MaxSteps,   // step budget exhausted while events were still occurring
    Invalid     // invalid model / options (see diagnostic)
};

// A hinge formed by the driver (stage 4b): member, released local dof (4/5/10/11), and the
// SIGNED residual Mp it carries from formation on (local end-force convention).
struct CollapseHingeEvent {
    MemberId member = 0;
    int      dof    = 0;
    real     Mp     = 0;
};

// One driver step. The event lists are what was APPLIED at the START of this step (step 0
// carries initialRemovals and may be empty; every later step removes exactly one element).
// The state fields describe the model AFTER removal -> fragment cleanup -> re-solve.
struct CollapseStep {
    int step = 0;

    std::vector<MemberId> removedMembers;   // applied at the start of this step
    std::vector<int>      removedShells;    // shell facets removed (mode == ShellVonMises)
    std::vector<CollapseHingeEvent> formedHinges;   // hinges formed (mode == Bending, 4b)
    FailMode mode = FailMode::None;         // governing mode that selected this step's event
                                            // (None for the scenario-imposed step 0)
    real triggerRatio = 0;                  // the ratio that condemned it: D/C for removals,
                                            // |M|/Mp for hinges (0 at step 0)

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
// extents, literature places LSP at roughly +/-30%). Members are screened by the section
// screen (checkSection), shells by the surface von Mises screen (checkShellSurface, stage
// 3d) -- both elastic screening grades, not code checks. Removal order is deterministic:
// worst D/C, ties broken member-before-shell, then smallest id.
FRAMECORE_API CollapseHistory runProgressiveCollapse(const FrameModel& model,
                                                     const CollapseOptions& opts = {});

} // namespace frame
