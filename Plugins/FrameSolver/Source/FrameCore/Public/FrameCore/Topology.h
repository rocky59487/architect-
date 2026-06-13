// Topology.h — S7 evolutionary hard-kill topology optimization (BESO) + N2 collapse-robustness.
//
// runBESO iteratively DEACTIVATES the lowest-sensitivity members (Member::active, the S1 removal
// hook — zero-intrusion) until a target volume fraction is reached, ranking members by their elastic
// STRAIN ENERGY (the standard BESO sensitivity alpha_e = 1/2 u_e^T K_e u_e, reconstructed from the
// recovered local end forces * end displacements; bit-exact strain energy for UDL-free members, so
// sum(alpha_e) == 1/2 * F.u to round-off — the energy-balance oracle).
//
// Three must-haves the WS-I prototype exposed, plus the literature's standard fix:
//   (1) two-step sensitivity history averaging (Huang & Xie)        -> damps 0/1 oscillation
//   (2) compliance-jump stop  (C_k > tol * C_{k-1})                  -> the prototype's tail blew up 52x
//   (3) mechanism guard via the N1 ReSolveSession capacitance test  -> removal that makes a mechanism
//                                                                       is rolled back (no fresh factor)
//   (4) compliance-BEST rollback -> output the pre-blow-up topology, not the final one (WS-I remedy A)
//
// N2 (opt-in collapse robustness): every k steps, probe single-member-removal scenarios with the
// progressive-collapse driver; if any scenario Collapses, roll back this step's removal and LOCK
// those members (protected, never removable). Prior art: fail-safe truss TO (Stolpe 2019, Zhu 2023)
// optimizes discrete members against single-member-removal (Jansen 2014 is the SIMP original). Here
// the evaluator is an LSP progressive-collapse driver (not LP/SDP), so this is NOT a global-optimum
// guarantee and must NOT be called "fail-safe".
//
// Honest scope: a heuristic (no global optimum; trusses' global optimum is LP/SDP ground structure);
// linear-elastic sensitivity; FIRST RELEASE IS ONE-WAY hard-kill (no bi-directional re-admission --
// a removed member has no internal forces to score); shell BESO not yet (ShellQuad::active is wired,
// but the Michell oracle is a truss). NOTE: TargetReached/converged means the VOLUME target was hit,
// NOT the Huang & Xie 2N-window compliance-stability criterion (not implemented -- inspect
// complianceHistory to judge stability). Pure POD/std public API (no Eigen, no UE).
#pragma once
#include "FrameCore/FrameModel.h"
#include "FrameCore/SolveResult.h"
#include "FrameCore/SolveOptions.h"
#include "FrameCore/Collapse.h"   // CollapseOptions for the N2 redundancy probe
#include <vector>

namespace frame {

// Element strain energy alpha_e = 1/2 u_e^T K_e u_e of member `memberIndex`, reconstructed from the
// recovered LOCAL end forces dotted with the LOCAL end displacements -- the BESO sensitivity number.
// For a UDL-free member this is the EXACT element strain energy, so the sum over all active members
// equals 1/2 * (external work F.u) to round-off (the energy-balance identity used as the F45 oracle).
// Weights split the four energy components (axial / bending / shear / torsion); all 1 = full energy.
// Also useful as an energy / efficiency field for visualization (C7). 0 for an out-of-range index.
//
// CAVEAT (UDL): for a member carrying a uniformly-distributed load, recover() returns end forces that
// include the fixed-end forces (Q = k_local u_local + Qf), so this computes 1/2 Q.u, NOT 1/2 u^T K u
// -- the energy-balance identity breaks and the value is an APPROXIMATE sensitivity. The BESO oracles
// and ground structures are nodal-load / UDL-free; on UDL models read this as a screen, not an energy.
FRAMECORE_API real memberStrainEnergy(const FrameModel& model, const SolveResult& result, int memberIndex,
                                      real wAxial = 1, real wBending = 1, real wShear = 1, real wTorsion = 1);

struct BESOOptions {
    real targetVolFrac = 0.5;        // target retained volume / initial volume (sum A_e*L_e, design set)
    real evolRate      = 0.02;       // per-step removal quota = evolRate * current volume (0.02~0.05)
    int  maxIter       = 200;        // iteration cap

    bool sensHistory          = true;  // (1) two-step sensitivity history average
    real complianceJumpTol    = 2.0;   // (2) stop when C_k > tol * C_{k-1} (single-step blow-up)
    bool complianceBestRollback = true; // (4) report bestActive (pre-blow-up) instead of finalActive

    // Per-component sensitivity weights (Karamba WTension/WCompr/WShear/WMoment analogue). All 1 =
    // full element strain energy. NOTE: only with all weights == 1 does sum(alpha) == 1/2 C hold
    // (the energy-balance oracle uses 1,1,1,1).
    real wAxial = 1, wBending = 1, wShear = 1, wTorsion = 1;

    // N2 collapse-robustness constraint (0 = OFF = plain BESO).
    int  redundancyCheckEvery = 0;     // probe every N steps; 0 disables
    int  redundancySamples    = 0;     // scenarios = the m highest-D/C active design members; 0 = all
    CollapseOptions redundancy;        // collapse options for the probe (dlf / removeThreshold passthrough)

    SolveOptions solve;                // threaded into baseline assembleAndFactor / ReSolve / collapse
};

enum class BESOStop {
    TargetReached,   // reached targetVolFrac cleanly (converged)
    ComplianceJump,  // single-step compliance blow-up tripped the guard (2)
    Stalled,         // no removable member left, or the N2 rollback state cycled (FNV-1a guard)
    Mechanism,       // a re-solve went singular unexpectedly
    MaxIter,         // iteration budget exhausted
    Invalid          // invalid model / options
};

struct BESOResult {
    std::vector<char>     finalActive;      // member active at stop (model.members order; 1/0)
    std::vector<char>     bestActive;       // compliance-best rollback topology (== finalActive if no blow-up)
    std::vector<real>     volFracHistory;   // retained volume fraction per iteration
    std::vector<real>     complianceHistory;// C = sum over nodal loads of F.u (external work) per iteration
    std::vector<MemberId> protectedMembers; // N2-locked members (ascending)
    int      bestIter   = -1;               // 1-based iteration of bestActive
    int      iterations = 0;
    BESOStop reason     = BESOStop::Invalid;
    bool     converged  = false;            // true iff reason == TargetReached
};

// Evolutionary hard-kill topology optimization. `designMembers` = member INDICES eligible for removal
// (empty = every active screenable member); members outside the set are kept (support links, etc.).
// The caller's model is NEVER mutated (internal working copy). See the header note for honest scope /
// the N2 fail-safe positioning.
FRAMECORE_API BESOResult runBESO(const FrameModel& model, const BESOOptions& opts,
                                 const std::vector<int>& designMembers = {});

}  // namespace frame
