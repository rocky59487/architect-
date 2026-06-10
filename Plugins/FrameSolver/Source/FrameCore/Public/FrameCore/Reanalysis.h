#pragma once
#include "FrameCore/FrameSolver.h"   // FRAMECORE_API, FrameModel, PreparedSystem, SolveOptions, SolveResult, MemberId, real
#include <memory>
#include <string>

namespace frame {

// Options for ReSolveSession (S1 ReSolve ladder). POD boundary; SolveOptions threads through to
// the baseline assembleAndFactor / rebaseline (e.g. enableReleases, useTimoshenko, pivotTol).
struct ReanalysisOptions {
    int  maxRank      = 96;      // Tier-1 cumulative rank cap (~16 beams); above -> Tier-3 rebaseline
    real pcgTol       = 1e-10;   // (Tier-2, reserved) relative residual
    int  pcgMaxIter   = 500;     // (Tier-2, reserved)
    bool allowTier2   = true;    // (Tier-2, reserved) when false, above maxRank goes straight to Tier-3
    real mechPivotTol = 1e-10;   // Tier-1 capacitance |piv|min/|piv|max below this -> mechanism
    SolveOptions solve;          // threaded into the baseline assembleAndFactor / rebaseline
};

// Per-solve diagnostics.
struct ReanalysisStats {
    int  tier        = 0;        // tier taken: 0 = no increment, 1 = Woodbury, 3 = rebaseline
    int  rank        = 0;        // current accumulated low-rank update rank
    int  pcgIters    = 0;        // (Tier-2, reserved)
    real relResidual = 0;        // (Tier-2, reserved)
    bool refactored  = false;    // this solve triggered a rebaseline (Tier-3)
    bool mechanism   = false;    // Tier-1 capacitance singular (the removed set forms a mechanism)
};

// Incremental re-analysis after element (member/shell) deactivate/restore, REUSING the baseline
// LDLᵀ factorization instead of a fresh assembleAndFactor — the interactive / collapse-driver path.
//
//   Tier-0  no change        -> baseline back-substitution.
//   Tier-1  rank <= maxRank   -> EXACT Woodbury low-rank update on the baseline factor; the
//                                capacitance matrix is singular  <=>  K' is singular  <=>  the
//                                removed set made a mechanism (detected from the factor, not topology).
//   Tier-3  rank >  maxRank    -> rebaseline (fresh assembleAndFactor on the current active set),
//                                which is always correct.
//
// (Tier-2 stale-LDLT preconditioned CG — the middle regime — is a follow-up commit; until it lands,
// `allowTier2`/`pcg*` are reserved and the ladder steps straight from Tier-1 to Tier-3 rebaseline.)
//
// The caller's model is never mutated (an internal working copy carries the active flags). Honest
// scope: SAME-TOPOLOGY increments only — the node set, support flags, and material/section VALUES
// must be unchanged (those need a fresh assembleAndFactor). Tier-1 is exact: vs a fresh
// assembleAndFactor+solveLoad it agrees to factorization round-off (~1e-13 on the test models).
class FRAMECORE_API ReSolveSession {
public:
    explicit ReSolveSession(const FrameModel& base, const ReanalysisOptions& opts = {});
    ~ReSolveSession();
    ReSolveSession(ReSolveSession&&) noexcept;
    ReSolveSession& operator=(ReSolveSession&&) noexcept;
    ReSolveSession(const ReSolveSession&) = delete;
    ReSolveSession& operator=(const ReSolveSession&) = delete;

    bool valid() const;                          // baseline validate() passed AND non-singular
    const std::string& diagnostic() const;       // reason when !valid()

    // Toggle an element's activity. Returns false for an unknown id; a no-op (same as current
    // state) returns true. Each state change appends one signed low-rank update to the ladder.
    bool setMemberActive(MemberId id, bool active);
    bool setShellActive(int shellId, bool active);

    // Solve the CURRENT active set. The result (u / reactions / member+shell forces) matches a fresh
    // assembleAndFactor+solveLoad to factorization round-off; on a mechanism, SolveResult.singular.
    SolveResult solve(ReanalysisStats* stats = nullptr);

    // Force Tier-3: rebuild the baseline factorization on the current active set and clear the ladder.
    void rebaseline();

private:
    struct Impl;
    std::unique_ptr<Impl> p_;
};

} // namespace frame
