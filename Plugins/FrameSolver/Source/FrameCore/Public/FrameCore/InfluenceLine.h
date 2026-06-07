#pragma once
#include "FrameCore/FrameSolver.h"
#include <vector>

namespace frame {

// Reaction INFLUENCE LINE: a unit downward (global -Z) load is marched over `loadNodes`,
// REUSING the factorization (one cheap back-substitution per position — the canonical
// same-K-many-RHS / moving-load application). Returns the reaction at (reactNode, reactDof)
// for each load position. The caller's `model` is taken by const-ref and is NEVER mutated —
// a local working copy carries the marched nodal load — so this is safe to call concurrently
// on one shared model + PreparedSystem (only the fingerprint-excluded nodal load varies).
FRAMECORE_API std::vector<real> reactionInfluenceLine(
    const PreparedSystem& prepared, const FrameModel& model,
    const std::vector<NodeId>& loadNodes, NodeId reactNode, int reactDof);

}  // namespace frame
