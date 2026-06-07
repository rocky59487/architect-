#pragma once
#include "FrameCore/FrameSolver.h"
#include <vector>

namespace frame {

// Reaction INFLUENCE LINE: a unit downward (global -Z) load is marched over `loadNodes`,
// REUSING the factorization (one cheap back-substitution per position — the canonical
// same-K-many-RHS / moving-load application). Returns the reaction at (reactNode, reactDof)
// for each load position. `model` is mutated transiently (its nodal loads) and restored.
FRAMECORE_API std::vector<real> reactionInfluenceLine(
    const PreparedSystem& prepared, FrameModel& model,
    const std::vector<NodeId>& loadNodes, NodeId reactNode, int reactDof);

}  // namespace frame
