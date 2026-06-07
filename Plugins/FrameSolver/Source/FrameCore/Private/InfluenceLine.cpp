#include "FrameCore/InfluenceLine.h"

namespace frame {

std::vector<real> reactionInfluenceLine(const PreparedSystem& prepared, FrameModel& model,
                                        const std::vector<NodeId>& loadNodes,
                                        NodeId reactNode, int reactDof) {
    std::vector<real> il;
    il.reserve(loadNodes.size());
    const std::vector<NodalLoad> saved = model.nodalLoads;   // preserve caller's loads
    const int rIdx = model.nodeIndex(reactNode);
    for (NodeId ln : loadNodes) {
        model.nodalLoads.clear();
        NodalLoad nl; nl.node = ln; nl.comp[Uz] = -1.0;       // unit downward
        model.nodalLoads.push_back(nl);
        const SolveResult r = solveLoad(prepared, model);     // reuse the factorization
        il.push_back(rIdx >= 0 ? r.reaction(rIdx, reactDof) : 0.0);
    }
    model.nodalLoads = saved;
    return il;
}

}  // namespace frame
