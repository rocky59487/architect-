#include "FrameCore/InfluenceLine.h"

namespace frame {

std::vector<real> reactionInfluenceLine(const PreparedSystem& prepared, const FrameModel& model,
                                        const std::vector<NodeId>& loadNodes,
                                        NodeId reactNode, int reactDof) {
    std::vector<real> il;
    il.reserve(loadNodes.size());
    // Work on a local copy so the caller's model is never mutated (thread-safe: many threads
    // may march influence lines over one shared model + PreparedSystem). The copy keeps the
    // same structural fingerprint, so solveLoad's reuse guard still accepts it; only the
    // (fingerprint-excluded) nodal load varies per position.
    FrameModel work = model;
    const int rIdx = work.nodeIndex(reactNode);
    for (NodeId ln : loadNodes) {
        work.nodalLoads.clear();
        NodalLoad nl; nl.node = ln; nl.comp[Uz] = -1.0;       // unit downward
        work.nodalLoads.push_back(nl);
        const SolveResult r = solveLoad(prepared, work);      // reuse the factorization
        il.push_back(rIdx >= 0 ? r.reaction(rIdx, reactDof) : 0.0);
    }
    return il;
}

}  // namespace frame
