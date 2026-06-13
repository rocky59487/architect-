#pragma once
//
// Shared collapse-driver helpers used by BOTH the static (Collapse.cpp) and the dynamic
// (DynamicCollapse.cpp) progressive-collapse drivers. Header-only inline so the two TUs share
// ONE definition — each previously kept a file-local copy, which collided under a unity build
// (see FrameCore.Build.cs). POD/std only; no Eigen.
//
#include "FrameCore/FrameModel.h"
#include <algorithm>

namespace frame {

// Temporarily ground a node that no longer carries any active element (a debris node or a bare
// free node). Fixing every DOF at zero is mathematically inert -- nothing couples to it any more,
// so it adds no stiffness and draws no reaction -- but it removes the zero-pivot rows that would
// otherwise read as a spurious mechanism. Its nodal loads leave with it (the load belonged to
// what fell; leaking it into the grounded remainder would be wrong).
inline void pinNode(FrameModel& work, NodeId nid) {
    const int k = work.nodeIndex(nid);
    if (k < 0) return;
    for (int d = 0; d < DOF_PER_NODE; ++d) {
        work.nodes[(size_t)k].fixed[d] = true;
        work.nodes[(size_t)k].prescribed[d] = 0;
    }
    work.nodalLoads.erase(std::remove_if(work.nodalLoads.begin(), work.nodalLoads.end(),
                                         [nid](const NodalLoad& nl) { return nl.node == nid; }),
                          work.nodalLoads.end());
}

// True while any member or shell is still active (the structure has not fully detached).
inline bool anyActive(const FrameModel& work) {
    for (const Member& m : work.members)   if (m.active) return true;
    for (const ShellQuad& s : work.shells) if (s.active) return true;
    return false;
}

}  // namespace frame
