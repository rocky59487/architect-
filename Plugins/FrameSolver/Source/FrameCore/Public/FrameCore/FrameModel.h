#pragma once
#include "FrameCore/Node.h"
#include "FrameCore/Member.h"
#include "FrameCore/Material.h"
#include "FrameCore/Section.h"
#include "FrameCore/Load.h"
#include <vector>
#include <string>

namespace frame {

// Owns the whole structural model. Materials/Sections are referenced by pointer
// from Members, so the caller must keep them alive for the model's lifetime
// (the fixtures store them inside the model via these vectors).
struct FrameModel {
    std::vector<Node>            nodes;
    std::vector<Member>          members;
    std::vector<Material>        materials;   // optional storage to keep refs alive
    std::vector<Section>         sections;    // optional storage to keep refs alive
    std::vector<NodalLoad>       nodalLoads;
    std::vector<MemberUDL>       memberUDLs;

    int dofCount() const { return DOF_PER_NODE * static_cast<int>(nodes.size()); }

    // Index into nodes[] for a given NodeId, or -1 if not found.
    int nodeIndex(NodeId id) const;

    // Structural sanity (ids resolve, non-null mat/sec, positive A/I/E/L). On
    // failure returns false and fills `why`.
    FRAMECORE_API bool validate(std::string& why) const;
};

} // namespace frame
