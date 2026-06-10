#pragma once
#include "FrameCore/Node.h"
#include "FrameCore/Member.h"
#include "FrameCore/Material.h"
#include "FrameCore/Section.h"
#include "FrameCore/Load.h"
#include "FrameCore/Shell.h"
#include "FrameCore/Hinge.h"
#include <vector>
#include <string>

namespace frame {

// Owns the whole structural model. Materials/Sections are referenced by INDEX
// (matIdx/secIdx) from Members and ShellQuads — never by pointer — so growing any
// vector can never dangle a reference; validate() range-checks the indices.
struct FrameModel {
    std::vector<Node>            nodes;
    std::vector<Member>          members;
    std::vector<ShellQuad>       shells;        // MITC4 flat-shell facets (parallel to members)
    std::vector<Material>        materials;   // optional storage to keep refs alive
    std::vector<Section>         sections;    // optional storage to keep refs alive
    std::vector<NodalLoad>       nodalLoads;
    std::vector<MemberUDL>       memberUDLs;
    std::vector<ShellPressure>   shellPressures;
    std::vector<PlasticHinge>    hinges;        // formed plastic hinges (stage 4a) -- model
                                                // STATE, honoured regardless of enableReleases;
                                                // a hinge on an INACTIVE member is inert

    int dofCount() const { return DOF_PER_NODE * static_cast<int>(nodes.size()); }

    // Index into nodes[] for a given NodeId, or -1 if not found.
    int nodeIndex(NodeId id) const;

    // Structural sanity (ids resolve, non-null mat/sec, positive A/I/E/L). On
    // failure returns false and fills `why`.
    FRAMECORE_API bool validate(std::string& why) const;
};

} // namespace frame
