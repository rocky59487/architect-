#include "FrameCore/FrameModel.h"

namespace frame {

int FrameModel::nodeIndex(NodeId id) const {
    for (size_t k = 0; k < nodes.size(); ++k)
        if (nodes[k].id == id) return static_cast<int>(k);
    return -1;
}

bool FrameModel::validate(std::string& why) const {
    if (nodes.empty())   { why = "no nodes";   return false; }
    if (members.empty()) { why = "no members"; return false; }
    for (const auto& m : members) {
        const int ni = nodeIndex(m.i), nj = nodeIndex(m.j);
        if (ni < 0 || nj < 0)            { why = "member references missing node"; return false; }
        if (m.i == m.j)                  { why = "member endpoints identical (i == j)"; return false; }
        if (!m.mat || !m.sec)            { why = "member missing material/section"; return false; }
        if (m.mat->E <= 0 || m.mat->G <= 0) { why = "non-positive E or G"; return false; }
        if (m.sec->A <= 0 || m.sec->Iy <= 0 || m.sec->Iz <= 0 || m.sec->J <= 0)
                                         { why = "non-positive section property"; return false; }
        if (norm(nodes[nj].pos - nodes[ni].pos) <= 0) { why = "coincident member endpoints"; return false; }
    }
    // Loads must reference existing nodes/members, else the solver would silently drop them
    // (nodal loads are matched by node id, member UDLs by member id) and report a model that
    // is quietly missing load -- a dangerous "successful" solve.
    for (const auto& nl : nodalLoads) {
        if (nodeIndex(nl.node) < 0) { why = "nodal load references missing node"; return false; }
    }
    for (const auto& u : memberUDLs) {
        bool found = false;
        for (const auto& m : members) { if (m.id == u.member) { found = true; break; } }
        if (!found) { why = "member UDL references missing member"; return false; }
    }
    return true;
}

} // namespace frame
