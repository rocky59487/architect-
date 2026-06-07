#include "FrameCore/FrameModel.h"

namespace frame {

int FrameModel::nodeIndex(NodeId id) const {
    for (size_t k = 0; k < nodes.size(); ++k)
        if (nodes[k].id == id) return static_cast<int>(k);
    return -1;
}

bool FrameModel::validate(std::string& why) const {
    if (nodes.empty())                  { why = "no nodes"; return false; }
    if (members.empty() && shells.empty()) { why = "no members or shells"; return false; }
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
    // Shell facets: 4 resolvable distinct nodes, valid material (needs nu for the
    // plane-stress/bending constitutive), positive thickness, non-degenerate quad.
    for (const auto& s : shells) {
        int idx[4];
        for (int k = 0; k < 4; ++k) {
            idx[k] = nodeIndex(s.n[k]);
            if (idx[k] < 0) { why = "shell references missing node"; return false; }
        }
        for (int a = 0; a < 4; ++a)
            for (int b = a + 1; b < 4; ++b)
                if (s.n[a] == s.n[b]) { why = "shell has duplicate corner nodes"; return false; }
        if (!s.mat)                          { why = "shell missing material"; return false; }
        if (s.mat->E <= 0 || s.mat->G <= 0)  { why = "shell non-positive E or G"; return false; }
        if (s.mat->nu < 0 || s.mat->nu >= 0.5) { why = "shell Poisson ratio out of [0,0.5)"; return false; }
        if (s.t <= 0)                        { why = "shell non-positive thickness"; return false; }
        const Vec3 nrm = cross(nodes[idx[2]].pos - nodes[idx[0]].pos,
                               nodes[idx[3]].pos - nodes[idx[1]].pos);
        if (norm(nrm) <= 0) { why = "degenerate (zero-area / collinear) shell quad"; return false; }
    }
    for (const auto& sp : shellPressures) {
        bool found = false;
        for (const auto& s : shells) { if (s.id == sp.shell) { found = true; break; } }
        if (!found) { why = "shell pressure references missing shell"; return false; }
    }
    return true;
}

} // namespace frame
