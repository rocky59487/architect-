#include "FrameCore/FrameModel.h"

#include <algorithm>
#include <cmath>

namespace frame {

int FrameModel::nodeIndex(NodeId id) const {
    for (size_t k = 0; k < nodes.size(); ++k)
        if (nodes[k].id == id) return static_cast<int>(k);
    return -1;
}

bool FrameModel::validate(std::string& why) const {
    if (nodes.empty())                  { why = "no nodes"; return false; }
    if (members.empty() && shells.empty()) { why = "no members or shells"; return false; }
    auto finiteVec3 = [](const Vec3& v) {
        return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
    };
    auto finiteDof6 = [](const std::array<real, 6>& v) {
        for (real x : v) if (!std::isfinite(x)) return false;
        return true;
    };
    for (const auto& n : nodes) {
        if (!finiteVec3(n.pos)) { why = "node has non-finite coordinate"; return false; }
        if (!finiteDof6(n.prescribed)) { why = "node has non-finite prescribed displacement"; return false; }
    }
    for (size_t a = 0; a < nodes.size(); ++a)
        for (size_t b = a + 1; b < nodes.size(); ++b)
            if (nodes[a].id == nodes[b].id) { why = "duplicate node id"; return false; }
    for (size_t a = 0; a < members.size(); ++a)
        for (size_t b = a + 1; b < members.size(); ++b)
            if (members[a].id == members[b].id) { why = "duplicate member id"; return false; }
    for (size_t a = 0; a < shells.size(); ++a)
        for (size_t b = a + 1; b < shells.size(); ++b)
            if (shells[a].id == shells[b].id) { why = "duplicate shell id"; return false; }
    for (const auto& m : members) {
        const int ni = nodeIndex(m.i), nj = nodeIndex(m.j);
        if (ni < 0 || nj < 0)            { why = "member references missing node"; return false; }
        if (m.i == m.j)                  { why = "member endpoints identical (i == j)"; return false; }
        if (!m.mat || !m.sec)            { why = "member missing material/section"; return false; }
        if (!finiteVec3(m.refVec))       { why = "member has non-finite reference vector"; return false; }
        if (!std::isfinite(m.mat->E) || !std::isfinite(m.mat->G) || !std::isfinite(m.mat->nu))
                                         { why = "member material has non-finite property"; return false; }
        if (m.mat->E <= 0 || m.mat->G <= 0) { why = "non-positive E or G"; return false; }
        if (!std::isfinite(m.sec->A) || !std::isfinite(m.sec->Iy) || !std::isfinite(m.sec->Iz) ||
            !std::isfinite(m.sec->J) || !std::isfinite(m.sec->cy) || !std::isfinite(m.sec->cz) ||
            !std::isfinite(m.sec->Asy) || !std::isfinite(m.sec->Asz))
                                         { why = "member section has non-finite property"; return false; }
        if (m.sec->A <= 0 || m.sec->Iy <= 0 || m.sec->Iz <= 0 || m.sec->J <= 0)
                                         { why = "non-positive section property"; return false; }
        if (norm(nodes[nj].pos - nodes[ni].pos) <= 0) { why = "coincident member endpoints"; return false; }
    }
    // Loads must reference existing nodes/members, else the solver would silently drop them
    // (nodal loads are matched by node id, member UDLs by member id) and report a model that
    // is quietly missing load -- a dangerous "successful" solve.
    for (const auto& nl : nodalLoads) {
        if (nodeIndex(nl.node) < 0) { why = "nodal load references missing node"; return false; }
        if (!finiteDof6(nl.comp)) { why = "nodal load has non-finite component"; return false; }
    }
    for (const auto& u : memberUDLs) {
        bool found = false;
        for (const auto& m : members) { if (m.id == u.member) { found = true; break; } }
        if (!found) { why = "member UDL references missing member"; return false; }
        if (!finiteVec3(u.w_local)) { why = "member UDL has non-finite component"; return false; }
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
        if (!std::isfinite(s.mat->E) || !std::isfinite(s.mat->G) || !std::isfinite(s.mat->nu))
                                             { why = "shell material has non-finite property"; return false; }
        if (s.mat->E <= 0 || s.mat->G <= 0)  { why = "shell non-positive E or G"; return false; }
        if (s.mat->nu < 0 || s.mat->nu >= 0.5) { why = "shell Poisson ratio out of [0,0.5)"; return false; }
        if (!std::isfinite(s.t))             { why = "shell non-finite thickness"; return false; }
        if (s.t <= 0)                        { why = "shell non-positive thickness"; return false; }
        const Vec3 nrm = cross(nodes[idx[2]].pos - nodes[idx[0]].pos,
                               nodes[idx[3]].pos - nodes[idx[1]].pos);
        const real area2 = norm(nrm);
        if (area2 <= 0) { why = "degenerate (zero-area / collinear) shell quad"; return false; }
        const Vec3 nhat = nrm * (1.0 / area2);
        real maxEdge = 0;
        for (int a = 0; a < 4; ++a)
            for (int b = a + 1; b < 4; ++b)
                maxEdge = std::max(maxEdge, norm(nodes[idx[b]].pos - nodes[idx[a]].pos));
        real maxWarp = 0;
        const Vec3 p0 = nodes[idx[0]].pos;
        for (int k = 0; k < 4; ++k)
            maxWarp = std::max(maxWarp, std::abs(dot(nodes[idx[k]].pos - p0, nhat)));
        if (maxWarp > 1.0e-6 * std::max(real(1), maxEdge)) {
            why = "non-coplanar shell quad (MITC4 flat facet requires planar corners)";
            return false;
        }
        Vec3 e1 = nodes[idx[1]].pos - nodes[idx[0]].pos;
        e1 = e1 - nhat * dot(e1, nhat);
        const real e1len = norm(e1);
        if (e1len <= 0) { why = "degenerate shell quad (edge parallel to normal)"; return false; }
        e1 = e1 * (1.0 / e1len);
        const Vec3 e2 = cross(nhat, e1);
        real x[4], y[4];
        for (int k = 0; k < 4; ++k) {
            const Vec3 r = nodes[idx[k]].pos - nodes[idx[0]].pos;
            x[k] = dot(r, e1);
            y[k] = dot(r, e2);
        }
        const real shapeTol = 1.0e-12 * std::max(real(1), maxEdge * maxEdge);
        for (int k = 0; k < 4; ++k) {
            const int a = k, b = (k + 1) % 4, c = (k + 2) % 4;
            const real ex1 = x[b] - x[a], ey1 = y[b] - y[a];
            const real ex2 = x[c] - x[b], ey2 = y[c] - y[b];
            const real turn = ex1 * ey2 - ey1 * ex2;
            if (turn <= shapeTol) {
                why = "invalid shell quad (concave, inverted, or inconsistent corner order)";
                return false;
            }
        }
        const real gp = 0.5773502691896258;
        const real detTol = shapeTol;
        for (const real xi : { -gp, gp })
            for (const real eta : { -gp, gp }) {
                const real dNxi[4]  = { -0.25 * (1 - eta),  0.25 * (1 - eta),
                                         0.25 * (1 + eta), -0.25 * (1 + eta) };
                const real dNeta[4] = { -0.25 * (1 - xi),  -0.25 * (1 + xi),
                                         0.25 * (1 + xi),   0.25 * (1 - xi) };
                real x_xi = 0, y_xi = 0, x_eta = 0, y_eta = 0;
                for (int k = 0; k < 4; ++k) {
                    x_xi  += dNxi[k]  * x[k]; y_xi  += dNxi[k]  * y[k];
                    x_eta += dNeta[k] * x[k]; y_eta += dNeta[k] * y[k];
                }
                const real detJ = x_xi * y_eta - y_xi * x_eta;
                if (detJ <= detTol) {
                    why = "invalid shell quad (non-positive Jacobian / concave or inverted)";
                    return false;
                }
            }
    }
    for (const auto& sp : shellPressures) {
        bool found = false;
        for (const auto& s : shells) { if (s.id == sp.shell) { found = true; break; } }
        if (!found) { why = "shell pressure references missing shell"; return false; }
        if (!std::isfinite(sp.p)) { why = "shell pressure has non-finite value"; return false; }
    }
    return true;
}

} // namespace frame
