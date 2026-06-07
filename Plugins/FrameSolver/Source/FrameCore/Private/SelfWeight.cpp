#include "FrameCore/SelfWeight.h"
#include "ElementStiffness.h"   // localAxes
#include "FrameEigen.h"

namespace frame {

// kg/m^3 -> tonne/mm^3 (see SelfWeight.h unit bridge).
static constexpr real kRhoConv = 1.0e-12;

void addSelfWeight(FrameModel& m, real g) {
    // ---- beams: gravity as a LOCAL uniformly distributed load (so sloped/vertical
    // members are correct). w_global = (0,0,-rho*A*g*1e-12); rotate into local axes. ----
    for (const auto& mem : m.members) {
        if (mem.matIdx < 0 || mem.matIdx >= static_cast<int>(m.materials.size()) ||
            mem.secIdx < 0 || mem.secIdx >= static_cast<int>(m.sections.size())) continue;
        const real wg = m.materials[static_cast<size_t>(mem.matIdx)].rho
                      * m.sections[static_cast<size_t>(mem.secIdx)].A * g * kRhoConv;   // N/mm (downward)
        if (wg == 0.0) continue;
        const int ni = m.nodeIndex(mem.i), nj = m.nodeIndex(mem.j);
        if (ni < 0 || nj < 0) continue;
        const Mat3 R = localAxes(m.nodes[ni].pos, m.nodes[nj].pos, mem.refVec);  // rows = local axes
        const Vec3e wLocal = R * Vec3e(0.0, 0.0, -wg);              // R * w_global -> local
        MemberUDL u;
        u.member  = mem.id;
        u.w_local = Vec3(wLocal(0), wLocal(1), wLocal(2));
        m.memberUDLs.push_back(u);
    }

    // ---- shells: body load lumped to the four corner nodes as a global -Z nodal load.
    // (On a regular mesh the A/4 lumping equals the consistent integral of the bilinear
    // shape functions, so it matches a transverse ShellPressure of rho*t*g.) ----
    for (const auto& sh : m.shells) {
        if (sh.matIdx < 0 || sh.matIdx >= static_cast<int>(m.materials.size())) continue;
        int idx[4];
        bool ok = true;
        for (int k = 0; k < 4; ++k) { idx[k] = m.nodeIndex(sh.n[k]); if (idx[k] < 0) ok = false; }
        if (!ok) continue;
        const Vec3 p0 = m.nodes[idx[0]].pos, p1 = m.nodes[idx[1]].pos;
        const Vec3 p2 = m.nodes[idx[2]].pos, p3 = m.nodes[idx[3]].pos;
        const real area = 0.5 * norm(cross(p2 - p0, p3 - p1));
        const real W    = m.materials[static_cast<size_t>(sh.matIdx)].rho * sh.t * g * kRhoConv * area;  // total facet weight (N)
        if (W == 0.0) continue;
        const real perNode = -0.25 * W;                              // global -Z, split 4 ways
        for (int k = 0; k < 4; ++k) {
            NodalLoad nl;
            nl.node = sh.n[k];
            nl.comp[Uz] = perNode;
            m.nodalLoads.push_back(nl);
        }
    }
}

}  // namespace frame
