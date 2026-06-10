// exp_beso_truss.cpp — WS-I / N2 ground-structure BESO experiment (scratch, NOT engine).
//
// 12x6 grid ground structure (hor/vert/both diagonals), left edge clamped, mid-right
// point load. Hard-kill BESO on the engine's Member.active flag:
//   sensitivity alpha_e = N^2 L / (2EA)   (axial strain energy; slender members)
//   per iteration remove the lowest-alpha members up to an evolution-rate length quota,
//   guarded by analyzeConnectivity (the loaded node must stay grounded) and a singularity
//   rollback. Tracks compliance; dumps the final topology for qualitative Michell reading.

#include "research_common.h"
#include "FrameCore/Connectivity.h"
#include <algorithm>
#include <cstdio>

using namespace research;

namespace {

constexpr int NX = 12, NZ = 6;
constexpr real SP = 500.0;

int nid(int i, int k) { return k * (NX + 1) + i; }

FrameModel makeGround() {
    FrameModel m;
    m.materials.emplace_back(200000.0, 76923.07692307692, 7850.0);
    m.sections.push_back(squareSection(40.0));
    for (int k = 0; k <= NZ; ++k)
        for (int i = 0; i <= NX; ++i) {
            Node n(nid(i, k), i * SP, 0.0, k * SP);
            if (i == 0) n.fixAll();
            else n.fixed[Uy] = true;
            m.nodes.push_back(n);
        }
    auto add = [&](int a, int b) {
        Member mem(static_cast<int>(m.members.size()), a, b, 0, 0);
        mem.refVec = refVecFor(m.nodes[static_cast<size_t>(a)].pos, m.nodes[static_cast<size_t>(b)].pos);
        m.members.push_back(mem);
    };
    for (int k = 0; k <= NZ; ++k)
        for (int i = 0; i < NX; ++i) add(nid(i, k), nid(i + 1, k));
    for (int k = 0; k < NZ; ++k)
        for (int i = 0; i <= NX; ++i) add(nid(i, k), nid(i, k + 1));
    for (int k = 0; k < NZ; ++k)
        for (int i = 0; i < NX; ++i) { add(nid(i, k), nid(i + 1, k + 1)); add(nid(i + 1, k), nid(i, k + 1)); }
    NodalLoad l; l.node = nid(NX, NZ / 2); l.comp[Uz] = -1.0e5;
    m.nodalLoads.push_back(l);
    return m;
}

} // namespace

int main() {
    FrameModel work = makeGround();
    const int loadNode = nid(NX, NZ / 2);
    const real E = work.materials[0].E, A = work.sections[0].A;

    auto memberLen = [&](const Member& mem) {
        return norm(work.nodes[static_cast<size_t>(work.nodeIndex(mem.j))].pos -
                    work.nodes[static_cast<size_t>(work.nodeIndex(mem.i))].pos);
    };
    real totalLen0 = 0;
    for (const Member& mem : work.members) totalLen0 += memberLen(mem);

    const real targetFrac = 0.30, er = 0.04;
    for (int it = 1; it <= 80; ++it) {
        PreparedSystem ps = assembleAndFactor(work);
        if (ps.impl->singular) { std::printf("[beso] FATAL singular at iter %d\n", it); return 2; }
        SolveResult r = solveLoad(ps, work);
        const real compliance = -1.0e5 * r.u[static_cast<size_t>(gdof(loadNode, Uz))];

        real curLen = 0;
        std::vector<std::pair<real, int>> alpha;   // (sensitivity, member index), active only
        for (size_t e = 0; e < work.members.size(); ++e) {
            const Member& mem = work.members[e];
            if (!mem.active) continue;
            const real L = memberLen(mem);
            curLen += L;
            const real N = r.memberForces[e].endI.N;
            alpha.emplace_back(N * N * L / (2.0 * E * A), static_cast<int>(e));
        }
        const real volFrac = curLen / totalLen0;
        std::printf("[beso] iter=%d volFrac=%.4f compliance=%.6g active=%zu\n",
                    it, static_cast<double>(volFrac), static_cast<double>(compliance), alpha.size());
        if (volFrac <= targetFrac) break;

        std::sort(alpha.begin(), alpha.end());
        const real quota = std::min(er * curLen, curLen - targetFrac * totalLen0);
        real removed = 0;
        int removedCount = 0;
        std::vector<int> batch;
        for (const auto& [a, e] : alpha) {
            if (removed >= quota) break;
            Member& mem = work.members[static_cast<size_t>(e)];
            mem.active = false;
            const ConnectivityResult cr = analyzeConnectivity(work);
            bool bad = !cr.valid;
            for (NodeId ln : cr.looseNodes) if (ln == loadNode) bad = true;
            for (const auto& fc : cr.detached)
                for (NodeId n : fc.nodes) if (n == loadNode) bad = true;
            if (!bad) {
                // connectivity cannot see bending mechanisms — per-member factor check
                // (model is small; ~ms per check; the engine-level version of this guard
                //  would ride the N1 ladder's capacitance singularity test instead)
                PreparedSystem chk = assembleAndFactor(work);
                bad = chk.impl->singular;
            }
            if (bad) { mem.active = true; continue; }
            removed += memberLen(mem);
            ++removedCount;
            batch.push_back(e);
        }
        if (removedCount == 0) { std::printf("[beso] stalled (no removable member) at iter %d\n", it); break; }
    }

    // dump final topology for plotting
    if (FILE* f = std::fopen("Research\\out\\beso_topology.txt", "w")) {
        for (const Member& mem : work.members) {
            if (!mem.active) continue;
            const Vec3 a = work.nodes[static_cast<size_t>(work.nodeIndex(mem.i))].pos;
            const Vec3 b = work.nodes[static_cast<size_t>(work.nodeIndex(mem.j))].pos;
            std::fprintf(f, "%.1f %.1f %.1f %.1f\n", static_cast<double>(a.x), static_cast<double>(a.z),
                         static_cast<double>(b.x), static_cast<double>(b.z));
        }
        std::fclose(f);
        std::printf("[beso] topology -> Research\\out\\beso_topology.txt\n");
    }
    std::printf("[done]\n");
    return 0;
}
