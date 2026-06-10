// exp_tension_only.cpp — WS-D tension-only iteration experiment (scratch, NOT engine).
//
// X-braced portal (3D, out-of-plane DOFs pinned at the free nodes). Designated
// tension-only members: the two diagonals. Iteration (Jacobi-style sweep):
//   active   brace with axial N > tolN (compression-positive)  -> deactivate
//   inactive brace with axial elongation (u_j - u_i) . xhat > 0 -> reactivate
// until the active set repeats (state-hash) or is stable.
//
// Case A (lateral): classic one-tension/one-compression. Oracle: converged displacements
//   must equal a model with the compressed brace OMITTED from the member list entirely
//   (engine stage-1 property: active=false === omission).
// Case B (vertical+lateral sweep): hunts a flip-flop window (both braces slack under
//   vertical squash; sway re-tensions one) to exercise the cycle guard, then re-runs the
//   cycling load with the monotone no-reactivation policy.

#include "research_common.h"
#include "FrameCore/MemberGeometry.h"
#include <set>
#include <cstring>

using namespace research;

namespace {

constexpr int BR_A = 3;   // member index: diagonal 0->3
constexpr int BR_B = 4;   // member index: diagonal 1->2

FrameModel makePortal(real H, real V, bool withBraceA = true, bool withBraceB = true) {
    FrameModel p;
    p.materials.emplace_back(200000.0, 76923.07692307692, 7850.0);
    p.sections.push_back(squareSection(250.0));   // columns/beam
    p.sections.push_back(squareSection(60.0));    // slender braces
    Node n0(0, 0, 0, 0);       n0.fixAll();
    Node n1(1, 6000, 0, 0);    n1.fixAll();
    Node n2(2, 0, 0, 3000);    n2.fixed[Uy] = true;
    Node n3(3, 6000, 0, 3000); n3.fixed[Uy] = true;
    p.nodes = { n0, n1, n2, n3 };
    auto add = [&](int id, int i, int j, int sec) {
        Member m(id, i, j, 0, sec);
        m.refVec = refVecFor(p.nodes[static_cast<size_t>(i)].pos, p.nodes[static_cast<size_t>(j)].pos);
        p.members.push_back(m);
    };
    add(0, 0, 2, 0);            // column L
    add(1, 1, 3, 0);            // column R
    add(2, 2, 3, 0);            // beam
    add(3, 0, 3, 1);            // brace A (0->3), index BR_A; placeholder if disabled
    add(4, 1, 2, 1);            // brace B (1->2), index BR_B
    if (!withBraceA) p.members[BR_A].active = false;
    if (!withBraceB) p.members[BR_B].active = false;
    NodalLoad l2; l2.node = 2; l2.comp[Ux] = H; l2.comp[Uz] = -V;
    NodalLoad l3; l3.node = 3; l3.comp[Uz] = -V;
    p.nodalLoads.push_back(l2);
    p.nodalLoads.push_back(l3);
    return p;
}

real braceElongation(const FrameModel& m, const SolveResult& r, int e) {
    const Member& mem = m.members[static_cast<size_t>(e)];
    const int ni = m.nodeIndex(mem.i), nj = m.nodeIndex(mem.j);
    Vec3 ax, ay, az;
    memberLocalAxes(m.nodes[static_cast<size_t>(ni)].pos, m.nodes[static_cast<size_t>(nj)].pos, mem.refVec, ax, ay, az);
    const Vec3 dui(r.u[static_cast<size_t>(gdof(ni, Ux))], r.u[static_cast<size_t>(gdof(ni, Uy))], r.u[static_cast<size_t>(gdof(ni, Uz))]);
    const Vec3 duj(r.u[static_cast<size_t>(gdof(nj, Ux))], r.u[static_cast<size_t>(gdof(nj, Uy))], r.u[static_cast<size_t>(gdof(nj, Uz))]);
    return dot(duj - dui, ax);
}

struct ToResult {
    bool converged = false, cycled = false, singular = false;
    int  iters = 0;
    bool activeA = true, activeB = true;
    SolveResult last;
};

ToResult iterate(FrameModel work, bool allowReactivation, int maxIter = 30) {
    ToResult R;
    const real tolN = 1e-9;     // N (force units); braces carry O(1e4) N
    std::set<unsigned> seen;
    for (int it = 1; it <= maxIter; ++it) {
        R.iters = it;
        const unsigned stateBits = (work.members[BR_A].active ? 1u : 0u) | (work.members[BR_B].active ? 2u : 0u);
        PreparedSystem ps = assembleAndFactor(work);
        if (ps.impl->singular) { R.singular = true; return R; }
        SolveResult r = solveLoad(ps, work);
        if (r.singular) { R.singular = true; return R; }
        R.last = r;

        bool changed = false;
        for (int e : { BR_A, BR_B }) {
            Member& mem = work.members[static_cast<size_t>(e)];
            if (mem.active) {
                const real N = r.memberForces[static_cast<size_t>(e)].endI.N;   // compression+
                if (N > tolN) { mem.active = false; changed = true; }
            } else if (allowReactivation) {
                if (braceElongation(work, r, e) > 0) { mem.active = true; changed = true; }
            }
        }
        if (!changed) {
            R.converged = true;
            R.activeA = work.members[BR_A].active;
            R.activeB = work.members[BR_B].active;
            return R;
        }
        const unsigned nextBits = (work.members[BR_A].active ? 1u : 0u) | (work.members[BR_B].active ? 2u : 0u);
        if (!seen.insert(stateBits * 4u + nextBits).second) {   // transition repeated => cycle
            R.cycled = true;
            R.activeA = work.members[BR_A].active;
            R.activeB = work.members[BR_B].active;
            return R;
        }
    }
    return R;
}

} // namespace

int main() {
    // ---------- case A: lateral load, classic eliminate-one ----------
    {
        FrameModel m = makePortal(5.0e4, 0.0);
        ToResult R = iterate(m, /*allowReactivation=*/true);

        // oracle: model with the compressed brace REMOVED from the member list entirely
        FrameModel ref = makePortal(5.0e4, 0.0);
        ref.members.erase(ref.members.begin() + (R.activeA ? BR_B : BR_A));
        SolveResult rr = solve(ref);
        const double relErr = relMaxDiffV(R.last.u, rr.u);

        const real Nten = R.last.memberForces[static_cast<size_t>(R.activeA ? BR_A : BR_B)].endI.N;
        std::printf("[caseA] converged=%d iters=%d activeA=%d activeB=%d tensionN=%.6g relErrVsOmitted=%.3e\n",
                    R.converged ? 1 : 0, R.iters, R.activeA ? 1 : 0, R.activeB ? 1 : 0,
                    static_cast<double>(Nten), relErr);
    }

    // ---------- case B: vertical squash + lateral sweep, hunt the flip-flop window ----------
    bool foundCycle = false;
    for (real Hl : { 1.0e2, 3.0e2, 1.0e3, 2.0e3, 3.0e3, 5.0e3, 1.0e4, 3.0e4 }) {
        FrameModel m = makePortal(Hl, 2.0e5);
        ToResult R = iterate(m, true);
        std::printf("[caseB] H=%.0f V=2e5 converged=%d cycled=%d iters=%d activeA=%d activeB=%d\n",
                    static_cast<double>(Hl), R.converged ? 1 : 0, R.cycled ? 1 : 0, R.iters,
                    R.activeA ? 1 : 0, R.activeB ? 1 : 0);
        if (R.cycled && !foundCycle) {
            foundCycle = true;
            ToResult Rm = iterate(makePortal(Hl, 2.0e5), /*allowReactivation=*/false);
            std::printf("[caseB-monotone] H=%.0f converged=%d iters=%d activeA=%d activeB=%d (no-reactivation policy)\n",
                        static_cast<double>(Hl), Rm.converged ? 1 : 0, Rm.iters, Rm.activeA ? 1 : 0, Rm.activeB ? 1 : 0);
        }
    }
    if (!foundCycle)
        std::printf("[caseB-note] no flip-flop window found in this sweep (honest negative result)\n");

    std::printf("[done]\n");
    return 0;
}
