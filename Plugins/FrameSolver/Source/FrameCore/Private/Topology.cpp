#include "FrameCore/Topology.h"
#include "FrameCore/FrameSolver.h"
#include "FrameCore/Reanalysis.h"
#include "FrameCore/Connectivity.h"
#include "FrameCore/ElasticAllowable.h"
#include "FrameCore/MemberGeometry.h"
#include "MemberQuery.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <set>
#include <vector>

namespace frame {

// (public; declared in Topology.h) Element strain energy alpha_e = 1/2 u_e^T K_e u_e, reconstructed
// from the recovered LOCAL end forces dotted with the LOCAL end displacements. Q = k_local u_local
// (UDL-free); recover() stores endI = Q[0..5], endJ = (-Q[6], Q[7..11]) -- the compression-positive
// flip is ONLY on endJ.N -- so the end-J axial term takes -endJ.N. UDL-free members: exact, so the
// sum over active members == 1/2 * F.u (energy balance). Weights split the four energy components.
real memberStrainEnergy(const FrameModel& m, const SolveResult& r, int e,
                        real wAxial, real wBending, real wShear, real wTorsion) {
    if (e < 0 || e >= (int)m.members.size() || (size_t)e >= r.memberForces.size()) return 0;
    const Member& mem = m.members[(size_t)e];
    const int ni = m.nodeIndex(mem.i), nj = m.nodeIndex(mem.j);
    if (ni < 0 || nj < 0) return 0;
    Vec3 ax, ay, az;
    memberLocalAxes(m.nodes[(size_t)ni].pos, m.nodes[(size_t)nj].pos, mem.refVec, ax, ay, az);
    auto loc = [&](const Vec3& v) { return Vec3(dot(ax, v), dot(ay, v), dot(az, v)); };
    const Vec3 dI = loc(Vec3(r.u[(size_t)gdof(ni, Ux)], r.u[(size_t)gdof(ni, Uy)], r.u[(size_t)gdof(ni, Uz)]));
    const Vec3 rI = loc(Vec3(r.u[(size_t)gdof(ni, Rx)], r.u[(size_t)gdof(ni, Ry)], r.u[(size_t)gdof(ni, Rz)]));
    const Vec3 dJ = loc(Vec3(r.u[(size_t)gdof(nj, Ux)], r.u[(size_t)gdof(nj, Uy)], r.u[(size_t)gdof(nj, Uz)]));
    const Vec3 rJ = loc(Vec3(r.u[(size_t)gdof(nj, Rx)], r.u[(size_t)gdof(nj, Ry)], r.u[(size_t)gdof(nj, Rz)]));
    const MemberEndForces& fi = r.memberForces[(size_t)e].endI;
    const MemberEndForces& fj = r.memberForces[(size_t)e].endJ;
    const real axial = wAxial   * (fi.N * dI.x + (-fj.N) * dJ.x);
    const real shear = wShear   * (fi.Vy * dI.y + fi.Vz * dI.z + fj.Vy * dJ.y + fj.Vz * dJ.z);
    const real tors  = wTorsion * (fi.T * rI.x + fj.T * rJ.x);
    const real bend  = wBending * (fi.My * rI.y + fi.Mz * rI.z + fj.My * rJ.y + fj.Mz * rJ.z);
    return real(0.5) * (axial + shear + tors + bend);
}

namespace {

// screenable / memberLen / memberDC are shared from MemberQuery.h.

// External work C = sum over nodal loads of F . u (the compliance). For a linear-elastic model with
// only nodal loads (no prescribed settlement), C == 2 * total strain energy == 2 * sum(alpha_e).
real complianceOf(const FrameModel& m, const SolveResult& r) {
    real C = 0;
    for (const NodalLoad& l : m.nodalLoads) {
        const int ni = m.nodeIndex(l.node);
        if (ni < 0) continue;
        for (int d = 0; d < 6; ++d) C += l.comp[d] * r.u[(size_t)gdof(ni, d)];
    }
    return C;
}

// Active-state cycle detection uses fnv1aHashBytes (MemberQuery.h): the N2 rollback+lock can
// revisit an active set; a repeat means a cycle, which the guard turns into finite termination.

}  // namespace

BESOResult runBESO(const FrameModel& model, const BESOOptions& opts,
                   const std::vector<int>& designMembers) {
    BESOResult R;
    std::string why;
    if (!model.validate(why)) { R.reason = BESOStop::Invalid; return R; }

    FrameModel work = model;                       // caller's model is never mutated
    const size_t nMem = work.members.size();

    // Resolve the design (removable) set -> member INDICES. Empty request = every screenable active
    // member. Members outside the set are never removed (support links, N2-locked bars).
    std::vector<char> isDesign(nMem, 0);
    if (designMembers.empty()) {
        for (size_t e = 0; e < nMem; ++e) if (screenable(work, e)) isDesign[e] = 1;
    } else {
        for (int e : designMembers)
            if (e >= 0 && e < (int)nMem && screenable(work, (size_t)e)) isDesign[(size_t)e] = 1;
    }
    std::vector<char> isProtected(nMem, 0);        // N2-locked (never removable)

    // Nodes that carry load must stay grounded (cheap connectivity guard, mirrors the prototype).
    std::set<NodeId> loadedNodes;
    for (const NodalLoad& l : work.nodalLoads) loadedNodes.insert(l.node);

    // Initial design volume (sum A_e * L_e over active design members).
    real totalVol0 = 0;
    for (size_t e = 0; e < nMem; ++e)
        if (isDesign[e] && work.members[e].active)
            totalVol0 += work.sections[(size_t)work.members[e].secIdx].A * memberLen(work, (int)e);
    if (totalVol0 <= 0) { R.reason = BESOStop::Invalid; return R; }

    ReanalysisOptions ro; ro.solve = opts.solve;
    ReSolveSession sess(work, ro);
    if (!sess.valid()) { R.reason = BESOStop::Mechanism; return R; }

    const ElasticAllowable screen;

    // bad(): a tentative removal is rejected if it ungrounds a loaded node (connectivity) OR makes a
    // mechanism. The N1 ReSolveSession capacitance test is the FAST primary guard (catches bending
    // mechanisms connectivity cannot see, no fresh factor); a fresh solve CONFIRMS the survivors,
    // because ReSolve's capacitance tolerance (mechPivotTol) and the fresh LDLT pivotTol differ, so
    // a planar edge mechanism near that pivot threshold can slip past the capacitance test alone.
    // sess reflects work's active flags (kept in lock-step below).
    auto bad = [&]() -> bool {
        const ConnectivityResult cr = analyzeConnectivity(work);
        for (NodeId ln : cr.looseNodes) if (loadedNodes.count(ln)) return true;
        for (const FragmentCluster& fc : cr.detached)
            for (NodeId n : fc.nodes) if (loadedNodes.count(n)) return true;
        if (sess.solve().singular) return true;            // N1 capacitance (fast primary)
        return solve(work, opts.solve).singular;           // fresh confirm (edge mechanisms)
    };

    auto curVol = [&]() -> real {
        real v = 0;
        for (size_t e = 0; e < nMem; ++e)
            if (isDesign[e] && work.members[e].active)
                v += work.sections[(size_t)work.members[e].secIdx].A * memberLen(work, (int)e);
        return v;
    };

    std::vector<real> prevAlpha(nMem, 0);
    bool havePrev = false;
    real Cprev = 0;
    std::set<uint64_t> seenStates;

    for (int it = 1; it <= opts.maxIter; ++it) {
        R.iterations = it;
        const SolveResult r = sess.solve();
        if (r.singular) { R.reason = BESOStop::Mechanism; break; }

        const real C = complianceOf(work, r);
        const real volFrac = curVol() / totalVol0;
        R.volFracHistory.push_back(volFrac);
        R.complianceHistory.push_back(C);

        // (2) single-step compliance blow-up -> stop, keep bestActive at the previous feasible step.
        if (it > 1 && std::isfinite(C) && std::isfinite(Cprev) && C > opts.complianceJumpTol * Cprev) {
            R.reason = BESOStop::ComplianceJump;
            break;
        }
        // this step is feasible -> it becomes the compliance-best topology
        R.bestActive.assign(nMem, 0);
        for (size_t e = 0; e < nMem; ++e) R.bestActive[e] = work.members[e].active ? 1 : 0;
        R.bestIter = it;

        if (volFrac <= opts.targetVolFrac) { R.reason = BESOStop::TargetReached; R.converged = true; break; }

        // (1) sensitivity = element strain energy, two-step history average.
        std::vector<real> alpha(nMem, 0);
        std::vector<std::pair<real, int>> order;   // (sensitivity, member index), active design only
        for (size_t e = 0; e < nMem; ++e) {
            if (!isDesign[e] || !work.members[e].active || isProtected[e]) continue;
            real a = memberStrainEnergy(work, r, (int)e, opts.wAxial, opts.wBending, opts.wShear, opts.wTorsion);
            if (opts.sensHistory && havePrev) a = real(0.5) * (a + prevAlpha[e]);
            alpha[e] = a;
            order.emplace_back(a, (int)e);
        }
        prevAlpha = alpha;
        havePrev = true;

        // (2-remove) remove the lowest-sensitivity bars up to the volume quota; guard each removal.
        std::sort(order.begin(), order.end());
        const real vol = curVol();
        const real quota = std::min(opts.evolRate * vol, vol - opts.targetVolFrac * totalVol0);
        real removedVol = 0;
        std::vector<int> batch;
        for (const auto& [a, e] : order) {
            (void)a;
            if (removedVol >= quota) break;
            const MemberId id = work.members[(size_t)e].id;
            work.members[(size_t)e].active = false;
            sess.setMemberActive(id, false);
            if (bad()) {                                   // (3) mechanism / ungrounded -> roll back
                work.members[(size_t)e].active = true;
                sess.setMemberActive(id, true);
                continue;
            }
            removedVol += work.sections[(size_t)work.members[(size_t)e].secIdx].A * memberLen(work, e);
            batch.push_back(e);
        }
        if (batch.empty()) { R.reason = BESOStop::Stalled; break; }

        // N2: every redundancyCheckEvery steps, probe single-member-removal scenarios. If ANY makes
        // the (already-trimmed) structure Collapse, this step removed too much -> roll the batch back
        // and LOCK those bars (protected). LSP-grade robustness, not fail-safe.
        if (opts.redundancyCheckEvery > 0 && (it % opts.redundancyCheckEvery) == 0) {
            // scenario set: active design bars (optionally the highest-D/C `redundancySamples`).
            std::vector<int> scen;
            for (size_t e = 0; e < nMem; ++e)
                if (isDesign[e] && work.members[e].active) scen.push_back((int)e);
            if (opts.redundancySamples > 0 && (int)scen.size() > opts.redundancySamples) {
                std::vector<std::pair<real, int>> byDC;
                for (int e : scen) byDC.emplace_back(memberDC(screen, work, r, (size_t)e), e);
                std::sort(byDC.begin(), byDC.end(), [](auto& A, auto& B) { return A.first > B.first; });
                scen.clear();
                for (int k = 0; k < opts.redundancySamples; ++k) scen.push_back(byDC[(size_t)k].second);
            }

            bool fragile = false;
            for (int es : scen) {
                CollapseOptions co = opts.redundancy;
                co.solve = opts.solve;
                co.initialRemovals = { work.members[(size_t)es].id };
                co.initialShellRemovals.clear();
                if (runProgressiveCollapse(work, co).outcome == CollapseOutcome::Collapsed) {
                    fragile = true;
                    break;
                }
            }
            if (fragile) {
                for (int e : batch) {
                    const MemberId id = work.members[(size_t)e].id;
                    work.members[(size_t)e].active = true;
                    sess.setMemberActive(id, true);
                    isProtected[(size_t)e] = 1;
                }
                // FNV-1a cycle guard on the restored active set (rollback+lock could loop).
                if (!seenStates.insert(fnv1aHashBytes([&] {
                        std::vector<char> a(nMem);
                        for (size_t e = 0; e < nMem; ++e) a[e] = work.members[e].active ? 1 : 0;
                        return a; }())).second) {
                    R.reason = BESOStop::Stalled;
                    break;
                }
            }
        }
        Cprev = C;
    }

    if (R.reason == BESOStop::Invalid && R.iterations >= opts.maxIter) R.reason = BESOStop::MaxIter;

    R.finalActive.assign(nMem, 0);
    for (size_t e = 0; e < nMem; ++e) R.finalActive[e] = work.members[e].active ? 1 : 0;
    if (!opts.complianceBestRollback || R.bestActive.empty()) R.bestActive = R.finalActive;
    for (size_t e = 0; e < nMem; ++e) if (isProtected[e]) R.protectedMembers.push_back(work.members[e].id);
    std::sort(R.protectedMembers.begin(), R.protectedMembers.end());
    return R;
}

}  // namespace frame
