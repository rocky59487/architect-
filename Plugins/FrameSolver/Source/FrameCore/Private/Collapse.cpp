// C2 progressive-collapse driver (collapse stage 3c): sequential linear analysis over the
// existing factorize-once machinery. Each step is an ordinary assembleAndFactor + solveLoad
// on a private working copy, so the solveLoad reuse fingerprint is honoured by construction
// (every structural change gets a fresh factorization).
#include "FrameCore/Collapse.h"
#include "FrameCore/FrameSolver.h"
#include "FrameCore/ElasticAllowable.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace frame {
namespace {

int memberIndexById(const FrameModel& m, MemberId id) {
    for (size_t e = 0; e < m.members.size(); ++e)
        if (m.members[e].id == id) return (int)e;
    return -1;
}

int shellIndexById(const FrameModel& m, int id) {
    for (size_t s = 0; s < m.shells.size(); ++s)
        if (m.shells[s].id == id) return (int)s;
    return -1;
}

// Temporarily ground a node that no longer carries any active element (a debris node or a
// bare free node). Fixing every DOF at zero is mathematically inert -- nothing couples to it
// any more, so it adds no stiffness and draws no reaction -- but it removes the zero-pivot
// rows that would otherwise read as a spurious mechanism. Its nodal loads leave with it
// (the load belonged to what fell; leaking it into the grounded remainder would be wrong).
void pinNode(FrameModel& work, NodeId nid) {
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

}  // namespace

CollapseHistory runProgressiveCollapse(const FrameModel& model, const CollapseOptions& opts) {
    CollapseHistory H;   // outcome defaults to Invalid until a terminal state is reached

    if (opts.maxSteps < 1) { H.diagnostic = "CollapseOptions.maxSteps must be >= 1"; return H; }
    if (!std::isfinite(opts.dlf) || opts.dlf <= 0) { H.diagnostic = "CollapseOptions.dlf must be finite and > 0"; return H; }
    if (!std::isfinite(opts.removeThreshold) || opts.removeThreshold < 0) {
        H.diagnostic = "CollapseOptions.removeThreshold must be finite and >= 0"; return H;
    }
    std::string why;
    if (!model.validate(why)) { H.diagnostic = "invalid model: " + why; return H; }

    FrameModel work = model;   // the caller's model is never mutated

    // Scale FORCE loads once by the dynamic load factor. Prescribed displacement values are
    // kinematic boundary data, not forces -- they are deliberately left untouched.
    for (auto& nl : work.nodalLoads)
        for (int d = 0; d < DOF_PER_NODE; ++d) nl.comp[d] *= opts.dlf;
    for (auto& u : work.memberUDLs) {
        u.w_local.x *= opts.dlf; u.w_local.y *= opts.dlf; u.w_local.z *= opts.dlf;
    }
    for (auto& sp : work.shellPressures) sp.p *= opts.dlf;

    for (MemberId id : opts.initialRemovals)
        if (memberIndexById(work, id) < 0) {
            H.diagnostic = "initialRemovals references missing member id " + std::to_string(id);
            return H;
        }

    const ElasticAllowable screen;
    std::vector<MemberId> pending = opts.initialRemovals;
    FailMode pendingMode = FailMode::None;   // step 0 is scenario-imposed, not D/C-selected
    real pendingRatio = 0;

    for (int step = 0;; ++step) {
        H.steps.emplace_back();
        CollapseStep& S = H.steps.back();
        S.step = step;
        S.mode = pendingMode;
        S.triggerRatio = pendingRatio;

        // 1) apply the pending removal(s)
        for (MemberId id : pending) {
            work.members[(size_t)memberIndexById(work, id)].active = false;
            S.removedMembers.push_back(id);
        }
        pending.clear();

        // 2) fragment cleanup BEFORE solving: anything no longer connected to a support is
        //    debris -- deactivate it wholesale, pin its nodes, drop its loads. Solving first
        //    would misread the floating piece as a mechanism of the whole model.
        const ConnectivityResult conn = analyzeConnectivity(work);
        if (!conn.valid) {   // cannot happen after the validate() above; defensive
            H.diagnostic = "connectivity analysis failed on the working model";
            return H;
        }
        for (const FragmentCluster& fc : conn.detached) {
            for (MemberId id : fc.members) work.members[(size_t)memberIndexById(work, id)].active = false;
            for (int sid : fc.shells)      work.shells[(size_t)shellIndexById(work, sid)].active = false;
            for (NodeId nid : fc.nodes)    pinNode(work, nid);
            S.detached.push_back(fc);
        }
        for (NodeId nid : conn.looseNodes) pinNode(work, nid);

        // 3) anything left to carry load?
        bool anyActive = false;
        for (const Member& mem : work.members)  anyActive = anyActive || mem.active;
        for (const ShellQuad& sh : work.shells) anyActive = anyActive || sh.active;
        if (!anyActive) {
            S.u.assign(work.nodes.size() * DOF_PER_NODE, 0.0);   // everything detached reads 0
            H.outcome = CollapseOutcome::Collapsed;
            H.diagnostic = "no active element remains grounded";
            return H;
        }

        // 4) fresh factorization + solve of the grounded remainder
        const PreparedSystem ps = assembleAndFactor(work, opts.solve);
        const SolveResult r = solveLoad(ps, work);
        if (r.singular) {
            S.pivotMargin = r.pivotMargin;
            H.outcome = CollapseOutcome::Collapsed;
            H.diagnostic = "mechanism in the grounded remainder: " + r.diagnostic;
            return H;
        }
        S.solved = true;
        S.u = r.u;
        S.pivotMargin = r.pivotMargin;

        // 5) D/C screen over active members (same math as worstUtilization, plus the explicit
        //    deterministic tie-break: worst ratio first, then smallest member id).
        bool any = false;
        real maxDC = 0;
        MemberId gov = 0;
        FailMode mode = FailMode::None;
        const size_t nM = std::min(work.members.size(), r.memberForces.size());
        for (size_t e = 0; e < nM; ++e) {
            const Member& mem = work.members[e];
            if (!mem.active) continue;
            if (mem.matIdx < 0 || mem.matIdx >= (int)work.materials.size()) continue;
            if (mem.secIdx < 0 || mem.secIdx >= (int)work.sections.size())  continue;
            const Section&  sec = work.sections[(size_t)mem.secIdx];
            const Capacity& cap = work.materials[(size_t)mem.matIdx].cap;
            const DemandResult di = screen.checkSection(r.memberForces[e].endI, sec, cap);
            const DemandResult dj = screen.checkSection(r.memberForces[e].endJ, sec, cap);
            const DemandResult& d = (di.risk >= dj.risk) ? di : dj;
            if (!any || d.risk > maxDC || (d.risk == maxDC && mem.id < gov)) {
                maxDC = d.risk; gov = mem.id; mode = d.mode;
            }
            any = true;
        }
        S.maxDC = any ? maxDC : 0;
        S.safetyFactor = any ? ((maxDC > 0) ? real(1) / maxDC : std::numeric_limits<real>::infinity())
                             : real(0);

        // 6) terminal checks, then queue the governing member for the next step
        if (!any || maxDC <= opts.removeThreshold) {
            H.outcome = CollapseOutcome::Stable;
            return H;
        }
        if ((int)H.steps.size() >= opts.maxSteps) {
            H.outcome = CollapseOutcome::MaxSteps;
            H.diagnostic = "step budget exhausted with D/C still above removeThreshold";
            return H;
        }
        pending = { gov };
        pendingMode = mode;
        pendingRatio = maxDC;
    }
}

}  // namespace frame
