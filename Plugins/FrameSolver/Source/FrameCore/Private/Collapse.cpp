// C2 progressive-collapse driver (collapse stage 3c): sequential linear analysis over the
// existing factorize-once machinery. Each step is an ordinary assembleAndFactor + solveLoad
// on a private working copy, so the solveLoad reuse fingerprint is honoured by construction
// (every structural change gets a fresh factorization).
#include "FrameCore/Collapse.h"
#include "FrameCore/FrameSolver.h"
#include "FrameCore/ElasticAllowable.h"
#include "FrameCore/MemberGeometry.h"
#include "FrameCore/NMInteraction.h"
#include "CollapseSupport.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace frame {
namespace {

// FrameModel::memberIndex / shellIndex (FrameModel.h) and pinNode / anyActive (CollapseSupport.h)
// are shared with the dynamic driver.

// The hinge's NODE-side channel (see Hinge.h): the joint moment -Mp * e_axis the condensed
// element can no longer deliver. `sign = -1` retracts a previously added moment (used when a
// driver-formed hinge's member is later removed BRITTLY while its node stays grounded --
// detached nodes need no retraction, pinNode() erases their loads wholesale).
void addHingeNodeMoment(FrameModel& work, const CollapseHingeEvent& h, real sign) {
    const int e = work.memberIndex(h.member);
    if (e < 0) return;
    const Member& mem = work.members[(size_t)e];
    const Vec3 pi = work.nodes[(size_t)work.nodeIndex(mem.i)].pos;
    const Vec3 pj = work.nodes[(size_t)work.nodeIndex(mem.j)].pos;
    Vec3 ax, ay, az;
    memberLocalAxes(pi, pj, mem.refVec, ax, ay, az);
    const Vec3 axis = (h.dof == 4 || h.dof == 10) ? ay : az;
    NodalLoad nl;
    nl.node = (h.dof < 6) ? mem.i : mem.j;
    nl.comp[Rx] = -sign * h.Mp * axis.x;
    nl.comp[Ry] = -sign * h.Mp * axis.y;
    nl.comp[Rz] = -sign * h.Mp * axis.z;
    work.nodalLoads.push_back(nl);
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
        if (work.memberIndex(id) < 0) {
            H.diagnostic = "initialRemovals references missing member id " + std::to_string(id);
            return H;
        }
    for (int sid : opts.initialShellRemovals)
        if (work.shellIndex(sid) < 0) {
            H.diagnostic = "initialShellRemovals references missing shell id " + std::to_string(sid);
            return H;
        }

    const ElasticAllowable screen;
    std::vector<MemberId> pending = opts.initialRemovals;
    std::vector<int> pendingShells = opts.initialShellRemovals;
    std::vector<CollapseHingeEvent> pendingHinges;     // 4b: at most one per step
    std::vector<CollapseHingeEvent> driverFormed;      // node-side moments WE added (for retraction)
    FailMode pendingMode = FailMode::None;   // step 0 is scenario-imposed, not D/C-selected
    real pendingRatio = 0;

    for (int step = 0;; ++step) {
        H.steps.emplace_back();
        CollapseStep& S = H.steps.back();
        S.step = step;
        S.mode = pendingMode;
        S.triggerRatio = pendingRatio;

        // 1) apply the pending event(s)
        for (MemberId id : pending) {
            work.members[(size_t)work.memberIndex(id)].active = false;
            S.removedMembers.push_back(id);
            // a brittle removal of a member the driver hinged earlier: its node-side moments
            // must leave with it (the node may stay grounded, so pinNode never sees them)
            for (const CollapseHingeEvent& h : driverFormed)
                if (h.member == id) addHingeNodeMoment(work, h, real(-1));
        }
        pending.clear();
        for (int sid : pendingShells) {
            work.shells[(size_t)work.shellIndex(sid)].active = false;
            S.removedShells.push_back(sid);
        }
        pendingShells.clear();
        for (const CollapseHingeEvent& h : pendingHinges) {
            work.hinges.push_back(PlasticHinge{ h.member, h.dof, h.Mp });
            addHingeNodeMoment(work, h, real(1));
            driverFormed.push_back(h);
            S.formedHinges.push_back(h);
        }
        pendingHinges.clear();

        // 2) fragment cleanup BEFORE solving: anything no longer connected to a support is
        //    debris -- deactivate it wholesale, pin its nodes, drop its loads. Solving first
        //    would misread the floating piece as a mechanism of the whole model.
        const ConnectivityResult conn = analyzeConnectivity(work);
        if (!conn.valid) {   // cannot happen after the validate() above; defensive
            H.diagnostic = "connectivity analysis failed on the working model";
            return H;
        }
        for (const FragmentCluster& fc : conn.detached) {
            for (MemberId id : fc.members) work.members[(size_t)work.memberIndex(id)].active = false;
            for (int sid : fc.shells)      work.shells[(size_t)work.shellIndex(sid)].active = false;
            for (NodeId nid : fc.nodes)    pinNode(work, nid);
            S.detached.push_back(fc);
        }
        for (NodeId nid : conn.looseNodes) pinNode(work, nid);

        // 3) anything left to carry load?
        if (!anyActive(work)) {
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

        // 5) D/C screen + event selection. The REPORTED maxDC is always the full allowable
        //    screen (worstUtilization / worstShellUtilization math). The EVENT is:
        //      - brittle mode (default): remove the worst element while D/C > removeThreshold;
        //      - hinge mode (4b): a hinge-capable member's bending is ductile -- it forms a
        //        hinge when |M| >= Mp = fy*Z instead of being removed; its brittle screen
        //        shrinks to the separable shear/torsion ratios; shells stay brittle.
        //    Deterministic tie-break: worst ratio, then removal-before-hinge (member before
        //    shell), then smallest id, then smallest dof.
        struct Event {
            real ratio = 0;
            int kind = -1;          // 0 = member removal, 1 = shell removal, 2 = hinge
            int id = 0;             // member or shell id
            int dof = 0;            // hinge dof
            real mpSigned = 0;      // hinge residual moment
            FailMode mode = FailMode::None;
            bool eligible = false;
        } best;
        auto consider = [&best](const Event& c) {
            if (!c.eligible) return;
            if (!best.eligible ||
                c.ratio > best.ratio ||
                (c.ratio == best.ratio && (c.kind < best.kind ||
                 (c.kind == best.kind && (c.id < best.id ||
                  (c.id == best.id && c.dof < best.dof))))))
                best = c;
        };
        auto ratioOf = [](real demand, real cap) -> real {
            if (cap > 0) return demand / cap;
            return demand > 0 ? std::numeric_limits<real>::infinity() : real(0);
        };

        bool any = false;
        real maxDC = 0;
        const size_t nM = std::min(work.members.size(), r.memberForces.size());
        for (size_t e = 0; e < nM; ++e) {
            const Member& mem = work.members[e];
            if (!mem.active) continue;
            if (mem.matIdx < 0 || mem.matIdx >= (int)work.materials.size()) continue;
            if (mem.secIdx < 0 || mem.secIdx >= (int)work.sections.size())  continue;
            const Section&  sec = work.sections[(size_t)mem.secIdx];
            const Material& mat = work.materials[(size_t)mem.matIdx];
            const MemberForcePair& mf = r.memberForces[e];
            const DemandResult di = screen.checkSection(mf.endI, sec, mat.cap);
            const DemandResult dj = screen.checkSection(mf.endJ, sec, mat.cap);
            const DemandResult& d = (di.risk >= dj.risk) ? di : dj;
            maxDC = std::max(maxDC, d.risk);
            any = true;

            const bool hingeCapable = opts.plasticHinges && mat.fy > 0 && sec.Zy > 0 && sec.Zz > 0;
            if (!hingeCapable) {
                consider(Event{ d.risk, 0, mem.id, 0, 0, d.mode, d.risk > opts.removeThreshold });
            } else {
                // ductile bending: hinge candidates per end/axis at |M|/Mp >= 1. With N-M
                // interaction (S10) the threshold AND the residual use the axially-reduced
                // Mp_eff(N) = Mp*(1-(N/Ny)^2) evaluated at THIS step's end axial force; once a
                // hinge actually forms its residual is stored in PlasticHinge.Mp and frozen
                // (the released end then recovers M = 0, so it is never re-evaluated). Without
                // the flag N is fed as 0 -> Mp_eff == Mp, bit-for-bit the stage-4b behaviour.
                const real MpY = mat.fy * sec.Zy, MpZ = mat.fy * sec.Zz;
                const real Ny  = mat.fy * sec.A;
                const real Ni  = opts.nmInteraction ? mf.endI.N : real(0);
                const real Nj  = opts.nmInteraction ? mf.endJ.N : real(0);
                const real Ms[4]   = { mf.endI.My, mf.endI.Mz, mf.endJ.My, mf.endJ.Mz };
                const int  dofs[4]  = { 4, 5, 10, 11 };
                const real Mp0[4]   = { MpY, MpZ, MpY, MpZ };
                const real Naxe[4]  = { Ni, Ni, Nj, Nj };   // end-matched axial force per dof
                for (int k = 0; k < 4; ++k) {
                    const real Mpk = reducedPlasticMoment(Mp0[k], Naxe[k], Ny);
                    if (!(Mpk > 0)) continue;   // axial-squashed: no ductile bending capacity
                                                // left -- the brittle axial screen below governs
                    const real ratio = std::fabs(Ms[k]) / Mpk;
                    consider(Event{ ratio, 2, mem.id, dofs[k],
                                    (Ms[k] >= 0 ? Mpk : -Mpk), FailMode::Bending,
                                    ratio >= 1.0 });
                }
                // brittle remainder: the SEPARABLE ratios -- pure axial (N/A, free of the
                // bending stress that is ductile here), shear, torsion. The combined-fibre
                // sM +/- sN allowable screen stays reporting-only for hinge-capable members.
                const real A = std::max(sec.A, real(1e-12));
                const real sNi = mf.endI.N / A, sNj = mf.endJ.N / A;   // compression-positive
                const real rC = ratioOf(std::max({ sNi, sNj, real(0) }), mat.cap.comp);
                const real rX = ratioOf(std::max({ -sNi, -sNj, real(0) }), mat.cap.tens);
                const real tau  = std::max(di.tau, dj.tau);
                const real sTor = std::max(di.sTor, dj.sTor);
                const real rS = ratioOf(tau, mat.cap.shear), rT = ratioOf(sTor, mat.cap.tors);
                real rB = rC; FailMode mB = FailMode::Crush;
                if (rX > rB) { rB = rX; mB = FailMode::Tension; }
                if (rS > rB) { rB = rS; mB = FailMode::Shear; }
                if (rT > rB) { rB = rT; mB = FailMode::Torsion; }
                consider(Event{ rB, 0, mem.id, 0, 0, mB, rB > opts.removeThreshold });
            }
        }
        const size_t nS = std::min(work.shells.size(), r.shellForces.size());
        for (size_t s = 0; s < nS; ++s) {
            const ShellQuad& sh = work.shells[s];
            if (!sh.active) continue;
            if (sh.matIdx < 0 || sh.matIdx >= (int)work.materials.size()) continue;
            const ShellDemandResult d = checkShellSurface(r.shellForces[s], sh.t,
                                                          work.materials[(size_t)sh.matIdx].cap);
            maxDC = std::max(maxDC, d.risk);
            any = true;
            consider(Event{ d.risk, 1, sh.id, 0, 0, FailMode::ShellVonMises,
                            d.risk > opts.removeThreshold });
        }

        S.maxDC = any ? maxDC : 0;
        S.safetyFactor = any ? ((maxDC > 0) ? real(1) / maxDC : std::numeric_limits<real>::infinity())
                             : real(0);

        // 6) terminal checks, then queue the winning event for the next step
        if (!best.eligible) {
            H.outcome = CollapseOutcome::Stable;
            return H;
        }
        if ((int)H.steps.size() >= opts.maxSteps) {
            H.outcome = CollapseOutcome::MaxSteps;
            H.diagnostic = "step budget exhausted with events still occurring";
            return H;
        }
        switch (best.kind) {
            case 0: pending = { best.id }; break;
            case 1: pendingShells = { best.id }; break;
            default: pendingHinges = { CollapseHingeEvent{ best.id, best.dof, best.mpSigned } }; break;
        }
        pendingMode = best.mode;
        pendingRatio = best.ratio;
    }
}

}  // namespace frame
