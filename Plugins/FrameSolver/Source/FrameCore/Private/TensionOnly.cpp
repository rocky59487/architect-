#include "FrameCore/TensionOnly.h"
#include "FrameCore/Reanalysis.h"
#include "FrameCore/MemberGeometry.h"
#include "MemberQuery.h"

#include <set>
#include <vector>
#include <cstdint>

namespace frame {

namespace {

// Axial elongation of member e under the current displacement field: (u_j - u_i) . xhat, where
// xhat is the member's local axis. Positive = the ends pulled apart, so a currently-slack
// tension-only member should re-activate. Valid for an INACTIVE member too — its end nodes still
// carry displacements driven by the rest of the structure.
real memberElongation(const FrameModel& m, const SolveResult& r, int e) {
    const Member& mem = m.members[(size_t)e];
    const int ni = m.nodeIndex(mem.i), nj = m.nodeIndex(mem.j);
    if (ni < 0 || nj < 0) return 0;
    Vec3 ax, ay, az;
    memberLocalAxes(m.nodes[(size_t)ni].pos, m.nodes[(size_t)nj].pos, mem.refVec, ax, ay, az);
    const Vec3 dui(r.u[(size_t)gdof(ni, Ux)], r.u[(size_t)gdof(ni, Uy)], r.u[(size_t)gdof(ni, Uz)]);
    const Vec3 duj(r.u[(size_t)gdof(nj, Ux)], r.u[(size_t)gdof(nj, Uy)], r.u[(size_t)gdof(nj, Uz)]);
    return dot(duj - dui, ax);
}

// State hashing for cycle detection uses fnv1aHashBytes (MemberQuery.h): the iteration is
// deterministic, so a repeated state means a flip-flop cycle; a finite state space (2^n_TO)
// guarantees a repeat within finitely many steps, which the guard turns into termination.

// One active-set fixed-point run on a fresh ReSolveSession (each flip = exact rank-6 Woodbury
// update). allowReact=false -> monotone (deactivate-only) policy. Fills converged / cycled /
// iterations / finalState / slack.
TensionOnlyResult iterateTO(const FrameModel& model, const std::vector<int>& toIdx,
                            const TensionOnlyOptions& opts, bool allowReact) {
    TensionOnlyResult R;
    ReanalysisOptions ro; ro.solve = opts.solve;
    ReSolveSession sess(model, ro);
    if (!sess.valid()) {
        R.finalState.singular = true;
        R.finalState.diagnostic = sess.diagnostic();
        return R;
    }

    // current active state of the TO members (start from the model's own active flags, mirrored
    // into the session baseline)
    std::vector<char> act(toIdx.size());
    for (size_t k = 0; k < toIdx.size(); ++k) act[k] = model.members[(size_t)toIdx[k]].active ? 1 : 0;

    std::set<uint64_t> seen;
    seen.insert(fnv1aHashBytes(act));   // the starting state counts as visited

    SolveResult sr;
    for (int it = 1; it <= opts.maxIter; ++it) {
        R.iterations = it;
        sr = sess.solve();
        if (sr.singular) {   // the current active set formed a mechanism
            R.finalState = sr;
            return R;
        }

        int flips = 0;
        for (size_t k = 0; k < toIdx.size(); ++k) {
            const int e = toIdx[k];
            const MemberId id = model.members[(size_t)e].id;
            if (act[k]) {
                const real N = sr.memberForces[(size_t)e].endI.N;   // compression-positive
                if (N > opts.axialTol) { sess.setMemberActive(id, false); act[k] = 0; ++flips; }
            } else if (allowReact) {
                if (memberElongation(model, sr, e) > 0) { sess.setMemberActive(id, true); act[k] = 1; ++flips; }
            }
        }

        if (flips == 0) {   // a self-consistent active set: every TO member is in tension or slack
            R.converged = true;
            R.finalState = sr;
            for (size_t k = 0; k < toIdx.size(); ++k)
                if (!act[k]) R.slack.push_back(model.members[(size_t)toIdx[k]].id);
            return R;
        }

        // cycle guard: the NEW active state was seen before -> a flip-flop, bail to the fallback
        if (!seen.insert(fnv1aHashBytes(act)).second) {
            R.cycled = true;
            R.finalState = sr;
            return R;
        }
    }

    // maxIter exhausted without a verdict: report the last solved state (not flagged converged)
    R.finalState = sr;
    for (size_t k = 0; k < toIdx.size(); ++k)
        if (!act[k]) R.slack.push_back(model.members[(size_t)toIdx[k]].id);
    return R;
}

}  // namespace

TensionOnlyResult runTensionOnly(const FrameModel& model, const TensionOnlyOptions& opts) {
    std::vector<int> toIdx;
    for (size_t e = 0; e < model.members.size(); ++e)
        if (model.members[e].tensionOnly) toIdx.push_back((int)e);

    TensionOnlyResult R = iterateTO(model, toIdx, opts, opts.allowReactivation);

    // a flip-flop cycle -> auto-downgrade to the monotone (deactivate-only) policy, which is
    // guaranteed to terminate in <= n_TO steps; keep `cycled` set so the caller knows the result
    // is the order-dependent-conservative monotone equilibrium.
    if (R.cycled && opts.allowReactivation) {
        TensionOnlyResult Rm = iterateTO(model, toIdx, opts, /*allowReact=*/false);
        Rm.cycled = true;
        return Rm;
    }
    return R;
}

}  // namespace frame
