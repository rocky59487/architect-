#pragma once
//
// Shared member-screening / state-hashing helpers used by the optimization passes (SizeOpt.cpp,
// Topology.cpp) and the tension-only active set (TensionOnly.cpp). Header-only inline so the
// helpers share ONE definition — each TU used to keep a file-local copy (a unity-build collision,
// see FrameCore.Build.cs). POD/std + ElasticAllowable only; no Eigen.
//
#include "FrameCore/FrameModel.h"
#include "FrameCore/SolveResult.h"
#include "FrameCore/ElasticAllowable.h"
#include <algorithm>
#include <cstdint>
#include <vector>

namespace frame {

// A member is screenable when it is active and its material/section indices are in range,
// so ElasticAllowable::checkSection can be called on it safely.
inline bool screenable(const FrameModel& m, size_t e) {
    const Member& mem = m.members[e];
    return mem.active
        && mem.matIdx >= 0 && mem.matIdx < (int)m.materials.size()
        && mem.secIdx >= 0 && mem.secIdx < (int)m.sections.size();
}

// Length of member e (distance between its end nodes; 0 if an end node is missing).
inline real memberLen(const FrameModel& m, int e) {
    const Member& mem = m.members[(size_t)e];
    const int ni = m.nodeIndex(mem.i), nj = m.nodeIndex(mem.j);
    if (ni < 0 || nj < 0) return 0;
    return norm(m.nodes[(size_t)nj].pos - m.nodes[(size_t)ni].pos);
}

// Worst (max) demand/capacity ratio of member e across its two ends, via the elastic screen.
inline real memberDC(const ElasticAllowable& screen, const FrameModel& m, const SolveResult& r, size_t e) {
    const Member& mem = m.members[e];
    const Section&  s = m.sections[(size_t)mem.secIdx];
    const Capacity& c = m.materials[(size_t)mem.matIdx].cap;
    const DemandResult di = screen.checkSection(r.memberForces[e].endI, s, c);
    const DemandResult dj = screen.checkSection(r.memberForces[e].endJ, s, c);
    return std::max(di.risk, dj.risk);
}

// FNV-1a hash of a 0/1 active-state vector, for cycle detection in active-set iterations
// (each entry contributes its boolean truth, so equal active sets hash equal).
inline uint64_t fnv1aHashBytes(const std::vector<char>& v) {
    uint64_t h = 1469598103934665603ull;
    for (char a : v) { h ^= (uint64_t)(a ? 1 : 0); h *= 1099511628211ull; }
    return h;
}

}  // namespace frame
