#pragma once
#include "FrameCore/FrameTypes.h"
#include "FrameCore/ISectionStrength.h"
#include <vector>
#include <string>

namespace frame {

struct MemberForcePair {
    MemberId        member = 0;
    MemberEndForces endI;   // local end forces at node i
    MemberEndForces endJ;   // local end forces at node j
};

// Stress resultants for a shell facet, in the facet's LOCAL frame. The scalar fields are at
// the element CENTRE (the average for a linearly-varying field). Bending/twist moments and
// transverse shears are per-unit-width (N*mm/mm = N); membrane forces are per-unit-width
// (N/mm). Per-CORNER bending moments are also provided for design PEAK recovery (see below).
struct ShellElementForces {
    int  shell = 0;                       // ShellQuad::id
    real Mxx = 0, Myy = 0, Mxy = 0;       // bending + twisting moments / width (centre)
    real Qx  = 0, Qy  = 0;                // transverse shears / width (centre)
    real Nxx = 0, Nyy = 0, Nxy = 0;       // membrane forces / width (centre)
    // Per-CORNER bending moments (corner order matches ShellQuad::n). These are LINEAR
    // quantities, so they combine()/envelope() correctly. The design PEAK moment is a
    // post-process: the max over corners of |Mxx|/|Myy|/|Mxy| — do NOT read the centre value
    // as a peak. For a constant-moment field every corner equals the centre.
    real MxxC[4] = { 0, 0, 0, 0 };
    real MyyC[4] = { 0, 0, 0, 0 };
    real MxyC[4] = { 0, 0, 0, 0 };
};

struct SolveResult {
    bool                            singular = false;   // mechanism / instability detected
    std::string                     diagnostic;
    std::vector<real>               u;                  // 6N global DOF displacements
    std::vector<real>               reactions;          // 6N global (nonzero at constrained DOF)
    std::vector<MemberForcePair>    memberForces;
    std::vector<ShellElementForces> shellForces;        // parallel to model.shells

    real disp(int nodeIndex, int dof) const { return u[static_cast<size_t>(gdof(nodeIndex, dof))]; }
    real reaction(int nodeIndex, int dof) const { return reactions[static_cast<size_t>(gdof(nodeIndex, dof))]; }
};

} // namespace frame
