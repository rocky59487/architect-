#pragma once
#include "FrameCore/FrameTypes.h"
#include "FrameCore/Material.h"
#include "FrameCore/Section.h"
#include <array>

namespace frame {

// Common member-end release patterns. Local DOF order per end is [Ux Uy Uz Rx Ry Rz].
enum class ReleasePreset { Rigid, TrussPin, HingeI, HingeJ };

// Build a release[12] mask for a preset:
//   Rigid    = no releases (moment-resisting frame member).
//   TrussPin = all 6 rotations (both ends) released -> a 2-force / axial member.
//   HingeI/J = the two bending rotations (Ry,Rz) at end i / j -> a moment hinge.
FRAMECORE_API std::array<bool, 12> makeRelease(ReleasePreset p);

// A prismatic 3D beam-column connecting node i -> node j.
struct Member {
    MemberId id  = 0;
    NodeId   i   = 0;
    NodeId   j   = 0;
    // Material / Section are referenced by INDEX into FrameModel::materials / ::sections,
    // NOT by raw pointer: pushing more nodes/members/materials can never dangle them, so
    // there is no "reserve before capturing pointers" rule. validate() range-checks them.
    int      matIdx = -1;
    int      secIdx = -1;

    // Reference vector defining the local x-y plane. Default global +Z; the
    // solver falls back to +Y automatically when refVec is parallel to the axis.
    Vec3 refVec { 0, 0, 1 };

    // Per-DOF end release (truss / hinge mode). MVP keeps all false; this is the
    // hook for the truss idealization (release rotational DOFs) added later.
    std::array<bool, 12> release { {} };

    // When false, the member is excluded from assembly: it contributes nothing to the
    // global stiffness K, adds no equivalent nodal loads, and its end forces stay zero.
    // This is the element-removal hook for progressive-collapse driving (disable the
    // governing member, re-factor, re-solve). Toggling it is a STRUCTURAL change, so it
    // is part of the solveLoad reuse fingerprint (a flipped flag rejects a stale factor).
    bool active = true;

    Member() = default;
    Member(MemberId id_, NodeId i_, NodeId j_, int matIdx_, int secIdx_)
        : id(id_), i(i_), j(j_), matIdx(matIdx_), secIdx(secIdx_) {}
};

} // namespace frame
