#pragma once
#include "FrameCore/FrameTypes.h"
#include "FrameCore/Material.h"

namespace frame {

// A 4-node MITC4 flat-shell facet (membrane + Reissner-Mindlin bending + drilling).
// The four corner nodes are ordered counter-clockwise when viewed from the +normal
// side; the element builds its own local frame from the corner geometry, so the
// caller only supplies connectivity + thickness + material (E, nu, G).
//
// Shells and members are PARALLEL element sources on the same node pool: each node
// still carries 6 DOF [Ux,Uy,Uz,Rx,Ry,Rz], so a shell facet contributes a 24-DOF
// (4 nodes x 6) local stiffness. A model may contain members only, shells only, or
// both.
struct ShellQuad {
    int    id  = 0;
    NodeId n[4]{ 0, 0, 0, 0 };   // corner node ids, CCW about the +normal
    int    matIdx = -1;          // index into FrameModel::materials (needs E, nu, G = E/[2(1+nu)])
    real   t   = 0;              // thickness (mm)

    ShellQuad() = default;
    ShellQuad(int id_, NodeId a, NodeId b, NodeId c, NodeId d, int matIdx_, real t_)
        : id(id_), n{ a, b, c, d }, matIdx(matIdx_), t(t_) {}
};

// Uniform transverse pressure on a shell facet (force / area), applied along the
// facet's +local-z (its normal). Lumped to consistent nodal loads by the element.
// In-plane / concentrated loads reuse the existing NodalLoad path.
struct ShellPressure {
    int  shell = 0;   // ShellQuad::id
    real p     = 0;   // pressure (N/mm^2), positive along the facet normal
};

} // namespace frame
