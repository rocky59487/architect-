#pragma once
#include "FrameCore/FrameModel.h"
#include <string>

namespace frame { namespace grillage {

// Grillage (woven beam-grid) preprocessing. An isotropic rectangular plate is idealized as a grid of longitudinal +
// transverse beams (Hambly, "Bridge Deck Behaviour"). Each beam represents a slab strip
// of tributary width w and thickness t. The strip rigidities are nu-inflated so the
// equivalent orthotropic plate matches the isotropic plate's Dx = Dy = 2H = D exactly
// (see Grillage.cpp stripSection for the derivation):
//   bending inertia  I = w * t^3 / [12 (1 - nu^2)]   (so E*I per width = D)
//   torsion constant J = w * t^3 / [6  (1 - nu )]    (so G*J per width = D, with the
//                                                     elastic G = E/[2(1+nu)])
// The plain Hambly recipe (I=w t^3/12, J=w t^3/6) instead leaves Dx=Dy=(1-nu^2)D and
// 2H=(1-nu)D — too flexible — which for nu=0.3 measures as a ~26% center-deflection
// over-estimate; the nu-inflations above remove it.
//
// Only the OUT-OF-PLANE action (vertical deflection Uz + bending rotations Rx,Ry +
// member torsion) is physical, so the in-plane DOFs (Ux,Uy,Rz) are restrained at every
// node — the classic grillage idealization. Simple supports fix Uz on the boundary.
//
// Loading: the uniform pressure q is lumped to consistent nodal loads (q * tributary
// area), so the total applied load equals q*a*b exactly (load-conservation oracle).
//
// Accuracy: with the matched rigidities the center deflection lands within ~2% of
// Kirchhoff/Timoshenko plate theory and is mesh-STABLE — it settles to a value a couple
// of percent off (the residual is load-lumping + the moment-distribution mismatch), it
// does NOT converge to the exact plate value. Transverse bending moments are known to be
// OVER-estimated by the grillage analogy (the Poisson cross-moment is traded for extra
// twisting moment), so this is an engineering idealization, not an exact plate solver.
struct PlateSpec {
    real a  = 4000.0;   // plate span along x (mm)
    real b  = 4000.0;   // plate span along y (mm)
    real t  = 250.0;    // thickness (mm)
    int  nx = 8;        // mesh cells along x (use even -> a center node exists)
    int  ny = 8;        // mesh cells along y
    real q  = 0.025;    // uniform downward pressure (N/mm^2); applied as -q in global z
    bool simplySupported = true;   // fix Uz on all 4 boundary edges (else only Uz corners)
};

// Fills the caller-owned `m` with the grillage idealization of `spec`, using material
// `mat` (its G must satisfy G = E/[2(1+nu)] for the plate's Poisson ratio so the torsion
// recipe matches isotropic plate theory). Returns false + sets `why` on bad input.
// Node index convention: node(i,j) = j*(nx+1) + i, i in [0,nx], j in [0,ny].
FRAMECORE_API bool buildGrillage(FrameModel& m, const PlateSpec& spec,
                                 const Material& mat, std::string& why);

// Convenience: index of node (i,j) under the convention above.
inline int gridNode(int i, int j, int nx) { return j * (nx + 1) + i; }

}} // namespace frame::grillage
