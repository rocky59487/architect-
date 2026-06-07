#pragma once
#include "FrameCore/FrameModel.h"

namespace frame {

// Append self-weight to the model as equivalent loads, derived from Material.rho and
// geometry. Gravity acts along global -Z with magnitude g (mm/s^2; default 9810 = 9.81 m/s^2).
//
// UNIT BRIDGE (critical): the engine is a consistent N-mm-tonne-s system, but rho is given
// in kg/m^3. mass[tonne] = mass[kg]/1000 and 1 m^3 = 1e9 mm^3, so
//     rho[tonne/mm^3] = rho[kg/m^3] * 1e-12.
// A beam's weight per length  w = rho*A*g*1e-12  [N/mm];
// a shell's weight per area    p = rho*t*g*1e-12  [MPa]  (verified: tonne*mm/s^2 = N).
//
// Beams receive a MemberUDL with the global gravity rotated into the member's local axes
// (so a sloped/vertical member is correct). Shells receive the body load lumped to their
// four corner nodes as global -Z NodalLoads. Call this BEFORE solve(); it only appends.
FRAMECORE_API void addSelfWeight(FrameModel& m, real g = 9810.0);

}  // namespace frame
