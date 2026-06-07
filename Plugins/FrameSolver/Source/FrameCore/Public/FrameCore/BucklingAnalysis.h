#pragma once
#include "FrameCore/FrameSolver.h"
#include "FrameCore/BucklingResult.h"

namespace frame {

// Linear (eigenvalue) BUCKLING. Applies the model's reference load (via solveLoad), builds the
// geometric stiffness Kg from the resulting member axial forces, and finds the smallest factor
// lambda such that (K + lambda*Kg) is singular. P_cr = lambda * (applied reference load). For
// a slender column this reproduces the Euler load pi^2 EI/(KL)^2. Beam-column geometric
// stiffness only (shell geometric stiffness is a future addition). Linear buckling = the onset
// eigenvalue, NOT a nonlinear post-buckling path.
FRAMECORE_API BucklingResult solveBuckling(const PreparedSystem& prepared, const FrameModel& model);

}  // namespace frame
