#pragma once
#include "FrameCore/FrameSolver.h"
#include "FrameCore/ModalResult.h"

namespace frame {

// Free-vibration MODAL analysis: solves the generalized eigenproblem  K phi = omega^2 M phi
// on the free DOFs and returns the lowest `opts.numModes` modes (ascending frequency). Reuses
// the prepared system's assembled stiffness + free-DOF map and assembles the consistent mass
// from the elements. Dense generalized eigensolver — intended for the modest models of an
// interactive sim; a sparse Lanczos path is the future scale-up. The mass requires
// Material.rho > 0.
FRAMECORE_API ModalResult solveModal(const PreparedSystem& prepared, const ModalOptions& opts = {});

}  // namespace frame
