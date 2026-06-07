#pragma once
//
// Private definition of PreparedSystem::Impl (Eigen-carrying). Shared by FrameSolver.cpp
// (assembleAndFactor / solveLoad) and the analysis modules (ModalAnalysis, Buckling, ...)
// that need the assembled K, the free-DOF map and the prepared elements. Never included by
// a public header — the opaque PreparedSystem keeps the public boundary Eigen-free.
//
#include "FrameCore/FrameSolver.h"
#include "FrameEigen.h"
#include "IElement.h"

#include <vector>
#include <string>
#include <memory>

namespace frame {

struct PreparedSystem::Impl {
    int                                    N = 0;
    SpMat                                  K;          // full global K (reactions, K_ff, M reduction)
    std::vector<int>                       fmap;       // free-DOF map (-1 = constrained)
    int                                    nf = 0;
    LDLTSolver                             ldlt;       // factorization of K_ff
    std::vector<std::unique_ptr<IElement>> elems;      // prepared elements (stiffness/mass/recovery)
    bool                                   singular = false;
    std::string                            diagnostic;
};

}  // namespace frame
