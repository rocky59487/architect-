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
#include <cstdint>

namespace frame {

// Extract the free-free (nf x nf) submatrix of a full N x N sparse matrix via the free-DOF map
// (fmap[g] >= 0 selects a free DOF and gives its reduced index; < 0 = constrained). `scale`
// multiplies every entry (use -1 to build -Kg for the buckling pencil; default 1 is bit-identical
// to an unscaled copy). Single definition shared by FrameSolver / Buckling / PDelta / Reanalysis /
// DynamicCollapse — each used to keep its own file-local copy (a unity-build name collision).
inline SpMat reduceFF(const SpMat& Afull, const std::vector<int>& fmap, int nf, real scale = real(1)) {
    std::vector<Triplet> t;
    t.reserve(static_cast<size_t>(Afull.nonZeros()));
    for (int c = 0; c < Afull.outerSize(); ++c)
        for (SpMat::InnerIterator it(Afull, c); it; ++it) {
            const int r = it.row();
            if (fmap[(size_t)r] >= 0 && fmap[(size_t)c] >= 0)
                t.emplace_back(fmap[(size_t)r], fmap[(size_t)c], scale * it.value());
        }
    SpMat R(nf, nf);
    R.setFromTriplets(t.begin(), t.end());
    R.makeCompressed();
    return R;
}

struct PreparedSystem::Impl {
    int                                    N = 0;
    SpMat                                  K;          // full global K (reactions, K_ff, M reduction)
    std::vector<int>                       fmap;       // free-DOF map (-1 = constrained)
    int                                    nf = 0;
    real                                   pivotMargin = 0;   // min/max |LDLT pivot| of K_ff (C4 criticality margin)
    LDLTSolver                             ldlt;       // factorization of K_ff
    std::vector<std::unique_ptr<IElement>> elems;      // prepared elements (stiffness/mass/recovery)
    bool                                   singular = false;
    std::string                            diagnostic;
    uint64_t                               fingerprint = 0;   // structural/geometry/UDL hash baseline;
                                                              // solveLoad rejects a model whose
                                                              // fingerprint changed (see FrameSolver.cpp)
};

}  // namespace frame
