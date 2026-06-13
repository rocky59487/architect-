#pragma once
//
// The SINGLE choke point for all Eigen includes inside FrameCore.
// Every FrameCore translation unit that needs Eigen includes THIS header,
// never <Eigen/...> directly.
//
//   * Standalone build : FRAMECORE_UE undefined  -> plain include.
//   * UE module build  : FRAMECORE_UE=1 (from FrameCore.Build.cs) -> Eigen is
//     wrapped in UE's third-party guard macros (pattern copied from Chaos
//     BlockSparseLinearSystem.cpp). UE-only TUs must #include "CoreMinimal.h"
//     BEFORE this header so the guard macros exist. Pure-core TUs must NEVER
//     pull in CoreMinimal.h (it would break the standalone build).
//
#ifndef EIGEN_MPL2_ONLY
#define EIGEN_MPL2_ONLY              // also passed via /D in the standalone build
#endif

#if defined(FRAMECORE_UE)
  #include "CoreMinimal.h"            // guarantees THIRD_PARTY_INCLUDES_* / PRAGMA_* macros (UE branch only)
  PRAGMA_DEFAULT_VISIBILITY_START
  THIRD_PARTY_INCLUDES_START
  #include <Eigen/Core>
  #include <Eigen/Dense>
  #include <Eigen/Sparse>
  #include <Eigen/SparseCholesky>    // SimplicialLDLT
  #include <Eigen/IterativeLinearSolvers>   // ConjugateGradient (ReSolve Tier-2)
  THIRD_PARTY_INCLUDES_END
  PRAGMA_DEFAULT_VISIBILITY_END
#else
  #include <Eigen/Core>
  #include <Eigen/Dense>
  #include <Eigen/Sparse>
  #include <Eigen/SparseCholesky>
  #include <Eigen/IterativeLinearSolvers>
#endif

#include "FrameCore/FrameTypes.h"

#include <algorithm>
#include <cmath>

namespace frame {

using SpMat      = Eigen::SparseMatrix<real>;          // column-major (SimplicialLDLT-friendly)
using Triplet    = Eigen::Triplet<real>;
using VecX       = Eigen::Matrix<real, Eigen::Dynamic, 1>;
using MatX       = Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic>;
using Mat12      = Eigen::Matrix<real, 12, 12>;
using Vec12      = Eigen::Matrix<real, 12, 1>;
using Mat3       = Eigen::Matrix<real, 3, 3>;
using Vec3e      = Eigen::Matrix<real, 3, 1>;
using LDLTSolver = Eigen::SimplicialLDLT<SpMat>;

inline Vec3e toE(const Vec3& v) { return Vec3e(v.x, v.y, v.z); }

// kPi / twoPi live in FrameTypes.h (visible to every TU without pulling in Eigen).

// Positive-definiteness test of an LDLT-factored matrix from its diagonal D. SimplicialLDLT
// reports info()==Success even for an indefinite matrix (it does not pivot for stability), so a
// non-positive / near-zero entry in D is how a singular-or-indefinite factor is detected — mirrors
// the mechanism-detection test in assembleAndFactor. tol is relative to max|D| (scale-invariant).
inline bool ldltPositiveDefinite(const LDLTSolver& s, real relTol) {
    const VecX D = s.vectorD();
    real maxAbs = 0;
    for (int i = 0; i < D.size(); ++i) maxAbs = std::max(maxAbs, std::abs(D(i)));
    const real tol = relTol * std::max<real>(1, maxAbs);
    for (int i = 0; i < D.size(); ++i)
        if (!(D(i) > tol)) return false;
    return true;
}

} // namespace frame
