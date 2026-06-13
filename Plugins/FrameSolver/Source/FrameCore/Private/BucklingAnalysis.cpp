#include "FrameCore/BucklingAnalysis.h"
#include "PreparedSystemImpl.h"
#include "SparseEigsolver.h"

#include <algorithm>
#include <vector>

namespace frame {

namespace {

// reduceFF (free-free submatrix extraction) is shared from PreparedSystemImpl.h.

// Historical DENSE path: (-Kg_ff) phi = gamma K_ff phi, criticalFactor = 1/gamma_max. Kept
// bit-identical to the pre-S1 implementation (so F23 is unchanged); the only refactor is that
// Kg is built once by the caller and passed in.
BucklingResult bucklingDense(const PreparedSystem::Impl& S, const SpMat& Kg) {
    BucklingResult R;
    const int N = S.N, nf = S.nf;
    MatX Kff = MatX::Zero(nf, nf), negKg = MatX::Zero(nf, nf);
    for (int c = 0; c < N; ++c)
        for (SpMat::InnerIterator it(S.K, c); it; ++it) {
            const int r = it.row();
            if (S.fmap[r] >= 0 && S.fmap[c] >= 0) Kff(S.fmap[r], S.fmap[c]) += it.value();
        }
    for (int c = 0; c < N; ++c)
        for (SpMat::InnerIterator it(Kg, c); it; ++it) {
            const int r = it.row();
            if (S.fmap[r] >= 0 && S.fmap[c] >= 0) negKg(S.fmap[r], S.fmap[c]) += -it.value();
        }
    Eigen::GeneralizedSelfAdjointEigenSolver<MatX> ges(negKg, Kff);
    if (ges.info() != Eigen::Success) { R.singular = true; R.diagnostic = "buckling eigensolve failed"; return R; }
    const VecX gam = ges.eigenvalues();
    const MatX vec = ges.eigenvectors();

    int idx = -1; real gmax = 0;
    for (int i = 0; i < gam.size(); ++i) if (gam(i) > gmax) { gmax = gam(i); idx = i; }
    if (idx < 0 || gmax <= 0) { R.singular = true; R.diagnostic = "no positive buckling eigenvalue (no compression?)"; return R; }

    R.criticalFactor = 1.0 / gmax;
    R.mode.assign(static_cast<size_t>(N), 0.0);
    for (int g = 0; g < N; ++g) if (S.fmap[g] >= 0) R.mode[static_cast<size_t>(g)] = vec(S.fmap[g], idx);
    return R;
}

// SPARSE path: subspaceSmallest(A = K_ff, Ainv = existing LDLT, B = -Kg_ff) returns the smallest
// lambda of  K_ff x = lambda (-Kg_ff) x,  which IS the critical factor (= 1/gamma_max of the
// dense pencil). It reuses K_ff's factorization (no re-factor). Returns true + fills R on
// success; false -> the caller falls back to the always-correct dense path. Honest scope: the
// result is iterative / tolerance-level (not bit-identical to dense). A loose pencil-residual
// gate (research measured 1e-7..1e-10 on the test models, so 1e-3 never trips on a well-posed
// model) rejects a converged-but-inaccurate eigenpair to the dense fallback.
bool bucklingSparse(const PreparedSystem::Impl& S, const SpMat& Kg, const BucklingOptions& opts,
                    BucklingResult& R) {
    const int N = S.N, nf = S.nf;
    const SpMat Kff     = reduceFF(S.K, S.fmap, nf, real(1));
    const SpMat negKgff = reduceFF(Kg,  S.fmap, nf, real(-1));

    VecX lambda; MatX vec;
    const int nev = std::max(1, opts.nev);
    if (!subspaceSmallest(Kff, S.ldlt, negKgff, nev, lambda, vec, opts.maxIter, opts.tol))
        return false;
    if (!(lambda(0) > 0)) return false;   // no positive critical factor -> let dense decide

    // Pencil residual ||(-Kg)phi - gamma K phi|| / (gamma ||K phi||), gamma = 1/lambda.
    const VecX v     = vec.col(0);
    const real gamma = real(1) / lambda(0);
    const VecX Kv    = Kff * v;
    const VecX Bv    = negKgff * v;
    const real residual = (Bv - gamma * Kv).norm() / std::max<real>(real(1e-300), gamma * Kv.norm());
    if (!(residual < real(1e-3))) return false;   // distrust -> dense fallback (always correct)

    R.criticalFactor = lambda(0);
    R.mode.assign(static_cast<size_t>(N), 0.0);
    for (int g = 0; g < N; ++g) if (S.fmap[g] >= 0) R.mode[static_cast<size_t>(g)] = v(S.fmap[g]);
    return true;
}

}  // namespace

BucklingResult solveBuckling(const PreparedSystem& prepared, const FrameModel& model,
                             const BucklingOptions& opts) {
    BucklingResult R;
    const PreparedSystem::Impl& S = *prepared.impl;
    if (S.singular) { R.singular = true; R.diagnostic = S.diagnostic; return R; }
    const int N = S.N, nf = S.nf;
    if (nf == 0) { R.singular = true; R.diagnostic = "no free DOF for buckling"; return R; }

    // 1) reference linear solve -> member axial forces (compression-positive)
    const SolveResult lin = solveLoad(prepared, model);
    if (lin.singular) { R.singular = true; R.diagnostic = "reference linear solve singular"; return R; }
    std::vector<real> axial(model.members.size(), 0.0);
    for (size_t e = 0; e < lin.memberForces.size() && e < axial.size(); ++e)
        axial[e] = lin.memberForces[e].endI.N;

    // 2) geometric stiffness Kg from those axial forces
    std::vector<Triplet> gtrips;
    for (const auto& el : S.elems) el->assembleGeometric(gtrips, axial);
    SpMat Kg(N, N);
    Kg.setFromTriplets(gtrips.begin(), gtrips.end());
    Kg.makeCompressed();

    // 3) reduced eigenproblem. Sparse subspace iteration for large models (reuses K_ff's LDLT),
    //    dense GeneralizedSelfAdjointEigenSolver otherwise. All-tension (gtrips.empty()) flows to
    //    the dense path, which reports the historical "no compression" diagnostic unchanged. The
    //    sparse path falls back to dense on any non-convergence, so correctness is never at risk.
    const bool wantSparse = !gtrips.empty() && (opts.denseThreshold <= 0 || nf > opts.denseThreshold);
    if (wantSparse && bucklingSparse(S, Kg, opts, R)) return R;
    return bucklingDense(S, Kg);
}

BucklingResult solveBuckling(const PreparedSystem& prepared, const FrameModel& model) {
    return solveBuckling(prepared, model, BucklingOptions{});
}

}  // namespace frame
