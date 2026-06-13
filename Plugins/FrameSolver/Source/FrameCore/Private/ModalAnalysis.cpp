#include "FrameCore/ModalAnalysis.h"
#include "PreparedSystemImpl.h"
#include "SparseEigsolver.h"

#include <algorithm>
#include <cmath>

namespace frame {

ModalResult solveModal(const PreparedSystem& prepared, const ModalOptions& opts) {
    ModalResult R;
    const PreparedSystem::Impl& S = *prepared.impl;
    if (S.singular) { R.singular = true; R.diagnostic = S.diagnostic; return R; }
    const int N = S.N, nf = S.nf;
    if (nf == 0) { R.singular = true; R.diagnostic = "no free DOF for modal analysis"; return R; }

    // assemble global consistent mass M
    std::vector<Triplet> mtrips;
    for (const auto& el : S.elems) el->assembleMass(mtrips);
    SpMat M(N, N);
    M.setFromTriplets(mtrips.begin(), mtrips.end());
    M.makeCompressed();

    // ---- scale-up path: sparse subspace iteration reusing the LDLT (opt-in). Same modes as
    // the dense path; falls through to dense on non-convergence. ----
    const int kReq = std::min(opts.numModes, nf);
    if (opts.useSparseSolver && kReq > 0) {
        std::vector<Triplet> kt, mt;
        for (int c = 0; c < N; ++c)
            for (SpMat::InnerIterator it(S.K, c); it; ++it) {
                const int r = it.row();
                if (S.fmap[r] >= 0 && S.fmap[c] >= 0) kt.emplace_back(S.fmap[r], S.fmap[c], it.value());
            }
        for (int c = 0; c < N; ++c)
            for (SpMat::InnerIterator it(M, c); it; ++it) {
                const int r = it.row();
                if (S.fmap[r] >= 0 && S.fmap[c] >= 0) mt.emplace_back(S.fmap[r], S.fmap[c], it.value());
            }
        SpMat Kff_s(nf, nf); Kff_s.setFromTriplets(kt.begin(), kt.end()); Kff_s.makeCompressed();
        SpMat Mff_s(nf, nf); Mff_s.setFromTriplets(mt.begin(), mt.end()); Mff_s.makeCompressed();
        VecX lambda; MatX vecs;
        if (subspaceSmallest(Kff_s, S.ldlt, Mff_s, kReq, lambda, vecs)) {
            for (int i = 0; i < kReq; ++i) {
                VecX phi = vecs.col(i);
                const real mm = phi.dot(Mff_s * phi);
                if (mm > 0) phi /= std::sqrt(mm);        // M-normalize (match the dense path)
                real lam = lambda(i); if (lam < 0) lam = 0;
                ModeShape ms;
                ms.omega  = std::sqrt(lam);
                ms.freqHz = ms.omega / twoPi;
                ms.shape.assign(static_cast<size_t>(N), 0.0);
                for (int g = 0; g < N; ++g) if (S.fmap[g] >= 0) ms.shape[static_cast<size_t>(g)] = phi(S.fmap[g]);
                R.modes.push_back(ms);
            }
            return R;
        }
        // non-convergence -> fall through to the dense path below
    }

    // reduce K_ff and M_ff to dense on the free DOFs (modest interactive models)
    MatX Kff = MatX::Zero(nf, nf), Mff = MatX::Zero(nf, nf);
    for (int c = 0; c < N; ++c)
        for (SpMat::InnerIterator it(S.K, c); it; ++it) {
            const int r = it.row();
            if (S.fmap[r] >= 0 && S.fmap[c] >= 0) Kff(S.fmap[r], S.fmap[c]) += it.value();
        }
    for (int c = 0; c < N; ++c)
        for (SpMat::InnerIterator it(M, c); it; ++it) {
            const int r = it.row();
            if (S.fmap[r] >= 0 && S.fmap[c] >= 0) Mff(S.fmap[r], S.fmap[c]) += it.value();
        }

    real mtot = 0;
    for (int i = 0; i < nf; ++i) mtot += Mff(i, i);
    if (mtot <= 0) { R.singular = true; R.diagnostic = "zero mass (set Material.rho > 0 for modal)"; return R; }

    // generalized symmetric eigenproblem  K phi = lambda M phi  (lambda = omega^2, ascending)
    Eigen::GeneralizedSelfAdjointEigenSolver<MatX> ges(Kff, Mff);
    if (ges.info() != Eigen::Success) { R.singular = true; R.diagnostic = "generalized eigensolve failed"; return R; }
    const VecX evals = ges.eigenvalues();
    const MatX evecs = ges.eigenvectors();

    const int k = std::min(opts.numModes, static_cast<int>(evals.size()));
    for (int i = 0; i < k; ++i) {
        real lam = evals(i);
        if (lam < 0) lam = 0;                       // guard tiny negative (roundoff / rigid mode)
        ModeShape ms;
        ms.omega  = std::sqrt(lam);
        ms.freqHz = ms.omega / twoPi;
        ms.shape.assign(static_cast<size_t>(N), 0.0);
        for (int g = 0; g < N; ++g) if (S.fmap[g] >= 0) ms.shape[static_cast<size_t>(g)] = evecs(S.fmap[g], i);
        R.modes.push_back(ms);
    }
    return R;
}

}  // namespace frame
