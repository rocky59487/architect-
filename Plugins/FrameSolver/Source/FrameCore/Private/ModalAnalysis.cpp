#include "FrameCore/ModalAnalysis.h"
#include "PreparedSystemImpl.h"

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
    const real twoPi = 2.0 * 3.14159265358979323846;
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
