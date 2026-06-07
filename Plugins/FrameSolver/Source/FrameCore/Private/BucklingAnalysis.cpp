#include "FrameCore/BucklingAnalysis.h"
#include "PreparedSystemImpl.h"

#include <algorithm>
#include <cmath>

namespace frame {

BucklingResult solveBuckling(const PreparedSystem& prepared, const FrameModel& model) {
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

    // 3) reduced eigenproblem  (-Kg_ff) phi = gamma K_ff phi  (B = K_ff SPD; gamma = 1/lambda).
    //    The buckling mode is the LARGEST positive gamma -> smallest critical lambda.
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

}  // namespace frame
