// exp_sparse_buckling.cpp — WS-B sparse linear-buckling prototype (scratch, NOT engine).
//
// Question: can the existing zero-dependency subspace iteration (Private/SparseEigsolver.h,
// currently used opt-in by modal) replace the DENSE GeneralizedSelfAdjointEigenSolver in
// BucklingAnalysis.cpp?
//
// Mapping:  engine dense path solves   (-Kg_ff) phi = gamma K_ff phi,  criticalFactor = 1/gamma_max.
//           sparse path: subspaceSmallest(A = K_ff, Ainv = existing ldlt, B = -Kg_ff, nev)
//           returns the smallest lambda of  K_ff x = lambda (-Kg_ff) x,  i.e. lambda_min = 1/gamma_max
//           = criticalFactor directly. B is only PSD (zero rows on unloaded DOFs) — that is the
//           numerical risk under test.
//
// Cases: meshed cantilever column (Euler analytic), towers S/M (dense vs sparse), tower L
// (sparse-only + residual check; dense memory estimated, not run), tension-only degenerate
// (Kg empty -> guard behaviour on both paths).

#include "research_common.h"
#include "SparseEigsolver.h"
#include "FrameCore/BucklingAnalysis.h"
#include <cstring>

using namespace research;

namespace {

struct CaseResult {
    bool   haveDense = false;
    double lamDense = 0, tDense = 0;
    bool   sparseOk = false;
    double lamSparse = 0, tSparse = 0, residual = 0;
};

// Mirrors BucklingAnalysis.cpp steps 1-2, then runs both eigen paths.
void runCase(const char* name, const FrameModel& model, bool doDense, double lamAnalytic = -1) {
    PreparedSystem ps = assembleAndFactor(model);
    const PreparedSystem::Impl& S = *ps.impl;
    if (S.singular) { std::printf("[%s] FATAL singular baseline: %s\n", name, S.diagnostic.c_str()); return; }

    SolveResult lin = solveLoad(ps, model);
    if (lin.singular) { std::printf("[%s] FATAL linear solve singular\n", name); return; }
    std::vector<real> axial(model.members.size(), 0.0);
    for (size_t e = 0; e < lin.memberForces.size() && e < axial.size(); ++e)
        axial[e] = lin.memberForces[e].endI.N;

    std::vector<Triplet> gtrips;
    for (const auto& el : S.elems) el->assembleGeometric(gtrips, axial);

    if (gtrips.empty()) {
        // degenerate: no compression anywhere -> Kg = 0 -> lambda_cr = +inf. Engine dense path?
        BucklingResult R = solveBuckling(ps, model);
        std::printf("[%s] degenerate KgEmpty=1 guard=lambda_cr_inf engineSingular=%d engineDiag=\"%s\"\n",
                    name, R.singular ? 1 : 0, R.diagnostic.c_str());
        return;
    }

    SpMat Kg(S.N, S.N);
    Kg.setFromTriplets(gtrips.begin(), gtrips.end());
    Kg.makeCompressed();
    const SpMat negKgff = research::reduceFF(Kg, S.fmap, S.nf, -1.0);
    const SpMat Kff     = research::reduceFF(S.K, S.fmap, S.nf);

    CaseResult cr;

    if (doDense) {
        Timer td;
        BucklingResult R = solveBuckling(ps, model);   // engine dense path (includes its own re-solve)
        cr.tDense = td.ms();
        if (!R.singular) { cr.haveDense = true; cr.lamDense = R.criticalFactor; }
    }

    {
        Timer tsp;
        VecX lambda; MatX vec;
        cr.sparseOk = subspaceSmallest(Kff, S.ldlt, negKgff, 3, lambda, vec);
        cr.tSparse = tsp.ms();
        if (cr.sparseOk) {
            cr.lamSparse = lambda(0);
            // residual of the converged pencil:  || negKg v - gamma K v || / (gamma ||K v||)
            const VecX v = vec.col(0);
            const real gamma = real(1) / lambda(0);
            const VecX Kv = Kff * v;
            const VecX Bv = negKgff * v;
            cr.residual = (Bv - gamma * Kv).norm() / std::max<real>(1e-300, gamma * Kv.norm());
        }
    }

    const double denseMemMiB = static_cast<double>(S.nf) * S.nf * 8.0 * 2.0 / (1024.0 * 1024.0);
    std::printf("[%s] nf=%d denseRun=%d lamDense=%.9g tDenseMs=%.1f | sparseOk=%d lamSparse=%.9g tSparseMs=%.1f residual=%.3e | relErrSparseVsDense=%.3e relErrVsAnalytic=%.3e denseMemEstMiB=%.0f\n",
                name, S.nf, cr.haveDense ? 1 : 0, cr.lamDense, cr.tDense,
                cr.sparseOk ? 1 : 0, cr.lamSparse, cr.tSparse, cr.residual,
                (cr.haveDense && cr.sparseOk) ? std::abs(cr.lamSparse - cr.lamDense) / std::abs(cr.lamDense) : -1.0,
                (lamAnalytic > 0 && cr.sparseOk) ? std::abs(cr.lamSparse - lamAnalytic) / lamAnalytic : -1.0,
                denseMemMiB);
}

} // namespace

int main(int argc, char** argv) {
    bool runLarge = true;
    for (int i = 1; i < argc; ++i)
        if (!std::strcmp(argv[i], "--noLarge")) runLarge = false;

    // 1) Euler column: fixed-free, P_cr = pi^2 EI / (4 L^2). Reference axial load P_ref.
    {
        const real L = 6000.0, side = 200.0, Pref = 1.0e5;
        const real E = 200000.0;
        const real I = side * side * side * side / 12.0;
        const real Pcr = 3.14159265358979323846 * 3.14159265358979323846 * E * I / (4.0 * L * L);
        FrameModel col = makeColumn(8, L, side, Pref, 0.0);
        runCase("column8", col, /*doDense=*/true, Pcr / Pref);
    }

    // 2) towers: dense vs sparse at growing nf
    runCase("tower-S(3x2x4)", makeTower(3, 2, 4), true);
    runCase("tower-M(5x4x8)", makeTower(5, 4, 8), true);

    // 3) large tower: sparse only (dense memory printed as estimate)
    if (runLarge)
        runCase("tower-L(12x9x24)", makeTower(12, 9, 24), /*doDense=*/false);

    // 4) degenerate: tip TENSION -> no compressive member -> Kg empty
    {
        FrameModel ten = makeColumn(8, 6000.0, 200.0, /*axialP=*/-1.0e5, 0.0);   // negative P = pull UP
        runCase("tension-degenerate", ten, true);
    }

    std::printf("[done]\n");
    return 0;
}
