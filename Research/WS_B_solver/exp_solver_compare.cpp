// exp_solver_compare.cpp — WS-B Eigen built-in solver comparison on the SAME K_ff / F_ff
// (scratch, NOT engine). Direct: SimplicialLDLT (engine default), SimplicialLLT, SparseLU.
// Iterative: CG+Jacobi, CG+IncompleteCholesky. All zero-new-dependency options.

#include "research_common.h"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>
#include <cstring>

using namespace research;

namespace {

double relResidual(const SpMat& K, const VecX& x, const VecX& F) {
    return static_cast<double>((K * x - F).norm() / std::max<real>(1e-300, F.norm()));
}

void runPreset(const char* name, int nx, int ny, int st) {
    FrameModel model = makeTower(nx, ny, st);
    assertNodalOnly(model, "exp_solver_compare");
    PreparedSystem ps = assembleAndFactor(model);
    const PreparedSystem::Impl& S = *ps.impl;
    if (S.singular) { std::printf("[%s] FATAL singular\n", name); return; }
    const SpMat Kff = research::reduceFF(S.K, S.fmap, S.nf);
    const VecX  Fff = reduceVec(nodalLoadVector(model), S.fmap, S.nf);
    std::printf("[case] %s nf=%d nnz=%lld\n", name, S.nf, static_cast<long long>(Kff.nonZeros()));

    {   // SimplicialLDLT (engine default path)
        Timer tc; LDLTSolver s; s.compute(Kff); const double tC = tc.ms();
        Timer tsv; const VecX x = s.solve(Fff); const double tS = tsv.ms();
        std::printf("[solver] %s SimplicialLDLT factorMs=%.1f solveMs=%.2f res=%.2e\n", name, tC, tS, relResidual(Kff, x, Fff));
    }
    {   // SimplicialLLT
        Timer tc; Eigen::SimplicialLLT<SpMat> s; s.compute(Kff); const double tC = tc.ms();
        Timer tsv; const VecX x = s.solve(Fff); const double tS = tsv.ms();
        std::printf("[solver] %s SimplicialLLT  factorMs=%.1f solveMs=%.2f res=%.2e\n", name, tC, tS, relResidual(Kff, x, Fff));
    }
    {   // SparseLU
        Timer tc; Eigen::SparseLU<SpMat> s; s.analyzePattern(Kff); s.factorize(Kff); const double tC = tc.ms();
        Timer tsv; const VecX x = s.solve(Fff); const double tS = tsv.ms();
        std::printf("[solver] %s SparseLU       factorMs=%.1f solveMs=%.2f res=%.2e\n", name, tC, tS, relResidual(Kff, x, Fff));
    }
    {   // CG + Jacobi
        Timer tc;
        Eigen::ConjugateGradient<SpMat, Eigen::Lower | Eigen::Upper> cg;
        cg.setTolerance(1e-10); cg.setMaxIterations(20000); cg.compute(Kff);
        const double tC = tc.ms();
        Timer tsv; const VecX x = cg.solve(Fff); const double tS = tsv.ms();
        std::printf("[solver] %s CG+Jacobi      setupMs=%.1f solveMs=%.1f iters=%d cgErr=%.2e res=%.2e\n",
                    name, tC, tS, static_cast<int>(cg.iterations()), static_cast<double>(cg.error()), relResidual(Kff, x, Fff));
    }
    {   // CG + IncompleteCholesky
        Timer tc;
        Eigen::ConjugateGradient<SpMat, Eigen::Lower | Eigen::Upper,
                                 Eigen::IncompleteCholesky<real>> cg;
        cg.setTolerance(1e-10); cg.setMaxIterations(20000); cg.compute(Kff);
        const double tC = tc.ms();
        Timer tsv; const VecX x = cg.solve(Fff); const double tS = tsv.ms();
        std::printf("[solver] %s CG+IC          setupMs=%.1f solveMs=%.1f iters=%d cgErr=%.2e res=%.2e\n",
                    name, tC, tS, static_cast<int>(cg.iterations()), static_cast<double>(cg.error()), relResidual(Kff, x, Fff));
    }
}

} // namespace

int main(int argc, char** argv) {
    bool large = true;
    for (int i = 1; i < argc; ++i)
        if (!std::strcmp(argv[i], "--noLarge")) large = false;
    runPreset("tower-M(5x4x8)", 5, 4, 8);
    if (large) runPreset("tower-XXL(12x9x24)", 12, 9, 24);
    std::printf("[done]\n");
    return 0;
}
