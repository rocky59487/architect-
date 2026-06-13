// exp_matrix_free_pcg.cpp -- WS-B research-only matrix-free LinearOperator/PCG validation.
//
// This file intentionally stays in Research/. It uses FrameCore Private headers through
// research_common.h, keeps the public API Eigen-free, and validates against the already
// assembled K_ff oracle from assembleAndFactor().

#include "research_common.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

using namespace research;

namespace {

struct ApplyStats {
    double relL2 = 0.0;
    double relMax = 0.0;
};

double safeNormDenom(const VecX& v) {
    return std::max<double>(1e-300, static_cast<double>(v.norm()));
}

double relL2Diff(const VecX& a, const VecX& b) {
    return static_cast<double>((a - b).norm()) / safeNormDenom(b);
}

double relMaxDiffEigen(const VecX& a, const VecX& b) {
    double maxDiff = 0.0;
    double maxRef = 0.0;
    const int n = std::min<int>(a.size(), b.size());
    for (int i = 0; i < n; ++i) {
        maxDiff = std::max(maxDiff, std::abs(static_cast<double>(a(i) - b(i))));
        maxRef = std::max(maxRef, std::abs(static_cast<double>(b(i))));
    }
    return (maxRef > 0.0) ? (maxDiff / maxRef) : maxDiff;
}

VecX deterministicVec(int n, int seed) {
    VecX x(n);
    const double a = 0.113 + 0.019 * seed;
    const double b = 0.071 + 0.011 * seed;
    for (int i = 0; i < n; ++i) {
        const double k = static_cast<double>(i + 1);
        x(i) = static_cast<real>(std::sin(a * k) + 0.5 * std::cos(b * (k + 3.0)));
    }
    return x;
}

struct SparseKffOperator {
    const SpMat& K;
    VecX diag;

    explicit SparseKffOperator(const SpMat& Kff) : K(Kff), diag(VecX::Zero(Kff.rows())) {
        for (int c = 0; c < K.outerSize(); ++c) {
            for (SpMat::InnerIterator it(K, c); it; ++it) {
                if (it.row() == it.col()) diag(static_cast<int>(it.row())) += it.value();
            }
        }
    }

    int rows() const { return static_cast<int>(K.rows()); }

    void apply(const VecX& x, VecX& y) const {
        y = K * x;
    }

    const VecX& diagonal() const { return diag; }
};

struct ElementWiseKffOperator {
    struct Entry {
        int r = 0;
        int c = 0;
        real v = 0.0;
    };
    struct Block {
        std::vector<Entry> entries;
    };

    int n = 0;
    std::vector<Block> blocks;
    VecX diag;
    size_t entryCount = 0;

    explicit ElementWiseKffOperator(const PreparedSystem::Impl& S)
        : n(S.nf), diag(VecX::Zero(S.nf)) {
        blocks.reserve(S.elems.size());
        for (const auto& el : S.elems) {
            std::vector<Triplet> trips;
            const int ldof = std::max(0, el->localDof());
            trips.reserve(static_cast<size_t>(ldof) * static_cast<size_t>(ldof));
            el->assemble(trips);

            Block b;
            b.entries.reserve(trips.size());
            for (const Triplet& t : trips) {
                const int gr = static_cast<int>(t.row());
                const int gc = static_cast<int>(t.col());
                if (gr < 0 || gc < 0 ||
                    gr >= static_cast<int>(S.fmap.size()) ||
                    gc >= static_cast<int>(S.fmap.size())) {
                    continue;
                }
                const int rr = S.fmap[static_cast<size_t>(gr)];
                const int cc = S.fmap[static_cast<size_t>(gc)];
                if (rr < 0 || cc < 0) continue;
                const real value = t.value();
                b.entries.push_back({ rr, cc, value });
                if (rr == cc) diag(rr) += value;
            }
            if (!b.entries.empty()) {
                entryCount += b.entries.size();
                blocks.push_back(std::move(b));
            }
        }
    }

    int rows() const { return n; }

    void apply(const VecX& x, VecX& y) const {
        y.setZero(n);
        for (const Block& b : blocks) {
            for (const Entry& e : b.entries) {
                y(e.r) += e.v * x(e.c);
            }
        }
    }

    const VecX& diagonal() const { return diag; }
};

template <typename Operator>
ApplyStats validateApply(const char* caseName, const char* opName, const Operator& A, const SpMat& Kff) {
    ApplyStats stats;
    VecX yOp(A.rows()), yRef(A.rows());
    for (int seed = 1; seed <= 3; ++seed) {
        const VecX x = deterministicVec(A.rows(), seed);
        A.apply(x, yOp);
        yRef = Kff * x;
        stats.relL2 = std::max(stats.relL2, relL2Diff(yOp, yRef));
        stats.relMax = std::max(stats.relMax, relMaxDiffEigen(yOp, yRef));
    }
    std::printf("[apply] %s %-12s relL2=%.3e relMax=%.3e\n",
                caseName, opName, stats.relL2, stats.relMax);
    return stats;
}

struct PcgResult {
    VecX x;
    int iterations = 0;
    double relResidual = std::numeric_limits<double>::infinity();
    bool converged = false;
    bool breakdown = false;
};

bool validJacobi(const VecX& diag, std::string& why) {
    for (int i = 0; i < diag.size(); ++i) {
        if (!(diag(i) > real(0)) || !std::isfinite(static_cast<double>(diag(i)))) {
            why = "non-positive/non-finite diagonal at reduced dof " + std::to_string(i);
            return false;
        }
    }
    return true;
}

template <typename Operator>
PcgResult solvePcg(const Operator& A, const VecX& b, double tol, int maxIter) {
    PcgResult out;
    const int n = A.rows();
    out.x = VecX::Zero(n);

    VecX Ax(n);
    A.apply(out.x, Ax);
    VecX r = b - Ax;
    const real bNorm = std::max<real>(real(1e-300), b.norm());
    out.relResidual = static_cast<double>(r.norm() / bNorm);
    if (out.relResidual <= tol) {
        out.converged = true;
        return out;
    }

    const VecX& diag = A.diagonal();
    VecX z = r.cwiseQuotient(diag);
    VecX p = z;
    real rz = r.dot(z);
    if (!(rz > real(0)) || !std::isfinite(static_cast<double>(rz))) {
        out.breakdown = true;
        return out;
    }

    VecX Ap(n);
    for (int iter = 0; iter < maxIter; ++iter) {
        A.apply(p, Ap);
        const real denom = p.dot(Ap);
        if (!(denom > real(0)) || !std::isfinite(static_cast<double>(denom))) {
            out.breakdown = true;
            return out;
        }

        const real alpha = rz / denom;
        out.x += alpha * p;
        r -= alpha * Ap;
        out.iterations = iter + 1;
        out.relResidual = static_cast<double>(r.norm() / bNorm);
        if (out.relResidual <= tol) {
            out.converged = true;
            return out;
        }

        z = r.cwiseQuotient(diag);
        const real rzNext = r.dot(z);
        if (!(rzNext >= real(0)) || !std::isfinite(static_cast<double>(rzNext))) {
            out.breakdown = true;
            return out;
        }
        const real beta = rzNext / rz;
        p = z + beta * p;
        rz = rzNext;
    }
    return out;
}

template <typename Operator>
bool runPcgCheck(const char* caseName,
                 const char* opName,
                 const Operator& A,
                 const SpMat& Kff,
                 const VecX& Fff,
                 const VecX& xRef,
                 double tol,
                 int maxIter,
                 double solutionTol) {
    std::string why;
    if (!validJacobi(A.diagonal(), why)) {
        std::printf("[pcg]   %s %-12s FAIL setup: %s\n", caseName, opName, why.c_str());
        return false;
    }

    Timer t;
    PcgResult r = solvePcg(A, Fff, tol, maxIter);
    const double ms = t.ms();
    const VecX trueResidual = Kff * r.x - Fff;
    const double relTrueRes = static_cast<double>(trueResidual.norm() / std::max<real>(real(1e-300), Fff.norm()));
    const double relX = relL2Diff(r.x, xRef);
    const bool ok = r.converged && !r.breakdown && relTrueRes <= 10.0 * tol && relX <= solutionTol;
    std::printf("[pcg]   %s %-12s %s iters=%d ms=%.2f pcgRel=%.3e trueRel=%.3e relX=%.3e\n",
                caseName, opName, ok ? "PASS" : "FAIL",
                r.iterations, ms, r.relResidual, relTrueRes, relX);
    return ok;
}

bool runCase(const char* name,
             int nx,
             int ny,
             int stories,
             double tol,
             int maxIter,
             bool withPcg) {
    Timer buildTimer;
    FrameModel model = makeTower(nx, ny, stories);
    const double buildMs = buildTimer.ms();
    assertNodalOnly(model, "exp_matrix_free_pcg");

    Timer factorTimer;
    PreparedSystem ps = assembleAndFactor(model);
    const double factorMs = factorTimer.ms();
    const PreparedSystem::Impl& S = *ps.impl;
    if (S.singular) {
        std::printf("[case] %s FATAL singular: %s\n", name, S.diagnostic.c_str());
        return false;
    }

    Timer reduceTimer;
    const SpMat Kff = research::reduceFF(S.K, S.fmap, S.nf);
    const VecX Fff = reduceVec(nodalLoadVector(model), S.fmap, S.nf);
    const double reduceMs = reduceTimer.ms();

    Timer ewTimer;
    const ElementWiseKffOperator elemOp(S);
    const double ewSetupMs = ewTimer.ms();
    const SparseKffOperator sparseOp(Kff);

    std::printf("[case] %s nx=%d ny=%d stories=%d nodes=%zu elems=%zu nf=%d nnzKff=%lld elemBlocks=%zu elemEntries=%zu buildMs=%.2f factorMs=%.2f reduceMs=%.2f ewSetupMs=%.2f\n",
                name, nx, ny, stories, model.nodes.size(), S.elems.size(), S.nf,
                static_cast<long long>(Kff.nonZeros()), elemOp.blocks.size(), elemOp.entryCount,
                buildMs, factorMs, reduceMs, ewSetupMs);

    constexpr double applyTol = 5e-12;
    const ApplyStats sparseStats = validateApply(name, "assembled", sparseOp, Kff);
    const ApplyStats elemStats = validateApply(name, "elementwise", elemOp, Kff);
    bool ok = sparseStats.relL2 <= applyTol && sparseStats.relMax <= applyTol &&
              elemStats.relL2 <= applyTol && elemStats.relMax <= applyTol;

    if (withPcg) {
        Timer refTimer;
        const VecX xRef = S.ldlt.solve(Fff);
        const double refMs = refTimer.ms();
        const double refRel = static_cast<double>((Kff * xRef - Fff).norm() /
                                                  std::max<real>(real(1e-300), Fff.norm()));
        std::printf("[ref]   %s LDLT solveMs=%.2f rel=%.3e\n", name, refMs, refRel);

        constexpr double solutionTol = 1e-6;
        ok = runPcgCheck(name, "assembled", sparseOp, Kff, Fff, xRef, tol, maxIter, solutionTol) && ok;
        ok = runPcgCheck(name, "elementwise", elemOp, Kff, Fff, xRef, tol, maxIter, solutionTol) && ok;
    }

    std::printf("[case-done] %s %s\n", name, ok ? "PASS" : "FAIL");
    return ok;
}

void usage() {
    std::printf("usage: exp_matrix_free_pcg [--tol v] [--maxIter n] [--noPcg] [--medium] [--large] [--single nx ny stories]\n");
}

} // namespace

int main(int argc, char** argv) {
    double tol = 1e-10;
    int maxIter = 20000;
    bool withPcg = true;
    bool medium = false;
    bool large = false;
    bool single = false;
    int sx = 0, sy = 0, ss = 0;

    for (int i = 1; i < argc; ++i) {
        auto next = [&]() -> const char* { return (i + 1 < argc) ? argv[++i] : "0"; };
        if (!std::strcmp(argv[i], "--tol")) {
            tol = std::atof(next());
        } else if (!std::strcmp(argv[i], "--maxIter")) {
            maxIter = std::atoi(next());
        } else if (!std::strcmp(argv[i], "--noPcg")) {
            withPcg = false;
        } else if (!std::strcmp(argv[i], "--medium")) {
            medium = true;
        } else if (!std::strcmp(argv[i], "--large")) {
            large = true;
        } else if (!std::strcmp(argv[i], "--single")) {
            single = true;
            sx = std::atoi(next());
            sy = std::atoi(next());
            ss = std::atoi(next());
        } else if (!std::strcmp(argv[i], "--help") || !std::strcmp(argv[i], "-h")) {
            usage();
            return 0;
        } else {
            std::printf("[args] unknown argument: %s\n", argv[i]);
            usage();
            return 2;
        }
    }

    bool ok = true;
    if (single) {
        ok = runCase("tower-custom", sx, sy, ss, tol, maxIter, withPcg) && ok;
    } else {
        ok = runCase("tower-XS(2x1x3)", 2, 1, 3, tol, maxIter, withPcg) && ok;
        ok = runCase("tower-S(3x2x4)", 3, 2, 4, tol, maxIter, withPcg) && ok;
        if (medium || large) ok = runCase("tower-M(5x4x8)", 5, 4, 8, tol, maxIter, withPcg) && ok;
        if (large) ok = runCase("tower-L(8x6x12)", 8, 6, 12, tol, maxIter, withPcg) && ok;
    }

    std::printf("[done] %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
