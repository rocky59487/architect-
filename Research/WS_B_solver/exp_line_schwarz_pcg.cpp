// exp_line_schwarz_pcg.cpp -- research-only line/column additive Schwarz PCG.
//
// This experiment stays entirely in Research/. It builds a fixed-storage
// element-by-element member operator, compares it against assembled K_ff, then
// runs Jacobi PCG versus Jacobi + overlapping line Schwarz PCG.

#include "research_common.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

using namespace research;

namespace {

struct Args {
    std::string preset = "small"; // small | xxl
    int pcgMaxIter = 5000;
    double pcgTol = 1e-10;
    double localPivotTol = 1e-12;
    double schwarzScale = 1.0;
    bool floorLines = true;
};

struct CaseDef {
    int nx = 5;
    int ny = 4;
    int stories = 8;
};

Args parseArgs(int argc, char** argv) {
    Args a;
    for (int i = 1; i < argc; ++i) {
        auto next = [&]() -> const char* { return (i + 1 < argc) ? argv[++i] : ""; };
        if (!std::strcmp(argv[i], "--preset")) {
            a.preset = next();
        } else if (!std::strcmp(argv[i], "--pcgMaxIter")) {
            a.pcgMaxIter = std::max(1, std::atoi(next()));
        } else if (!std::strcmp(argv[i], "--pcgTol")) {
            a.pcgTol = std::atof(next());
        } else if (!std::strcmp(argv[i], "--localPivotTol")) {
            a.localPivotTol = std::atof(next());
        } else if (!std::strcmp(argv[i], "--schwarzScale")) {
            a.schwarzScale = std::atof(next());
        } else if (!std::strcmp(argv[i], "--columnsOnly")) {
            a.floorLines = false;
        } else if (!std::strcmp(argv[i], "--floorLines")) {
            a.floorLines = true;
        }
    }
    return a;
}

CaseDef caseForPreset(const std::string& preset) {
    if (preset == "xxl") return { 12, 9, 24 };
    return { 5, 4, 8 };
}

double safeNorm(const VecX& v) {
    return std::max<double>(1e-300, static_cast<double>(v.norm()));
}

double relNorm(const VecX& a, const VecX& b) {
    return static_cast<double>((a - b).norm()) / safeNorm(b);
}

double trueRelResidual(const SpMat& K, const VecX& x, const VecX& b) {
    return static_cast<double>((K * x - b).norm()) / safeNorm(b);
}

VecX deterministicVec(int n) {
    VecX x(n);
    for (int i = 0; i < n; ++i) {
        const double k = static_cast<double>(i + 1);
        x(i) = static_cast<real>(std::sin(0.017 * k) + 0.5 * std::cos(0.011 * (k + 7.0)));
    }
    return x;
}

bool denseLdltLooksSpd(const Eigen::LDLT<MatX>& ldlt, double relTol) {
    if (ldlt.info() != Eigen::Success) return false;
    const VecX D = ldlt.vectorD();
    if (D.size() == 0) return false;

    real maxAbs = 0;
    real minAbs = std::numeric_limits<real>::max();
    for (int i = 0; i < D.size(); ++i) {
        const real di = D(i);
        if (!std::isfinite(static_cast<double>(di)) || !(di > real(0))) return false;
        const real adi = std::abs(di);
        maxAbs = std::max(maxAbs, adi);
        minAbs = std::min(minAbs, adi);
    }
    if (!(maxAbs > real(0))) return false;
    return minAbs / std::max<real>(real(1), maxAbs) >= static_cast<real>(relTol);
}

struct ElementBlock {
    int n = 0;
    std::array<int, 12> rid{};
    std::array<real, 144> k{};
};

struct ElementOperator {
    int nf = 0;
    VecX diag;
    std::vector<ElementBlock> blocks;

    bool build(const FrameModel& model, const PreparedSystem::Impl& S, std::string& why) {
        nf = S.nf;
        diag = VecX::Zero(nf);
        blocks.clear();
        blocks.reserve(model.members.size());

        SolveOptions opts;
        for (size_t e = 0; e < model.members.size(); ++e) {
            if (!model.members[e].active) continue;

            std::vector<Triplet> trips;
            if (!memberGlobalK(model, static_cast<int>(e), opts, trips, why)) return false;

            int gd[12];
            memberGlobalDofs(model, static_cast<int>(e), gd);

            std::array<int, 12> fullToReducedLocal;
            fullToReducedLocal.fill(-1);

            ElementBlock b;
            b.n = 0;
            b.rid.fill(-1);
            b.k.fill(real(0));

            for (int i = 0; i < 12; ++i) {
                const int g = gd[i];
                if (g < 0 || g >= static_cast<int>(S.fmap.size())) continue;
                const int f = S.fmap[static_cast<size_t>(g)];
                if (f < 0) continue;
                fullToReducedLocal[static_cast<size_t>(i)] = b.n;
                b.rid[static_cast<size_t>(b.n)] = f;
                ++b.n;
            }
            if (b.n == 0) continue;

            auto localOfGlobal = [&](int g) -> int {
                for (int i = 0; i < 12; ++i) {
                    if (gd[i] == g) return i;
                }
                return -1;
            };

            for (const Triplet& t : trips) {
                const int lrFull = localOfGlobal(static_cast<int>(t.row()));
                const int lcFull = localOfGlobal(static_cast<int>(t.col()));
                if (lrFull < 0 || lcFull < 0) continue;
                const int lr = fullToReducedLocal[static_cast<size_t>(lrFull)];
                const int lc = fullToReducedLocal[static_cast<size_t>(lcFull)];
                if (lr < 0 || lc < 0) continue;
                b.k[static_cast<size_t>(lr * 12 + lc)] += t.value();
            }

            for (int i = 0; i < b.n; ++i) {
                diag(b.rid[static_cast<size_t>(i)]) += b.k[static_cast<size_t>(i * 12 + i)];
            }
            blocks.push_back(std::move(b));
        }
        return true;
    }

    void apply(const VecX& x, VecX& y) const {
        y.setZero(nf);
        for (const ElementBlock& b : blocks) {
            real xe[12];
            for (int i = 0; i < b.n; ++i) xe[i] = x(b.rid[static_cast<size_t>(i)]);

            for (int r = 0; r < b.n; ++r) {
                const real* row = &b.k[static_cast<size_t>(r * 12)];
                real acc = 0;
                for (int c = 0; c < b.n; ++c) acc += row[c] * xe[c];
                y(b.rid[static_cast<size_t>(r)]) += acc;
            }
        }
    }
};

bool validJacobi(const VecX& diag) {
    for (int i = 0; i < diag.size(); ++i) {
        if (!(diag(i) > real(0)) || !std::isfinite(static_cast<double>(diag(i)))) return false;
    }
    return true;
}

struct JacobiPreconditioner {
    const VecX* diag = nullptr;

    void apply(const VecX& r, VecX& z) const {
        z.resize(r.size());
        for (int i = 0; i < r.size(); ++i) z(i) = r(i) / (*diag)(i);
    }
};

void addNodeReducedDofs(const FrameModel& model,
                        const PreparedSystem::Impl& S,
                        int nodeId,
                        std::vector<int>& dofs) {
    const int ni = model.nodeIndex(nodeId);
    if (ni < 0) return;
    for (int d = 0; d < 6; ++d) {
        const int g = gdof(ni, d);
        if (g < 0 || g >= static_cast<int>(S.fmap.size())) continue;
        const int f = S.fmap[static_cast<size_t>(g)];
        if (f >= 0) dofs.push_back(f);
    }
}

void sortUnique(std::vector<int>& v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
}

std::vector<std::vector<int>> makeLineDofSets(const FrameModel& model,
                                              const PreparedSystem::Impl& S,
                                              const CaseDef& c,
                                              bool floorLines) {
    std::vector<std::vector<int>> lines;
    lines.reserve(static_cast<size_t>((c.nx + 1) * (c.ny + 1) +
                                      c.stories * ((c.ny + 1) + (c.nx + 1))));

    for (int j = 0; j <= c.ny; ++j) {
        for (int i = 0; i <= c.nx; ++i) {
            std::vector<int> dofs;
            for (int k = 0; k <= c.stories; ++k) {
                addNodeReducedDofs(model, S, towerNodeId(i, j, k, c.nx, c.ny), dofs);
            }
            sortUnique(dofs);
            lines.push_back(std::move(dofs));
        }
    }

    if (!floorLines) return lines;

    for (int k = 1; k <= c.stories; ++k) {
        for (int j = 0; j <= c.ny; ++j) {
            std::vector<int> dofs;
            for (int i = 0; i <= c.nx; ++i) {
                addNodeReducedDofs(model, S, towerNodeId(i, j, k, c.nx, c.ny), dofs);
            }
            sortUnique(dofs);
            lines.push_back(std::move(dofs));
        }
    }

    for (int k = 1; k <= c.stories; ++k) {
        for (int i = 0; i <= c.nx; ++i) {
            std::vector<int> dofs;
            for (int j = 0; j <= c.ny; ++j) {
                addNodeReducedDofs(model, S, towerNodeId(i, j, k, c.nx, c.ny), dofs);
            }
            sortUnique(dofs);
            lines.push_back(std::move(dofs));
        }
    }

    return lines;
}

struct SchwarzSubdomain {
    std::vector<int> dofs;
    Eigen::LDLT<MatX> ldlt;
    VecX sqrtInvCount;
};

struct SchwarzPreconditioner {
    const VecX* diag = nullptr;
    int nf = 0;
    int attempted = 0;
    int skipped = 0;
    double scale = 1.0;
    std::vector<SchwarzSubdomain> subdomains;

    bool build(const ElementOperator& A,
               const std::vector<std::vector<int>>& lineDofs,
               double pivotTol,
               double schwarzScale) {
        diag = &A.diag;
        nf = A.nf;
        attempted = static_cast<int>(lineDofs.size());
        skipped = 0;
        scale = schwarzScale;
        subdomains.clear();
        subdomains.reserve(lineDofs.size());

        std::vector<int> local(nf, -1);
        for (const std::vector<int>& dofs : lineDofs) {
            const int n = static_cast<int>(dofs.size());
            if (n <= 0) {
                ++skipped;
                continue;
            }

            for (int i = 0; i < n; ++i) local[static_cast<size_t>(dofs[static_cast<size_t>(i)])] = i;

            MatX K = MatX::Zero(n, n);
            for (const ElementBlock& b : A.blocks) {
                for (int r = 0; r < b.n; ++r) {
                    const int lr = local[static_cast<size_t>(b.rid[static_cast<size_t>(r)])];
                    if (lr < 0) continue;
                    const real* row = &b.k[static_cast<size_t>(r * 12)];
                    for (int c = 0; c < b.n; ++c) {
                        const int lc = local[static_cast<size_t>(b.rid[static_cast<size_t>(c)])];
                        if (lc < 0) continue;
                        K(lr, lc) += row[c];
                    }
                }
            }

            for (int i = 0; i < n; ++i) local[static_cast<size_t>(dofs[static_cast<size_t>(i)])] = -1;

            K = real(0.5) * (K + K.transpose()).eval();

            SchwarzSubdomain sd;
            sd.dofs = dofs;
            sd.ldlt.compute(K);
            if (!denseLdltLooksSpd(sd.ldlt, pivotTol)) {
                ++skipped;
                continue;
            }

            const VecX probe = VecX::Ones(n);
            const VecX solved = sd.ldlt.solve(probe);
            if (!solved.allFinite()) {
                ++skipped;
                continue;
            }
            subdomains.push_back(std::move(sd));
        }

        std::vector<int> count(static_cast<size_t>(nf), 0);
        for (const SchwarzSubdomain& sd : subdomains) {
            for (int d : sd.dofs) ++count[static_cast<size_t>(d)];
        }

        for (SchwarzSubdomain& sd : subdomains) {
            const int n = static_cast<int>(sd.dofs.size());
            sd.sqrtInvCount.resize(n);
            for (int i = 0; i < n; ++i) {
                const int d = sd.dofs[static_cast<size_t>(i)];
                const int c = std::max(1, count[static_cast<size_t>(d)]);
                sd.sqrtInvCount(i) = static_cast<real>(1.0 / std::sqrt(static_cast<double>(c)));
            }
        }
        return true;
    }

    void apply(const VecX& r, VecX& z) const {
        z.resize(nf);
        for (int i = 0; i < nf; ++i) z(i) = r(i) / (*diag)(i);

        VecX rs;
        VecX ys;
        for (const SchwarzSubdomain& sd : subdomains) {
            const int n = static_cast<int>(sd.dofs.size());
            rs.resize(n);
            for (int i = 0; i < n; ++i) {
                rs(i) = sd.sqrtInvCount(i) * r(sd.dofs[static_cast<size_t>(i)]);
            }
            ys = sd.ldlt.solve(rs);
            if (!ys.allFinite()) continue;
            for (int i = 0; i < n; ++i) {
                z(sd.dofs[static_cast<size_t>(i)]) +=
                    static_cast<real>(scale) * sd.sqrtInvCount(i) * ys(i);
            }
        }
    }
};

struct PcgStats {
    VecX x;
    int iters = 0;
    double rel = std::numeric_limits<double>::infinity();
    double trueRel = std::numeric_limits<double>::infinity();
    double ms = 0.0;
    bool converged = false;
    bool breakdown = false;
};

template <class Preconditioner>
PcgStats pcg(const ElementOperator& A,
             const SpMat& Kff,
             const VecX& b,
             const Preconditioner& M,
             int maxIter,
             double tol) {
    PcgStats st;
    st.x = VecX::Zero(A.nf);

    VecX Ax(A.nf);
    VecX Ap(A.nf);
    VecX z(A.nf);
    A.apply(st.x, Ax);
    VecX r = b - Ax;
    const real bnorm = std::max<real>(real(1e-300), b.norm());
    st.rel = static_cast<double>(r.norm() / bnorm);
    if (st.rel <= tol) {
        st.converged = true;
        st.trueRel = trueRelResidual(Kff, st.x, b);
        return st;
    }

    Timer timer;
    M.apply(r, z);
    VecX p = z;
    real rz = r.dot(z);
    if (!(rz > real(0)) || !std::isfinite(static_cast<double>(rz))) {
        st.breakdown = true;
        st.ms = timer.ms();
        st.trueRel = trueRelResidual(Kff, st.x, b);
        return st;
    }

    for (int it = 0; it < maxIter; ++it) {
        A.apply(p, Ap);
        const real denom = p.dot(Ap);
        if (!(denom > real(0)) || !std::isfinite(static_cast<double>(denom))) {
            st.breakdown = true;
            break;
        }

        const real alpha = rz / denom;
        st.x += alpha * p;
        r -= alpha * Ap;
        st.iters = it + 1;
        st.rel = static_cast<double>(r.norm() / bnorm);
        if (st.rel <= tol) {
            st.converged = true;
            break;
        }

        M.apply(r, z);
        const real rzNext = r.dot(z);
        if (!(rzNext > real(0)) || !std::isfinite(static_cast<double>(rzNext))) {
            st.breakdown = true;
            break;
        }
        const real beta = rzNext / rz;
        p = z + beta * p;
        rz = rzNext;
    }

    st.ms = timer.ms();
    st.trueRel = trueRelResidual(Kff, st.x, b);
    return st;
}

} // namespace

int main(int argc, char** argv) {
    const Args args = parseArgs(argc, argv);
    const CaseDef c = caseForPreset(args.preset);

    FrameModel model = makeTower(c.nx, c.ny, c.stories);
    assertNodalOnly(model, "exp_line_schwarz_pcg");

    PreparedSystem ps = assembleAndFactor(model);
    const PreparedSystem::Impl& S = *ps.impl;
    if (S.singular) {
        std::printf("[line_schwarz] preset=%s status=singular diagnostic=\"%s\"\n",
                    args.preset.c_str(), S.diagnostic.c_str());
        return 2;
    }

    const SpMat Kff = research::reduceFF(S.K, S.fmap, S.nf);
    const VecX Fff = reduceVec(nodalLoadVector(model), S.fmap, S.nf);

    Timer setupTimer;
    std::string why;
    ElementOperator A;
    if (!A.build(model, S, why)) {
        std::printf("[line_schwarz] preset=%s status=build_failed why=\"%s\"\n",
                    args.preset.c_str(), why.c_str());
        return 2;
    }
    if (!validJacobi(A.diag)) {
        std::printf("[line_schwarz] preset=%s status=bad_jacobi\n", args.preset.c_str());
        return 2;
    }

    const std::vector<std::vector<int>> lineDofs = makeLineDofSets(model, S, c, args.floorLines);
    SchwarzPreconditioner schwarz;
    schwarz.build(A, lineDofs, args.localPivotTol, args.schwarzScale);
    const double setupMs = setupTimer.ms();

    const VecX probe = deterministicVec(S.nf);
    VecX yElem(S.nf);
    A.apply(probe, yElem);
    const VecX yRef = Kff * probe;
    const double applyRel = relNorm(yElem, yRef);

    const VecX xLdlt = S.ldlt.solve(Fff);

    JacobiPreconditioner jacobi;
    jacobi.diag = &A.diag;
    const PcgStats base = pcg(A, Kff, Fff, jacobi, args.pcgMaxIter, args.pcgTol);
    const PcgStats sch = pcg(A, Kff, Fff, schwarz, args.pcgMaxIter, args.pcgTol);

    const double pcgVsLdlt = relNorm(sch.x, xLdlt);
    const double maxTrueRel = std::max(base.trueRel, sch.trueRel);
    const double speedupIters = static_cast<double>(base.iters) / std::max(1, sch.iters);
    const double speedupMs = base.ms / std::max(1e-300, sch.ms);

    std::printf("[line_schwarz] preset=%s floorLines=%d nf=%d blocks=%zu subdomains=%zu skipped=%d applyRel=%.3e baseIters=%d schwarzIters=%d baseMs=%.3f schwarzMs=%.3f speedupIters=%.3f speedupMs=%.3f pcgVsLdlt=%.3e maxTrueRel=%.3e setupMs=%.3f\n",
                args.preset.c_str(), args.floorLines ? 1 : 0, S.nf, A.blocks.size(), schwarz.subdomains.size(),
                schwarz.skipped, applyRel, base.iters, sch.iters, base.ms, sch.ms,
                speedupIters, speedupMs, pcgVsLdlt, maxTrueRel, setupMs);

    const bool ok = applyRel <= 1e-10 &&
                    base.converged && sch.converged &&
                    maxTrueRel <= 1e-8 &&
                    pcgVsLdlt <= 1e-6;
    return ok ? 0 : 1;
}
