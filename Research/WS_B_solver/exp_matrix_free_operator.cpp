// exp_matrix_free_operator.cpp -- research-only element-by-element operator path.
//
// Goal: prove a no-new-dependency LinearOperator shape before touching FrameCore.
// The operator stores per-element dense free-free blocks and applies K*x by
// gather -> dense matvec -> scatter. It is compared against the engine's
// assembled K_ff at machine precision, then used by a small custom PCG.

#include "research_common.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

using namespace research;

namespace {

enum PrecondKind {
    PRECOND_DIAG = 0,
    PRECOND_BLOCK6 = 1,
    PRECOND_DIAG_COARSE = 2,
    PRECOND_BLOCK6_COARSE = 3,
};

int parsePrecond(const char* raw) {
    const std::string s = raw ? raw : "";
    if (s == "diag") return PRECOND_DIAG;
    if (s == "block6") return PRECOND_BLOCK6;
    if (s == "diag_coarse") return PRECOND_DIAG_COARSE;
    if (s == "block6_coarse") return PRECOND_BLOCK6_COARSE;
    return PRECOND_BLOCK6_COARSE;
}

const char* precondName(int kind) {
    switch (kind) {
    case PRECOND_DIAG: return "diag";
    case PRECOND_BLOCK6: return "block6";
    case PRECOND_DIAG_COARSE: return "diag_coarse";
    default: return "block6_coarse";
    }
}

struct Args {
    std::string preset = "small"; // small | xxl | mega
    int repeat = 20;
    int pcgMaxIter = 5000;
    double pcgTol = 1e-10;
    int precond = PRECOND_BLOCK6_COARSE;
    int coarseBinsX = 1;
    int coarseBinsY = 1;
    int maxCoarseDofs = 4096;
};

Args parseArgs(int argc, char** argv) {
    Args a;
    for (int i = 1; i < argc; ++i) {
        auto next = [&]() -> const char* { return (i + 1 < argc) ? argv[++i] : ""; };
        if      (!std::strcmp(argv[i], "--preset")) a.preset = next();
        else if (!std::strcmp(argv[i], "--repeat")) a.repeat = std::max(1, std::atoi(next()));
        else if (!std::strcmp(argv[i], "--pcgMaxIter")) a.pcgMaxIter = std::max(1, std::atoi(next()));
        else if (!std::strcmp(argv[i], "--pcgTol")) a.pcgTol = std::atof(next());
        else if (!std::strcmp(argv[i], "--precond")) a.precond = parsePrecond(next());
        else if (!std::strcmp(argv[i], "--coarseBinsX")) a.coarseBinsX = std::max(1, std::atoi(next()));
        else if (!std::strcmp(argv[i], "--coarseBinsY")) a.coarseBinsY = std::max(1, std::atoi(next()));
        else if (!std::strcmp(argv[i], "--maxCoarseDofs")) a.maxCoarseDofs = std::max(6, std::atoi(next()));
    }
    return a;
}

FrameModel makeCase(const std::string& preset) {
    if (preset == "xxl")  return makeTower(12, 9, 24);
    if (preset == "mega") return makeTower(18, 14, 36);
    return makeTower(5, 4, 8);
}

double relNorm(const VecX& a, const VecX& b) {
    return static_cast<double>((a - b).norm() / std::max<real>(real(1e-300), b.norm()));
}

double checksum(const VecX& v) {
    double s = 0.0;
    const int step = std::max(1, static_cast<int>(v.size()) / 4096);
    for (int i = 0; i < v.size(); i += step) s += std::abs(static_cast<double>(v(i)));
    return s;
}

struct ElementBlock {
    int n = 0;
    std::array<int, 12> rid{};       // reduced free dofs touched by this element
    std::array<real, 144> k{};       // row-major n x n block, fixed storage for the hot loop
};

struct PreconditionerScratch {
    VecX rc;
    VecX zc;
};

struct ElementOperator {
    int nf = 0;
    std::vector<ElementBlock> blocks;
    VecX diag;
    std::vector<std::array<real, 36>> blockInv6;
    std::vector<int> reducedBlockToCoarseGroup;
    MatX coarseInv;
    int precond = PRECOND_BLOCK6_COARSE;
    int coarseBinsX = 1;
    int coarseBinsY = 1;
    int maxCoarseDofs = 4096;

    bool build(const FrameModel& model, const PreparedSystem::Impl& S, std::string& why) {
        nf = S.nf;
        diag = VecX::Zero(nf);
        std::vector<std::array<real, 36>> diag6;
        if (nf % 6 == 0) {
            diag6.resize(static_cast<size_t>(nf / 6));
            for (auto& b : diag6) b.fill(real(0));
        }
        blocks.clear();
        blocks.reserve(model.members.size());

        SolveOptions opts;
        for (size_t e = 0; e < model.members.size(); ++e) {
            if (!model.members[e].active) continue;
            std::vector<Triplet> trips;
            if (!memberGlobalK(model, static_cast<int>(e), opts, trips, why)) return false;

            int gd[12];
            memberGlobalDofs(model, static_cast<int>(e), gd);
            std::vector<int> keep;
            std::vector<int> rid;
            keep.reserve(12);
            rid.reserve(12);
            for (int i = 0; i < 12; ++i) {
                const int f = S.fmap[static_cast<size_t>(gd[i])];
                if (f >= 0) {
                    keep.push_back(i);
                    rid.push_back(f);
                }
            }
            const int n = static_cast<int>(rid.size());
            if (n == 0) continue;

            auto localIndex = [&](int g) -> int {
                for (int i = 0; i < 12; ++i) if (gd[i] == g) return i;
                return -1;
            };

            MatX k = MatX::Zero(n, n);
            for (const Triplet& t : trips) {
                const int lr = localIndex(static_cast<int>(t.row()));
                const int lc = localIndex(static_cast<int>(t.col()));
                if (lr < 0 || lc < 0) continue;
                const auto ir = std::find(keep.begin(), keep.end(), lr);
                const auto ic = std::find(keep.begin(), keep.end(), lc);
                if (ir == keep.end() || ic == keep.end()) continue;
                const int rr = static_cast<int>(ir - keep.begin());
                const int cc = static_cast<int>(ic - keep.begin());
                k(rr, cc) += t.value();
            }
            for (int i = 0; i < n; ++i) diag(rid[static_cast<size_t>(i)]) += k(i, i);
            if (!diag6.empty()) {
                for (int r = 0; r < n; ++r) {
                    const int rr = rid[static_cast<size_t>(r)];
                    const int br = rr / 6;
                    const int lr = rr - 6 * br;
                    for (int c = 0; c < n; ++c) {
                        const int cc = rid[static_cast<size_t>(c)];
                        const int bc = cc / 6;
                        if (br != bc) continue;
                        const int lc = cc - 6 * bc;
                        diag6[static_cast<size_t>(br)][static_cast<size_t>(lr * 6 + lc)] += k(r, c);
                    }
                }
            }
            ElementBlock b;
            b.n = n;
            for (int i = 0; i < n; ++i) b.rid[static_cast<size_t>(i)] = rid[static_cast<size_t>(i)];
            for (int r = 0; r < n; ++r)
                for (int c = 0; c < n; ++c)
                    b.k[static_cast<size_t>(r * 12 + c)] = k(r, c);
            blocks.push_back(std::move(b));
        }

        buildFloorCoarse(model, S);
        blockInv6.clear();
        if (!diag6.empty()) {
            blockInv6.resize(diag6.size());
            using Mat6 = Eigen::Matrix<real, 6, 6>;
            const Mat6 I = Mat6::Identity();
            bool ok = true;
            for (size_t bi = 0; bi < diag6.size(); ++bi) {
                Mat6 A;
                for (int r = 0; r < 6; ++r)
                    for (int c = 0; c < 6; ++c)
                        A(r, c) = diag6[bi][static_cast<size_t>(r * 6 + c)];
                Eigen::LDLT<Mat6> ldlt(A);
                if (ldlt.info() != Eigen::Success) { ok = false; break; }
                const Mat6 inv = ldlt.solve(I);
                if (!inv.allFinite()) { ok = false; break; }
                for (int r = 0; r < 6; ++r)
                    for (int c = 0; c < 6; ++c)
                        blockInv6[bi][static_cast<size_t>(r * 6 + c)] = inv(r, c);
            }
            if (!ok) blockInv6.clear();
        }
        return true;
    }

    void buildFloorCoarse(const FrameModel& model, const PreparedSystem::Impl& S) {
        reducedBlockToCoarseGroup.clear();
        coarseInv.resize(0, 0);
        if (nf % 6 != 0) return;
        const int nb = nf / 6;
        reducedBlockToCoarseGroup.assign(static_cast<size_t>(nb), -1);

        auto reducedBlockForNode = [&](size_t ni) -> int {
            int block = -1;
            for (int d = 0; d < 6; ++d) {
                const int g = gdof(static_cast<int>(ni), d);
                if (g < 0 || g >= static_cast<int>(S.fmap.size())) return -1;
                const int f = S.fmap[static_cast<size_t>(g)];
                if (f < 0) return -1;
                if (d == 0) {
                    if (f % 6 != 0) return -1;
                    block = f / 6;
                } else if (f != block * 6 + d) {
                    return -1;
                }
            }
            return block;
        };

        std::vector<real> levels;
        real minX = std::numeric_limits<real>::max();
        real maxX = -std::numeric_limits<real>::max();
        real minY = std::numeric_limits<real>::max();
        real maxY = -std::numeric_limits<real>::max();
        for (size_t ni = 0; ni < model.nodes.size(); ++ni) {
            if (reducedBlockForNode(ni) < 0) continue;
            minX = std::min(minX, model.nodes[ni].pos.x);
            maxX = std::max(maxX, model.nodes[ni].pos.x);
            minY = std::min(minY, model.nodes[ni].pos.y);
            maxY = std::max(maxY, model.nodes[ni].pos.y);
            const real z = model.nodes[ni].pos.z;
            bool found = false;
            for (real existing : levels) {
                if (std::abs(existing - z) <= real(1e-7) * std::max<real>(real(1), std::abs(z))) {
                    found = true;
                    break;
                }
            }
            if (!found) levels.push_back(z);
        }
        std::sort(levels.begin(), levels.end());
        if (levels.empty()) return;

        auto groupForZ = [&](real z) -> int {
            int best = 0;
            real bestDist = std::abs(z - levels[0]);
            for (int i = 1; i < static_cast<int>(levels.size()); ++i) {
                const real dist = std::abs(z - levels[static_cast<size_t>(i)]);
                if (dist < bestDist) { best = i; bestDist = dist; }
            }
            return best;
        };

        auto binFor = [](real value, real lo, real hi, int bins) -> int {
            if (bins <= 1 || !(hi > lo)) return 0;
            const real t = (value - lo) / (hi - lo);
            int b = static_cast<int>(std::floor(static_cast<double>(t * bins)));
            if (b < 0) b = 0;
            if (b >= bins) b = bins - 1;
            return b;
        };

        const int binsX = std::max(1, coarseBinsX);
        const int binsY = std::max(1, coarseBinsY);
        const int groupsPerLevel = binsX * binsY;
        const int nc = static_cast<int>(levels.size()) * groupsPerLevel * 6;
        if (nc > maxCoarseDofs) {
            reducedBlockToCoarseGroup.clear();
            return;
        }

        for (size_t ni = 0; ni < model.nodes.size(); ++ni) {
            const int block = reducedBlockForNode(ni);
            if (block >= 0 && block < nb) {
                const int floorGroup = groupForZ(model.nodes[ni].pos.z);
                const int xb = binFor(model.nodes[ni].pos.x, minX, maxX, binsX);
                const int yb = binFor(model.nodes[ni].pos.y, minY, maxY, binsY);
                reducedBlockToCoarseGroup[static_cast<size_t>(block)] =
                    floorGroup * groupsPerLevel + yb * binsX + xb;
            }
        }
        for (int bi = 0; bi < nb; ++bi) {
            if (reducedBlockToCoarseGroup[static_cast<size_t>(bi)] >= 0) continue;
            reducedBlockToCoarseGroup.clear();
            return;
        }

        MatX Kc = MatX::Zero(nc, nc);
        for (const ElementBlock& b : blocks) {
            for (int r = 0; r < b.n; ++r) {
                const int rr = b.rid[static_cast<size_t>(r)];
                const int gr = reducedBlockToCoarseGroup[static_cast<size_t>(rr / 6)];
                if (gr < 0) continue;
                const int cr = gr * 6 + (rr % 6);
                for (int c = 0; c < b.n; ++c) {
                    const int cc = b.rid[static_cast<size_t>(c)];
                    const int gc = reducedBlockToCoarseGroup[static_cast<size_t>(cc / 6)];
                    if (gc < 0) continue;
                    const int col = gc * 6 + (cc % 6);
                    Kc(cr, col) += b.k[static_cast<size_t>(r * 12 + c)];
                }
            }
        }
        Eigen::LDLT<MatX> ldlt(Kc);
        if (ldlt.info() != Eigen::Success) return;
        coarseInv = ldlt.solve(MatX::Identity(nc, nc));
        if (!coarseInv.allFinite()) coarseInv.resize(0, 0);
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

    void jacobi(const VecX& r, VecX& z, PreconditionerScratch& scratch) const {
        z.resize(nf);
        const bool useBlock = (precond == PRECOND_BLOCK6 || precond == PRECOND_BLOCK6_COARSE) && !blockInv6.empty();
        const bool useCoarse = (precond == PRECOND_DIAG_COARSE || precond == PRECOND_BLOCK6_COARSE)
                            && coarseInv.size() > 0 && !reducedBlockToCoarseGroup.empty();
        if (useBlock) {
            for (int bi = 0; bi < static_cast<int>(blockInv6.size()); ++bi) {
                const real* inv = blockInv6[static_cast<size_t>(bi)].data();
                const int off = bi * 6;
                for (int row = 0; row < 6; ++row) {
                    const real* invRow = &inv[row * 6];
                    z(off + row) = invRow[0] * r(off + 0) + invRow[1] * r(off + 1) + invRow[2] * r(off + 2)
                                 + invRow[3] * r(off + 3) + invRow[4] * r(off + 4) + invRow[5] * r(off + 5);
                }
            }
        } else {
            for (int i = 0; i < nf; ++i) z(i) = (diag(i) != real(0)) ? r(i) / diag(i) : r(i);
        }
        if (useCoarse) {
            const int nc = static_cast<int>(coarseInv.rows());
            scratch.rc.resize(nc);
            scratch.rc.setZero();
            for (int bi = 0; bi < static_cast<int>(reducedBlockToCoarseGroup.size()); ++bi) {
                const int g = reducedBlockToCoarseGroup[static_cast<size_t>(bi)];
                if (g < 0) continue;
                for (int d = 0; d < 6; ++d) scratch.rc(g * 6 + d) += r(bi * 6 + d);
            }
            scratch.zc.resize(nc);
            scratch.zc.noalias() = coarseInv * scratch.rc;
            for (int bi = 0; bi < static_cast<int>(reducedBlockToCoarseGroup.size()); ++bi) {
                const int g = reducedBlockToCoarseGroup[static_cast<size_t>(bi)];
                if (g < 0) continue;
                for (int d = 0; d < 6; ++d) z(bi * 6 + d) += scratch.zc(g * 6 + d);
            }
            return;
        }
    }

    void jacobi(const VecX& r, VecX& z) const {
        PreconditionerScratch scratch;
        jacobi(r, z, scratch);
    }
};

struct PcgStats {
    int iters = 0;
    double relResidual = 0.0;
    double ms = 0.0;
};

PcgStats pcg(const ElementOperator& A, const VecX& b, VecX& x, int maxIter, double tol) {
    PcgStats st;
    VecX Ax(b.size());
    A.apply(x, Ax);
    VecX r = b - Ax;
    VecX z;
    PreconditionerScratch precondScratch;
    A.jacobi(r, z, precondScratch);
    VecX p = z;
    real rzOld = r.dot(z);
    const real bnorm = std::max<real>(real(1e-300), b.norm());
    VecX Ap(b.size());
    Timer t;
    for (int it = 0; it < maxIter; ++it) {
        A.apply(p, Ap);
        const real denom = p.dot(Ap);
        if (!(denom > real(0)) || !std::isfinite(static_cast<double>(denom))) break;
        const real alpha = rzOld / denom;
        x += alpha * p;
        r -= alpha * Ap;
        st.relResidual = static_cast<double>(r.norm() / bnorm);
        st.iters = it + 1;
        if (st.relResidual <= tol) break;
        A.jacobi(r, z, precondScratch);
        const real rzNew = r.dot(z);
        if (!(rzNew > real(0)) || !std::isfinite(static_cast<double>(rzNew))) break;
        const real beta = rzNew / rzOld;
        p = z + beta * p;
        rzOld = rzNew;
    }
    st.ms = t.ms();
    return st;
}

} // namespace

int main(int argc, char** argv) {
    const Args args = parseArgs(argc, argv);
    FrameModel model = makeCase(args.preset);
    assertNodalOnly(model, "exp_matrix_free_operator");

    Timer tPrep;
    PreparedSystem ps = assembleAndFactor(model);
    const double prepMs = tPrep.ms();
    const PreparedSystem::Impl& S = *ps.impl;
    if (S.singular) {
        std::printf("[mfop] status=singular diagnostic=\"%s\"\n", S.diagnostic.c_str());
        return 2;
    }
    const SpMat Kff = research::reduceFF(S.K, S.fmap, S.nf);
    const VecX Fff = reduceVec(nodalLoadVector(model), S.fmap, S.nf);

    std::string why;
    ElementOperator A;
    A.precond = args.precond;
    A.coarseBinsX = args.coarseBinsX;
    A.coarseBinsY = args.coarseBinsY;
    A.maxCoarseDofs = args.maxCoarseDofs;
    Timer tBuild;
    if (!A.build(model, S, why)) {
        std::printf("[mfop] status=build_failed why=\"%s\"\n", why.c_str());
        return 2;
    }
    const double opBuildMs = tBuild.ms();

    VecX x = VecX::Zero(S.nf);
    for (int i = 0; i < x.size(); ++i) x(i) = std::sin(0.001 * static_cast<double>(i + 1));
    VecX yRef = Kff * x;
    VecX yMf(S.nf);
    A.apply(x, yMf);
    const double applyRel = relNorm(yMf, yRef);

    Timer tSparse;
    double csSparse = 0.0;
    for (int r = 0; r < args.repeat; ++r) {
        yRef.noalias() = Kff * x;
        csSparse += checksum(yRef);
    }
    const double sparseMs = tSparse.ms() / args.repeat;

    Timer tMf;
    double csMf = 0.0;
    for (int r = 0; r < args.repeat; ++r) {
        A.apply(x, yMf);
        csMf += checksum(yMf);
    }
    const double mfMs = tMf.ms() / args.repeat;

    Timer tLdlt;
    const VecX xLdlt = S.ldlt.solve(Fff);
    const double ldltSolveMs = tLdlt.ms();

    VecX xPcg = VecX::Zero(S.nf);
    const PcgStats st = pcg(A, Fff, xPcg, args.pcgMaxIter, args.pcgTol);
    const double pcgVsLdlt = relNorm(xPcg, xLdlt);

    std::printf("[mfop] preset=%s precond=%s coarseBins=%dx%d nodes=%zu members=%zu nf=%d nnz=%lld blocks=%zu coarseDofs=%lld prepMs=%.3f opBuildMs=%.3f applyRel=%.3e sparseApplyMs=%.6f mfApplyMs=%.6f speedup=%.3f checksumRel=%.3e ldltSolveMs=%.6f pcgMs=%.3f pcgIters=%d pcgRel=%.3e pcgVsLdlt=%.3e\n",
                args.preset.c_str(), precondName(args.precond), args.coarseBinsX, args.coarseBinsY,
                model.nodes.size(), model.members.size(), S.nf,
                static_cast<long long>(Kff.nonZeros()), A.blocks.size(),
                static_cast<long long>(A.coarseInv.rows()), prepMs, opBuildMs,
                applyRel, sparseMs, mfMs, sparseMs / std::max(1e-300, mfMs),
                std::abs(csMf - csSparse) / std::max(1e-300, std::abs(csSparse)),
                ldltSolveMs, st.ms, st.iters, st.relResidual, pcgVsLdlt);
    return (applyRel <= 1e-10 && pcgVsLdlt <= 1e-6) ? 0 : 1;
}
