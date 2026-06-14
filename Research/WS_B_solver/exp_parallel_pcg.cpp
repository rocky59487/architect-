// exp_parallel_pcg.cpp -- research-only persistent-thread element PCG.
//
// This isolates the hardware question: if matrix-free element apply is executed
// through a persistent thread pool with thread-local y buffers, how much of the
// apply microbenchmark speedup survives inside a full PCG loop?

#include "research_common.h"

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <array>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <limits>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

using namespace research;

namespace {

struct Args {
    std::string preset = "small"; // small | xxl
    int threads = 4;
    int pcgMaxIter = 5000;
    double pcgTol = 1e-10;
    std::string precond = "diag"; // diag | block6 | block6_coarse
    int coarseBinsX = 1;
    int coarseBinsY = 1;
    int maxCoarseDofs = 4096;
    int rhsCount = 1;
    int basisMax = 0;
    bool seedLoadBasis = false;
    bool parallelPrecond = false;     // route block6 precond map through the thread pool
    std::string coarseSolve = "dense"; // dense | banded (block-tridiagonal floor solve)
    int deflation = 0;                // k = number of harvested Ritz vectors (0 = off)
    int deflationWindow = 0;          // eigCG Lanczos window m (0 = auto = max(64, 6k))
    bool symApply = false;            // store only the lower triangle of each element block
    int towerNx = 0, towerNy = 0, towerStories = 0;  // override tower dims for scaling sweeps
    bool noSerialBaseline = false;    // skip the per-RHS serial PCG reference (large-batch sweeps)
    // Real game-engine load family (default off = original synthetic family). Validates the
    // seeded lane under gravity + moving sparse contact point loads instead of a smooth combo
    // of the same modes it seeds. See makeGameLoadFamily.
    bool gameLoad = false;
    int gameContactCandidates = 12;   // seedable contact nodes (the in-subspace set)
    int gameContactActive = 2;        // active contacts per frame
    double gameOutFraction = 0.0;     // [0,1] fraction of frames whose contacts leave the seed span
    bool gameAddHorizontal = false;   // add a small Ux/Uy component to contacts
    bool gameVerify = false;          // dense V (VtKV)^-1 Vt b cross-check of initialGuess (small cases)
};

Args parseArgs(int argc, char** argv) {
    Args a;
    for (int i = 1; i < argc; ++i) {
        auto next = [&]() -> const char* { return (i + 1 < argc) ? argv[++i] : ""; };
        if      (!std::strcmp(argv[i], "--preset")) a.preset = next();
        else if (!std::strcmp(argv[i], "--threads")) a.threads = std::max(1, std::atoi(next()));
        else if (!std::strcmp(argv[i], "--pcgMaxIter")) a.pcgMaxIter = std::max(1, std::atoi(next()));
        else if (!std::strcmp(argv[i], "--pcgTol")) a.pcgTol = std::atof(next());
        else if (!std::strcmp(argv[i], "--precond")) a.precond = next();
        else if (!std::strcmp(argv[i], "--coarseBinsX")) a.coarseBinsX = std::max(1, std::atoi(next()));
        else if (!std::strcmp(argv[i], "--coarseBinsY")) a.coarseBinsY = std::max(1, std::atoi(next()));
        else if (!std::strcmp(argv[i], "--maxCoarseDofs")) a.maxCoarseDofs = std::max(6, std::atoi(next()));
        else if (!std::strcmp(argv[i], "--rhs")) a.rhsCount = std::max(1, std::atoi(next()));
        else if (!std::strcmp(argv[i], "--basisMax")) a.basisMax = std::max(0, std::atoi(next()));
        else if (!std::strcmp(argv[i], "--seedLoadBasis")) a.seedLoadBasis = true;
        else if (!std::strcmp(argv[i], "--parallelPrecond")) a.parallelPrecond = true;
        else if (!std::strcmp(argv[i], "--coarseSolve")) a.coarseSolve = next();
        else if (!std::strcmp(argv[i], "--deflation")) a.deflation = std::max(0, std::atoi(next()));
        else if (!std::strcmp(argv[i], "--deflationWindow")) a.deflationWindow = std::max(0, std::atoi(next()));
        else if (!std::strcmp(argv[i], "--symApply")) a.symApply = true;
        else if (!std::strcmp(argv[i], "--towerNx")) a.towerNx = std::max(0, std::atoi(next()));
        else if (!std::strcmp(argv[i], "--towerNy")) a.towerNy = std::max(0, std::atoi(next()));
        else if (!std::strcmp(argv[i], "--towerStories")) a.towerStories = std::max(0, std::atoi(next()));
        else if (!std::strcmp(argv[i], "--noSerialBaseline")) a.noSerialBaseline = true;
        else if (!std::strcmp(argv[i], "--gameLoad")) a.gameLoad = true;
        else if (!std::strcmp(argv[i], "--gameContactCandidates")) a.gameContactCandidates = std::max(1, std::atoi(next()));
        else if (!std::strcmp(argv[i], "--gameContactActive")) a.gameContactActive = std::max(1, std::atoi(next()));
        else if (!std::strcmp(argv[i], "--gameOutFraction")) a.gameOutFraction = std::min(1.0, std::max(0.0, std::atof(next())));
        else if (!std::strcmp(argv[i], "--gameAddHorizontal")) a.gameAddHorizontal = true;
        else if (!std::strcmp(argv[i], "--gameVerify")) a.gameVerify = true;
    }
    if (a.preset != "xxl") a.preset = "small";
    return a;
}

FrameModel makeCase(const Args& a) {
    if (a.towerNx > 0 && a.towerNy > 0 && a.towerStories > 0)
        return makeTower(a.towerNx, a.towerNy, a.towerStories);
    if (a.preset == "xxl") return makeTower(12, 9, 24);
    return makeTower(5, 4, 8);
}

struct FreeDofMap {
    int N = 0;
    int nf = 0;
    std::vector<int> fmap;
    bool singular = false;
    std::string diagnostic;
};

FreeDofMap makeFreeDofMap(const FrameModel& model) {
    FreeDofMap out;
    std::string why;
    if (!model.validate(why)) {
        out.singular = true;
        out.diagnostic = "invalid model: " + why;
        out.N = std::max(0, model.dofCount());
        return out;
    }
    out.N = model.dofCount();
    out.fmap.assign(static_cast<size_t>(out.N), -1);
    for (size_t k = 0; k < model.nodes.size(); ++k) {
        for (int d = 0; d < 6; ++d) {
            const int g = gdof(static_cast<int>(k), d);
            if (!model.nodes[k].fixed[d]) out.fmap[static_cast<size_t>(g)] = out.nf++;
        }
    }
    if (out.nf == 0) {
        out.singular = true;
        out.diagnostic = "fully constrained (no free DOF)";
    }
    return out;
}

double relNorm(const VecX& a, const VecX& b) {
    return static_cast<double>((a - b).norm() / std::max<real>(real(1e-300), b.norm()));
}

struct ElementBlock12 {
    std::array<int, 12> rid{};
    std::array<real, 144> k{};
};

// Symmetric element block: only the lower triangle (incl. diagonal) of the 12x12 is
// stored (78 reals vs 144), halving the element-data streamed each matrix-free apply.
// klo[i*(i+1)/2 + j] holds K_e(i,j) for j <= i.
struct ElementBlockSym {
    std::array<int, 12> rid{};
    std::array<real, 78> klo{};
};

struct ElementOperator {
    int nf = 0;
    bool symmetric = false;
    std::vector<ElementBlock12> blocks;      // full storage (symmetric == false)
    std::vector<ElementBlockSym> blocksSym;  // lower-triangle storage (symmetric == true)
    VecX diag;
    std::vector<std::array<real, 36>> diag6;

    size_t numBlocks() const { return symmetric ? blocksSym.size() : blocks.size(); }
    const std::array<int, 12>& ridOf(size_t bi) const {
        return symmetric ? blocksSym[bi].rid : blocks[bi].rid;
    }
    // (r,c) entry of element block bi (reconstructed from the lower triangle when symmetric).
    real kEntry(size_t bi, int r, int c) const {
        if (symmetric) {
            const int i = r >= c ? r : c, j = r >= c ? c : r;
            return blocksSym[bi].klo[static_cast<size_t>(i * (i + 1) / 2 + j)];
        }
        return blocks[bi].k[static_cast<size_t>(r * 12 + c)];
    }

    bool build(const FrameModel& model, const std::vector<int>& fmap, int freeCount, bool sym,
               std::string& why) {
        nf = freeCount;
        symmetric = sym;
        diag = VecX::Zero(nf);
        diag6.clear();
        if (nf % 6 == 0) {
            diag6.resize(static_cast<size_t>(nf / 6));
            for (auto& b : diag6) b.fill(real(0));
        }
        blocks.clear();
        blocksSym.clear();
        blocks.reserve(model.members.size());
        if (symmetric) blocksSym.reserve(model.members.size());

        SolveOptions opts;
        for (size_t e = 0; e < model.members.size(); ++e) {
            if (!model.members[e].active) continue;
            std::vector<Triplet> trips;
            if (!memberGlobalK(model, static_cast<int>(e), opts, trips, why)) return false;

            int gd[12];
            memberGlobalDofs(model, static_cast<int>(e), gd);

            ElementBlock12 b;
            b.rid.fill(-1);
            b.k.fill(real(0));
            int freeDofs = 0;
            for (int i = 0; i < 12; ++i) {
                const int g = gd[i];
                if (g >= 0 && g < static_cast<int>(fmap.size())) {
                    b.rid[static_cast<size_t>(i)] = fmap[static_cast<size_t>(g)];
                    if (b.rid[static_cast<size_t>(i)] >= 0) ++freeDofs;
                }
            }
            if (freeDofs == 0) continue;

            auto localIndex = [&](int g) -> int {
                for (int i = 0; i < 12; ++i) if (gd[i] == g) return i;
                return -1;
            };

            for (const Triplet& t : trips) {
                const int lr = localIndex(static_cast<int>(t.row()));
                const int lc = localIndex(static_cast<int>(t.col()));
                if (lr >= 0 && lc >= 0)
                    b.k[static_cast<size_t>(lr * 12 + lc)] += t.value();
            }
            for (int r = 0; r < 12; ++r) {
                const int rr = b.rid[static_cast<size_t>(r)];
                if (rr >= 0) diag(rr) += b.k[static_cast<size_t>(r * 12 + r)];
                if (rr < 0 || diag6.empty()) continue;
                const int br = rr / 6;
                const int lr = rr - br * 6;
                for (int c = 0; c < 12; ++c) {
                    const int cc = b.rid[static_cast<size_t>(c)];
                    if (cc < 0 || cc / 6 != br) continue;
                    const int lc = cc - br * 6;
                    diag6[static_cast<size_t>(br)][static_cast<size_t>(lr * 6 + lc)] +=
                        b.k[static_cast<size_t>(r * 12 + c)];
                }
            }
            if (symmetric) {
                ElementBlockSym bs;
                bs.rid = b.rid;
                // symmetrize against tiny asymmetry, store lower triangle
                for (int i = 0; i < 12; ++i)
                    for (int j = 0; j <= i; ++j)
                        bs.klo[static_cast<size_t>(i * (i + 1) / 2 + j)] =
                            real(0.5) * (b.k[static_cast<size_t>(i * 12 + j)] +
                                         b.k[static_cast<size_t>(j * 12 + i)]);
                blocksSym.push_back(std::move(bs));
            } else {
                blocks.push_back(std::move(b));
            }
        }
        return true;
    }

    void accumulateRange(size_t begin, size_t end, const VecX& x, std::vector<real>& y) const {
        if (symmetric) {
            for (size_t bi = begin; bi < end; ++bi) {
                const ElementBlockSym& b = blocksSym[bi];
                real xe[12], ye[12];
                for (int i = 0; i < 12; ++i) {
                    const int id = b.rid[static_cast<size_t>(i)];
                    xe[i] = id >= 0 ? x(id) : real(0);
                    ye[i] = real(0);
                }
                for (int i = 0; i < 12; ++i) {
                    const real* krow = &b.klo[static_cast<size_t>(i * (i + 1) / 2)];
                    const real xi = xe[i];
                    real yi = krow[i] * xi;  // diagonal
                    for (int j = 0; j < i; ++j) {
                        const real kij = krow[j];
                        yi += kij * xe[j];
                        ye[j] += kij * xi;
                    }
                    ye[i] += yi;
                }
                for (int i = 0; i < 12; ++i) {
                    const int rr = b.rid[static_cast<size_t>(i)];
                    if (rr >= 0) y[static_cast<size_t>(rr)] += ye[i];
                }
            }
            return;
        }
        for (size_t bi = begin; bi < end; ++bi) {
            const ElementBlock12& b = blocks[bi];
            real xe[12];
            for (int i = 0; i < 12; ++i) {
                const int id = b.rid[static_cast<size_t>(i)];
                xe[i] = id >= 0 ? x(id) : real(0);
            }
            for (int r = 0; r < 12; ++r) {
                const int rr = b.rid[static_cast<size_t>(r)];
                if (rr < 0) continue;
                const real* row = &b.k[static_cast<size_t>(r * 12)];
                real acc = real(0);
                for (int c = 0; c < 12; ++c) acc += row[c] * xe[c];
                y[static_cast<size_t>(rr)] += acc;
            }
        }
    }

    void applySerial(const VecX& x, VecX& y) const {
        std::vector<real> yy(static_cast<size_t>(nf), real(0));
        accumulateRange(0, numBlocks(), x, yy);
        y.resize(nf);
        for (int i = 0; i < nf; ++i) y(i) = yy[static_cast<size_t>(i)];
    }
};

struct PreconditionerScratch {
    VecX rc;
    VecX zc;
    VecX cg;  // banded coarse: forward-sweep temp
    VecX ch;  // banded coarse: diagonal-solve temp
};

struct Preconditioner {
    const ElementOperator* A = nullptr;
    bool useBlock6 = false;
    bool useCoarse = false;
    std::vector<std::array<real, 36>> blockInv6;
    std::vector<int> reducedBlockToCoarseGroup;
    MatX coarseInv;
    int coarseDim = 0;

    // Optional parallel block6 map (set by main() once the thread pool exists).
    std::function<void(const VecX&, VecX&, const std::vector<std::array<real, 36>>&)> parallelBlock6;
    bool wantParallel = false;

    // Banded (block-tridiagonal) coarse solve: stores a block-Thomas LDL^T factor of
    // the floor-aggregated coarse matrix instead of a dense inverse. The coarse DOFs
    // are floor-contiguous (group = floor*groupsPerLevel + bin), and frame elements
    // only couple same/adjacent floors, so the coarse matrix is block-tridiagonal.
    bool bandedRequested = false;
    bool bandedOk = false;
    int coarseFb = 0;      // floor-block size = groupsPerLevel * 6
    int coarseFloors = 0;  // number of floor blocks
    double bandedResidual = 0;  // one-time self-check: ||Kc*bandedSolve(probe)-probe||/||probe||
    double bandedOffBandRel = 0;  // max off-tridiagonal block magnitude / max|Kc| (should be ~0)
    std::vector<Eigen::LDLT<MatX>> coarseDfac;  // diagonal floor-block factors D_i
    std::vector<MatX> coarseEsub;               // sub-diagonal multipliers E_i (i>=1)

    // Sparse coarse solve: factor the (grid-structured, sparse) coarse matrix with the
    // engine's own SimplicialLDLT. Exploits ALL sparsity, so fine coarse spaces (3x3,
    // 4x4, 8x8) cost far less than the dense-per-floor block-Thomas (~fb^2) -> enables a
    // much richer coarse space at low per-apply cost.
    bool sparseRequested = false;
    bool sparseOk = false;
    LDLTSolver coarseSparseFac;
    double coarseNnz = 0;  // nonzeros in the sparse coarse matrix (cache-footprint proxy)

    void build(const ElementOperator& op, const FrameModel& model, const std::vector<int>& fmap,
               const Args& args) {
        A = &op;
        useBlock6 = args.precond == "block6" || args.precond == "block6_coarse";
        useCoarse = args.precond == "block6_coarse";
        wantParallel = args.parallelPrecond;
        bandedRequested = args.coarseSolve == "banded";
        sparseRequested = args.coarseSolve == "sparse";
        buildBlock6();
        if (useCoarse) buildCoarse(model, fmap, args.coarseBinsX, args.coarseBinsY, args.maxCoarseDofs);
    }

    void apply(const VecX& r, VecX& z, PreconditionerScratch& scratch, bool parallel) const {
        z.resize(A->nf);
        if (useBlock6 && !blockInv6.empty()) {
            if (parallel && parallelBlock6) {
                parallelBlock6(r, z, blockInv6);
            } else {
                for (int bi = 0; bi < static_cast<int>(blockInv6.size()); ++bi) {
                    const real* inv = blockInv6[static_cast<size_t>(bi)].data();
                    const int off = bi * 6;
                    for (int row = 0; row < 6; ++row) {
                        const real* invRow = &inv[row * 6];
                        z(off + row) = invRow[0] * r(off + 0) + invRow[1] * r(off + 1) + invRow[2] * r(off + 2)
                                     + invRow[3] * r(off + 3) + invRow[4] * r(off + 4) + invRow[5] * r(off + 5);
                    }
                }
            }
        } else {
            for (int i = 0; i < A->nf; ++i) z(i) = (A->diag(i) != real(0)) ? r(i) / A->diag(i) : r(i);
        }
        if (useCoarse && coarseDim > 0 && !reducedBlockToCoarseGroup.empty()) {
            const int nc = coarseDim;
            scratch.rc.resize(nc);
            scratch.rc.setZero();
            for (int bi = 0; bi < static_cast<int>(reducedBlockToCoarseGroup.size()); ++bi) {
                const int g = reducedBlockToCoarseGroup[static_cast<size_t>(bi)];
                if (g < 0) continue;
                for (int d = 0; d < 6; ++d) scratch.rc(g * 6 + d) += r(bi * 6 + d);
            }
            scratch.zc.resize(nc);
            if (sparseOk) scratch.zc = coarseSparseFac.solve(scratch.rc);
            else if (bandedOk) bandedSolve(scratch);
            else scratch.zc.noalias() = coarseInv * scratch.rc;
            for (int bi = 0; bi < static_cast<int>(reducedBlockToCoarseGroup.size()); ++bi) {
                const int g = reducedBlockToCoarseGroup[static_cast<size_t>(bi)];
                if (g < 0) continue;
                for (int d = 0; d < 6; ++d) z(bi * 6 + d) += scratch.zc(g * 6 + d);
            }
        }
    }

    void apply(const VecX& r, VecX& z, PreconditionerScratch& scratch) const {
        apply(r, z, scratch, wantParallel);
    }

    void apply(const VecX& r, VecX& z) const {
        PreconditionerScratch scratch;
        apply(r, z, scratch, wantParallel);
    }

private:
    // Symmetric block-tridiagonal solve via the precomputed block-Thomas LDL^T factor.
    // Forward (L g = rc), diagonal (D h = g), back (L^T zc = h); all fb-block sweeps.
    void bandedApply(const VecX& rc, VecX& cg, VecX& ch, VecX& zc) const {
        const int fb = coarseFb, L = coarseFloors;
        cg.resize(L * fb);
        ch.resize(L * fb);
        zc.resize(L * fb);
        cg.segment(0, fb) = rc.segment(0, fb);
        for (int i = 1; i < L; ++i)
            cg.segment(i * fb, fb).noalias() =
                rc.segment(i * fb, fb) - coarseEsub[static_cast<size_t>(i)] * cg.segment((i - 1) * fb, fb);
        for (int i = 0; i < L; ++i)
            ch.segment(i * fb, fb) = coarseDfac[static_cast<size_t>(i)].solve(cg.segment(i * fb, fb));
        zc.segment((L - 1) * fb, fb) = ch.segment((L - 1) * fb, fb);
        for (int i = L - 2; i >= 0; --i)
            zc.segment(i * fb, fb).noalias() =
                ch.segment(i * fb, fb) -
                coarseEsub[static_cast<size_t>(i + 1)].transpose() * zc.segment((i + 1) * fb, fb);
    }

    void bandedSolve(PreconditionerScratch& s) const { bandedApply(s.rc, s.cg, s.ch, s.zc); }

    // Build the block-Thomas LDL^T factor from the assembled dense coarse matrix Kc.
    // Returns false (caller falls back to dense) if Kc is not block-tridiagonal in the
    // floor ordering or a diagonal floor block is not SPD.
    bool buildBandedFactor(const MatX& Kc, int nc, int groupsPerLevel) {
        const int fb = groupsPerLevel * 6;
        // On any failure clear partial factors so the dense fallback never inherits stale
        // floor blocks (memory hygiene; the tower is always block-tridiagonal so this path
        // is defensive only).
        auto bail = [this]() { coarseDfac.clear(); coarseEsub.clear(); return false; };
        if (fb <= 0 || nc % fb != 0) return bail();
        const int L = nc / fb;
        const real ref = Kc.cwiseAbs().maxCoeff();
        const real tol = real(1e-9) * std::max<real>(real(1), ref);
        real maxOffBand = 0;
        for (int i = 0; i < L; ++i)
            for (int j = 0; j < L; ++j)
                if (std::abs(i - j) > 1)
                    maxOffBand = std::max(maxOffBand, Kc.block(i * fb, j * fb, fb, fb).cwiseAbs().maxCoeff());
        bandedOffBandRel = static_cast<double>(maxOffBand / std::max<real>(real(1e-300), ref));
        if (maxOffBand > tol) return bail();  // not block-tridiagonal -> dense fallback
        coarseDfac.assign(static_cast<size_t>(L), Eigen::LDLT<MatX>());
        coarseEsub.assign(static_cast<size_t>(L), MatX::Zero(fb, fb));
        coarseDfac[0].compute(Kc.block(0, 0, fb, fb));
        if (coarseDfac[0].info() != Eigen::Success) return bail();
        for (int i = 1; i < L; ++i) {
            const MatX Bi = Kc.block(i * fb, (i - 1) * fb, fb, fb);  // A_{i,i-1}
            // E_i = B_i D_{i-1}^{-1};  E_i^T = D_{i-1}^{-1} B_i^T = Dfac[i-1].solve(B_i^T)
            coarseEsub[static_cast<size_t>(i)] = coarseDfac[static_cast<size_t>(i - 1)].solve(Bi.transpose()).transpose();
            const MatX Di = Kc.block(i * fb, i * fb, fb, fb) - coarseEsub[static_cast<size_t>(i)] * Bi.transpose();
            coarseDfac[static_cast<size_t>(i)].compute(Di);
            if (coarseDfac[static_cast<size_t>(i)].info() != Eigen::Success) return bail();
        }
        coarseFb = fb;
        coarseFloors = L;
        // One-time self-check: the block-Thomas factor must solve Kc exactly. A reproducible
        // probe is run through forward/diag/back and the true coarse residual is measured.
        VecX probe(nc);
        for (int i = 0; i < nc; ++i)
            probe(i) = std::sin(real(0.013) * static_cast<real>(i + 1)) +
                       real(0.4) * std::cos(real(0.207) * static_cast<real>(i + 3)) +
                       ((i & 1) ? real(-0.25) : real(0.25));
        VecX cg, ch, zc;
        bandedApply(probe, cg, ch, zc);
        bandedResidual = static_cast<double>((Kc * zc - probe).norm() / std::max<real>(real(1e-300), probe.norm()));
        if (bandedResidual > 1e-8) return false;  // factor not accurate -> dense fallback
        return true;
    }

    void buildBlock6() {
        blockInv6.clear();
        if (!useBlock6 || A->diag6.empty()) return;
        blockInv6.resize(A->diag6.size());
        using Mat6 = Eigen::Matrix<real, 6, 6>;
        const Mat6 I = Mat6::Identity();
        bool ok = true;
        for (size_t bi = 0; bi < A->diag6.size(); ++bi) {
            Mat6 M;
            for (int r = 0; r < 6; ++r)
                for (int c = 0; c < 6; ++c)
                    M(r, c) = A->diag6[bi][static_cast<size_t>(r * 6 + c)];
            Eigen::LDLT<Mat6> ldlt(M);
            if (ldlt.info() != Eigen::Success) { ok = false; break; }
            const Mat6 inv = ldlt.solve(I);
            if (!inv.allFinite()) { ok = false; break; }
            for (int r = 0; r < 6; ++r)
                for (int c = 0; c < 6; ++c)
                    blockInv6[bi][static_cast<size_t>(r * 6 + c)] = inv(r, c);
        }
        if (!ok) blockInv6.clear();
    }

    void buildCoarse(const FrameModel& model, const std::vector<int>& fmap,
                     int coarseBinsX, int coarseBinsY, int maxCoarseDofs) {
        reducedBlockToCoarseGroup.clear();
        coarseInv.resize(0, 0);
        if (A->nf % 6 != 0) return;
        const int nb = A->nf / 6;
        reducedBlockToCoarseGroup.assign(static_cast<size_t>(nb), -1);

        auto reducedBlockForNode = [&](size_t ni) -> int {
            int block = -1;
            for (int d = 0; d < 6; ++d) {
                const int g = gdof(static_cast<int>(ni), d);
                const int f = fmap[static_cast<size_t>(g)];
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
        real minX = std::numeric_limits<real>::max(), maxX = -std::numeric_limits<real>::max();
        real minY = std::numeric_limits<real>::max(), maxY = -std::numeric_limits<real>::max();
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
            int b = static_cast<int>(std::floor(static_cast<double>(((value - lo) / (hi - lo)) * bins)));
            if (b < 0) b = 0;
            if (b >= bins) b = bins - 1;
            return b;
        };

        const int bx = std::max(1, coarseBinsX);
        const int by = std::max(1, coarseBinsY);
        const int groupsPerLevel = bx * by;
        const int nc = static_cast<int>(levels.size()) * groupsPerLevel * 6;
        // The dense inverse is O(nc^2) memory; banded/sparse factors scale far better, so
        // only the dense path is bound by maxCoarseDofs.
        const bool denseMode = !sparseRequested && !bandedRequested;
        const int hardCap = denseMode ? maxCoarseDofs : std::max(maxCoarseDofs, A->nf);
        if (nc > hardCap) {
            reducedBlockToCoarseGroup.clear();
            return;
        }

        for (size_t ni = 0; ni < model.nodes.size(); ++ni) {
            const int block = reducedBlockForNode(ni);
            if (block < 0 || block >= nb) continue;
            const int floorGroup = groupForZ(model.nodes[ni].pos.z);
            const int xb = binFor(model.nodes[ni].pos.x, minX, maxX, bx);
            const int yb = binFor(model.nodes[ni].pos.y, minY, maxY, by);
            reducedBlockToCoarseGroup[static_cast<size_t>(block)] =
                floorGroup * groupsPerLevel + yb * bx + xb;
        }
        for (int bi = 0; bi < nb; ++bi) {
            if (reducedBlockToCoarseGroup[static_cast<size_t>(bi)] >= 0) continue;
            reducedBlockToCoarseGroup.clear();
            return;
        }

        // Assemble coarse contributions once as triplets; the sparse path consumes them
        // directly, the dense/banded paths materialize a dense Kc from them.
        std::vector<Triplet> ctrip;
        ctrip.reserve(A->numBlocks() * 36);
        for (size_t bi = 0; bi < A->numBlocks(); ++bi) {
            const std::array<int, 12>& rid = A->ridOf(bi);
            for (int r = 0; r < 12; ++r) {
                const int rr = rid[static_cast<size_t>(r)];
                if (rr < 0) continue;
                const int gr = reducedBlockToCoarseGroup[static_cast<size_t>(rr / 6)];
                if (gr < 0) continue;
                const int cr = gr * 6 + rr % 6;
                for (int c = 0; c < 12; ++c) {
                    const int cc = rid[static_cast<size_t>(c)];
                    if (cc < 0) continue;
                    const int gc = reducedBlockToCoarseGroup[static_cast<size_t>(cc / 6)];
                    if (gc < 0) continue;
                    ctrip.emplace_back(cr, gc * 6 + cc % 6, A->kEntry(bi, r, c));
                }
            }
        }
        coarseDim = nc;

        if (sparseRequested) {
            SpMat Kcs(nc, nc);
            Kcs.setFromTriplets(ctrip.begin(), ctrip.end());
            Kcs.makeCompressed();
            coarseNnz = static_cast<double>(Kcs.nonZeros());
            coarseSparseFac.compute(Kcs);
            if (coarseSparseFac.info() == Eigen::Success) {
                sparseOk = true;
                return;
            }
            reducedBlockToCoarseGroup.clear();  // factor failed -> coarse off
            coarseDim = 0;
            return;
        }

        MatX Kc = MatX::Zero(nc, nc);
        for (const Triplet& t : ctrip)
            Kc(static_cast<int>(t.row()), static_cast<int>(t.col())) += t.value();
        if (bandedRequested && buildBandedFactor(Kc, nc, groupsPerLevel)) {
            bandedOk = true;
            return;  // banded factor ready; skip the dense inverse entirely
        }
        Eigen::LDLT<MatX> ldlt(Kc);
        if (ldlt.info() != Eigen::Success) {
            reducedBlockToCoarseGroup.clear();
            coarseDim = 0;
            return;
        }
        coarseInv = ldlt.solve(MatX::Identity(nc, nc));
        if (!coarseInv.allFinite()) {
            coarseInv.resize(0, 0);
            reducedBlockToCoarseGroup.clear();
            coarseDim = 0;
        }
    }
};

// Persistent thread pool that executes two data-parallel jobs over the SAME warm
// worker set: (1) matrix-free element apply with per-thread touched-DOF reduction,
// and (2) a block6 Jacobi preconditioner map z[block] = inv6[block] * r[block].
// The block6 map writes disjoint output slices (one DOF-block per worker range), so
// it needs no reduction. Reusing the warm threads avoids a second spawn/teardown.
class ThreadApplyPool {
public:
    enum class Job { Apply, Precond };

    ThreadApplyPool(const ElementOperator& opIn, int requestedThreads)
        : op(opIn), nt(std::max(1, requestedThreads)) {
        const size_t nblk = op.numBlocks();
        if (nblk > 0) nt = std::min(nt, static_cast<int>(nblk));
        localY.assign(static_cast<size_t>(nt), std::vector<real>(static_cast<size_t>(op.nf), real(0)));
        const size_t chunk = (nblk + static_cast<size_t>(nt) - 1) / static_cast<size_t>(nt);
        touched.resize(static_cast<size_t>(nt));
        // Precond DOF-block partition: contiguous ranges over [0, nf/6) blocks.
        const int nb = (op.nf % 6 == 0) ? op.nf / 6 : 0;
        const int pchunk = (nb + nt - 1) / std::max(1, nt);
        pcRange.assign(static_cast<size_t>(nt), {0, 0});
        applyRange.assign(static_cast<size_t>(nt), {0, 0});
        workers.reserve(static_cast<size_t>(nt));
        for (int tid = 0; tid < nt; ++tid) {
            const size_t begin = std::min(nblk, static_cast<size_t>(tid) * chunk);
            const size_t end = std::min(nblk, begin + chunk);
            applyRange[static_cast<size_t>(tid)] = {begin, end};
            const int pb = std::min(nb, tid * pchunk);
            const int pe = std::min(nb, pb + pchunk);
            pcRange[static_cast<size_t>(tid)] = {pb, pe};
            buildTouched(tid, begin, end);
            workers.emplace_back([this, tid]() { workerLoop(tid); });
        }
    }

    ~ThreadApplyPool() {
        {
            std::lock_guard<std::mutex> lock(mu);
            stop = true;
            ++generation;
        }
        cvStart.notify_all();
        for (std::thread& w : workers) if (w.joinable()) w.join();
    }

    int threads() const { return nt; }

    void apply(const VecX& x, VecX& y) {
        y.setZero(op.nf);
        {
            std::lock_guard<std::mutex> lock(mu);
            job = Job::Apply;
            xPtr = &x;
            done = 0;
            ++generation;
        }
        cvStart.notify_all();
        {
            std::unique_lock<std::mutex> lock(mu);
            cvDone.wait(lock, [&]() { return done == nt; });
        }
        for (int tid = 0; tid < nt; ++tid) {
            const std::vector<real>& yy = localY[static_cast<size_t>(tid)];
            for (int id : touched[static_cast<size_t>(tid)]) y(id) += yy[static_cast<size_t>(id)];
        }
    }

    // Parallel block6 Jacobi map: z = blockInv6 (.) r, block-diagonal. z must be sized nf.
    void precondBlock6(const VecX& r, VecX& z, const std::vector<std::array<real, 36>>& inv6) {
        {
            std::lock_guard<std::mutex> lock(mu);
            job = Job::Precond;
            rPtr = &r;
            zPtr = &z;
            invPtr = &inv6;
            done = 0;
            ++generation;
        }
        cvStart.notify_all();
        std::unique_lock<std::mutex> lock(mu);
        cvDone.wait(lock, [&]() { return done == nt; });
    }

private:
    const ElementOperator& op;
    int nt = 1;
    std::vector<std::vector<real>> localY;
    std::vector<std::vector<int>> touched;
    std::vector<std::pair<size_t, size_t>> applyRange;
    std::vector<std::pair<int, int>> pcRange;
    std::vector<std::thread> workers;
    std::mutex mu;
    std::condition_variable cvStart;
    std::condition_variable cvDone;
    Job job = Job::Apply;
    const VecX* xPtr = nullptr;
    const VecX* rPtr = nullptr;
    VecX* zPtr = nullptr;
    const std::vector<std::array<real, 36>>* invPtr = nullptr;
    int generation = 0;
    int done = 0;
    bool stop = false;

    void buildTouched(int tid, size_t begin, size_t end) {
        std::vector<int> ids;
        ids.reserve((end - begin) * 12);
        for (size_t bi = begin; bi < end; ++bi) {
            const std::array<int, 12>& rid = op.ridOf(bi);
            for (int d = 0; d < 12; ++d) {
                const int id = rid[static_cast<size_t>(d)];
                if (id >= 0) ids.push_back(id);
            }
        }
        std::sort(ids.begin(), ids.end());
        ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
        touched[static_cast<size_t>(tid)] = std::move(ids);
    }

    void runApply(int tid) {
        const auto range = applyRange[static_cast<size_t>(tid)];
        std::vector<real>& yy = localY[static_cast<size_t>(tid)];
        for (int id : touched[static_cast<size_t>(tid)]) yy[static_cast<size_t>(id)] = real(0);
        op.accumulateRange(range.first, range.second, *xPtr, yy);
    }

    void runPrecond(int tid) {
        const auto range = pcRange[static_cast<size_t>(tid)];
        const VecX& r = *rPtr;
        VecX& z = *zPtr;
        const std::vector<std::array<real, 36>>& inv6 = *invPtr;
        for (int bi = range.first; bi < range.second; ++bi) {
            const real* inv = inv6[static_cast<size_t>(bi)].data();
            const int off = bi * 6;
            for (int row = 0; row < 6; ++row) {
                const real* ir = &inv[row * 6];
                z(off + row) = ir[0] * r(off + 0) + ir[1] * r(off + 1) + ir[2] * r(off + 2)
                             + ir[3] * r(off + 3) + ir[4] * r(off + 4) + ir[5] * r(off + 5);
            }
        }
    }

    void workerLoop(int tid) {
        int seen = 0;
        for (;;) {
            Job localJob;
            {
                std::unique_lock<std::mutex> lock(mu);
                cvStart.wait(lock, [&]() { return stop || generation != seen; });
                if (stop) return;
                seen = generation;
                localJob = job;
            }

            if (localJob == Job::Apply) runApply(tid);
            else runPrecond(tid);

            {
                std::lock_guard<std::mutex> lock(mu);
                ++done;
                if (done == nt) cvDone.notify_one();
            }
        }
    }
};

struct PcgStats {
    int iters = 0;
    double initialRelResidual = 0.0;
    double relResidual = 0.0;
    double trueRel = 0.0;
    double ms = 0.0;
    double initialApplyMs = 0.0;
    double initialPrecondMs = 0.0;
    double applyMs = 0.0;
    double precondMs = 0.0;
};

template <typename ApplyFn>
PcgStats pcg(const ElementOperator& A, const Preconditioner& P, const SpMat& Kff, const VecX& b, VecX& x,
             int maxIter, double tol, ApplyFn apply, bool parallelPrecond = false) {
    PcgStats st;
    Timer t;
    VecX Ax(b.size());
    {
        Timer tApply;
        apply(x, Ax);
        st.initialApplyMs += tApply.ms();
    }
    VecX r = b - Ax;
    const real bnorm = std::max<real>(real(1e-300), b.norm());
    st.relResidual = static_cast<double>(r.norm() / bnorm);
    st.initialRelResidual = st.relResidual;
    if (st.relResidual <= tol) {
        st.ms = t.ms();
        st.trueRel = relNorm(Kff * x, b);
        return st;
    }

    VecX z;
    PreconditionerScratch precondScratch;
    {
        Timer tPrecond;
        P.apply(r, z, precondScratch, parallelPrecond);
        st.initialPrecondMs += tPrecond.ms();
    }
    VecX p = z;
    real rzOld = r.dot(z);
    VecX Ap(b.size());

    for (int it = 0; it < maxIter; ++it) {
        {
            Timer tApply;
            apply(p, Ap);
            st.applyMs += tApply.ms();
        }
        const real denom = p.dot(Ap);
        if (!(denom > real(0)) || !std::isfinite(static_cast<double>(denom))) break;
        const real alpha = rzOld / denom;
        x += alpha * p;
        r -= alpha * Ap;
        st.iters = it + 1;
        st.relResidual = static_cast<double>(r.norm() / bnorm);
        if (st.relResidual <= tol) break;
        {
            Timer tPrecond;
            P.apply(r, z, precondScratch, parallelPrecond);
            st.precondMs += tPrecond.ms();
        }
        const real rzNew = r.dot(z);
        if (!(rzNew > real(0)) || !std::isfinite(static_cast<double>(rzNew))) break;
        p = z + (rzNew / rzOld) * p;
        rzOld = rzNew;
    }
    st.ms = t.ms();
    st.trueRel = relNorm(Kff * x, b);
    return st;
}

struct RhsFamily {
    std::vector<VecX> sequence;
    std::vector<VecX> loadModes;
    std::vector<char> isOutFrame;  // per-frame: contacts left the seed span (gameLoad only)
};

RhsFamily makeRhsFamily(const FrameModel& model, const PreparedSystem::Impl& S,
                        const VecX& base, int count) {
    RhsFamily out;
    out.sequence.reserve(static_cast<size_t>(count));

    real minX = std::numeric_limits<real>::max(), maxX = -std::numeric_limits<real>::max();
    real minY = std::numeric_limits<real>::max(), maxY = -std::numeric_limits<real>::max();
    real minZ = std::numeric_limits<real>::max(), maxZ = -std::numeric_limits<real>::max();
    for (const Node& n : model.nodes) {
        minX = std::min(minX, n.pos.x); maxX = std::max(maxX, n.pos.x);
        minY = std::min(minY, n.pos.y); maxY = std::max(maxY, n.pos.y);
        minZ = std::min(minZ, n.pos.z); maxZ = std::max(maxZ, n.pos.z);
    }
    const real spanX = std::max<real>(real(1), maxX - minX);
    const real spanY = std::max<real>(real(1), maxY - minY);
    const real spanZ = std::max<real>(real(1), maxZ - minZ);

    VecX windX = VecX::Zero(S.nf);
    VecX windY = VecX::Zero(S.nf);
    VecX liveZ = VecX::Zero(S.nf);
    VecX torsion = VecX::Zero(S.nf);
    for (size_t ni = 0; ni < model.nodes.size(); ++ni) {
        const Node& node = model.nodes[ni];
        const real xn = (node.pos.x - minX) / spanX - real(0.5);
        const real yn = (node.pos.y - minY) / spanY - real(0.5);
        const real zn = (node.pos.z - minZ) / spanZ;
        for (int d = 0; d < 6; ++d) {
            const int f = S.fmap[static_cast<size_t>(gdof(static_cast<int>(ni), d))];
            if (f < 0) continue;
            if (d == Ux) windX(f) = real(1600.0) * zn * (real(1) + real(0.15) * yn);
            if (d == Uy) windY(f) = real(1200.0) * zn * (real(1) - real(0.10) * xn);
            if (d == Uz) liveZ(f) = real(-4500.0) * (real(0.25) + zn) * (real(1) + real(0.10) * xn);
            if (d == Rz) torsion(f) = real(2.5e5) * zn * (xn - yn);
        }
    }

    out.loadModes.reserve(5);
    out.loadModes.push_back(base);
    out.loadModes.push_back(windX);
    out.loadModes.push_back(windY);
    out.loadModes.push_back(liveZ);
    out.loadModes.push_back(torsion);

    if (count <= 1) {
        out.sequence.push_back(base);
        return out;
    }
    for (int k = 0; k < count; ++k) {
        const double t = static_cast<double>(k);
        VecX b = (1.0 + 0.035 * std::sin(0.37 * t)) * base;
        b.noalias() += (0.65 + 0.35 * std::sin(0.23 * t + 0.4)) * windX;
        b.noalias() += (0.45 + 0.25 * std::cos(0.31 * t)) * windY;
        b.noalias() += (0.20 + 0.15 * std::sin(0.17 * t + 1.1)) * liveZ;
        b.noalias() += (0.25 * std::sin(0.41 * t)) * torsion;
        out.sequence.push_back(std::move(b));
    }
    return out;
}

// Real game-engine load family: constant gravity + a few MOVING sparse contact point
// loads. makeRhsFamily above is "honest cheating": every frame is a smooth sinusoidal
// combo of the SAME 5 modes it seeds, so every frame is in-subspace by construction.
// A real game re-solves a fixed structure each frame under gravity (constant) plus a
// handful of contact point loads whose LOCATIONS move. The seed basis here is
// { gravity } U { unit -Z response of each in-candidate contact node }. A frame is
// in-subspace iff all its contacts sit on seeded nodes (its solution is then a linear
// combo of seeded responses -> Galerkin-exact initial guess). `outFraction` of the
// frames deliberately place contacts on NON-seeded nodes, so the projection residual
// jumps -> exercises the gate. All randomness is a stateless Knuth hash (reproducible;
// no rand/time, so artifacts reproduce bit-for-bit).
RhsFamily makeGameLoadFamily(const FrameModel& model, const PreparedSystem::Impl& S,
                             const VecX& base, int count, int candidates, int active,
                             double outFraction, bool addHorizontal) {
    RhsFamily out;
    out.sequence.reserve(static_cast<size_t>(count));
    out.isOutFrame.reserve(static_cast<size_t>(count));

    auto hash32 = [](std::uint32_t x) -> std::uint32_t {
        x *= 2654435761u; x ^= x >> 16; x *= 0x45d9f3bu; x ^= x >> 16; return x;
    };
    auto frameHash = [&](int frame, int slot) -> std::uint32_t {
        return hash32(static_cast<std::uint32_t>(frame) * 2654435761u +
                      static_cast<std::uint32_t>(slot) * 40503u + 0x9e3779b9u);
    };

    // Nodes that own a FREE vertical DOF. The fixed base layer (k=0 fixAll) is excluded
    // automatically because its Uz maps to fmap < 0. node index == node id (makeTower
    // pushes nodes in towerNodeId order).
    std::vector<int> freeZNodes;
    freeZNodes.reserve(model.nodes.size());
    for (size_t ni = 0; ni < model.nodes.size(); ++ni)
        if (S.fmap[static_cast<size_t>(gdof(static_cast<int>(ni), Uz))] >= 0)
            freeZNodes.push_back(static_cast<int>(ni));

    // in-candidates (seeded) come from the mid-region of the free-Z nodes (a character's
    // activity floors, not the base); out-candidates are a disjoint free set used only by
    // out-of-subspace frames.
    const int nFree = static_cast<int>(freeZNodes.size());
    const int cand = std::max(1, std::min(candidates, std::max(1, nFree / 2)));
    const int inStart = std::min(std::max(0, nFree / 4), std::max(0, nFree - 2 * cand));
    std::vector<int> inNodes, outNodes;
    for (int t = 0; t < cand; ++t) inNodes.push_back(freeZNodes[static_cast<size_t>(inStart + t)]);
    for (int t = 0; t < cand; ++t) {
        const int idx = inStart + cand + t;
        outNodes.push_back(freeZNodes[static_cast<size_t>(idx < nFree ? idx : (idx % nFree))]);
    }
    // Out-candidates must stay disjoint from in-candidates, else "out-of-subspace" frames
    // would secretly land on seeded nodes and R3's residual would be falsely small. This
    // needs >= 2*cand distinct free-Z nodes; warn loudly rather than mislead silently.
    if (nFree < 2 * cand)
        std::fprintf(stderr, "[gameload] WARNING: only %d free-Z nodes for in+out = 2x%d; "
                     "out-candidates overlap in-candidates -> R3 out-of-subspace is INVALID\n",
                     nFree, cand);

    // Fixed, height-independent magnitude keeps every seed response of similar norm so
    // the K-orthonormal Gram-Schmidt stays well-conditioned.
    constexpr real kContactZ = real(5.0e4);
    constexpr real kContactH = real(5.0e3);

    auto addPoint = [&](VecX& b, int nodeId, real fz, real fx, real fy) {
        const int dz = S.fmap[static_cast<size_t>(gdof(nodeId, Uz))];
        if (dz >= 0) b(dz) += fz;
        if (addHorizontal) {
            const int dx = S.fmap[static_cast<size_t>(gdof(nodeId, Ux))];
            const int dy = S.fmap[static_cast<size_t>(gdof(nodeId, Uy))];
            if (dx >= 0) b(dx) += fx;
            if (dy >= 0) b(dy) += fy;
        }
    };

    // The seed basis spans gravity + the -Z unit response of each in-candidate. It does
    // NOT include horizontal directions, so with --gameAddHorizontal each frame's Ux/Uy
    // component lies OUTSIDE the seed span by design (a controlled perturbation): in that
    // mode even in-frames carry a small non-zero initialRel -- expected, not a bug. The
    // 0-iter result is the pure -Z case (the default, addHorizontal=false).
    out.loadModes.reserve(static_cast<size_t>(cand) + 1);
    out.loadModes.push_back(base);                       // gravity response is seeded
    for (int p = 0; p < cand; ++p) {                     // each in-candidate's unit -Z response
        VecX e = VecX::Zero(S.nf);
        addPoint(e, inNodes[static_cast<size_t>(p)], -kContactZ, real(0), real(0));
        out.loadModes.push_back(std::move(e));
    }

    const int na = std::max(1, active);
    for (int fr = 0; fr < count; ++fr) {
        const bool isOut = (frameHash(fr, 99) % 1000u) <
                           static_cast<std::uint32_t>(outFraction * 1000.0);
        const std::vector<int>& poolNodes = isOut ? outNodes : inNodes;
        VecX b = base;
        for (int s = 0; s < na; ++s) {
            const std::uint32_t h = frameHash(fr, s);
            const int node = poolNodes[static_cast<size_t>(h % poolNodes.size())];
            const real strength = real(0.5) +
                static_cast<real>((h >> 8) & 0xFFFFu) / real(65535.0);   // [0.5, 1.5]
            const real fx = addHorizontal
                ? kContactH * (real(-0.5) + static_cast<real>((h >> 12) & 0xFFu) / real(255.0)) : real(0);
            const real fy = addHorizontal
                ? kContactH * (real(-0.5) + static_cast<real>((h >> 20) & 0xFFu) / real(255.0)) : real(0);
            addPoint(b, node, -kContactZ * strength, fx, fy);
        }
        out.sequence.push_back(std::move(b));
        out.isOutFrame.push_back(isOut ? char(1) : char(0));
    }
    return out;
}

struct RecycleBasis {
    int maxSize = 0;
    std::vector<VecX> v;
    std::vector<VecX> av;
    double projectMs = 0.0;
    double addMs = 0.0;
    int accepted = 0;
    int rejected = 0;

    explicit RecycleBasis(int maxBasis) : maxSize(maxBasis) {}

    VecX initialGuess(int n, const VecX& b) {
        Timer t;
        VecX x = VecX::Zero(n);
        for (size_t i = 0; i < v.size(); ++i) x.noalias() += v[i].dot(b) * v[i];
        projectMs += t.ms();
        return x;
    }

    template <typename ApplyFn>
    void add(const VecX& x, ApplyFn apply) {
        if (maxSize <= 0) return;
        Timer t;
        VecX q = x;
        VecX Aq(x.size());
        apply(q, Aq);
        for (size_t i = 0; i < v.size(); ++i) {
            const real c = v[i].dot(Aq);
            q.noalias() -= c * v[i];
            Aq.noalias() -= c * av[i];
        }
        const real normA2 = q.dot(Aq);
        const double scaleRef = std::max(1e-300, static_cast<double>(x.norm()));
        if (!(normA2 > real(0)) || !std::isfinite(static_cast<double>(normA2)) ||
            static_cast<double>(q.norm()) <= 1e-12 * scaleRef) {
            ++rejected;
            addMs += t.ms();
            return;
        }
        const real inv = real(1) / std::sqrt(normA2);
        q *= inv;
        Aq *= inv;
        if (static_cast<int>(v.size()) >= maxSize) {
            v.erase(v.begin());
            av.erase(av.begin());
        }
        v.push_back(std::move(q));
        av.push_back(std::move(Aq));
        ++accepted;
        addMs += t.ms();
    }
};

// ---------------------------------------------------------------- spectral deflation
// RHS-independent acceleration: deflate the smallest eigenvectors of M^{-1}K (the
// slowest-converging modes). W (n x k) holds approximate Ritz vectors harvested via
// eigCG from the first solve; E = W^T K W and KW = K W are precomputed once and reused
// for every subsequent RHS (Type-1 Deflated CG, Saad/Yeung/Erhel/Guyomarc'h 2000).
struct DeflationSpace {
    MatX W;                    // n x k, Ritz vectors of M^{-1}K (smallest eigenvalues)
    MatX KW;                   // n x k = K * W (built via the matrix-free apply)
    Eigen::LDLT<MatX> Efac;    // factor of E = W^T K W (k x k)
    int k = 0;
    bool ready = false;
    double maxRitzResidual = 0.0;  // max_i ||K w_i - theta_i w_i|| / |theta_i|

    template <typename ApplyFn>
    void finalize(ApplyFn apply) {
        k = static_cast<int>(W.cols());
        if (k == 0) return;
        const int n = static_cast<int>(W.rows());
        KW.resize(n, k);
        VecX col(n), kcol(n);
        for (int j = 0; j < k; ++j) {
            col = W.col(j);
            apply(col, kcol);
            KW.col(j) = kcol;
        }
        MatX E = W.transpose() * KW;          // k x k
        E = real(0.5) * (E + E.transpose());  // symmetrize against roundoff
        Efac.compute(E);
        ready = (Efac.info() == Eigen::Success);
    }

    // Corrected initial guess x0 = W E^{-1} (W^T b); makes W^T (b - K x0) = 0.
    VecX initGuess(const VecX& b) const {
        return W * Efac.solve(W.transpose() * b);
    }

    // out = (I - W E^{-1} (KW)^T) z : the K-orthogonal-to-W component of z. Uses the
    // precomputed KW so no extra matrix-free apply is needed inside the PCG loop.
    void deflate(const VecX& z, VecX& out) const {
        out.noalias() = z - W * Efac.solve(KW.transpose() * z);
    }
};

// First-solve PCG that snapshots preconditioned residuals z_j = M^{-1} r_j into a
// subspace, then does a Rayleigh-Ritz on the K-inner-product projected pencil to
// extract the k smallest eigenvectors of M^{-1}K (the slowest CG modes) into Wout.
// This direct projected eigenproblem avoids the orthogonality loss of naive
// CG-Lanczos and gives reliable smallest-eigenvector approximations. Solves b too.
template <typename ApplyFn>
PcgStats harvestRR(const ElementOperator& A, const Preconditioner& P, const SpMat& Kff,
                   const VecX& b, VecX& x, int maxIter, double tol, ApplyFn apply,
                   bool parallelPrecond, int k, int window, MatX& Wout, double& ritzResOut) {
    PcgStats st;
    Timer t;
    const int n = static_cast<int>(b.size());
    VecX Ax(n);
    apply(x, Ax);
    VecX r = b - Ax;
    const real bnorm = std::max<real>(real(1e-300), b.norm());
    st.relResidual = static_cast<double>(r.norm() / bnorm);
    st.initialRelResidual = st.relResidual;

    VecX z;
    PreconditionerScratch ps;
    P.apply(r, z, ps, parallelPrecond);
    real rz = r.dot(z);

    // Sample preconditioned residuals across the whole convergence history (stride) so
    // the subspace captures the slow modes that linger in late residuals, not only the
    // early Krylov vectors. We do not know the iteration count ahead, so we keep the
    // most recent `window` snapshots taken every `stride` steps and grow stride adaptively.
    // Store the first `window` consecutive preconditioned residuals: a clean Krylov
    // basis K_window(M^{-1}K, z0) in which Rayleigh-Ritz resolves the extreme
    // eigenvectors. (Decimating across the whole history aliased the subspace.)
    std::vector<VecX> Z;
    Z.reserve(static_cast<size_t>(window) + 1);
    if (rz > real(0)) Z.push_back(z);

    VecX p = z;
    VecX Ap(n);
    for (int it = 0; it < maxIter; ++it) {
        apply(p, Ap);
        const real denom = p.dot(Ap);
        if (!(denom > real(0)) || !std::isfinite(static_cast<double>(denom))) break;
        const real alpha = rz / denom;
        x += alpha * p;
        r -= alpha * Ap;
        st.iters = it + 1;
        st.relResidual = static_cast<double>(r.norm() / bnorm);
        if (st.relResidual <= tol) break;
        P.apply(r, z, ps, parallelPrecond);
        const real rzNew = r.dot(z);
        if (!(rzNew > real(0)) || !std::isfinite(static_cast<double>(rzNew))) break;
        const real beta = rzNew / rz;
        if (static_cast<int>(Z.size()) < window) Z.push_back(z);
        p = z + beta * p;
        rz = rzNew;
    }
    st.ms = t.ms();
    st.trueRel = relNorm(Kff * x, b);
    ritzResOut = 0.0;

    const int m = static_cast<int>(Z.size());
    if (m <= 0 || k <= 0) { Wout.resize(n, 0); return st; }

    // Euclidean-orthonormalize the snapshot subspace -> Q (n x mq).
    MatX Zmat(n, m);
    for (int j = 0; j < m; ++j) Zmat.col(j) = Z[static_cast<size_t>(j)];
    Eigen::HouseholderQR<MatX> qr(Zmat);
    MatX Q = qr.householderQ() * MatX::Identity(n, m);

    // K-inner-product Rayleigh-Ritz for M^{-1}K:  H s = theta G s,
    // G = Q^T K Q,  H = (KQ)^T M^{-1} (KQ).  theta = eigenvalue of M^{-1}K.
    MatX KQ(n, m), MiKQ(n, m);
    VecX col(n), kcol(n), mcol(n);
    for (int j = 0; j < m; ++j) {
        col = Q.col(j);
        apply(col, kcol);
        KQ.col(j) = kcol;
        P.apply(kcol, mcol, ps, parallelPrecond);
        MiKQ.col(j) = mcol;
    }
    MatX G = Q.transpose() * KQ;
    MatX H = KQ.transpose() * MiKQ;
    G = real(0.5) * (G + G.transpose());
    H = real(0.5) * (H + H.transpose());
    Eigen::GeneralizedSelfAdjointEigenSolver<MatX> ges(H, G);
    if (ges.info() != Eigen::Success) { Wout.resize(n, 0); return st; }
    const int kk = std::min(k, m);
    const MatX S = ges.eigenvectors().leftCols(kk);   // ascending theta
    Wout.noalias() = Q * S;                            // n x kk Ritz vectors

    // Honest eigenvector residual ||M^{-1}K w - theta w|| / |theta| over the kept set.
    for (int j = 0; j < kk; ++j) {
        const real theta = ges.eigenvalues()(j);
        const VecX wv = Wout.col(j);
        const VecX mikw = MiKQ * S.col(j);            // M^{-1}K (Q s) = MiKQ s
        const double denom = std::abs(static_cast<double>(theta)) * static_cast<double>(wv.norm()) + 1e-300;
        ritzResOut = std::max(ritzResOut, static_cast<double>((mikw - theta * wv).norm()) / denom);
    }
    return st;
}

// Type-1 Deflated PCG for a subsequent RHS, reusing a frozen deflation space.
template <typename ApplyFn>
PcgStats deflatedPcg(const ElementOperator& A, const Preconditioner& P, const SpMat& Kff,
                     const VecX& b, VecX& x, int maxIter, double tol, ApplyFn apply,
                     bool parallelPrecond, const DeflationSpace& D) {
    PcgStats st;
    Timer t;
    const int n = static_cast<int>(b.size());
    x = D.initGuess(b);
    VecX Ax(n);
    apply(x, Ax);
    VecX r = b - Ax;
    const real bnorm = std::max<real>(real(1e-300), b.norm());
    st.relResidual = static_cast<double>(r.norm() / bnorm);
    st.initialRelResidual = st.relResidual;
    if (st.relResidual <= tol) {
        st.ms = t.ms();
        st.trueRel = relNorm(Kff * x, b);
        return st;
    }
    VecX z;
    PreconditionerScratch ps;
    {
        Timer tp;
        P.apply(r, z, ps, parallelPrecond);
        st.initialPrecondMs += tp.ms();
    }
    VecX p(n);
    D.deflate(z, p);  // p0 = (I - QK) z0, K-orthogonal to W
    real rz = r.dot(z);
    VecX Ap(n), wproj(n);
    for (int it = 0; it < maxIter; ++it) {
        {
            Timer ta;
            apply(p, Ap);
            st.applyMs += ta.ms();
        }
        const real denom = p.dot(Ap);
        if (!(denom > real(0)) || !std::isfinite(static_cast<double>(denom))) break;
        const real alpha = rz / denom;
        x += alpha * p;
        r -= alpha * Ap;
        st.iters = it + 1;
        st.relResidual = static_cast<double>(r.norm() / bnorm);
        if (st.relResidual <= tol) break;
        {
            Timer tp;
            P.apply(r, z, ps, parallelPrecond);
            st.precondMs += tp.ms();
        }
        const real rzNew = r.dot(z);
        if (!(rzNew > real(0)) || !std::isfinite(static_cast<double>(rzNew))) break;
        D.deflate(z, wproj);  // (I - QK) z
        const real beta = rzNew / rz;
        p = wproj + beta * p;
        rz = rzNew;
    }
    st.ms = t.ms();
    st.trueRel = relNorm(Kff * x, b);
    return st;
}

} // namespace

int main(int argc, char** argv) {
    const Args args = parseArgs(argc, argv);
    FrameModel model = makeCase(args);
    assertNodalOnly(model, "exp_parallel_pcg");

    Timer tHpSetup;
    const FreeDofMap freeMap = makeFreeDofMap(model);
    if (freeMap.singular) {
        std::printf("[parallel_pcg] status=singular diagnostic=\"%s\"\n", freeMap.diagnostic.c_str());
        return 2;
    }
    std::string why;
    ElementOperator A;
    Timer tBuild;
    if (!A.build(model, freeMap.fmap, freeMap.nf, args.symApply, why)) {
        std::printf("[parallel_pcg] status=build_failed why=\"%s\"\n", why.c_str());
        return 2;
    }
    const double opBuildMs = tBuild.ms();
    Preconditioner P;
    Timer tPrecond;
    P.build(A, model, freeMap.fmap, args);
    const double precondBuildMs = tPrecond.ms();
    ThreadApplyPool pool(A, args.threads);
    P.parallelBlock6 = [&pool](const VecX& r, VecX& z, const std::vector<std::array<real, 36>>& inv6) {
        pool.precondBlock6(r, z, inv6);
    };
    const bool usePoolPrecond = args.parallelPrecond;
    const double hpSetupMs = tHpSetup.ms();

    Timer tPrep;
    PreparedSystem ps = assembleAndFactor(model);
    const double prepMs = tPrep.ms();
    const PreparedSystem::Impl& S = *ps.impl;
    if (S.singular) {
        std::printf("[parallel_pcg] status=oracle_singular diagnostic=\"%s\"\n", S.diagnostic.c_str());
        return 2;
    }
    if (S.nf != freeMap.nf || S.fmap != freeMap.fmap) {
        std::printf("[parallel_pcg] status=map_mismatch hpNf=%d oracleNf=%d\n", freeMap.nf, S.nf);
        return 2;
    }

    const SpMat Kff = research::reduceFF(S.K, S.fmap, S.nf);
    const VecX Fff = reduceVec(nodalLoadVector(model), freeMap.fmap, freeMap.nf);

    VecX probe = VecX::Zero(S.nf);
    // Mixed low/high frequency plus a sign-alternating term so a single probe still
    // excites per-DOF coupling and sign errors, not just smooth global modes.
    for (int i = 0; i < S.nf; ++i)
        probe(i) = std::sin(0.0017 * static_cast<double>(i + 1)) +
                   real(0.5) * std::cos(0.131 * static_cast<double>(i + 7)) +
                   ((i & 1) ? real(-0.3) : real(0.3));
    VecX ySer(S.nf), yPar(S.nf);
    A.applySerial(probe, ySer);
    pool.apply(probe, yPar);
    const double applyRel = std::max(relNorm(ySer, Kff * probe), relNorm(yPar, ySer));

    if (args.deflation > 0) {
        const RhsFamily rhsFamily = makeRhsFamily(model, S, Fff, args.rhsCount);
        const std::vector<VecX>& rhs = rhsFamily.sequence;
        const int nr = static_cast<int>(rhs.size());
        const int window = args.deflationWindow > 0 ? args.deflationWindow
                                                    : std::max(64, 6 * args.deflation);
        auto poolApply = [&](const VecX& x, VecX& y) { pool.apply(x, y); };

        // Harvest the deflation space from the first RHS (this solve is the RHS0 cost).
        MatX W;
        double ritzRes = 0.0;
        VecX xHarvest = VecX::Zero(S.nf);
        Timer tHarvest;
        const PcgStats harv = harvestRR(A, P, Kff, rhs[0], xHarvest, args.pcgMaxIter, args.pcgTol,
                                        poolApply, usePoolPrecond, args.deflation, window, W, ritzRes);
        const double harvestMs = tHarvest.ms();
        DeflationSpace D;
        D.W = W;
        Timer tDefSetup;
        D.finalize(poolApply);
        D.maxRitzResidual = ritzRes;
        const double defSetupMs = tDefSetup.ms();
        const double maxHarvestVsLdlt =
            static_cast<double>(relNorm(xHarvest, S.ldlt.solve(rhs[0])));

        int plainIters = 0, deflIters = 0, reuseN = 0;
        double plainMs = 0.0, deflMs = 0.0, ldltSolveMs = 0.0;
        double deflApplyMs = 0.0, deflPrecondMs = 0.0, plainApplyMs = 0.0, plainPrecondMs = 0.0;
        double maxPlainTrueRel = 0.0, maxDeflTrueRel = 0.0, maxDeflVsLdlt = 0.0, maxDeflInitRel = 0.0;
        for (int i = 1; i < nr; ++i) {  // RHS0 was the harvest; reuse-solve the rest
            const VecX& b = rhs[static_cast<size_t>(i)];
            Timer tL;
            const VecX xLdlt = S.ldlt.solve(b);
            ldltSolveMs += tL.ms();

            VecX xPlain = VecX::Zero(S.nf);
            const PcgStats plain = pcg(A, P, Kff, b, xPlain, args.pcgMaxIter, args.pcgTol,
                                       poolApply, usePoolPrecond);
            plainIters += plain.iters;
            plainMs += plain.ms;
            plainApplyMs += plain.initialApplyMs + plain.applyMs;
            plainPrecondMs += plain.initialPrecondMs + plain.precondMs;
            maxPlainTrueRel = std::max(maxPlainTrueRel, plain.trueRel);

            VecX xDefl(S.nf);
            const PcgStats defl = deflatedPcg(A, P, Kff, b, xDefl, args.pcgMaxIter, args.pcgTol,
                                              poolApply, usePoolPrecond, D);
            deflIters += defl.iters;
            deflMs += defl.ms;
            deflApplyMs += defl.initialApplyMs + defl.applyMs;
            deflPrecondMs += defl.initialPrecondMs + defl.precondMs;
            maxDeflInitRel = std::max(maxDeflInitRel, defl.initialRelResidual);
            maxDeflTrueRel = std::max(maxDeflTrueRel, defl.trueRel);
            maxDeflVsLdlt = std::max(maxDeflVsLdlt, relNorm(xDefl, xLdlt));
            ++reuseN;
        }
        const double rn = static_cast<double>(std::max(1, reuseN));
        const double plainItersAvg = plainIters / rn;
        const double deflItersAvg = deflIters / rn;
        const double plainMsAvg = plainMs / rn;
        const double deflMsAvg = deflMs / rn;
        const double deflApplyMsAvg = deflApplyMs / rn;
        const double deflPrecondMsAvg = deflPrecondMs / rn;
        const double deflOtherMsAvg = deflMsAvg - deflApplyMsAvg - deflPrecondMsAvg;
        const double ldltSolveMsAvg = ldltSolveMs / rn;
        // Honest batch: LDLT factors once then solves all nr; HP harvests RHS0 then
        // deflate-solves RHS1..nr-1. Both amortize the same nr loads.
        const double directBatchMs = prepMs + ldltSolveMs + static_cast<double>(0);  // rhs0 ldlt below
        Timer tL0;
        const VecX xLdlt0 = S.ldlt.solve(rhs[0]);
        const double directBatchMsFull = directBatchMs + tL0.ms();
        const double hpBatchMs = hpSetupMs + harvestMs + defSetupMs + deflMs;
        const double factorBypassBatchSpeedup = directBatchMsFull / std::max(1e-300, hpBatchMs);
        const double speedupIters = plainItersAvg / std::max(1e-300, deflItersAvg);
        const double speedupMs = plainMsAvg / std::max(1e-300, deflMsAvg);
        const bool ok = applyRel <= 1e-10 && D.ready && maxDeflTrueRel <= 1e-8 &&
                        maxDeflVsLdlt <= 1e-6 && maxHarvestVsLdlt <= 1e-6;

        std::printf("[deflation_hpfem] preset=%s parallelPrecond=%d symApply=%d coarseSolve=%s bandedOk=%d sparseOk=%d coarseNnz=%.0f threads=%d nf=%d blocks=%zu deflK=%d kept=%d window=%d harvestIters=%d harvestMs=%.3f defSetupMs=%.3f maxRitzResidual=%.3e reuseN=%d prepMs=%.3f hpSetupMs=%.3f ldltSolveMsAvg=%.3f directBatchMs=%.3f hpBatchMs=%.3f factorBypassBatchSpeedup=%.3f applyRel=%.3e plainItersAvg=%.2f deflItersAvg=%.2f speedupIters=%.3f plainMsAvg=%.3f deflMsAvg=%.3f deflApplyMsAvg=%.3f deflPrecondMsAvg=%.3f deflOtherMsAvg=%.3f speedupMs=%.3f maxDeflInitRel=%.3e maxHarvestVsLdlt=%.3e maxPlainTrueRel=%.3e maxDeflTrueRel=%.3e maxDeflVsLdlt=%.3e\n",
                    args.preset.c_str(), args.parallelPrecond ? 1 : 0, args.symApply ? 1 : 0, args.coarseSolve.c_str(),
                    P.bandedOk ? 1 : 0, P.sparseOk ? 1 : 0, P.coarseNnz, pool.threads(), S.nf, A.numBlocks(), args.deflation,
                    D.k, window, harv.iters, harvestMs, defSetupMs, D.maxRitzResidual, reuseN,
                    prepMs, hpSetupMs, ldltSolveMsAvg, directBatchMsFull, hpBatchMs,
                    factorBypassBatchSpeedup, applyRel, plainItersAvg, deflItersAvg, speedupIters,
                    plainMsAvg, deflMsAvg, deflApplyMsAvg, deflPrecondMsAvg, deflOtherMsAvg, speedupMs,
                    maxDeflInitRel, maxHarvestVsLdlt, maxPlainTrueRel, maxDeflTrueRel, maxDeflVsLdlt);
        return ok ? 0 : 1;
    }

    if (args.rhsCount > 1 || args.basisMax > 0) {
        const RhsFamily rhsFamily = args.gameLoad
            ? makeGameLoadFamily(model, S, Fff, args.rhsCount, args.gameContactCandidates,
                                 args.gameContactActive, args.gameOutFraction, args.gameAddHorizontal)
            : makeRhsFamily(model, S, Fff, args.rhsCount);
        const std::vector<VecX>& rhs = rhsFamily.sequence;
        // Seeded game loads must keep EVERY seed response (gravity + each contact unit
        // response). RecycleBasis evicts FIFO, so a too-small basisMax would silently drop
        // early seeds and break the in-subspace guarantee -> auto-grow it (and report it).
        int effBasisMax = args.basisMax;
        if (args.gameLoad && args.seedLoadBasis &&
            effBasisMax < static_cast<int>(rhsFamily.loadModes.size())) {
            std::fprintf(stderr, "[gameload] note: raising basisMax %d -> %zu to fit all seed modes\n",
                         effBasisMax, rhsFamily.loadModes.size());
            effBasisMax = static_cast<int>(rhsFamily.loadModes.size());
        }
        RecycleBasis basis(effBasisMax);
        int seedIters = 0;
        int seedCount = 0;
        double seedMs = 0.0;
        double seedApplyMs = 0.0;
        double seedPrecondMs = 0.0;
        double maxSeedTrueRel = 0.0;
        double verifyOffDiag = -1.0, verifyCrossRel = -1.0;  // gameVerify dense cross-check
        if (args.seedLoadBasis && effBasisMax > 0) {
            for (const VecX& seedRhs : rhsFamily.loadModes) {
                if (seedRhs.norm() <= real(0)) continue;
                // Recycle the accumulating A-orthogonal basis as the seed's initial guess so
                // later load-mode seeds start closer (first seed gets a zero guess). PCG still
                // converges to tol from any start, so this only trims seed iterations.
                VecX xSeed = basis.initialGuess(S.nf, seedRhs);
                const PcgStats seed = pcg(A, P, Kff, seedRhs, xSeed, args.pcgMaxIter, args.pcgTol,
                                          [&](const VecX& x, VecX& y) { pool.apply(x, y); }, usePoolPrecond);
                seedIters += seed.iters;
                seedMs += seed.ms;
                seedApplyMs += seed.initialApplyMs + seed.applyMs;
                seedPrecondMs += seed.initialPrecondMs + seed.precondMs;
                maxSeedTrueRel = std::max(maxSeedTrueRel, seed.trueRel);
                ++seedCount;
                basis.add(xSeed, [&](const VecX& x, VecX& y) { pool.apply(x, y); });
            }
        }

        // Dense second-path cross-check of the cheap Euclidean initialGuess: build V from
        // the seeded basis, form VtKV via the matrix-free apply, and compare the
        // Galerkin-optimal V (VtKV)^-1 Vt b against initialGuess(b). verifyOffDiag is how
        // far VtKV is from I (the precondition that makes the two equal). Small cases only.
        if (args.gameLoad && args.gameVerify && !basis.v.empty() && !rhs.empty()) {
            const int kdim = static_cast<int>(basis.v.size());
            MatX V(S.nf, kdim);
            for (int j = 0; j < kdim; ++j) V.col(j) = basis.v[static_cast<size_t>(j)];
            MatX KV(S.nf, kdim);
            for (int j = 0; j < kdim; ++j) {
                VecX vj = V.col(j), kv(S.nf);
                pool.apply(vj, kv);
                KV.col(j) = kv;
            }
            const MatX VtKV = V.transpose() * KV;
            verifyOffDiag = (VtKV - MatX::Identity(kdim, kdim)).norm();
            const VecX& bP = rhs[0];
            Eigen::LDLT<MatX> Ef(VtKV);
            const VecX x0dense = V * Ef.solve(V.transpose() * bP);
            const VecX x0fast = basis.initialGuess(S.nf, bP);
            verifyCrossRel = (x0dense - x0fast).norm() /
                             std::max(1e-300, static_cast<double>(x0dense.norm()));
        }

        int serialIters = 0;
        int combinedIters = 0;
        int combinedItersSkip1 = 0;
        double serialMs = 0.0;
        double combinedMs = 0.0;
        double serialApplyMs = 0.0;
        double serialPrecondMs = 0.0;
        double combinedApplyMs = 0.0;
        double combinedPrecondMs = 0.0;
        double ldltSolveMs = 0.0;
        double maxSerialTrueRel = 0.0;
        double maxCombinedTrueRel = 0.0;
        double maxCombinedInitialRel = 0.0;
        double maxCombinedVsLdlt = 0.0;
        int skipN = 0;
        int inFrames = 0, outFrames = 0;
        long long inIters = 0, outIters = 0;
        double inMaxInitialRel = 0.0, outMaxInitialRel = 0.0;

        for (int i = 0; i < static_cast<int>(rhs.size()); ++i) {
            const VecX& b = rhs[static_cast<size_t>(i)];

            Timer tLdltSolve;
            const VecX xLdlt = S.ldlt.solve(b);
            ldltSolveMs += tLdltSolve.ms();

            if (!args.noSerialBaseline) {
                VecX xSerial = VecX::Zero(S.nf);
                const PcgStats serial = pcg(A, P, Kff, b, xSerial, args.pcgMaxIter, args.pcgTol,
                                            [&](const VecX& x, VecX& y) { A.applySerial(x, y); });
                serialIters += serial.iters;
                serialMs += serial.ms;
                serialApplyMs += serial.initialApplyMs + serial.applyMs;
                serialPrecondMs += serial.initialPrecondMs + serial.precondMs;
                maxSerialTrueRel = std::max(maxSerialTrueRel, serial.trueRel);
            }

            VecX xCombined = basis.initialGuess(S.nf, b);
            const PcgStats combined = pcg(A, P, Kff, b, xCombined, args.pcgMaxIter, args.pcgTol,
                                          [&](const VecX& x, VecX& y) { pool.apply(x, y); }, usePoolPrecond);
            combinedIters += combined.iters;
            combinedMs += combined.ms;
            combinedApplyMs += combined.initialApplyMs + combined.applyMs;
            combinedPrecondMs += combined.initialPrecondMs + combined.precondMs;
            maxCombinedInitialRel = std::max(maxCombinedInitialRel, combined.initialRelResidual);
            maxCombinedTrueRel = std::max(maxCombinedTrueRel, combined.trueRel);
            maxCombinedVsLdlt = std::max(maxCombinedVsLdlt, relNorm(xCombined, xLdlt));
            if (args.gameLoad && i < static_cast<int>(rhsFamily.isOutFrame.size())) {
                if (rhsFamily.isOutFrame[static_cast<size_t>(i)]) {
                    ++outFrames; outIters += combined.iters;
                    outMaxInitialRel = std::max(outMaxInitialRel, combined.initialRelResidual);
                } else {
                    ++inFrames; inIters += combined.iters;
                    inMaxInitialRel = std::max(inMaxInitialRel, combined.initialRelResidual);
                }
            }
            if (i > 0) {
                combinedItersSkip1 += combined.iters;
                ++skipN;
            }
            if (!args.seedLoadBasis)
                basis.add(xCombined, [&](const VecX& x, VecX& y) { pool.apply(x, y); });
        }

        const double n = static_cast<double>(rhs.size());
        const double serialItersAvg = static_cast<double>(serialIters) / n;
        const double combinedItersAvg = static_cast<double>(combinedIters) / n;
        const double combinedItersSkip1Avg = static_cast<double>(combinedItersSkip1) /
                                             static_cast<double>(std::max(1, skipN));
        const double serialMsAvg = serialMs / n;
        const double combinedMsAvg = combinedMs / n;
        const double serialApplyMsAvg = serialApplyMs / n;
        const double serialPrecondMsAvg = serialPrecondMs / n;
        const double serialOtherMsAvg = serialMsAvg - serialApplyMsAvg - serialPrecondMsAvg;
        const double combinedApplyMsAvg = combinedApplyMs / n;
        const double combinedPrecondMsAvg = combinedPrecondMs / n;
        const double combinedOtherMsAvg = combinedMsAvg - combinedApplyMsAvg - combinedPrecondMsAvg;
        const double ldltSolveMsAvg = ldltSolveMs / n;
        const double directBatchMs = prepMs + ldltSolveMs;
        const double hpBatchMs = hpSetupMs + seedMs + combinedMs;
        const double factorBypassSetupSpeedup = prepMs / std::max(1e-300, hpSetupMs);
        const double factorBypassFirstSolveSpeedup =
            (prepMs + ldltSolveMsAvg) / std::max(1e-300, hpSetupMs + seedMs + combinedMsAvg);
        const double factorBypassBatchSpeedup = directBatchMs / std::max(1e-300, hpBatchMs);
        const bool ok = applyRel <= 1e-10 &&
                        maxSeedTrueRel <= 1e-8 &&
                        maxCombinedTrueRel <= 1e-8 &&
                        maxCombinedVsLdlt <= 1e-6;

        std::printf("[combined_hpfem] preset=%s precond=%s parallelPrecond=%d symApply=%d coarseSolve=%s bandedOk=%d sparseOk=%d coarseNnz=%.0f bandedResidual=%.3e bandedOffBandRel=%.3e coarseBins=%dx%d coarseDofs=%lld threads=%d rhs=%d basisMax=%d seedLoadBasis=%d serialBaseline=%d seedCount=%d seedIters=%d seedMs=%.3f seedApplyMs=%.3f seedPrecondMs=%.3f basisAccepted=%d basisRejected=%d nf=%d blocks=%zu prepMs=%.3f hpSetupMs=%.3f opBuildMs=%.3f precondBuildMs=%.3f ldltSolveMsAvg=%.3f factorBypassSetupSpeedup=%.3f factorBypassFirstSolveSpeedup=%.3f factorBypassBatchSpeedup=%.3f directBatchMs=%.3f hpBatchMs=%.3f applyRel=%.3e serialItersAvg=%.2f combinedItersAvg=%.2f combinedItersSkip1Avg=%.2f speedupIters=%.3f serialMsAvg=%.3f serialApplyMsAvg=%.3f serialPrecondMsAvg=%.3f serialOtherMsAvg=%.3f combinedMsAvg=%.3f combinedApplyMsAvg=%.3f combinedPrecondMsAvg=%.3f combinedOtherMsAvg=%.3f speedupMs=%.3f projectMsAvg=%.6f basisAddMs=%.3f maxCombinedInitialRel=%.3e maxSeedTrueRel=%.3e maxSerialTrueRel=%.3e maxCombinedTrueRel=%.3e maxCombinedVsLdlt=%.3e\n",
                    args.preset.c_str(), args.precond.c_str(), args.parallelPrecond ? 1 : 0,
                    args.symApply ? 1 : 0, args.coarseSolve.c_str(), P.bandedOk ? 1 : 0,
                    P.sparseOk ? 1 : 0, P.coarseNnz, P.bandedResidual, P.bandedOffBandRel,
                    args.coarseBinsX, args.coarseBinsY,
                    static_cast<long long>(P.coarseDim), pool.threads(), static_cast<int>(rhs.size()),
                    effBasisMax, args.seedLoadBasis ? 1 : 0, args.noSerialBaseline ? 0 : 1,
                    seedCount, seedIters, seedMs, seedApplyMs,
                    seedPrecondMs, basis.accepted, basis.rejected, S.nf, A.numBlocks(), prepMs, hpSetupMs,
                    opBuildMs, precondBuildMs, ldltSolveMsAvg,
                    factorBypassSetupSpeedup, factorBypassFirstSolveSpeedup, factorBypassBatchSpeedup,
                    directBatchMs, hpBatchMs, applyRel, serialItersAvg, combinedItersAvg, combinedItersSkip1Avg,
                    combinedItersAvg > 0.0 ? serialItersAvg / combinedItersAvg : -1.0,
                    serialMsAvg, serialApplyMsAvg, serialPrecondMsAvg, serialOtherMsAvg,
                    combinedMsAvg, combinedApplyMsAvg, combinedPrecondMsAvg, combinedOtherMsAvg,
                    serialMsAvg / std::max(1e-300, combinedMsAvg),
                    basis.projectMs / n, basis.addMs, maxCombinedInitialRel, maxSeedTrueRel, maxSerialTrueRel,
                    maxCombinedTrueRel, maxCombinedVsLdlt);
        if (args.gameLoad) {
            // HONEST per-frame cost INCLUDES the Galerkin projection (basis.projectMs),
            // which pcg's combinedMs does NOT time. For the seeded game lane the k=
            // seedModes projection is a real per-frame cost (it grows with the contact-
            // candidate count), so the meaningful speedup vs reused LDLT is computed
            // against combinedMsAvg + projectMsAvg, not combinedMsAvg alone.
            const double perFrameMsWithProj = combinedMsAvg + basis.projectMs / n;
            std::printf("[gameload_meta] gameLoad=1 candidates=%d active=%d outFraction=%.3f "
                        "addHorizontal=%d seedModes=%zu inFrames=%d outFrames=%d "
                        "inMaxInitialRel=%.3e outMaxInitialRel=%.3e inAvgIters=%.2f outAvgIters=%.2f "
                        "basisAccepted=%d basisRejected=%d effBasisMax=%d "
                        "verifyOffDiagVtKV=%.3e verifyCrossRel=%.3e "
                        "perFrameMsWithProj=%.4f perFrameSpeedupWithProj=%.2f\n",
                        args.gameContactCandidates, args.gameContactActive, args.gameOutFraction,
                        args.gameAddHorizontal ? 1 : 0, rhsFamily.loadModes.size(),
                        inFrames, outFrames, inMaxInitialRel, outMaxInitialRel,
                        inFrames ? static_cast<double>(inIters) / inFrames : 0.0,
                        outFrames ? static_cast<double>(outIters) / outFrames : 0.0,
                        basis.accepted, basis.rejected, effBasisMax,
                        verifyOffDiag, verifyCrossRel,
                        perFrameMsWithProj, ldltSolveMsAvg / std::max(1e-300, perFrameMsWithProj));
        }
        return ok ? 0 : 1;
    }

    VecX xSerial = VecX::Zero(S.nf);
    const PcgStats serial = pcg(A, P, Kff, Fff, xSerial, args.pcgMaxIter, args.pcgTol,
                                [&](const VecX& x, VecX& y) { A.applySerial(x, y); });

    VecX xParallel = VecX::Zero(S.nf);
    const PcgStats parallel = pcg(A, P, Kff, Fff, xParallel, args.pcgMaxIter, args.pcgTol,
                                  [&](const VecX& x, VecX& y) { pool.apply(x, y); }, usePoolPrecond);

    Timer tLdltSolve;
    const VecX xLdlt = S.ldlt.solve(Fff);
    const double ldltSolveMs = tLdltSolve.ms();
    const double pcgVsLdlt = relNorm(xParallel, xLdlt);
    const double factorBypassSetupSpeedup = prepMs / std::max(1e-300, hpSetupMs);
    const double directFirstMs = prepMs + ldltSolveMs;
    const double hpFirstMs = hpSetupMs + parallel.ms;
    const double factorBypassFirstSolveSpeedup = directFirstMs / std::max(1e-300, hpFirstMs);
    const double serialApplyMs = serial.initialApplyMs + serial.applyMs;
    const double serialPrecondMs = serial.initialPrecondMs + serial.precondMs;
    const double parallelApplyMs = parallel.initialApplyMs + parallel.applyMs;
    const double parallelPrecondMs = parallel.initialPrecondMs + parallel.precondMs;
    const double serialOtherMs = serial.ms - serialApplyMs - serialPrecondMs;
    const double parallelOtherMs = parallel.ms - parallelApplyMs - parallelPrecondMs;
    const bool ok = applyRel <= 1e-10 && parallel.trueRel <= 1e-8 && pcgVsLdlt <= 1e-6;

    std::printf("[parallel_pcg] preset=%s precond=%s parallelPrecond=%d symApply=%d coarseSolve=%s bandedOk=%d sparseOk=%d coarseNnz=%.0f bandedResidual=%.3e bandedOffBandRel=%.3e coarseBins=%dx%d coarseDofs=%lld threads=%d nf=%d blocks=%zu prepMs=%.3f hpSetupMs=%.3f opBuildMs=%.3f precondBuildMs=%.3f ldltSolveMs=%.3f factorBypassSetupSpeedup=%.3f factorBypassFirstSolveSpeedup=%.3f directFirstMs=%.3f hpFirstMs=%.3f applyRel=%.3e serialIters=%d parallelIters=%d serialMs=%.3f serialApplyMs=%.3f serialPrecondMs=%.3f serialOtherMs=%.3f parallelMs=%.3f parallelApplyMs=%.3f parallelPrecondMs=%.3f parallelOtherMs=%.3f speedupMs=%.3f serialTrueRel=%.3e parallelTrueRel=%.3e pcgVsLdlt=%.3e\n",
                args.preset.c_str(), args.precond.c_str(), args.parallelPrecond ? 1 : 0,
                args.symApply ? 1 : 0, args.coarseSolve.c_str(), P.bandedOk ? 1 : 0, P.sparseOk ? 1 : 0, P.coarseNnz, P.bandedResidual, P.bandedOffBandRel,
                args.coarseBinsX, args.coarseBinsY,
                static_cast<long long>(P.coarseDim), pool.threads(), S.nf, A.numBlocks(), prepMs,
                hpSetupMs, opBuildMs, precondBuildMs, ldltSolveMs,
                factorBypassSetupSpeedup, factorBypassFirstSolveSpeedup, directFirstMs, hpFirstMs,
                applyRel, serial.iters, parallel.iters, serial.ms, serialApplyMs, serialPrecondMs, serialOtherMs,
                parallel.ms, parallelApplyMs, parallelPrecondMs, parallelOtherMs,
                serial.ms / std::max(1e-300, parallel.ms), serial.trueRel, parallel.trueRel, pcgVsLdlt);
    return ok ? 0 : 1;
}
