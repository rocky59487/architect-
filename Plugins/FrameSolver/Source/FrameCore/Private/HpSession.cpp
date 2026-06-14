//
// Seeded HP-FEM solve session — see Public/FrameCore/HpSession.h.
//
// A stateful accelerator for a FIXED factorized structure solved repeatedly under a low-dimensional
// load family (the game-engine niche: gravity + a few contacts). solveFrame reuses the SAME RHS
// assembly, prescribed reduction, scatter, reactions (K*u - F) and element recovery as solveLoad —
// copied verbatim below and marked "sync if changed" — and replaces ONLY the linear solve with:
//
//   1. a seeded, A-orthonormal (K-inner-product) basis of load responses, whose Galerkin initial
//      guess is exact when the load lies in the seeded subspace (=> ~0 PCG iterations);
//   2. a matrix-free preconditioned conjugate gradient that refines an out-of-subspace guess —
//      per-element dense apply, a block6 (or scalar) Jacobi, and an optional floor-aggregated
//      coarse-grid correction; the apply and the block6 map parallelize across a persistent worker
//      pool when opts.threads > 1;
//   3. the PreparedSystem's LDLT factorization as the always-correct oracle and fallback.
//
// The matrix-free operator is built from the prepared elements (S.elems), so it equals K_ff and the
// accelerated result matches the direct solve to tolerance (the serial apply is bit-consistent with
// K_ff; the threaded reduction is only tolerance-equal, not bit-identical). Eigen and the
// threading are confined to this Private .cpp; the public header is POD. The seeded basis and the
// LDLT stay single-threaded (Eigen is not thread-safe) — only the element apply + block6 map run on
// the pool. The matrix-free apply, worker pool, block6 and coarse-grid components are adapted from
// the validated solver-research prototype in Research/WS_B_solver; the net-neutral / negative
// experiments there (symmetric apply, finer / sparse coarse, spectral deflation) are deliberately
// not carried over.
//
#include "FrameCore/HpSession.h"
#include "PreparedSystemImpl.h"   // PreparedSystem::Impl (K/fmap/nf/ldlt/elems/...) + FrameEigen.h
#include "IElement.h"             // el->assemble / addEquivalentNodalLoads / recover / localDof
#include "ModelHash.h"            // modelFingerprint (the shared reuse-validity guard)

#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <array>
#include <utility>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>

namespace frame {

namespace {

// ---- A-orthonormal seeded basis (K-inner-product Gram-Schmidt) ----------------------------------
// Invariant: v_i^T K v_j = delta_ij (cross-checked ||VtKV - I|| ~ 3.4e-15), so initialGuess(b) =
// sum_i (v_i . b) v_i is the EXACT Galerkin solution whenever K^{-1} b lies in span{v_i} (=> ~0 PCG
// iterations in-subspace).
struct RecycleBasis {
    int maxSize = 0;
    std::vector<VecX> v;    // K-orthonormal responses
    std::vector<VecX> av;   // av[i] = K * v[i]

    RecycleBasis() = default;
    explicit RecycleBasis(int maxBasis) : maxSize(maxBasis) {}

    // Galerkin projection x0 = sum_i (v_i . b) v_i, together with its image Kx0 = K x0 =
    // sum_i (v_i . b) av_i — assembled for FREE from the stored av_i = K v_i, so the caller gets
    // the initial residual (b - Kx0) WITHOUT a matrix-free apply. The per-frame apply is the
    // bottleneck, so this makes an in-subspace (converged) frame cost ZERO applies.
    void initialGuess(int n, const VecX& b, VecX& x0, VecX& Kx0) const {
        x0  = VecX::Zero(n);
        Kx0 = VecX::Zero(n);
        for (size_t i = 0; i < v.size(); ++i) {
            const real c = v[i].dot(b);
            x0.noalias()  += c * v[i];
            Kx0.noalias() += c * av[i];
        }
    }

    // Append a response x (= K^{-1} b for some seed load b), A-orthonormalizing it against the
    // existing basis. `apply(q, Aq)` must set Aq = K * q. A linearly dependent / null seed is
    // rejected. Capacity is sized up-front to hold every seed, so FIFO eviction never fires.
    template <typename ApplyFn>
    void add(const VecX& x, ApplyFn apply) {
        if (maxSize <= 0) return;
        VecX q = x;
        VecX Aq(x.size());
        apply(q, Aq);
        for (size_t i = 0; i < v.size(); ++i) {
            const real c = v[i].dot(Aq);
            q.noalias()  -= c * v[i];
            Aq.noalias() -= c * av[i];
        }
        const real normA2 = q.dot(Aq);
        const double scaleRef = std::max(1e-300, static_cast<double>(x.norm()));
        if (!(normA2 > real(0)) || !std::isfinite(static_cast<double>(normA2)) ||
            static_cast<double>(q.norm()) <= 1e-12 * scaleRef) {
            return;   // linearly dependent / null seed -> reject
        }
        const real inv = real(1) / std::sqrt(normA2);
        q  *= inv;
        Aq *= inv;
        if (static_cast<int>(v.size()) >= maxSize) {
            v.erase(v.begin());
            av.erase(av.begin());
        }
        v.push_back(std::move(q));
        av.push_back(std::move(Aq));
    }
};

// ---- per-element dense operator (matrix-free K_ff) ----------------------------------------------
// A dense 12x12 per beam-column element with its 12 reduced-DOF indices (rid[i] < 0 = constrained /
// unused). Summing the element-local dense matvecs equals K_ff * x by construction: it is built from
// the prepared elements' assemble() triplets (the exact triplets reduceFF builds K_ff from), so the
// apply matches the LDLT operator to the last bit. accumulateRange uses an element-local dense matvec
// to avoid the random-access scatter of a COO apply.
struct ElementBlock12 {
    std::array<int, 12>   rid{};   // reduced-DOF index per local slot (-1 = constrained / unused)
    std::array<real, 144> k{};     // dense 12x12, row-major
};

struct ElementOperator {
    int nf = 0;
    std::vector<ElementBlock12>       blocks;
    VecX                              diag;    // scalar Jacobi diagonal
    std::vector<std::array<real, 36>> diag6;   // 6x6 diagonal blocks (block6); empty if nf%6 != 0

    size_t numBlocks() const { return blocks.size(); }
    const std::array<int, 12>& ridOf(size_t bi) const { return blocks[bi].rid; }
    real kEntry(size_t bi, int r, int c) const { return blocks[bi].k[static_cast<size_t>(r * 12 + c)]; }

    // Build from the prepared elements. Frame-only: returns false on a non-12-DOF element (shell).
    bool build(const PreparedSystem::Impl& S) {
        nf = S.nf;
        diag = VecX::Zero(nf);
        diag6.clear();
        if (nf > 0 && nf % 6 == 0) {
            diag6.resize(static_cast<size_t>(nf / 6));
            for (auto& b : diag6) b.fill(real(0));
        }
        blocks.clear();
        blocks.reserve(S.elems.size());
        for (const auto& el : S.elems) {
            if (el->localDof() != 12) return false;   // shell present -> frame-only, caller routes to LDLT
            std::vector<Triplet> trips;
            el->assemble(trips);
            // distinct global DOFs this element touches (12 for a beam; its diagonal stiffness is
            // positive so every touched DOF appears as a row) -> a compact local index 0..ne-1.
            std::vector<int> gd;
            gd.reserve(trips.size());
            for (const auto& t : trips) gd.push_back(static_cast<int>(t.row()));
            std::sort(gd.begin(), gd.end());
            gd.erase(std::unique(gd.begin(), gd.end()), gd.end());
            if (gd.size() > 12) return false;   // not a 12-DOF beam -> frame-only, caller routes to LDLT
            auto lof = [&](int g) -> int {
                const auto it = std::lower_bound(gd.begin(), gd.end(), g);
                return (it != gd.end() && *it == g) ? static_cast<int>(it - gd.begin()) : -1;
            };
            ElementBlock12 b; b.rid.fill(-1); b.k.fill(real(0));
            for (const auto& t : trips) {
                const int lr = lof(static_cast<int>(t.row()));
                const int lc = lof(static_cast<int>(t.col()));
                if (lr >= 0 && lc >= 0) b.k[static_cast<size_t>(lr * 12 + lc)] += t.value();
            }
            int freeDofs = 0;
            for (int i = 0; i < static_cast<int>(gd.size()); ++i) {
                b.rid[static_cast<size_t>(i)] = S.fmap[static_cast<size_t>(gd[static_cast<size_t>(i)])];
                if (b.rid[static_cast<size_t>(i)] >= 0) ++freeDofs;
            }
            if (freeDofs == 0) continue;   // fully-constrained element: no free contribution
            for (int r = 0; r < 12; ++r) {
                const int rr = b.rid[static_cast<size_t>(r)];
                if (rr >= 0) diag(rr) += b.k[static_cast<size_t>(r * 12 + r)];
                if (rr < 0 || diag6.empty()) continue;
                const int br = rr / 6, lr = rr - br * 6;
                for (int c = 0; c < 12; ++c) {
                    const int cc = b.rid[static_cast<size_t>(c)];
                    if (cc < 0 || cc / 6 != br) continue;
                    const int lc = cc - br * 6;
                    diag6[static_cast<size_t>(br)][static_cast<size_t>(lr * 6 + lc)] +=
                        b.k[static_cast<size_t>(r * 12 + c)];
                }
            }
            blocks.push_back(std::move(b));
        }
        return true;
    }

    // y[begin..end blocks] += K_block * x  (writes only touched reduced DOFs; caller pre-zeros them).
    void accumulateRange(size_t begin, size_t end, const VecX& x, std::vector<real>& y) const {
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

// Invert each 6x6 diagonal block (block6 Jacobi). Returns false (and clears inv6) if any block is
// singular -> caller falls back to scalar Jacobi.
bool buildBlock6(const std::vector<std::array<real, 36>>& diag6,
                 std::vector<std::array<real, 36>>& inv6) {
    inv6.clear();
    if (diag6.empty()) return false;
    inv6.resize(diag6.size());
    using Mat6 = Eigen::Matrix<real, 6, 6>;
    const Mat6 I = Mat6::Identity();
    for (size_t bi = 0; bi < diag6.size(); ++bi) {
        Mat6 M;
        for (int r = 0; r < 6; ++r)
            for (int c = 0; c < 6; ++c) M(r, c) = diag6[bi][static_cast<size_t>(r * 6 + c)];
        Eigen::LDLT<Mat6> ldlt(M);
        if (ldlt.info() != Eigen::Success) { inv6.clear(); return false; }
        const Mat6 inv = ldlt.solve(I);
        if (!inv.allFinite()) { inv6.clear(); return false; }
        for (int r = 0; r < 6; ++r)
            for (int c = 0; c < 6; ++c) inv6[bi][static_cast<size_t>(r * 6 + c)] = inv(r, c);
    }
    return true;
}

// ---- coarse-grid correction (floor-aggregated multi-level preconditioner) -----------------------
// Floor-aggregation by node z (one 6-DOF coarse group per z-level) + a block-Thomas LDL^T (banded)
// factor of the aggregated coarse matrix when it is block-tridiagonal in the floor ordering (a
// tower), else a dense LDL^T. correct(r, z) adds z += P (Kc^{-1} (P^T r)), P = floor prolongation.
// A non-tower model yields floors==1 (or a non-block-tridiagonal Kc -> dense), so the correction
// degrades to a near-no-op and the block6 / scalar Jacobi carries the solve. Serial (main thread).
struct CoarseOperator {
    int  dim    = 0;       // coarse DOF count (0 = inactive)
    bool banded = false;
    int  floors = 0;
    int  fb     = 0;       // floor-block size (= 6; one 6-DOF coarse group per level)
    std::vector<int>               blockToGroup;   // reduced 6-block -> coarse group (-1 = unmapped)
    std::vector<Eigen::LDLT<MatX>> Dfac;           // banded: diagonal floor-block factors
    std::vector<MatX>              Esub;           // banded: sub-diagonal multipliers
    MatX                           denseInv;       // dense fallback inverse
    mutable VecX rc, zc, cg, ch;                   // coarse-solve scratch (main thread only)

    bool active() const { return dim > 0 && !blockToGroup.empty(); }

    void reset() {
        dim = 0; banded = false; floors = 0; fb = 0;
        blockToGroup.clear(); Dfac.clear(); Esub.clear(); denseInv.resize(0, 0);
    }

    // Build the floor-aggregated coarse operator from the element blocks + node z coordinates.
    // Returns true if a usable correction is active (banded for a tower, else dense).
    bool build(const ElementOperator& op, const FrameModel& model, const std::vector<int>& fmap,
               int maxDenseDofs) {
        reset();
        const int nf = op.nf;
        if (nf <= 0 || nf % 6 != 0) return false;
        const int nb = nf / 6;
        blockToGroup.assign(static_cast<size_t>(nb), -1);

        // A node maps to a reduced 6-block iff all 6 of its DOFs are free and contiguous.
        auto blockForNode = [&](size_t ni) -> int {
            int block = -1;
            for (int d = 0; d < 6; ++d) {
                const int g = gdof(static_cast<int>(ni), d);
                const int f = fmap[static_cast<size_t>(g)];
                if (f < 0) return -1;
                if (d == 0) { if (f % 6 != 0) return -1; block = f / 6; }
                else if (f != block * 6 + d) return -1;
            }
            return block;
        };

        // distinct z-levels (sorted) over the mappable nodes
        std::vector<real> levels;
        for (size_t ni = 0; ni < model.nodes.size(); ++ni) {
            if (blockForNode(ni) < 0) continue;
            const real z = model.nodes[ni].pos.z;
            bool found = false;
            for (real e : levels)
                if (std::abs(e - z) <= real(1e-7) * std::max<real>(real(1), std::abs(z))) { found = true; break; }
            if (!found) levels.push_back(z);
        }
        std::sort(levels.begin(), levels.end());
        if (levels.empty()) return false;
        auto groupForZ = [&](real z) -> int {
            int best = 0; real bd = std::abs(z - levels[0]);
            for (int i = 1; i < static_cast<int>(levels.size()); ++i) {
                const real d = std::abs(z - levels[static_cast<size_t>(i)]);
                if (d < bd) { best = i; bd = d; }
            }
            return best;
        };
        floors = static_cast<int>(levels.size());
        fb = 6;                                       // one 6-DOF coarse group per level
        const int nc = floors * fb;
        // Defensive cap on the coarse size (Kc is dense nc x nc): floor-only aggregation gives
        // nc = floors*6 <= nf, so this never fires for a normal tower — it only guards a pathological
        // every-node-its-own-level model, where coarse turns off and the block6 Jacobi carries on.
        if (nc > maxDenseDofs) { reset(); return false; }

        for (size_t ni = 0; ni < model.nodes.size(); ++ni) {
            const int block = blockForNode(ni);
            if (block < 0 || block >= nb) continue;
            blockToGroup[static_cast<size_t>(block)] = groupForZ(model.nodes[ni].pos.z);
        }
        for (int bi = 0; bi < nb; ++bi)
            if (blockToGroup[static_cast<size_t>(bi)] < 0) { reset(); return false; }

        // Kc = P^T K P: aggregate every element's reduced entries into their floor groups.
        MatX Kc = MatX::Zero(nc, nc);
        for (size_t bi = 0; bi < op.numBlocks(); ++bi) {
            const std::array<int, 12>& rid = op.ridOf(bi);
            for (int r = 0; r < 12; ++r) {
                const int rr = rid[static_cast<size_t>(r)]; if (rr < 0) continue;
                const int gr = blockToGroup[static_cast<size_t>(rr / 6)]; if (gr < 0) continue;
                const int cr = gr * 6 + rr % 6;
                for (int c = 0; c < 12; ++c) {
                    const int cc = rid[static_cast<size_t>(c)]; if (cc < 0) continue;
                    const int gc = blockToGroup[static_cast<size_t>(cc / 6)]; if (gc < 0) continue;
                    Kc(cr, gc * 6 + cc % 6) += op.kEntry(bi, r, c);
                }
            }
        }
        dim = nc;
        if (buildBanded(Kc, nc)) { banded = true; return true; }   // tower: block-tridiagonal
        Eigen::LDLT<MatX> ldlt(Kc);                                // else dense fallback
        if (ldlt.info() != Eigen::Success) { reset(); return false; }
        denseInv = ldlt.solve(MatX::Identity(nc, nc));
        if (!denseInv.allFinite()) { reset(); return false; }
        return true;
    }

    // block-Thomas LDL^T factor; false (clears factors) if Kc is not block-tridiagonal or a block
    // is not SPD -> caller uses the dense fallback.
    bool buildBanded(const MatX& Kc, int nc) {
        auto bail = [this]() { Dfac.clear(); Esub.clear(); return false; };
        if (fb <= 0 || nc % fb != 0) return bail();
        const int L = nc / fb;
        const real ref = Kc.cwiseAbs().maxCoeff();
        const real tol = real(1e-9) * std::max<real>(real(1), ref);
        real maxOff = 0;
        for (int i = 0; i < L; ++i)
            for (int j = 0; j < L; ++j)
                if (std::abs(i - j) > 1)
                    maxOff = std::max(maxOff, Kc.block(i * fb, j * fb, fb, fb).cwiseAbs().maxCoeff());
        if (maxOff > tol) return bail();              // not block-tridiagonal -> dense fallback
        Dfac.assign(static_cast<size_t>(L), Eigen::LDLT<MatX>());
        Esub.assign(static_cast<size_t>(L), MatX::Zero(fb, fb));
        Dfac[0].compute(Kc.block(0, 0, fb, fb));
        if (Dfac[0].info() != Eigen::Success) return bail();
        for (int i = 1; i < L; ++i) {
            const MatX Bi = Kc.block(i * fb, (i - 1) * fb, fb, fb);     // A_{i,i-1}
            Esub[static_cast<size_t>(i)] = Dfac[static_cast<size_t>(i - 1)].solve(Bi.transpose()).transpose();
            const MatX Di = Kc.block(i * fb, i * fb, fb, fb) - Esub[static_cast<size_t>(i)] * Bi.transpose();
            Dfac[static_cast<size_t>(i)].compute(Di);
            if (Dfac[static_cast<size_t>(i)].info() != Eigen::Success) return bail();
        }
        // one-time self-check: the block-Thomas factor must solve Kc accurately (reproducible probe).
        VecX probe(nc);
        for (int i = 0; i < nc; ++i)
            probe(i) = std::sin(real(0.013) * static_cast<real>(i + 1)) +
                       real(0.4) * std::cos(real(0.207) * static_cast<real>(i + 3)) +
                       ((i & 1) ? real(-0.25) : real(0.25));
        rc = probe; bandedSolve();
        const double resid = static_cast<double>((Kc * zc - probe).norm() /
                                                 std::max<real>(real(1e-300), probe.norm()));
        if (resid > 1e-8) return bail();
        return true;
    }

    // forward (L g = rc), diagonal (D h = g), back (L^T zc = h); fills zc from rc via cg/ch.
    void bandedSolve() const {
        const int L = floors;
        cg.resize(dim); ch.resize(dim); zc.resize(dim);
        cg.segment(0, fb) = rc.segment(0, fb);
        for (int i = 1; i < L; ++i)
            cg.segment(i * fb, fb).noalias() =
                rc.segment(i * fb, fb) - Esub[static_cast<size_t>(i)] * cg.segment((i - 1) * fb, fb);
        for (int i = 0; i < L; ++i)
            ch.segment(i * fb, fb) = Dfac[static_cast<size_t>(i)].solve(cg.segment(i * fb, fb));
        zc.segment((L - 1) * fb, fb) = ch.segment((L - 1) * fb, fb);
        for (int i = L - 2; i >= 0; --i)
            zc.segment(i * fb, fb).noalias() =
                ch.segment(i * fb, fb) - Esub[static_cast<size_t>(i + 1)].transpose() * zc.segment((i + 1) * fb, fb);
    }

    // z += P (Kc^{-1} (P^T r)) — restrict the residual to floor groups, coarse-solve, prolong back.
    void correct(const VecX& r, VecX& z) const {
        rc.setZero(dim);
        for (int bi = 0; bi < static_cast<int>(blockToGroup.size()); ++bi) {
            const int g = blockToGroup[static_cast<size_t>(bi)]; if (g < 0) continue;
            for (int d = 0; d < 6; ++d) rc(g * 6 + d) += r(bi * 6 + d);
        }
        if (banded) bandedSolve();
        else        zc.noalias() = denseInv * rc;
        for (int bi = 0; bi < static_cast<int>(blockToGroup.size()); ++bi) {
            const int g = blockToGroup[static_cast<size_t>(bi)]; if (g < 0) continue;
            for (int d = 0; d < 6; ++d) z(bi * 6 + d) += zc(g * 6 + d);
        }
    }
};

// ---- persistent worker pool (data-parallel element apply + block6 map) --------------------------
// Warm threads serve two jobs over the element operator: (1) a matrix-free apply (each worker
// accumulates its block range into a private localY, the main thread reduces the disjoint touched-
// DOF sets), and (2) a block6 Jacobi map (each worker writes a disjoint output slice, no reduction).
// A condition_variable generation protocol wakes the workers; reusing them avoids per-apply spawn.
// Eigen is not touched off the main thread — workers only read the const operator and write into
// plain std::vector<real> / VecX slices.
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
        const int nb = (op.nf % 6 == 0) ? op.nf / 6 : 0;
        const int pchunk = (nb + nt - 1) / std::max(1, nt);
        pcRange.assign(static_cast<size_t>(nt), {0, 0});
        applyRange.assign(static_cast<size_t>(nt), {0, 0});
        workers.reserve(static_cast<size_t>(nt));
        try {
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
        } catch (...) {
            // A worker failed to spawn (e.g. thread-resource exhaustion). Wake and join the ones
            // already started, then rethrow — otherwise their still-joinable std::thread dtors would
            // call std::terminate. The HpSession ctor catches this and falls back to the serial path.
            {
                std::lock_guard<std::mutex> lock(mu);
                stop = true;
                ++generation;
            }
            cvStart.notify_all();
            for (std::thread& w : workers) if (w.joinable()) w.join();
            throw;
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

    ThreadApplyPool(const ThreadApplyPool&) = delete;
    ThreadApplyPool& operator=(const ThreadApplyPool&) = delete;

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
            else                        runPrecond(tid);
            {
                std::lock_guard<std::mutex> lock(mu);
                ++done;
                if (done == nt) cvDone.notify_one();
            }
        }
    }
};

}  // namespace

// ============================================================================
// PIMPL: non-owning view of the factorization + the seeded basis + apply/precond (serial or pooled).
// ============================================================================
struct HpSession::Impl {
    const PreparedSystem::Impl*       ps_raw = nullptr;   // NON-owning; caller guarantees it outlives us
    HpSessionOptions                  opts;
    bool                              baseValid = false;   // HP path ready (non-null, non-singular, frame-only)
    std::string                       diag;
    RecycleBasis                      basis;
    int                               effectiveBasisMax = 0;
    int                               threads = 1;
    ElementOperator                   op;                  // per-element dense operator (frame-only)
    VecX                              cachedFequiv;        // structure-fixed distributed-load equiv. forces
    std::vector<std::array<real, 36>> blockInv6;           // block6 Jacobi (empty -> scalar Jacobi)
    bool                              useBlock6 = false;
    CoarseOperator                    coarse;               // A2b coarse correction (inactive until prepareCoarse)
    // Declared AFTER op so it is DESTROYED FIRST: ~ThreadApplyPool joins the workers (which hold a
    // reference to op) before op itself is torn down.
    std::unique_ptr<ThreadApplyPool>  pool;

    Impl(const PreparedSystem& prepared, const HpSessionOptions& o)
        : ps_raw(prepared.impl.get()), opts(o),
          basis(std::max(0, o.basisMax)), effectiveBasisMax(std::max(0, o.basisMax)),
          threads(std::max(1, o.threads)) {
        if (!ps_raw) { diag = "HpSession: null PreparedSystem"; return; }
        const PreparedSystem::Impl& S = *ps_raw;
        if (S.singular) { diag = "HpSession: PreparedSystem is singular: " + S.diagnostic; return; }

        // Distributed-load (UDL / pressure / self-weight) equivalent nodal forces are structure-fixed
        // — the fingerprint guards their inputs — so cache them once instead of re-scattering every
        // frame. Built for both frame and shell models (the shell path still needs them for LDLT).
        cachedFequiv = VecX::Zero(S.N);
        for (const auto& el : S.elems) el->addEquivalentNodalLoads(cachedFequiv);

        if (!op.build(S)) {
            // a non-12-DOF element (shell): frame-only -> HP path not ready, solveFrame uses LDLT.
            baseValid = false;
            diag = "HpSession: shell element present (frame-only; solveFrame uses LDLT)";
            return;
        }
        baseValid = true;
        useBlock6 = buildBlock6(op.diag6, blockInv6);   // false (empty inv6) -> scalar Jacobi
        if (threads > 1 && op.numBlocks() > 0) {
            try {
                pool = std::make_unique<ThreadApplyPool>(op, threads);
            } catch (...) {
                pool.reset();   // thread-spawn failure -> fall back to serial (still fully functional)
            }
        }
        const int actual = pool ? pool->threads() : 1;
        diag = "HpSession: ready (frame-only seeded lane, " +
               std::string(actual > 1 ? "parallel x" + std::to_string(actual) : "serial") + ", " +
               std::string(useBlock6 ? "block6" : "scalar") + " precond)";
    }

    // y = K_ff * x (matrix-free; pooled when a thread pool exists, else serial). Bit-consistent with
    // the LDLT K_ff up to the reduction order (threaded reduction is NOT bit-identical to serial).
    void applyKff(const VecX& x, VecX& y) {
        if (pool) pool->apply(x, y);
        else      op.applySerial(x, y);
    }

    // z = M^{-1} r: block6 Jacobi when available (pooled or serial), else scalar Jacobi.
    void precond(const VecX& rr, VecX& zz) {
        zz.resize(op.nf);
        if (useBlock6 && !blockInv6.empty()) {
            if (pool) {
                pool->precondBlock6(rr, zz, blockInv6);
            } else {
                for (int bi = 0; bi < static_cast<int>(blockInv6.size()); ++bi) {
                    const real* inv = blockInv6[static_cast<size_t>(bi)].data();
                    const int off = bi * 6;
                    for (int row = 0; row < 6; ++row) {
                        const real* ir = &inv[row * 6];
                        zz(off + row) = ir[0] * rr(off + 0) + ir[1] * rr(off + 1) + ir[2] * rr(off + 2)
                                      + ir[3] * rr(off + 3) + ir[4] * rr(off + 4) + ir[5] * rr(off + 5);
                    }
                }
            }
        } else {
            for (int i = 0; i < op.nf; ++i) zz(i) = (op.diag(i) != real(0)) ? rr(i) / op.diag(i) : rr(i);
        }
        if (coarse.active()) coarse.correct(rr, zz);   // A2b coarse-grid correction (out-of-subspace relief)
    }
};

// ============================================================================
// Public surface.
// ============================================================================
HpSession::HpSession(const PreparedSystem& prepared, const HpSessionOptions& opts)
    : p_(std::make_unique<Impl>(prepared, opts)) {}
HpSession::~HpSession() = default;
HpSession::HpSession(HpSession&&) noexcept = default;
HpSession& HpSession::operator=(HpSession&&) noexcept = default;

bool HpSession::valid() const { return p_->baseValid; }
const std::string& HpSession::diagnostic() const { return p_->diag; }

bool HpSession::setLoadBasis(const std::vector<std::vector<real>>& loadVectors) {
    Impl& P = *p_;
    if (!P.baseValid) return false;
    const PreparedSystem::Impl& S = *P.ps_raw;
    const int N = S.N;

    // FIFO eviction would drop early seeds and break the in-subspace guarantee, so size the cap to
    // hold every seed (A1/B hit this blind spot). Rebuild the basis from scratch.
    P.effectiveBasisMax = std::max(P.opts.basisMax, (int)loadVectors.size());
    P.basis = RecycleBasis(P.effectiveBasisMax);

    for (const auto& lv : loadVectors) {
        if ((int)lv.size() != N) {
            P.diag = "HpSession::setLoadBasis: load vector size " + std::to_string(lv.size()) +
                     " != 6N (" + std::to_string(N) + ")";
            return false;
        }
        // The seed is a PURE nodal load (prescribed = 0): reduce the 6N global load to the free DOFs,
        // solve via the already-factored LDLT, and A-orthonormalize the response into the basis.
        VecX Ff = VecX::Zero(S.nf);
        for (int g = 0; g < N; ++g)
            if (S.fmap[g] >= 0) Ff(S.fmap[g]) += lv[(size_t)g];
        const VecX resp = S.ldlt.solve(Ff);
        if (S.ldlt.info() != Eigen::Success || !resp.allFinite()) continue;   // skip a degenerate seed
        P.basis.add(resp, [&](const VecX& q, VecX& Aq) { P.applyKff(q, Aq); });
    }
    return true;
}

HpCoarseInfo HpSession::prepareCoarse(const FrameModel& model) {
    Impl& P = *p_;
    HpCoarseInfo info;
    if (!P.baseValid) return info;
    const PreparedSystem::Impl& S = *P.ps_raw;
    // The model's node coordinates define the floor aggregation; its fingerprint must still match the
    // prepared system or the coarse map would not line up with K_ff.
    if (modelFingerprint(model) != S.fingerprint) {
        P.diag = "HpSession::prepareCoarse: model changed since assembleAndFactor (fingerprint mismatch)";
        return info;
    }
    const bool ok = P.coarse.build(P.op, model, S.fmap, /*maxDenseDofs=*/P.op.nf);
    info.built  = ok;
    info.banded = P.coarse.banded;
    info.floors = P.coarse.floors;
    info.dim    = P.coarse.dim;
    P.diag = "HpSession: coarse " +
             std::string(!ok ? "unavailable (block6/scalar carries the solve)"
                        : P.coarse.banded ? "banded block-Thomas (floors=" + std::to_string(P.coarse.floors) + ")"
                                          : "dense (dim=" + std::to_string(P.coarse.dim) + ")");
    return info;
}

SolveResult HpSession::solveFrame(const FrameModel& model, HpSessionStats* stats) {
    Impl& P = *p_;
    HpSessionStats st;
    st.basisSize = (int)P.basis.v.size();
    SolveResult R;

    if (!P.ps_raw) {
        R.singular = true; R.diagnostic = "HpSession::solveFrame: null PreparedSystem";
        if (stats) *stats = st;
        return R;
    }
    const PreparedSystem::Impl& S = *P.ps_raw;
    const int N = S.N;
    R.u.assign((size_t)std::max(0, N), 0.0);
    R.reactions.assign((size_t)std::max(0, N), 0.0);
    if (S.singular) {
        R.singular = true; R.diagnostic = S.diagnostic;
        if (stats) *stats = st;
        return R;
    }

    // Reuse-validity guard — identical to solveLoad / solveLoadHP (shared modelFingerprint).
    if (modelFingerprint(model) != S.fingerprint) {
        R.singular = true;
        R.diagnostic = "HpSession::solveFrame: model changed since assembleAndFactor (geometry/topology/"
                       "support flags/distributed loads). Re-run assembleAndFactor.";
        if (stats) *stats = st;
        return R;
    }

    // ---- RHS assembly: per-frame nodal loads + the cached structure-fixed distributed-load forces
    // (solveLoad re-scatters addEquivalentNodalLoads every call; the result is identical because
    // those equivalent forces only depend on fingerprint-guarded geometry/loads). ---------------
    VecX F = P.cachedFequiv;
    for (const auto& nl : model.nodalLoads) {
        const int ni = model.nodeIndex(nl.node);
        if (ni < 0) continue;
        for (int d = 0; d < 6; ++d) F(gdof(ni, d)) += nl.comp[d];
    }

    std::vector<real> presc((size_t)N, 0.0);
    bool anyPresc = false;
    for (size_t k = 0; k < model.nodes.size(); ++k)
        for (int d = 0; d < 6; ++d)
            if (model.nodes[k].fixed[d]) {
                presc[(size_t)gdof((int)k, d)] = model.nodes[k].prescribed[d];
                if (model.nodes[k].prescribed[d] != real(0)) anyPresc = true;
            }

    VecX Ff = VecX::Zero(S.nf);
    // The prescribed (support-settlement) reduce is an O(nnz) sweep of K; skip it entirely when no
    // support carries a nonzero prescribed value — the common per-frame case (gravity + contacts are
    // nodal loads). An all-zero presc contributes nothing, so the result is identical to solveLoad.
    if (anyPresc)
        for (int c = 0; c < N; ++c)
            for (SpMat::InnerIterator it(S.K, c); it; ++it) {
                const int r = it.row();
                if (S.fmap[r] < 0) continue;
                if (S.fmap[c] < 0 && presc[(size_t)c] != 0.0) Ff(S.fmap[r]) -= it.value() * presc[(size_t)c];
            }
    for (int g = 0; g < N; ++g) if (S.fmap[g] >= 0) Ff(S.fmap[g]) += F(g);

    // ---- the ONLY departure: seeded projection / warm PCG / full PCG / LDLT safety net ----
    const bool canHp = P.opts.enabled && P.baseValid && S.nf > 0 && !P.basis.v.empty();
    const real bnorm = Ff.norm();
    VecX uf;
    bool solved = false;
    std::string note;

    if (canHp && !(bnorm > real(0))) {
        // A zero RHS has the exact solution uf = 0 (same hardening as HpSolver.cpp). Trivially in-subspace.
        uf = VecX::Zero(S.nf);
        solved = true;
        st.usedProjection = true;
        st.initialRel = 0;
        note = "[HpSession] zero RHS -> uf=0";
    } else if (canHp) {
        VecX x0, Kx0;
        P.basis.initialGuess(S.nf, Ff, x0, Kx0);                 // Kx0 = K x0 assembled from av_i: no apply
        const real initRel = (Ff - Kx0).norm() / bnorm;          // bnorm > 0 here -> well-posed
        st.initialRel = initRel;
        const bool inSub = (static_cast<double>(initRel) < P.opts.projGateTol);
        const int  maxIt = inSub ? P.opts.pcgWarmIter : P.opts.pcgMaxIter;

        // Warm-start PCG from x0 (a residual safety net atop the projection, never raw x0). Count a
        // step only AFTER x advances and test convergence right after, so pcgIters is honest.
        VecX x = x0;
        VecX r = Ff - Kx0;
        bool conv = static_cast<double>(r.norm() / bnorm) <= P.opts.pcgTol;
        int  iters = 0;
        if (!conv) {
            st.usedCoarse = P.coarse.active();   // the (coarse-augmented) preconditioner is exercised here
            VecX z(S.nf), p(S.nf), Ap(S.nf);
            P.precond(r, z);
            p = z;
            real rz = r.dot(z);
            while (!conv && iters < maxIt) {
                P.applyKff(p, Ap);
                const real pAp = p.dot(Ap);
                if (!(pAp > real(0)) || !std::isfinite(static_cast<double>(pAp))) break;   // indefinite -> LDLT
                const real alpha = rz / pAp;
                x.noalias() += alpha * p;
                r.noalias() -= alpha * Ap;
                ++iters;
                if (static_cast<double>(r.norm() / bnorm) <= P.opts.pcgTol) { conv = true; break; }
                P.precond(r, z);
                const real rzNew = r.dot(z);
                if (!(rzNew > real(0)) || !std::isfinite(static_cast<double>(rzNew))) break;
                p = z + (rzNew / rz) * p;
                rz = rzNew;
            }
        }
        st.pcgIters = iters;
        if (conv && x.allFinite()) {
            uf = x;
            solved = true;
            if (inSub) { st.usedProjection = true; note = "[HpSession] seeded projection + warm PCG (" + std::to_string(iters) + " iter)"; }
            else       { st.usedPcg        = true; note = "[HpSession] full PCG (" + std::to_string(iters) + " iter)"; }
        }
        // not converged (or indefinite) -> fall through to the LDLT safety net below
    }

    if (!solved) {
        if (canHp && !P.opts.fallbackOnFail) {
            R.singular = true;
            R.diagnostic = "HpSession::solveFrame: HP path did not converge; fallback disabled";
            if (stats) *stats = st;
            return R;
        }
        uf = S.ldlt.solve(Ff);                      // the oracle / safety net
        if (S.ldlt.info() != Eigen::Success || !uf.allFinite()) {
            R.singular = true;
            R.diagnostic = "HpSession::solveFrame: LDLT fallback produced non-finite displacements (mechanism)";
            if (stats) *stats = st;
            return R;
        }
        st.usedLdlt = true;
        if (note.empty()) note = "[HpSession] LDLT (disabled / un-seeded / shell / non-converge)";
    }

    // ---- scatter (always); reactions + recover (full result, or skipped on the fast path) ----
    VecX u = VecX::Zero(N);
    for (int g = 0; g < N; ++g) u(g) = (S.fmap[g] >= 0) ? uf(S.fmap[g]) : presc[(size_t)g];
    for (int g = 0; g < N; ++g) R.u[(size_t)g] = u(g);

    if (!P.opts.displacementsOnly) {
        // reactions / recover: verbatim from solveLoad (sync if changed). On the displacements-only
        // fast path these O(nnz)/O(elements) sweeps are skipped: reactions stay all-zero (assigned
        // above) and memberForces/shellForces stay empty.
        const VecX Rv = S.K * u - F;
        for (int g = 0; g < N; ++g) R.reactions[(size_t)g] = Rv(g);

        R.memberForces.resize(model.members.size());
        R.shellForces.resize(model.shells.size());
        for (size_t e = 0; e < model.members.size(); ++e) R.memberForces[e].member = model.members[e].id;
        for (size_t s = 0; s < model.shells.size(); ++s)  R.shellForces[s].shell   = model.shells[s].id;
        for (const auto& el : S.elems) el->recover(u, R);
    }
    R.pivotMargin = S.pivotMargin;
    R.diagnostic = note;   // which path ran; does not affect the numeric result
    if (stats) *stats = st;
    return R;
}

}  // namespace frame
