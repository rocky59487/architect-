#include "FrameCore/Reanalysis.h"
#include "PreparedSystemImpl.h"
#include "IElement.h"
#include "BeamColumnElement.h"
#include "MITC4ShellElement.h"
#include "ElementFactory.h"

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

namespace frame {

namespace {
// makeMemberElem / makeShellElem are shared from ElementFactory.h.

// Eigen preconditioner whose solve() forwards to the STALE baseline LDLT — turns Eigen's CG into a
// stale-factor-preconditioned solve of the modified K' (Tier-2). Mirrors the research prototype.
struct StalePrecond {
    const LDLTSolver* ldlt = nullptr;
    StalePrecond() = default;
    template <class M> explicit StalePrecond(const M&) {}
    template <class M> StalePrecond& analyzePattern(const M&) { return *this; }
    template <class M> StalePrecond& factorize(const M&) { return *this; }
    template <class M> StalePrecond& compute(const M&) { return *this; }
    template <class Rhs> VecX solve(const Rhs& b) const { return ldlt->solve(b); }
    Eigen::ComputationInfo info() const { return Eigen::Success; }
};

// reduceFF (free-free submatrix extraction) is shared from PreparedSystemImpl.h.
}  // namespace

// ============================================================================
// PIMPL: baseline factorization + the Woodbury ladder.
// ============================================================================
struct ReSolveSession::Impl {
    FrameModel        base;          // immutable reference state (caller's model is copied in)
    ReanalysisOptions opts;

    PreparedSystem    ps;            // current baseline factorization (rebuilt on ctor + rebaseline)
    bool              baseValid = false;
    std::string       diag;

    std::vector<char> memberActive;  // current active flags relative to base.members
    std::vector<char> shellActive;   // current active flags relative to base.shells

    // K_cur_ff = K0_ff + W diag(s) Wᵀ  (the accumulated change vs the baseline active set).
    MatX W, Z, Mc;
    VecX s;
    int  R = 0;

    Impl(const FrameModel& b, const ReanalysisOptions& o) : base(b), opts(o) {
        memberActive.resize(base.members.size());
        shellActive.resize(base.shells.size());
        for (size_t i = 0; i < base.members.size(); ++i) memberActive[i] = base.members[i].active ? 1 : 0;
        for (size_t i = 0; i < base.shells.size(); ++i)  shellActive[i]  = base.shells[i].active  ? 1 : 0;
        buildBaseline();
    }

    const PreparedSystem::Impl& S() const { return *ps.impl; }

    FrameModel workModel() const {
        FrameModel m = base;
        for (size_t i = 0; i < m.members.size(); ++i) m.members[i].active = memberActive[i] != 0;
        for (size_t i = 0; i < m.shells.size(); ++i)  m.shells[i].active  = shellActive[i]  != 0;
        return m;
    }

    void clearLadder() { W.resize(0, 0); Z.resize(0, 0); Mc.resize(0, 0); s.resize(0); R = 0; }

    // Rebuild the baseline on the current active set and clear the ladder (= Tier-3, always correct).
    void buildBaseline() {
        ps = assembleAndFactor(workModel(), opts.solve);
        baseValid = !S().singular;
        diag      = S().diagnostic;
        clearLadder();
    }

    // Append one toggled element's signed low-rank contribution to the ladder.
    //   kind: 0 = member, 1 = shell;  idx = element index;  sign = +1 activate / -1 deactivate.
    // Returns false (sets diag) on a prepare error. A fully-fixed element adds rank 0 (returns true).
    bool addToLadder(int kind, int idx, int sign) {
        const PreparedSystem::Impl& Si = S();
        std::unique_ptr<IElement> el = (kind == 0) ? makeMemberElem(idx) : makeShellElem(idx);
        std::string why;
        if (!el->prepare(base, opts.solve, why)) { diag = why; return false; }

        std::vector<Triplet> trips;
        el->assemble(trips);   // global-DOF triplets (bit-consistent with assembleAndFactor)

        // distinct global DOFs the element touches (its diagonal stiffness is positive, so every
        // touched DOF appears as a row) -> a compact local index 0..ne-1.
        std::vector<int> dofs;
        dofs.reserve(24);
        for (const auto& t : trips) dofs.push_back((int)t.row());
        std::sort(dofs.begin(), dofs.end());
        dofs.erase(std::unique(dofs.begin(), dofs.end()), dofs.end());
        const int ne = (int)dofs.size();
        auto lof = [&](int g) -> int {
            const auto it = std::lower_bound(dofs.begin(), dofs.end(), g);
            return (it != dofs.end() && *it == g) ? (int)(it - dofs.begin()) : -1;
        };
        MatX Ke = MatX::Zero(ne, ne);
        for (const auto& t : trips) {
            const int lr = lof((int)t.row()), lc = lof((int)t.col());
            if (lr >= 0 && lc >= 0) Ke(lr, lc) += t.value();
        }

        // restrict to free DOFs (constrained DOFs drop out of the reduced Woodbury update)
        std::vector<int> keep, rid;
        for (int k = 0; k < ne; ++k) {
            const int f = Si.fmap[(size_t)dofs[k]];
            if (f >= 0) { keep.push_back(k); rid.push_back(f); }
        }
        const int nred = (int)keep.size();
        if (nred == 0) return true;   // fully-fixed element: no free contribution

        MatX Ks(nred, nred);
        for (int a = 0; a < nred; ++a)
            for (int b = 0; b < nred; ++b) Ks(a, b) = Ke(keep[(size_t)a], keep[(size_t)b]);
        Eigen::SelfAdjointEigenSolver<MatX> es(Ks);
        const VecX lam  = es.eigenvalues();
        const real lmax = lam.cwiseAbs().maxCoeff();
        std::vector<int> kept;
        for (int i = 0; i < nred; ++i) if (lam(i) > real(1e-9) * lmax) kept.push_back(i);  // PSD: keep positive modes
        const int r = (int)kept.size();
        if (r == 0) return true;

        MatX Wn = MatX::Zero(Si.nf, r);
        VecX sn(r);
        for (int q = 0; q < r; ++q) {
            for (int a = 0; a < nred; ++a) Wn(rid[(size_t)a], q) = es.eigenvectors()(a, kept[(size_t)q]);
            sn(q) = sign * lam(kept[(size_t)q]);
        }
        const MatX Zn = Si.ldlt.solve(Wn);

        MatX McNew(R + r, R + r);
        if (R > 0) {
            const MatX B1 = Wn.transpose() * Z;     // r x R  (= (W_old^T Zn)^T by symmetry of K0^-1)
            McNew.topLeftCorner(R, R)    = Mc;
            McNew.topRightCorner(R, r)   = B1.transpose();
            McNew.bottomLeftCorner(r, R) = B1;
        }
        McNew.bottomRightCorner(r, r) = Wn.transpose() * Zn;

        W.conservativeResize(Si.nf, R + r); W.rightCols(r) = Wn;
        Z.conservativeResize(Si.nf, R + r); Z.rightCols(r) = Zn;
        s.conservativeResize(R + r);        s.tail(r) = sn;
        Mc = McNew;
        R += r;
        return true;
    }

    // Woodbury correction: given u0ff = K0^-1 Ff, return K_cur^-1 Ff; set mechanism / pivot ratio.
    VecX woodbury(const VecX& u0ff, bool& mech, real& pivRatio) const {
        if (R == 0) { mech = false; pivRatio = 1; return u0ff; }
        MatX C = Mc;
        for (int i = 0; i < R; ++i) C(i, i) += real(1) / s(i);   // C = diag(1/s) + Wᵀ K0^-1 W
        Eigen::FullPivLU<MatX> lu(C);
        const VecX d = lu.matrixLU().diagonal().cwiseAbs();
        const real dmax = d.maxCoeff(), dmin = d.minCoeff();
        pivRatio = (dmax > 0) ? dmin / dmax : 0;
        mech = (pivRatio < opts.mechPivotTol);                   // C singular <=> K' singular <=> mechanism
        const VecX q = lu.solve(W.transpose() * u0ff);
        return u0ff - Z * q;
    }
};

// ============================================================================
// Public surface.
// ============================================================================
ReSolveSession::ReSolveSession(const FrameModel& base, const ReanalysisOptions& opts)
    : p_(std::make_unique<Impl>(base, opts)) {}
ReSolveSession::~ReSolveSession() = default;
ReSolveSession::ReSolveSession(ReSolveSession&&) noexcept = default;
ReSolveSession& ReSolveSession::operator=(ReSolveSession&&) noexcept = default;

bool ReSolveSession::valid() const { return p_->baseValid; }
const std::string& ReSolveSession::diagnostic() const { return p_->diag; }
void ReSolveSession::rebaseline() { p_->buildBaseline(); }

bool ReSolveSession::setMemberActive(MemberId id, bool active) {
    Impl& P = *p_;
    int idx = -1;
    for (size_t i = 0; i < P.base.members.size(); ++i)
        if (P.base.members[i].id == id) { idx = (int)i; break; }
    if (idx < 0) return false;
    if ((P.memberActive[(size_t)idx] != 0) == active) return true;   // no-op
    P.memberActive[(size_t)idx] = active ? 1 : 0;
    if (P.baseValid && !P.addToLadder(0, idx, active ? +1 : -1))
        P.buildBaseline();   // prepare failure -> fall back to the always-correct rebaseline
    return true;
}

bool ReSolveSession::setShellActive(int shellId, bool active) {
    Impl& P = *p_;
    int idx = -1;
    for (size_t i = 0; i < P.base.shells.size(); ++i)
        if (P.base.shells[i].id == shellId) { idx = (int)i; break; }
    if (idx < 0) return false;
    if ((P.shellActive[(size_t)idx] != 0) == active) return true;    // no-op
    P.shellActive[(size_t)idx] = active ? 1 : 0;
    if (P.baseValid && !P.addToLadder(1, idx, active ? +1 : -1))
        P.buildBaseline();
    return true;
}

SolveResult ReSolveSession::solve(ReanalysisStats* stats) {
    Impl& P = *p_;
    SolveResult R;
    ReanalysisStats st;
    if (!P.baseValid) {
        R.singular = true; R.diagnostic = P.diag;
        if (stats) *stats = st;
        return R;
    }

    // fmap / N / nf are support-determined and invariant under element activity, so they survive a
    // rebaseline; copy them so the references stay valid even if ps is rebuilt mid-solve.
    const std::vector<int> fmap = P.S().fmap;
    const int N  = P.S().N;
    const int nf = P.S().nf;
    const FrameModel work = P.workModel();

    // --- current active elements (transient, bit-consistent with assembleAndFactor) ---
    std::vector<std::unique_ptr<IElement>> elems;
    elems.reserve(work.members.size() + work.shells.size());
    for (size_t e = 0; e < work.members.size(); ++e)
        if (work.members[e].active) elems.push_back(makeMemberElem((int)e));
    for (size_t sh = 0; sh < work.shells.size(); ++sh)
        if (work.shells[sh].active) elems.push_back(makeShellElem((int)sh));
    for (auto& el : elems) {
        std::string why;
        if (!el->prepare(work, P.opts.solve, why)) {
            R.singular = true; R.diagnostic = "ReSolve: element prepare failed: " + why;
            if (stats) *stats = st;
            return R;
        }
    }

    // --- K_cur (assembled, NOT factored) for the prescribed term + reactions ---
    std::vector<Triplet> trips;
    for (auto& el : elems) el->assemble(trips);
    SpMat Kcur(N, N);
    Kcur.setFromTriplets(trips.begin(), trips.end());
    Kcur.makeCompressed();

    // --- F_cur(N) = nodal loads + active elements' equivalent (UDL/pressure) loads ---
    //     This is the F-increment: a deactivated member/shell drops out of the active set, so its
    //     equivalent load leaves F automatically; a restored one re-enters.
    VecX F = VecX::Zero(N);
    for (const auto& nl : work.nodalLoads) {
        const int ni = work.nodeIndex(nl.node);
        if (ni < 0) continue;
        for (int d = 0; d < 6; ++d) F(gdof(ni, d)) += nl.comp[d];
    }
    for (auto& el : elems) el->addEquivalentNodalLoads(F);

    // --- prescribed (support-displacement) values ---
    std::vector<real> presc((size_t)N, 0.0);
    for (size_t k = 0; k < work.nodes.size(); ++k)
        for (int d = 0; d < 6; ++d)
            if (work.nodes[k].fixed[d]) presc[(size_t)gdof((int)k, d)] = work.nodes[k].prescribed[d];

    // --- reduced RHS  Ff = F_f - K_cur_fc u_c ---
    VecX Ff = VecX::Zero(nf);
    for (int c = 0; c < N; ++c)
        for (SpMat::InnerIterator it(Kcur, c); it; ++it) {
            const int rr = it.row();
            if (fmap[(size_t)rr] < 0) continue;
            if (fmap[(size_t)c] < 0 && presc[(size_t)c] != 0.0) Ff(fmap[(size_t)rr]) -= it.value() * presc[(size_t)c];
        }
    for (int g = 0; g < N; ++g) if (fmap[(size_t)g] >= 0) Ff(fmap[(size_t)g]) += F(g);

    // --- displacement: reuse the baseline factor (Tier-0/1), else rebaseline (Tier-3) ---
    bool mech = false; real pivRatio = 1; VecX uf;
    if (P.R == 0) {
        st.tier = 0;
        uf = P.S().ldlt.solve(Ff);
    } else if (P.R <= P.opts.maxRank) {
        st.tier = 1;
        const VecX u0ff = P.S().ldlt.solve(Ff);
        uf = P.woodbury(u0ff, mech, pivRatio);
    } else if (P.opts.allowTier2 && P.R <= 2 * P.opts.maxRank) {
        // Tier-2: stale-LDLT preconditioned CG on the (assembled, not factored) K_cur_ff. Reuses the
        // baseline factor as the preconditioner + the baseline solve as the warm-start guess.
        // Tolerance-level (NOT bit-identical); on non-convergence fall through to Tier-3 (always correct).
        const VecX  u0ff = P.S().ldlt.solve(Ff);
        const SpMat Kff  = reduceFF(Kcur, fmap, nf);
        Eigen::ConjugateGradient<SpMat, Eigen::Lower | Eigen::Upper, StalePrecond> cg;
        cg.preconditioner().ldlt = &P.S().ldlt;
        cg.setTolerance(P.opts.pcgTol);
        cg.setMaxIterations(P.opts.pcgMaxIter);
        cg.compute(Kff);
        uf = cg.solveWithGuess(Ff, u0ff);
        st.pcgIters    = (int)cg.iterations();
        st.relResidual = cg.error();
        if (cg.info() == Eigen::Success && uf.allFinite()) {
            st.tier = 2;
        } else {
            st.tier = 3; st.refactored = true;               // CG stalled -> rebaseline
            P.buildBaseline();
            if (!P.baseValid) { R.singular = true; R.diagnostic = P.diag; if (stats) *stats = st; return R; }
            uf = P.S().ldlt.solve(Ff);
        }
    } else {
        st.tier = 3; st.refactored = true;
        P.buildBaseline();                                   // rebaseline on the current set
        if (!P.baseValid) {
            R.singular = true; R.diagnostic = P.diag;
            if (stats) *stats = st;
            return R;
        }
        uf = P.S().ldlt.solve(Ff);                           // Ff still valid (fmap invariant)
    }
    st.rank      = P.R;
    st.mechanism = mech;
    const real pivMargin = P.S().pivotMargin;

    if (mech) {
        R.singular   = true;
        R.diagnostic = "ReSolve: capacitance singular (removed set forms a mechanism)";
        R.pivotMargin = pivMargin;
        if (stats) *stats = st;
        return R;
    }
    if (!uf.allFinite()) {
        R.singular   = true;
        R.diagnostic = "ReSolve: non-finite displacements (mechanism)";
        if (stats) *stats = st;
        return R;
    }

    // --- assemble the full result (mirror solveLoad) ---
    R.u.assign((size_t)N, 0.0);
    R.reactions.assign((size_t)N, 0.0);
    VecX u = VecX::Zero(N);
    for (int g = 0; g < N; ++g) u(g) = (fmap[(size_t)g] >= 0) ? uf(fmap[(size_t)g]) : presc[(size_t)g];
    for (int g = 0; g < N; ++g) R.u[(size_t)g] = u(g);
    const VecX Rv = Kcur * u - F;
    for (int g = 0; g < N; ++g) R.reactions[(size_t)g] = Rv(g);

    R.memberForces.resize(work.members.size());
    R.shellForces.resize(work.shells.size());
    for (size_t e = 0; e < work.members.size(); ++e) R.memberForces[e].member = work.members[e].id;
    for (size_t sh = 0; sh < work.shells.size(); ++sh) R.shellForces[sh].shell = work.shells[sh].id;
    for (auto& el : elems) el->recover(u, R);
    R.pivotMargin = pivMargin;
    if (stats) *stats = st;
    return R;
}

}  // namespace frame
