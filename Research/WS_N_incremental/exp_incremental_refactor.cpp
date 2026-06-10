// exp_incremental_refactor.cpp — N1 "ReSolve ladder" research experiment (scratch, NOT engine).
//
// Question: after element removal/restore, can we solve the modified system EXACTLY on top
// of the baseline SimplicialLDLT instead of a fresh assembleAndFactor?
//
//   Tier-1  exact Woodbury rank-(<=6 per member) updates:
//             K' = K0 + W diag(s) W^T,
//             K'^-1 b = K0^-1 b - Z C^-1 (W^T K0^-1 b),   Z = K0^-1 W,
//             C = diag(1/s) + W^T Z  (capacitance; SINGULAR  <=>  K' singular  <=> mechanism)
//   Tier-2  stale-LDLT preconditioned CG for batch changes too large for tier-1
//           (preconditioner = the baseline factorization; iterations/time vs fresh refactor
//            and vs Jacobi-PCG as contrast).
//
// Measures: accuracy vs fresh refactor, drift over long remove+restore sequences, speedup,
// mechanism detection, tier-2 iteration counts. One structured stdout line per measurement.
//
// Honest scope: ladder timings are u-only (no member-force recovery); the fresh path includes
// the engine's recover loop. The pure back-substitution time is printed for a fair reading.

#include "research_common.h"
#include <Eigen/IterativeLinearSolvers>
#include <algorithm>
#include <cstring>
#include <string>
#include <vector>

using namespace research;

namespace {

struct Args {
    std::string preset = "xxl";   // small | xxl | mega
    int  steps    = 50;
    int  refEvery = 1;            // fresh reference every k-th tier-1 step
    bool doRestore = true;
    bool doTier2   = true;
    bool doMech    = true;
};

Args parseArgs(int argc, char** argv) {
    Args a;
    for (int i = 1; i < argc; ++i) {
        auto next = [&]() -> const char* { return (i + 1 < argc) ? argv[++i] : ""; };
        if      (!std::strcmp(argv[i], "--preset"))   a.preset = next();
        else if (!std::strcmp(argv[i], "--steps"))    a.steps = std::atoi(next());
        else if (!std::strcmp(argv[i], "--refEvery")) a.refEvery = std::max(1, std::atoi(next()));
        else if (!std::strcmp(argv[i], "--noRestore")) a.doRestore = false;
        else if (!std::strcmp(argv[i], "--noTier2"))   a.doTier2 = false;
        else if (!std::strcmp(argv[i], "--noMech"))    a.doMech = false;
    }
    return a;
}

// ------------------------------------------------------------------ tier-1 ladder
struct Ladder {
    const PreparedSystem::Impl* S = nullptr;
    VecX Fff, u0ff;
    MatX W;          // nf x R
    VecX s;          // R   (signed: -lambda removal, +lambda restore)
    MatX Z;          // nf x R = K0^-1 W
    MatX Mc;         // R x R = W^T Z
    int  R = 0;

    void init(const PreparedSystem::Impl* impl, const VecX& Fglobal) {
        S = impl;
        Fff = reduceVec(Fglobal, S->fmap, S->nf);
        u0ff = S->ldlt.solve(Fff);
        W.resize(S->nf, 0); Z.resize(S->nf, 0); s.resize(0); Mc.resize(0, 0); R = 0;
    }

    // sign = -1 remove, +1 restore. Returns rank added (0 = fully-fixed element), -1 = error.
    int addMember(const FrameModel& m, int memberIdx, int sign, std::string& why) {
        std::vector<Triplet> trips;
        if (!memberGlobalK(m, memberIdx, SolveOptions{}, trips, why)) return -1;
        int dofs[12];
        memberGlobalDofs(m, memberIdx, dofs);
        Mat12 Ke = Mat12::Zero();
        auto lof = [&](int g) -> int { for (int k = 0; k < 12; ++k) if (dofs[k] == g) return k; return -1; };
        for (const auto& t : trips) {
            const int lr = lof(static_cast<int>(t.row())), lc = lof(static_cast<int>(t.col()));
            if (lr >= 0 && lc >= 0) Ke(lr, lc) += t.value();
        }
        std::vector<int> keep, rid;
        for (int k = 0; k < 12; ++k) {
            const int f = S->fmap[static_cast<size_t>(dofs[k])];
            if (f >= 0) { keep.push_back(k); rid.push_back(f); }
        }
        const int nred = static_cast<int>(keep.size());
        if (nred == 0) return 0;
        MatX Ks(nred, nred);
        for (int a = 0; a < nred; ++a)
            for (int b = 0; b < nred; ++b) Ks(a, b) = Ke(keep[static_cast<size_t>(a)], keep[static_cast<size_t>(b)]);
        Eigen::SelfAdjointEigenSolver<MatX> es(Ks);
        const VecX lam = es.eigenvalues();
        const real lmax = lam.cwiseAbs().maxCoeff();
        std::vector<int> kept;
        for (int i = 0; i < nred; ++i)
            if (lam(i) > 1e-9 * lmax) kept.push_back(i);   // element K is PSD: keep the positive modes
        const int r = static_cast<int>(kept.size());
        if (r == 0) return 0;

        MatX Wn = MatX::Zero(S->nf, r);
        VecX sn(r);
        for (int q = 0; q < r; ++q) {
            for (int a = 0; a < nred; ++a) Wn(rid[static_cast<size_t>(a)], q) = es.eigenvectors()(a, kept[static_cast<size_t>(q)]);
            sn(q) = sign * lam(kept[static_cast<size_t>(q)]);
        }
        MatX Zn = S->ldlt.solve(Wn);

        MatX McNew(R + r, R + r);
        if (R > 0) {
            const MatX B1 = Wn.transpose() * Z;       // r x R   (= (W_old^T Zn)^T by symmetry of K0^-1)
            McNew.topLeftCorner(R, R)     = Mc;
            McNew.topRightCorner(R, r)    = B1.transpose();
            McNew.bottomLeftCorner(r, R)  = B1;
        }
        McNew.bottomRightCorner(r, r) = Wn.transpose() * Zn;

        W.conservativeResize(S->nf, R + r);  W.rightCols(r) = Wn;
        Z.conservativeResize(S->nf, R + r);  Z.rightCols(r) = Zn;
        s.conservativeResize(R + r);         s.tail(r) = sn;
        Mc = McNew;
        R += r;
        return r;
    }

    VecX solve(bool& mech, real& pivRatio) const {
        if (R == 0) { mech = false; pivRatio = 1; return u0ff; }
        MatX C = Mc;
        for (int i = 0; i < R; ++i) C(i, i) += real(1) / s(i);
        Eigen::FullPivLU<MatX> lu(C);
        const VecX d = lu.matrixLU().diagonal().cwiseAbs();
        const real dmax = d.maxCoeff(), dmin = d.minCoeff();
        pivRatio = (dmax > 0) ? dmin / dmax : 0;
        mech = (pivRatio < 1e-10);
        const VecX t = W.transpose() * u0ff;
        const VecX q = lu.solve(t);
        return u0ff - Z * q;
    }
};

// ------------------------------------------------------------------ tier-2 stale preconditioner
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

// Sparse W / diag(s) for a batch of removals (for assembling K' = Kff + W S W^T).
void batchUpdateSparse(const FrameModel& m, const PreparedSystem::Impl& S,
                       const std::vector<int>& memberIdxs, SpMat& Wsp, VecX& sv) {
    std::vector<Triplet> wtr;
    std::vector<real> svec;
    int col = 0;
    std::string why;
    for (int e : memberIdxs) {
        std::vector<Triplet> trips;
        if (!memberGlobalK(m, e, SolveOptions{}, trips, why)) continue;
        int dofs[12];
        memberGlobalDofs(m, e, dofs);
        Mat12 Ke = Mat12::Zero();
        auto lof = [&](int g) -> int { for (int k = 0; k < 12; ++k) if (dofs[k] == g) return k; return -1; };
        for (const auto& t : trips) {
            const int lr = lof(static_cast<int>(t.row())), lc = lof(static_cast<int>(t.col()));
            if (lr >= 0 && lc >= 0) Ke(lr, lc) += t.value();
        }
        std::vector<int> keep, rid;
        for (int k = 0; k < 12; ++k) {
            const int f = S.fmap[static_cast<size_t>(dofs[k])];
            if (f >= 0) { keep.push_back(k); rid.push_back(f); }
        }
        const int nred = static_cast<int>(keep.size());
        if (nred == 0) continue;
        MatX Ks(nred, nred);
        for (int a = 0; a < nred; ++a)
            for (int b = 0; b < nred; ++b) Ks(a, b) = Ke(keep[static_cast<size_t>(a)], keep[static_cast<size_t>(b)]);
        Eigen::SelfAdjointEigenSolver<MatX> es(Ks);
        const VecX lam = es.eigenvalues();
        const real lmax = lam.cwiseAbs().maxCoeff();
        for (int i = 0; i < nred; ++i) {
            if (lam(i) <= 1e-9 * lmax) continue;
            for (int a = 0; a < nred; ++a) {
                const real v = es.eigenvectors()(a, i);
                if (v != 0) wtr.emplace_back(rid[static_cast<size_t>(a)], col, v);
            }
            svec.push_back(-lam(i));   // removal
            ++col;
        }
    }
    Wsp.resize(S.nf, col);
    Wsp.setFromTriplets(wtr.begin(), wtr.end());
    Wsp.makeCompressed();
    sv.resize(col);
    for (int i = 0; i < col; ++i) sv(i) = svec[static_cast<size_t>(i)];
}

// ------------------------------------------------------------------ mechanism probe
void mechProbe() {
    {   // case A: 2-element cantilever chain, remove BASE element -> floating top = mechanism
        FrameModel chain = makeColumn(2, 2000.0, 100.0, 0.0, 1000.0);
        assertNodalOnly(chain, "mechProbe");
        PreparedSystem ps = assembleAndFactor(chain);
        Ladder L;
        L.init(ps.impl.get(), nodalLoadVector(chain));
        std::string why;
        L.addMember(chain, 0, -1, why);
        bool mech = false; real piv = 0;
        (void)L.solve(mech, piv);
        FrameModel work = chain;
        work.members[0].active = false;
        PreparedSystem ps2 = assembleAndFactor(work);
        std::printf("[mech] case=chain ladder_mech=%d piv=%.3e fresh_singular=%d\n",
                    mech ? 1 : 0, static_cast<double>(piv), ps2.impl->singular ? 1 : 0);
    }
    {   // case B: portal frame, remove the BEAM -> two cantilever columns = stable
        FrameModel p;
        p.materials.emplace_back(200000.0, 76923.07692307692, 7850.0);
        p.sections.push_back(squareSection(300.0));
        Node n0(0, 0, 0, 0);      n0.fixAll();
        Node n1(1, 6000, 0, 0);   n1.fixAll();
        Node n2(2, 0, 0, 3000);
        Node n3(3, 6000, 0, 3000);
        p.nodes = { n0, n1, n2, n3 };
        Member c0(0, 0, 2, 0, 0); c0.refVec = Vec3(1, 0, 0);
        Member c1(1, 1, 3, 0, 0); c1.refVec = Vec3(1, 0, 0);
        Member bm(2, 2, 3, 0, 0); bm.refVec = Vec3(0, 0, 1);
        p.members = { c0, c1, bm };
        NodalLoad l; l.node = 2; l.comp[Ux] = 5.0e4;
        p.nodalLoads.push_back(l);

        PreparedSystem ps = assembleAndFactor(p);
        Ladder L;
        L.init(ps.impl.get(), nodalLoadVector(p));
        std::string why;
        L.addMember(p, 2, -1, why);
        bool mech = false; real piv = 0;
        const VecX uf = L.solve(mech, piv);
        FrameModel work = p;
        work.members[2].active = false;
        PreparedSystem ps2 = assembleAndFactor(work);
        SolveResult r2 = solveLoad(ps2, work);
        const double relErr = relMaxDiff(scatterVec(uf, ps.impl->fmap, ps.impl->N), r2.u);
        std::printf("[mech] case=portal ladder_mech=%d piv=%.3e fresh_singular=%d relErr=%.3e\n",
                    mech ? 1 : 0, static_cast<double>(piv), ps2.impl->singular ? 1 : 0, relErr);
    }
}

} // namespace

int main(int argc, char** argv) {
    const Args a = parseArgs(argc, argv);

    int nx = 12, ny = 9, st = 24;
    if (a.preset == "small")     { nx = 3;  ny = 2;  st = 4; }
    else if (a.preset == "mega") { nx = 18; ny = 14; st = 36; }

    FrameModel model = makeTower(nx, ny, st);
    assertNodalOnly(model, "exp_incremental_refactor");

    Timer tf;
    PreparedSystem ps = assembleAndFactor(model);
    const double factorMs = tf.ms();
    const PreparedSystem::Impl& S = *ps.impl;
    if (S.singular) { std::printf("FATAL baseline singular: %s\n", S.diagnostic.c_str()); return 2; }

    const VecX Fglobal = nodalLoadVector(model);

    Ladder L;
    L.init(&S, Fglobal);

    // engine sanity + pure back-substitution timing
    Timer ts0;
    SolveResult r0 = solveLoad(ps, model);
    const double engineSolveMs = ts0.ms();
    Timer tbs;
    const VecX ubs = S.ldlt.solve(L.Fff);
    const double backsubMs = tbs.ms();
    const double sanity = relMaxDiff(scatterVec(ubs, S.fmap, S.N), r0.u);
    std::printf("[baseline] preset=%s nodes=%zu members=%zu nf=%d factorMs=%.1f engineSolveMs=%.1f backsubMs=%.1f sanityRel=%.3e\n",
                a.preset.c_str(), model.nodes.size(), model.members.size(), S.nf,
                factorMs, engineSolveMs, backsubMs, sanity);

    // removal candidates: braces (secIdx == 2), evenly spread
    std::vector<int> braces;
    for (size_t e = 0; e < model.members.size(); ++e)
        if (model.members[e].secIdx == 2) braces.push_back(static_cast<int>(e));
    const int steps = std::min<int>(a.steps, static_cast<int>(braces.size()));
    std::vector<int> cand;
    for (int i = 0; i < steps; ++i)
        cand.push_back(braces[static_cast<size_t>(i) * braces.size() / static_cast<size_t>(steps)]);

    // ---------------- tier-1 sequential removal ----------------
    FrameModel work = model;
    std::string why;
    for (int k = 0; k < steps; ++k) {
        const int e = cand[static_cast<size_t>(k)];
        Timer tu;
        const int radd = L.addMember(model, e, -1, why);
        bool mech = false; real piv = 0;
        const VecX uf = L.solve(mech, piv);
        const double tLadder = tu.ms();
        if (radd < 0) { std::printf("FATAL addMember: %s\n", why.c_str()); return 2; }

        work.members[static_cast<size_t>(e)].active = false;
        double tFresh = -1, relErr = -1;
        int freshSingular = 0;
        if ((k % a.refEvery) == 0 || k == steps - 1) {
            Timer tr;
            PreparedSystem ps2 = assembleAndFactor(work);
            SolveResult r2 = solveLoad(ps2, work);
            tFresh = tr.ms();
            freshSingular = r2.singular ? 1 : 0;
            relErr = relMaxDiff(scatterVec(uf, S.fmap, S.N), r2.u);
        }
        std::printf("[tier1] step=%d member=%d rankAdd=%d R=%d ladderMs=%.2f freshMs=%.1f speedup=%.1f relErr=%.3e mech=%d piv=%.2e freshSing=%d\n",
                    k + 1, model.members[static_cast<size_t>(e)].id, radd, L.R, tLadder, tFresh,
                    (tFresh > 0 && tLadder > 0) ? tFresh / tLadder : -1.0, relErr, mech ? 1 : 0,
                    static_cast<double>(piv), freshSingular);
    }

    // ---------------- restore drift ----------------
    if (a.doRestore && a.preset != "mega") {
        for (int k = steps - 1; k >= 0; --k)
            L.addMember(model, cand[static_cast<size_t>(k)], +1, why);
        bool mech = false; real piv = 0;
        Timer tsv;
        const VecX uf = L.solve(mech, piv);
        const double tSolve = tsv.ms();
        const double drift = relMaxDiff(scatterVec(uf, S.fmap, S.N), r0.u);
        std::printf("[restore] restored=%d R=%d solveMs=%.2f relErrVsBaseline=%.3e mech=%d\n",
                    steps, L.R, tSolve, drift, mech ? 1 : 0);
    }

    // ---------------- tier-2 batches ----------------
    if (a.doTier2) {
        const SpMat Kff = reduceFF(S.K, S.fmap, S.nf);
        for (int nb : { 16, 40, 80, 160 }) {
            if (nb > static_cast<int>(braces.size())) break;
            std::vector<int> batch;
            for (int i = 0; i < nb; ++i)
                batch.push_back(braces[static_cast<size_t>(i) * braces.size() / static_cast<size_t>(nb)]);

            // fresh reference
            FrameModel wk = model;
            for (int e : batch) wk.members[static_cast<size_t>(e)].active = false;
            Timer trf;
            PreparedSystem ps2 = assembleAndFactor(wk);
            SolveResult r2 = solveLoad(ps2, wk);
            const double tFresh = trf.ms();

            // K' = Kff + W S W^T (explicit assembly, timed)
            Timer tas;
            SpMat Wsp; VecX sv;
            batchUpdateSparse(model, S, batch, Wsp, sv);
            SpMat D(static_cast<int>(sv.size()), static_cast<int>(sv.size()));
            {
                std::vector<Triplet> dt;
                for (int i = 0; i < sv.size(); ++i) dt.emplace_back(i, i, sv(i));
                D.setFromTriplets(dt.begin(), dt.end());
            }
            SpMat Kp = Kff + SpMat(Wsp * D * Wsp.transpose());
            Kp.makeCompressed();
            const double tAssemble = tas.ms();

            // stale-LDLT PCG
            Timer tpc;
            Eigen::ConjugateGradient<SpMat, Eigen::Lower | Eigen::Upper, StalePrecond> cg;
            cg.preconditioner().ldlt = &S.ldlt;
            cg.setTolerance(1e-10);
            cg.setMaxIterations(500);
            cg.compute(Kp);
            const VecX x = cg.solveWithGuess(L.Fff, L.u0ff);
            const double tPCG = tpc.ms();
            const double relErr = relMaxDiff(scatterVec(x, S.fmap, S.N), r2.u);

            // Jacobi-PCG contrast
            Timer tj;
            Eigen::ConjugateGradient<SpMat, Eigen::Lower | Eigen::Upper> cgj;
            cgj.setTolerance(1e-10);
            cgj.setMaxIterations(2000);
            cgj.compute(Kp);
            const VecX xj = cgj.solveWithGuess(L.Fff, L.u0ff);
            const double tJac = tj.ms();

            std::printf("[tier2] nb=%d rank=%d staleIters=%d staleMs=%.1f assembleMs=%.1f totalMs=%.1f freshMs=%.1f speedup=%.1f relErr=%.3e cgErr=%.2e | jacobiIters=%d jacobiMs=%.1f jacobiErr=%.2e\n",
                        nb, static_cast<int>(sv.size()), static_cast<int>(cg.iterations()), tPCG, tAssemble,
                        tPCG + tAssemble, tFresh, tFresh / (tPCG + tAssemble), relErr,
                        static_cast<double>(cg.error()), static_cast<int>(cgj.iterations()), tJac,
                        static_cast<double>(cgj.error()));
        }
    }

    // ---------------- mechanism probe ----------------
    if (a.doMech) mechProbe();

    std::printf("[done]\n");
    return 0;
}
