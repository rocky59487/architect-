// exp_dynamic_inherit.cpp — N4 "momentum-inheriting dynamic collapse" experiment (scratch, NOT engine).
//
// Part 1  EQUIVALENCE: small tower under a suddenly-applied (step) load, one brace removed
//         mid-vibration. Reference = full-system Newmark (avg acceleration, C=0) with the
//         state (u,v) carried across the event. Candidate = modal superposition with the
//         SAME integrator + M'-orthonormal projection of (u,v) onto the post-event modes.
//         With the FULL modal basis the two are mathematically identical -> machine-level
//         agreement expected. Truncated bases (m modes) quantify the inheritance error,
//         with and without mode-acceleration static correction.
//
// Part 2  MOMENTUM IDENTITY: a chain whose upper part detaches. Fragment linear/angular
//         momentum extracted from the FE consistent mass (rigid-mode projection) is compared
//         against the engine's FragmentCluster closed-form mass/inertia (analyzeConnectivity)
//         under imposed rigid test motions. p is exact; L differs only by section-level
//         rotary inertia (slender-rod model) — the data that decides the S2 handoff spec.
//
// Part 3  WARM-START: subspace iteration (research copy of SparseEigsolver.h with an X0
//         parameter) seeded with the PRE-event modes vs the deterministic random start.

#include "research_common.h"
#include "FrameCore/Connectivity.h"
#include <random>
#include <cstring>

using namespace research;

namespace {

// ---------------------------------------------------------------- Newmark (beta=1/4, gamma=1/2, C=0)
struct Newmark {
    real dt = 0, a0 = 0, a2 = 0;
    LDLTSolver ldltHat;
    const SpMat* K = nullptr;
    const SpMat* M = nullptr;
    bool ok = false;

    void setup(const SpMat& Kff, const SpMat& Mff, real dt_) {
        K = &Kff; M = &Mff; dt = dt_;
        a0 = 4.0 / (dt * dt); a2 = 4.0 / dt;
        SpMat Khat = Kff + a0 * Mff;
        ldltHat.compute(Khat);
        ok = (ldltHat.info() == Eigen::Success);
    }
    // one step in place
    void step(const VecX& F, VecX& u, VecX& v, VecX& a) const {
        const VecX Fhat = F + (*M) * (a0 * u + a2 * v + a);          // a3 = 1
        const VecX un = ldltHat.solve(Fhat);
        const VecX an = a0 * (un - u) - a2 * v - a;
        v = v + (dt / 2.0) * (a + an);
        u = un; a = an;
    }
};

real energy(const SpMat& K, const SpMat& M, const VecX& F, const VecX& u, const VecX& v) {
    return 0.5 * v.dot(M * v) + 0.5 * u.dot(K * u) - F.dot(u);
}

// ---------------------------------------------------------------- modal pieces
struct ModalBasis {
    VecX w2;    // eigenvalues omega^2 (ascending)
    MatX Phi;   // M-orthonormal columns
};

ModalBasis denseModes(const SpMat& Kff, const SpMat& Mff) {
    ModalBasis B;
    const MatX Kd = MatX(Kff), Md = MatX(Mff);
    Eigen::GeneralizedSelfAdjointEigenSolver<MatX> ges(Kd, Md);   // K phi = w2 M phi, phi^T M phi = I
    B.w2 = ges.eigenvalues();
    B.Phi = ges.eigenvectors();
    return B;
}

// per-mode Newmark on  q'' + w2 q = f  (modal mass 1), same constants as the physical path
struct ModalNewmark {
    real dt = 0, a0 = 0, a2 = 0;
    VecX w2, khat;
    void setup(const VecX& w2_, real dt_) {
        dt = dt_; a0 = 4.0 / (dt * dt); a2 = 4.0 / dt;
        w2 = w2_; khat = w2.array() + a0;
    }
    void step(const VecX& f, VecX& q, VecX& qd, VecX& qdd) const {
        const VecX qn = (f + a0 * q + a2 * qd + qdd).cwiseQuotient(khat);
        const VecX qddn = a0 * (qn - q) - a2 * qd - qdd;
        qd = qd + (dt / 2.0) * (qdd + qddn);
        q = qn; qdd = qddn;
    }
};

// ---------------------------------------------------------------- research copy of SparseEigsolver.h
// (subspaceSmallest with an explicit start block X0 and an iteration-count output; the
//  engine header keeps its deterministic random start — copy attributed, engine untouched.)
bool subspaceSmallestX0(const SpMat& A, const LDLTSolver& Ainv, const SpMat& B,
                        int nev, VecX& lambda, MatX& vec, int& itersOut,
                        const MatX* X0 = nullptr, int maxIter = 300, real tol = 1e-11) {
    const int n = static_cast<int>(A.rows());
    if (nev <= 0 || n <= 0 || nev > n) return false;
    const int p = std::min(n, std::max(nev + 4, 2 * nev));

    std::mt19937 rng(20260607u);
    std::uniform_real_distribution<real> dist(real(-1), real(1));
    MatX X(n, p);
    for (int j = 0; j < p; ++j)
        for (int i = 0; i < n; ++i) X(i, j) = dist(rng);
    if (X0) {
        const int c = std::min<int>(p, static_cast<int>(X0->cols()));
        X.leftCols(c) = X0->leftCols(c);
    }

    VecX muPrev = VecX::Constant(p, real(1e30));
    VecX mu = VecX::Zero(p);
    MatX Xconv = X;
    bool converged = false;
    int iter = 0;
    for (; iter < maxIter; ++iter) {
        const MatX Y = B * X;
        const MatX Xbar = Ainv.solve(Y);
        if (Ainv.info() != Eigen::Success) return false;
        const MatX Ar = Xbar.transpose() * (A * Xbar);
        const MatX Br = Xbar.transpose() * (B * Xbar);
        Eigen::GeneralizedSelfAdjointEigenSolver<MatX> es(Br, Ar);
        if (es.info() != Eigen::Success) return false;
        mu = es.eigenvalues();
        X = Xbar * es.eigenvectors();
        const MatX AX = A * X;
        for (int j = 0; j < p; ++j) {
            const real q = X.col(j).dot(AX.col(j));
            const real nrm = (q > 0) ? std::sqrt(q) : real(0);
            if (nrm > 0) X.col(j) /= nrm;
        }
        Xconv = X;
        real maxRel = 0;
        for (int k = 0; k < nev; ++k) {
            const int idx = p - 1 - k;
            const real d = std::abs(mu(idx) - muPrev(idx)) / std::max<real>(1e-30, std::abs(mu(idx)));
            maxRel = std::max(maxRel, d);
        }
        if (maxRel < tol) { converged = true; ++iter; break; }
        muPrev = mu;
    }
    itersOut = iter;
    if (!converged) return false;
    lambda.resize(nev);
    vec.resize(n, nev);
    for (int k = 0; k < nev; ++k) {
        const int idx = p - 1 - k;
        const real m = mu(idx);
        lambda(k) = (m > 0) ? real(1) / m : real(0);
        vec.col(k) = Xconv.col(idx);
    }
    return true;
}

SpMat massFF(const PreparedSystem::Impl& S) {
    std::vector<Triplet> mt;
    for (const auto& el : S.elems) el->assembleMass(mt);
    SpMat M(S.N, S.N);
    M.setFromTriplets(mt.begin(), mt.end());
    return reduceFF(M, S.fmap, S.nf);
}

// ---------------------------------------------------------------- part 1
void partEquivalence() {
    FrameModel model = makeTower(2, 1, 3);
    assertNodalOnly(model, "dyn-equiv");
    PreparedSystem ps1 = assembleAndFactor(model);
    const PreparedSystem::Impl& S1 = *ps1.impl;
    const SpMat K1 = reduceFF(S1.K, S1.fmap, S1.nf);
    const SpMat M1 = massFF(S1);
    const VecX  F  = reduceVec(nodalLoadVector(model), S1.fmap, S1.nf);

    // pick a mid brace (secIdx == 2), removal keeps the moment frame connected
    int braceIdx = -1;
    for (size_t e = 0; e < model.members.size(); ++e)
        if (model.members[e].secIdx == 2) { braceIdx = static_cast<int>(e); }
    FrameModel work = model;
    work.members[static_cast<size_t>(braceIdx)].active = false;
    PreparedSystem ps2 = assembleAndFactor(work);
    const PreparedSystem::Impl& S2 = *ps2.impl;
    const SpMat K2 = reduceFF(S2.K, S2.fmap, S2.nf);
    const SpMat M2 = massFF(S2);

    const ModalBasis B1 = denseModes(K1, M1);
    const ModalBasis B2 = denseModes(K2, M2);
    const real w1 = std::sqrt(B1.w2(0));
    const real T1 = 2.0 * 3.14159265358979323846 / w1;
    const real dt = T1 / 200.0;
    const int  nSteps = 600, eventStep = 240;

    // ---- path A: physical Newmark
    Newmark nm1; nm1.setup(K1, M1, dt);
    Newmark nm2; nm2.setup(K2, M2, dt);
    if (!nm1.ok || !nm2.ok) { std::printf("[equiv] FATAL Newmark setup\n"); return; }
    LDLTSolver Mldlt1; Mldlt1.compute(M1);
    LDLTSolver Mldlt2; Mldlt2.compute(M2);

    VecX u = VecX::Zero(S1.nf), v = VecX::Zero(S1.nf);
    VecX a = Mldlt1.solve(F);
    std::vector<VecX> traceA;
    traceA.reserve(static_cast<size_t>(nSteps));
    real Epre = 0, Epost = 0;
    for (int s = 0; s < nSteps; ++s) {
        if (s == eventStep) {
            Epre = energy(K1, M1, F, u, v);
            a = Mldlt2.solve(F - K2 * u);            // state (u,v) carried; a recomputed
            Epost = energy(K2, M2, F, u, v);
        }
        (s < eventStep ? nm1 : nm2).step(F, u, v, a);
        traceA.push_back(u);
    }

    // ---- path B: modal superposition (basis size m), optional mode-acceleration correction
    auto runModal = [&](int m, bool modeAccel) -> double {
        const int m1 = std::min<int>(m, static_cast<int>(B1.w2.size()));
        const MatX P1 = B1.Phi.leftCols(m1);  const VecX W1 = B1.w2.head(m1);
        const MatX P2 = B2.Phi.leftCols(m1);  const VecX W2 = B2.w2.head(m1);
        ModalNewmark mn1; mn1.setup(W1, dt);
        ModalNewmark mn2; mn2.setup(W2, dt);
        const VecX f1 = P1.transpose() * F, f2 = P2.transpose() * F;

        // static correction vectors (constant F): us - Phi (f ./ w2)
        LDLTSolver k1l, k2l;
        VecX sc1 = VecX::Zero(S1.nf), sc2 = VecX::Zero(S1.nf);
        if (modeAccel) {
            k1l.compute(K1); k2l.compute(K2);
            sc1 = k1l.solve(F) - P1 * f1.cwiseQuotient(W1);
            sc2 = k2l.solve(F) - P2 * f2.cwiseQuotient(W2);
        }

        VecX q = VecX::Zero(m1), qd = VecX::Zero(m1), qdd = f1;
        double maxErr = 0, maxRef = 0;
        for (int s = 0; s < nSteps; ++s) {
            if (s == eventStep) {
                const VecX ue = P1 * q + (modeAccel ? sc1 : VecX::Zero(S1.nf));
                const VecX ve = P1 * qd;
                q  = P2.transpose() * (M2 * ue);     // M'-orthonormal projection
                qd = P2.transpose() * (M2 * ve);
                qdd = f2 - W2.cwiseProduct(q);
            }
            (s < eventStep ? mn1 : mn2).step(s < eventStep ? f1 : f2, q, qd, qdd);
            const VecX ub = (s < eventStep ? P1 : P2) * q + (modeAccel ? (s < eventStep ? sc1 : sc2) : VecX::Zero(S1.nf));
            const VecX d = ub - traceA[static_cast<size_t>(s)];
            maxErr = std::max(maxErr, static_cast<double>(d.cwiseAbs().maxCoeff()));
            maxRef = std::max(maxRef, static_cast<double>(traceA[static_cast<size_t>(s)].cwiseAbs().maxCoeff()));
        }
        return maxErr / std::max(1e-300, maxRef);
    };

    std::printf("[equiv] nf=%d brace=%d dt=%.4e steps=%d eventStep=%d fullBasisRelErr=%.3e Epre=%.6g Epost=%.6g dE=%.4g\n",
                S1.nf, braceIdx, static_cast<double>(dt), nSteps, eventStep,
                runModal(S1.nf, false), static_cast<double>(Epre), static_cast<double>(Epost),
                static_cast<double>(Epost - Epre));
    for (int m : { 5, 10, 20, 40 })
        std::printf("[trunc] m=%d errPlain=%.3e errModeAccel=%.3e\n", m, runModal(m, false), runModal(m, true));
}

// ---------------------------------------------------------------- part 2
void partMomentum() {
    FrameModel chain;
    chain.materials.emplace_back(200000.0, 76923.07692307692, 7850.0);
    chain.sections.push_back(squareSection(150.0));
    for (int k = 0; k <= 3; ++k) {
        Node n(k, 0, 0, 1000.0 * k);
        if (k == 0) n.fixAll();
        chain.nodes.push_back(n);
    }
    for (int k = 0; k < 3; ++k) {
        Member m(k, k, k + 1, 0, 0);
        m.refVec = Vec3(1, 0, 0);
        chain.members.push_back(m);
    }
    NodalLoad l; l.node = 3; l.comp[Ux] = 1.0;   // load content irrelevant here
    chain.nodalLoads.push_back(l);

    // engine closed-form fragment (rod model): remove member 1 -> nodes {2,3} detach
    FrameModel work = chain;
    work.members[1].active = false;
    ConnectivityResult cr = analyzeConnectivity(work);
    if (!cr.valid || cr.detached.empty()) { std::printf("[momentum] FATAL no fragment\n"); return; }
    const FragmentCluster& fc = cr.detached[0];

    // FE consistent mass of the fragment's single member (2-3), 12x12 on nodes {2,3}
    BeamColumnElement el(2);
    std::string why;
    if (!el.prepare(chain, SolveOptions{}, why)) { std::printf("[momentum] FATAL %s\n", why.c_str()); return; }
    std::vector<Triplet> mt;
    el.assembleMass(mt);
    int dofs[12];
    memberGlobalDofs(chain, 2, dofs);
    auto lof = [&](int g) -> int { for (int k = 0; k < 12; ++k) if (dofs[k] == g) return k; return -1; };
    MatX M12 = MatX::Zero(12, 12);
    for (const auto& t : mt) {
        const int lr = lof(static_cast<int>(t.row())), lc = lof(static_cast<int>(t.col()));
        if (lr >= 0 && lc >= 0) M12(lr, lc) += t.value();
    }

    const Vec3 rn[2] = { chain.nodes[2].pos, chain.nodes[3].pos };
    auto trans = [&](int axis) { VecX T = VecX::Zero(12); T(axis) = 1; T(6 + axis) = 1; return T; };
    auto rot = [&](const Vec3& ek) {
        VecX T = VecX::Zero(12);
        for (int n = 0; n < 2; ++n) {
            const Vec3 rel = rn[n] - fc.com;
            const Vec3 tr = cross(ek, rel);
            T(6 * n + 0) = tr.x; T(6 * n + 1) = tr.y; T(6 * n + 2) = tr.z;
            T(6 * n + 3) = ek.x; T(6 * n + 4) = ek.y; T(6 * n + 5) = ek.z;
        }
        return T;
    };

    const VecX Tx = trans(0), Ty = trans(1), Tz = trans(2);
    const VecX Rx = rot(Vec3(1, 0, 0)), Ry = rot(Vec3(0, 1, 0)), Rz = rot(Vec3(0, 0, 1));

    // total mass identity
    const real mFE = Tx.dot(M12 * Tx);
    std::printf("[momentum] mFE=%.12g mCluster=%.12g relDiff=%.3e com=(%.6g,%.6g,%.6g)\n",
                static_cast<double>(mFE), static_cast<double>(fc.mass),
                std::abs(static_cast<double>(mFE - fc.mass)) / static_cast<double>(fc.mass),
                static_cast<double>(fc.com.x), static_cast<double>(fc.com.y), static_cast<double>(fc.com.z));

    // rigid translation: p = m v0 exactly, L_com = 0 exactly
    {
        const real v0 = 7.3;
        const VecX v = v0 * Tx;
        const real px = Tx.dot(M12 * v), py = Ty.dot(M12 * v), pz = Tz.dot(M12 * v);
        const real Ly = Ry.dot(M12 * v);
        std::printf("[momentum-trans] px=%.12g expected=%.12g rel=%.3e py=%.3e pz=%.3e Lcom_y=%.3e (exact zeros)\n",
                    static_cast<double>(px), static_cast<double>(mFE * v0),
                    std::abs(static_cast<double>(px - mFE * v0)) / static_cast<double>(mFE * v0),
                    static_cast<double>(py), static_cast<double>(pz), static_cast<double>(Ly));
    }
    // rigid rotation about com, transverse axis y: p ~ 0, L_y = I_yy w0 (FE includes section
    // rotary inertia; cluster = slender rod) — the relative gap is the documented model diff
    {
        const real w0 = 0.31;
        const VecX v = w0 * Ry;
        const real px = Tx.dot(M12 * v);
        const real LyFE = Ry.dot(M12 * v);
        const real IyyCl = fc.inertia[1];
        std::printf("[momentum-rot] px=%.3e (exact zero) LyFE=%.10g IyyCluster*w=%.10g relDiff=%.3e\n",
                    static_cast<double>(px), static_cast<double>(LyFE),
                    static_cast<double>(IyyCl * w0),
                    std::abs(static_cast<double>(LyFE - IyyCl * w0)) / std::abs(static_cast<double>(IyyCl * w0)));
    }
    // own-axis rotation z: cluster rod model carries ZERO own-axis inertia; FE carries the
    // section polar term — documented edge case for the handoff spec
    {
        const real w0 = 0.31;
        const VecX v = w0 * Rz;
        const real LzFE = Rz.dot(M12 * v);
        std::printf("[momentum-ownaxis] LzFE=%.6g IzzCluster=%.6g (rod model: 0 own-axis) note=section-polar-term\n",
                    static_cast<double>(LzFE), static_cast<double>(fc.inertia[2] * w0));
    }
}

// ---------------------------------------------------------------- part 3
void partWarmStart() {
    FrameModel model = makeTower(3, 2, 4);
    PreparedSystem ps1 = assembleAndFactor(model);
    const PreparedSystem::Impl& S1 = *ps1.impl;
    const SpMat K1 = reduceFF(S1.K, S1.fmap, S1.nf);
    const SpMat M1 = massFF(S1);

    int braceIdx = -1;
    for (size_t e = 0; e < model.members.size(); ++e)
        if (model.members[e].secIdx == 2) { braceIdx = static_cast<int>(e); break; }
    FrameModel work = model;
    work.members[static_cast<size_t>(braceIdx)].active = false;
    PreparedSystem ps2 = assembleAndFactor(work);
    const PreparedSystem::Impl& S2 = *ps2.impl;
    const SpMat K2 = reduceFF(S2.K, S2.fmap, S2.nf);
    const SpMat M2 = massFF(S2);

    const int nev = 10;
    VecX lam1, lam2c, lam2w;
    MatX v1, v2c, v2w;
    int itBase = 0, itCold = 0, itWarm = 0;
    const bool ok1 = subspaceSmallestX0(K1, S1.ldlt, M1, nev, lam1, v1, itBase);
    const bool ok2 = subspaceSmallestX0(K2, S2.ldlt, M2, nev, lam2c, v2c, itCold);
    const bool ok3 = subspaceSmallestX0(K2, S2.ldlt, M2, nev, lam2w, v2w, itWarm, &v1);
    const double dLam = (ok2 && ok3)
        ? static_cast<double>((lam2w - lam2c).cwiseAbs().maxCoeff() / lam2c.cwiseAbs().maxCoeff()) : -1.0;
    std::printf("[warmstart] nf=%d nev=%d ok=%d%d%d itersBase=%d itersCold=%d itersWarm=%d sameLamRel=%.3e\n",
                S1.nf, nev, ok1 ? 1 : 0, ok2 ? 1 : 0, ok3 ? 1 : 0, itBase, itCold, itWarm, dLam);
}

} // namespace

int main() {
    partEquivalence();
    partMomentum();
    partWarmStart();
    std::printf("[done]\n");
    return 0;
}
