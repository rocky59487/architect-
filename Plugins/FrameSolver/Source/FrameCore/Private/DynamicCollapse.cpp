// S2 N4 dynamic-collapse driver: modal-space Newmark over a sequence of brittle removal events,
// with cross-event state inheritance and momentum-preserving debris handoff. Each event re-factors
// FRESH (assembleAndFactor) -- the connectivity cleanup pins debris nodes (a support change beyond
// the ReSolve same-topology regime), and rebuilding the modal/Ritz basis needs the new K'_ff factor.
// Mirrors the static runProgressiveCollapse working-copy / cleanup contract.
#include "FrameCore/DynamicCollapse.h"
#include "FrameCore/FrameSolver.h"
#include "FrameCore/ElasticAllowable.h"
#include "PreparedSystemImpl.h"
#include "IElement.h"
#include "BeamColumnElement.h"
#include "MITC4ShellElement.h"
#include "FragmentMomentum.h"
#include "FrameEigen.h"
#include "CollapseSupport.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <string>
#include <vector>

namespace frame {
namespace {

// kPi (FrameEigen.h), reduceFF (PreparedSystemImpl.h) and FrameModel::memberIndex/shellIndex
// are shared — this TU no longer keeps file-local copies.

// Global consistent mass of the prepared (active) elements, reduced to free DOFs.
SpMat massFF(const PreparedSystem::Impl& S) {
    std::vector<Triplet> mt;
    for (const auto& el : S.elems) el->assembleMass(mt);
    SpMat M(S.N, S.N); M.setFromTriplets(mt.begin(), mt.end());
    return reduceFF(M, S.fmap, S.nf);
}

// Reduced load vector Ff = F_f - K_fc u_c (nodal + active-element equivalent loads + prescribed
// support term), mirroring Reanalysis.cpp's F-increment assembly.
VecX reducedLoad(const FrameModel& work, const PreparedSystem::Impl& S) {
    const int N = S.N, nf = S.nf; const auto& fmap = S.fmap;
    VecX F = VecX::Zero(N);
    for (const auto& nl : work.nodalLoads) {
        const int ni = work.nodeIndex(nl.node);
        if (ni < 0) continue;
        for (int d = 0; d < 6; ++d) F(gdof(ni, d)) += nl.comp[d];
    }
    for (const auto& el : S.elems) el->addEquivalentNodalLoads(F);
    std::vector<real> presc((size_t)N, 0.0);
    for (size_t k = 0; k < work.nodes.size(); ++k)
        for (int d = 0; d < 6; ++d)
            if (work.nodes[k].fixed[d]) presc[(size_t)gdof((int)k, d)] = work.nodes[k].prescribed[d];
    VecX Ff = VecX::Zero(nf);
    for (int c = 0; c < N; ++c)
        for (SpMat::InnerIterator it(S.K, c); it; ++it) {
            const int rr = it.row();
            if (fmap[(size_t)rr] < 0) continue;
            if (fmap[(size_t)c] < 0 && presc[(size_t)c] != 0.0) Ff(fmap[(size_t)rr]) -= it.value() * presc[(size_t)c];
        }
    for (int g = 0; g < N; ++g) if (fmap[(size_t)g] >= 0) Ff(fmap[(size_t)g]) += F(g);
    return Ff;
}

// Scatter a free-DOF vector to the full global vector. usePrescribed=true fills constrained DOFs
// with their prescribed displacement (A3); false fills 0 (velocity of a prescribed support is 0).
VecX scatterToGlobal(const VecX& uf, const FrameModel& work, const std::vector<int>& fmap,
                     int N, bool usePrescribed) {
    VecX u = VecX::Zero(N);
    for (int g = 0; g < N; ++g) {
        if (fmap[(size_t)g] >= 0) { u(g) = uf(fmap[(size_t)g]); continue; }
        if (usePrescribed) {
            const int ni = g / DOF_PER_NODE, d = g % DOF_PER_NODE;
            if (work.nodes[(size_t)ni].fixed[d]) u(g) = work.nodes[(size_t)ni].prescribed[d];
        }
    }
    return u;
}

// Reduce a full global vector onto the (post-event) free space. New free DOFs are a subset of the
// old free DOFs (pinning only adds constraints), so pinned-fragment DOFs simply fall away.
VecX reduceToFree(const VecX& u_N, const std::vector<int>& fmap, int nf) {
    VecX uf = VecX::Zero(nf);
    for (int g = 0; g < (int)fmap.size(); ++g) if (fmap[(size_t)g] >= 0) uf(fmap[(size_t)g]) = u_N(g);
    return uf;
}

// pinNode / anyActive are shared with the static driver (CollapseSupport.h).

// ---------------------------------------------------------------- basis generation
// Pure eigenmodes (dense GES): K phi = w2 M phi, phi M-orthonormal, ascending w2. Small models.
void denseModes(const SpMat& Kff, const SpMat& Mff, int m, MatX& Phi, VecX& W2) {
    const MatX Kd = MatX(Kff), Md = MatX(Mff);
    Eigen::GeneralizedSelfAdjointEigenSolver<MatX> ges(Kd, Md);
    const VecX ev = ges.eigenvalues();
    const MatX V  = ges.eigenvectors();
    const int mm = std::min(m, (int)ev.size());
    W2  = ev.head(mm);
    Phi = V.leftCols(mm);
}

// Load-dependent Ritz vectors (Wilson 1985). Seeds with `g`, builds an
// M-orthonormal Krylov-like block via K^{-1}, then a Rayleigh-Ritz projection. Full basis (m=nf)
// spans the whole space -> identical to denseModes there (F37b gate). Degenerate vectors are
// replaced by a deterministic random restart (reproducible builds).
void ritzBasis(const SpMat& Kff, const LDLTSolver& ldlt, const SpMat& Mff, const VecX& g,
               int m, MatX& Phi, VecX& W2) {
    const int nf = (int)Kff.rows();
    const int mm = std::min(m, nf);
    MatX X(nf, mm);
    std::mt19937 rng(20260607u);
    std::uniform_real_distribution<real> dist(real(-1), real(1));
    auto mNorm = [&](const VecX& x) { VecX mx = Mff * x; return std::sqrt(std::max<real>(0, x.dot(mx))); };

    // x1 = K^{-1} g, M-normalized (random restart if g is in the null content of M^{1/2}K^{-1}).
    VecX x = ldlt.solve(g);
    real nrm = mNorm(x);
    if (!(nrm > real(1e-300))) { for (int i = 0; i < nf; ++i) x(i) = dist(rng); x = ldlt.solve(Mff * x); nrm = mNorm(x); }
    X.col(0) = x / std::max<real>(nrm, real(1e-300));

    for (int i = 1; i < mm; ++i) {
        VecX xi = ldlt.solve(Mff * X.col(i - 1));
        for (int pass = 0; pass < 2; ++pass)                       // two-pass Gram-Schmidt (M-inner)
            for (int j = 0; j < i; ++j) xi -= (X.col(j).dot(Mff * xi)) * X.col(j);
        real beta = mNorm(xi);
        for (int tries = 0; tries < 3 && !(beta > real(1e-8)); ++tries) {   // degenerate -> random restart
            for (int k = 0; k < nf; ++k) xi(k) = dist(rng);
            for (int pass = 0; pass < 2; ++pass)
                for (int j = 0; j < i; ++j) xi -= (X.col(j).dot(Mff * xi)) * X.col(j);
            beta = mNorm(xi);
        }
        X.col(i) = xi / std::max<real>(beta, real(1e-300));
    }

    const MatX Kr = X.transpose() * (Kff * X);
    const MatX Mr = X.transpose() * (Mff * X);                     // ~ I; solved as a GES for safety
    Eigen::GeneralizedSelfAdjointEigenSolver<MatX> ges(Kr, Mr);
    W2  = ges.eigenvalues();
    Phi = X * ges.eigenvectors();                                  // M-orthonormal, ascending W2
}

// ---------------------------------------------------------------- per-mode Newmark (avg accel)
// beta=1/4, gamma=1/2; modal damping c_i = alpha + beta_r * w2_i (= Phi^T(aM+bK)Phi diagonal,
// no 1/omega -> no overflow at a rigid/near-zero mode). c=0 reproduces the prototype bit-for-bit.
void modalNewmarkStep(VecX& q, VecX& qd, VecX& qdd, const VecX& f, const VecX& W2,
                      real dt, real ralpha, real rbeta) {
    const real beta = real(0.25), gamma = real(0.5);
    const real b1 = real(1) / (beta * dt * dt), b2 = real(1) / (beta * dt), b3 = real(1) / (2 * beta) - 1;
    const real g1 = gamma / (beta * dt), g2 = gamma / beta - 1, g3 = dt * (gamma / (2 * beta) - 1);
    for (int i = 0; i < (int)q.size(); ++i) {
        const real k = W2(i);
        const real c = ralpha + rbeta * k;
        const real keff = b1 + g1 * c + k;
        const real rhs  = f(i) + b1 * q(i) + b2 * qd(i) + b3 * qdd(i)
                        + c * (g1 * q(i) + g2 * qd(i) + g3 * qdd(i));
        const real qn   = rhs / keff;
        const real qddn = b1 * (qn - q(i)) - b2 * qd(i) - b3 * qdd(i);
        const real qdn  = qd(i) + dt * ((1 - gamma) * qdd(i) + gamma * qddn);
        q(i) = qn; qd(i) = qdn; qdd(i) = qddn;
    }
}

// Fragment momentum handoff (fillFragmentVelocity / fragmentKE) lives in FragmentMomentum.h so the
// deep audit can exercise the identical extraction against an imposed rigid-body velocity field.

// ---------------------------------------------------------------- configuration system
struct ConfigSystem {
    PreparedSystem ps;                 // holds the fresh factorization + active elements (for recover)
    int N = 0, nf = 0;
    std::vector<int> fmap;
    SpMat Kff, Mff;
    VecX  Fff, u0ff, f, W2;            // f = Phi^T Fff (modal force), W2 = modal omega^2
    MatX  Phi;
    bool ok = false;
    std::string diag;
};

// Fresh assembleAndFactor + reduce + static equilibrium + modal/Ritz basis. inheritU_N (optional)
// supplies the inherited global displacement so the Ritz seed is the post-event residual r = F'-K'u'
// (B1: falls back to F' when ||r|| is tiny, then random when ||F'|| is tiny too).
ConfigSystem buildConfig(const FrameModel& work, const DynCollapseOptions& opts, const VecX* inheritU_N) {
    ConfigSystem cfg;
    PreparedSystem ps = assembleAndFactor(work, opts.solve);
    const PreparedSystem::Impl& S = *ps.impl;
    if (S.singular) { cfg.diag = S.diagnostic.empty() ? "singular configuration" : S.diagnostic; cfg.ps = std::move(ps); return cfg; }

    cfg.N = S.N; cfg.nf = S.nf; cfg.fmap = S.fmap;
    cfg.Kff = reduceFF(S.K, S.fmap, S.nf);
    cfg.Mff = massFF(S);
    cfg.Fff = reducedLoad(work, S);
    cfg.u0ff = S.ldlt.solve(cfg.Fff);

    const int m = std::min(opts.basisSize, cfg.nf);
    if (opts.useRitzVectors) {
        VecX seed = cfg.Fff;
        if (inheritU_N) {
            const VecX uIn = reduceToFree(*inheritU_N, cfg.fmap, cfg.nf);
            const VecX r = cfg.Fff - cfg.Kff * uIn;
            if (r.norm() > real(1e-6) * std::max<real>(cfg.Fff.norm(), real(1e-300))) seed = r;
        }
        ritzBasis(cfg.Kff, S.ldlt, cfg.Mff, seed, m, cfg.Phi, cfg.W2);
    } else {
        denseModes(cfg.Kff, cfg.Mff, m, cfg.Phi, cfg.W2);
    }
    cfg.f = cfg.Phi.transpose() * cfg.Fff;
    cfg.ok = true;
    cfg.ps = std::move(ps);
    return cfg;
}

real configEnergy(const ConfigSystem& cfg, const VecX& uff, const VecX& vff) {
    return real(0.5) * vff.dot(cfg.Mff * vff) + real(0.5) * uff.dot(cfg.Kff * uff) - cfg.Fff.dot(uff);
}

// Recover member/shell forces of the current active configuration into a SolveResult (mirrors
// solveLoad's id-stamping; inactive rows keep id + zero force).
SolveResult recoverForces(const ConfigSystem& cfg, const FrameModel& work, const VecX& u_N) {
    SolveResult r;
    r.memberForces.resize(work.members.size());
    r.shellForces.resize(work.shells.size());
    for (size_t e = 0; e < work.members.size(); ++e) r.memberForces[e].member = work.members[e].id;
    for (size_t s = 0; s < work.shells.size(); ++s) r.shellForces[s].shell = work.shells[s].id;
    for (const auto& el : cfg.ps.impl->elems) el->recover(u_N, r);
    return r;
}

// Worst brittle screening event (member section / shell von Mises), deterministic tie-break:
// worst ratio, member-before-shell, then smallest id. (S2 has no ductile hinge path -- reserved.)
struct ScreenEvent { real ratio = 0; int kind = -1; int id = 0; FailMode mode = FailMode::None; bool eligible = false; };
ScreenEvent screenWorst(const FrameModel& work, const SolveResult& r, real threshold) {
    const ElasticAllowable screen;
    ScreenEvent best;
    auto consider = [&best](const ScreenEvent& c) {
        if (!c.eligible) return;
        if (!best.eligible || c.ratio > best.ratio ||
            (c.ratio == best.ratio && (c.kind < best.kind || (c.kind == best.kind && c.id < best.id))))
            best = c;
    };
    const size_t nM = std::min(work.members.size(), r.memberForces.size());
    for (size_t e = 0; e < nM; ++e) {
        const Member& mem = work.members[e];
        if (!mem.active) continue;
        if (mem.matIdx < 0 || mem.matIdx >= (int)work.materials.size()) continue;
        if (mem.secIdx < 0 || mem.secIdx >= (int)work.sections.size()) continue;
        const Section&  sec = work.sections[(size_t)mem.secIdx];
        const Material& mat = work.materials[(size_t)mem.matIdx];
        const MemberForcePair& mf = r.memberForces[e];
        const DemandResult di = screen.checkSection(mf.endI, sec, mat.cap);
        const DemandResult dj = screen.checkSection(mf.endJ, sec, mat.cap);
        const DemandResult& d = (di.risk >= dj.risk) ? di : dj;
        consider(ScreenEvent{ d.risk, 0, mem.id, d.mode, d.risk > threshold });
    }
    const size_t nS = std::min(work.shells.size(), r.shellForces.size());
    for (size_t s = 0; s < nS; ++s) {
        const ShellQuad& sh = work.shells[s];
        if (!sh.active) continue;
        if (sh.matIdx < 0 || sh.matIdx >= (int)work.materials.size()) continue;
        const ShellDemandResult d = checkShellSurface(r.shellForces[s], sh.t, work.materials[(size_t)sh.matIdx].cap);
        consider(ScreenEvent{ d.risk, 1, sh.id, FailMode::ShellVonMises, d.risk > threshold });
    }
    return best;
}

// Connectivity cleanup: identify detached fragments, extract their momentum (BEFORE pinning, while
// v_N is intact = A1/P0), then deactivate + pin + drop loads. v_N==nullptr -> zero handoff velocity.
void cleanupFragments(FrameModel& work, DynCollapseEvent& ev, const VecX* v_N) {
    const ConnectivityResult conn = analyzeConnectivity(work);
    if (!conn.valid) return;
    for (FragmentCluster fc : conn.detached) {
        if (v_N) fillFragmentVelocity(fc, work, *v_N);
        for (MemberId id : fc.members) work.members[(size_t)work.memberIndex(id)].active = false;
        for (int sid : fc.shells)      work.shells[(size_t)work.shellIndex(sid)].active = false;
        for (NodeId nid : fc.nodes)    pinNode(work, nid);
        ev.detached.push_back(fc);
    }
    for (NodeId nid : conn.looseNodes) pinNode(work, nid);
}

}  // namespace

DynCollapseHistory runDynamicCollapse(const FrameModel& model, const DynCollapseOptions& opts) {
    DynCollapseHistory H;
    if (!(opts.dt > 0) || !std::isfinite(opts.dt)) { H.diagnostic = "DynCollapseOptions.dt must be finite and > 0"; return H; }
    if (!(opts.maxTime > 0)) { H.diagnostic = "DynCollapseOptions.maxTime must be > 0"; return H; }
    if (opts.basisSize < 1) { H.diagnostic = "DynCollapseOptions.basisSize must be >= 1"; return H; }
    if (opts.screenEvery < 1 || opts.frameStride < 1) { H.diagnostic = "screenEvery / frameStride must be >= 1"; return H; }
    if (opts.maxEvents < 0) { H.diagnostic = "maxEvents must be >= 0"; return H; }
    std::string why;
    if (!model.validate(why)) { H.diagnostic = "invalid model: " + why; return H; }

    FrameModel work = model;   // the caller's model is never mutated
    for (MemberId id : opts.initialRemovals)
        if (work.memberIndex(id) < 0) { H.diagnostic = "initialRemovals references missing member id " + std::to_string(id); return H; }
    for (int sid : opts.initialShellRemovals)
        if (work.shellIndex(sid) < 0) { H.diagnostic = "initialShellRemovals references missing shell id " + std::to_string(sid); return H; }

    // ---- t=0 scenario: apply initial removals + fragment cleanup (zero handoff velocity, static start)
    DynCollapseEvent ev0; ev0.t = 0; ev0.mode = FailMode::None;
    for (MemberId id : opts.initialRemovals) { work.members[(size_t)work.memberIndex(id)].active = false; ev0.removedMembers.push_back(id); }
    for (int sid : opts.initialShellRemovals) { work.shells[(size_t)work.shellIndex(sid)].active = false; ev0.removedShells.push_back(sid); }
    cleanupFragments(work, ev0, nullptr);
    const bool hasInitialEvent = !ev0.removedMembers.empty() || !ev0.removedShells.empty() || !ev0.detached.empty();

    if (!anyActive(work)) {
        if (hasInitialEvent) H.events.push_back(ev0);
        H.outcome = CollapseOutcome::Collapsed; H.diagnostic = "no active element remains grounded"; return H;
    }
    ConfigSystem cfg = buildConfig(work, opts, nullptr);
    if (!cfg.ok) {
        if (hasInitialEvent) H.events.push_back(ev0);
        H.outcome = CollapseOutcome::Collapsed; H.diagnostic = "mechanism at t=0: " + cfg.diag; return H;
    }
    if (hasInitialEvent) {
        ev0.energyAfter = configEnergy(cfg, cfg.u0ff, VecX::Zero(cfg.nf));   // static elastic energy; fragments at rest
        ev0.energyBefore = ev0.energyAfter;
        H.events.push_back(ev0);
    }

    // ---- modal state from the static equilibrium (q from u0 projection, qd=0)
    VecX q   = cfg.Phi.transpose() * (cfg.Mff * cfg.u0ff);
    VecX qd  = VecX::Zero(cfg.Phi.cols());
    VecX qdd = cfg.f - cfg.W2.cwiseProduct(q);
    real T1 = (cfg.W2.size() > 0 && cfg.W2(0) > 0) ? real(2) * kPi / std::sqrt(cfg.W2(0)) : opts.maxTime;

    auto storeFrame = [&](real t, const VecX& uN, const VecX& vN) {
        DynCollapseFrame fr; fr.t = t;
        fr.u.assign(uN.data(), uN.data() + uN.size());
        fr.v.assign(vN.data(), vN.data() + vN.size());
        H.frames.push_back(std::move(fr));
    };
    storeFrame(0, scatterToGlobal(cfg.Phi * q, work, cfg.fmap, cfg.N, true), VecX::Zero(cfg.N));

    real t = 0; long step = 0; real maxKE = 0; real quietTime = 0; int events = 0;
    const real tinyKE = real(1e-300);

    while (t < opts.maxTime) {
        modalNewmarkStep(q, qd, qdd, cfg.f, cfg.W2, opts.dt, opts.rayleighAlpha, opts.rayleighBeta);
        t += opts.dt; ++step;
        const real ke = real(0.5) * qd.squaredNorm();
        maxKE = std::max(maxKE, ke);

        if (step % opts.frameStride == 0)
            storeFrame(t, scatterToGlobal(cfg.Phi * q, work, cfg.fmap, cfg.N, true),
                          scatterToGlobal(cfg.Phi * qd, work, cfg.fmap, cfg.N, false));

        if (step % opts.screenEvery != 0) continue;

        const VecX uff = cfg.Phi * q;
        const VecX uN  = scatterToGlobal(uff, work, cfg.fmap, cfg.N, true);
        const SolveResult r = recoverForces(cfg, work, uN);
        const ScreenEvent best = screenWorst(work, r, opts.removeThreshold);

        if (!best.eligible) {                                   // no element over threshold -> stability check
            if (events == 0) { H.outcome = CollapseOutcome::Stable; return H; }   // undisturbed static equilibrium
            if (maxKE <= tinyKE || ke < opts.quietKineticRatio * maxKE) {
                quietTime += opts.screenEvery * opts.dt;
                if (quietTime >= T1) { H.outcome = CollapseOutcome::Stable; return H; }
            } else quietTime = 0;
            continue;
        }

        // ---- EVENT (strict order: energy/freeze -> condemn -> momentum (pre-pin) -> cleanup -> re-factor -> inherit)
        DynCollapseEvent ev; ev.t = t; ev.mode = best.mode;
        const VecX vff = cfg.Phi * qd;
        ev.energyBefore = configEnergy(cfg, uff, vff);
        const VecX uN_f = uN;
        const VecX vN_f = scatterToGlobal(vff, work, cfg.fmap, cfg.N, false);

        if (best.kind == 0) { work.members[(size_t)work.memberIndex(best.id)].active = false; ev.removedMembers.push_back(best.id); }
        else                { work.shells[(size_t)work.shellIndex(best.id)].active = false;  ev.removedShells.push_back(best.id); }
        cleanupFragments(work, ev, &vN_f);

        real fragKE = 0; for (const FragmentCluster& fc : ev.detached) fragKE += fragmentKE(fc);

        if (!anyActive(work)) {
            ev.energyAfter = fragKE; H.events.push_back(ev);
            H.outcome = CollapseOutcome::Collapsed; H.diagnostic = "no active element remains grounded"; return H;
        }
        ConfigSystem cfg2 = buildConfig(work, opts, &uN_f);
        if (!cfg2.ok) {
            ev.energyAfter = fragKE; H.events.push_back(ev);
            H.outcome = CollapseOutcome::Collapsed; H.diagnostic = "mechanism in the grounded remainder: " + cfg2.diag; return H;
        }

        const VecX upff = reduceToFree(uN_f, cfg2.fmap, cfg2.nf);
        const VecX vpff = reduceToFree(vN_f, cfg2.fmap, cfg2.nf);
        VecX q2   = cfg2.Phi.transpose() * (cfg2.Mff * upff);
        VecX qd2  = cfg2.Phi.transpose() * (cfg2.Mff * vpff);
        VecX qdd2 = cfg2.f - cfg2.W2.cwiseProduct(q2);
        const VecX proj = cfg2.Phi * q2;
        ev.truncationResidual = (upff - proj).norm() / std::max<real>(upff.norm(), real(1e-300));
        ev.energyAfter = configEnergy(cfg2, upff, vpff) + fragKE;
        H.events.push_back(ev);

        cfg = std::move(cfg2); q = std::move(q2); qd = std::move(qd2); qdd = std::move(qdd2);
        T1 = (cfg.W2.size() > 0 && cfg.W2(0) > 0) ? real(2) * kPi / std::sqrt(cfg.W2(0)) : opts.maxTime;
        storeFrame(t, scatterToGlobal(cfg.Phi * q, work, cfg.fmap, cfg.N, true),     // post-inheritance snapshot
                      scatterToGlobal(cfg.Phi * qd, work, cfg.fmap, cfg.N, false));   // (new config, inherited q/qd)
        ++events; quietTime = 0;
        if (events >= opts.maxEvents) { H.outcome = CollapseOutcome::MaxSteps; H.diagnostic = "event budget exhausted with events still occurring"; return H; }
    }

    H.outcome = CollapseOutcome::MaxSteps; H.diagnostic = "time horizon reached with motion still active";
    return H;
}

}  // namespace frame
