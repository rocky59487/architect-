#include "FrameCore/FrameSolver.h"
#include "PreparedSystemImpl.h"
#include "BeamColumnElement.h"
#include "MITC4ShellElement.h"

#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <memory>

namespace frame {

PreparedSystem::PreparedSystem() : impl(std::make_unique<Impl>()) {}
PreparedSystem::~PreparedSystem() = default;
PreparedSystem::PreparedSystem(PreparedSystem&&) noexcept = default;
PreparedSystem& PreparedSystem::operator=(PreparedSystem&&) noexcept = default;

// ============================================================================
// Phase 1: assemble + factorize. Everything that depends on geometry / topology /
// support flags / distributed loads is done ONCE here. The result is reused by any
// number of solveLoad() calls that only vary nodal loads + prescribed values.
// ============================================================================
PreparedSystem assembleAndFactor(const FrameModel& model, const SolveOptions& opts) {
    PreparedSystem ps;
    PreparedSystem::Impl& S = *ps.impl;

    std::string why;
    if (!model.validate(why)) {
        S.singular = true;
        S.diagnostic = "invalid model: " + why;
        S.N = std::max(0, model.dofCount());
        return ps;
    }
    S.N = model.dofCount();
    const int N = S.N;

    // build the element list (beams + shells), prepare geometry + baked distributed loads.
    S.elems.reserve(model.members.size() + model.shells.size());
    for (size_t e = 0; e < model.members.size(); ++e)
        S.elems.push_back(std::make_unique<BeamColumnElement>((int)e));
    for (size_t s = 0; s < model.shells.size(); ++s)
        S.elems.push_back(std::make_unique<MITC4ShellElement>((int)s));

    for (auto& el : S.elems) {
        std::string ewhy;
        if (!el->prepare(model, opts, ewhy)) {
            S.singular = true;
            S.diagnostic = ewhy;
        }
    }
    if (S.singular) return ps;

    // assemble global K
    std::vector<Triplet> trips;
    size_t cap = 0;
    for (const auto& el : S.elems) { const size_t d = (size_t)el->localDof(); cap += d * d; }
    trips.reserve(cap);
    for (const auto& el : S.elems) el->assemble(trips);
    SpMat K(N, N);
    K.setFromTriplets(trips.begin(), trips.end());
    K.makeCompressed();
    S.K = std::move(K);

    // free-DOF map from support FLAGS (prescribed VALUES are applied later, in solveLoad)
    S.fmap.assign(N, -1);
    S.nf = 0;
    for (size_t k = 0; k < model.nodes.size(); ++k)
        for (int d = 0; d < 6; ++d) {
            const int g = gdof((int)k, d);
            if (!model.nodes[k].fixed[d]) S.fmap[g] = S.nf++;
        }
    if (S.nf == 0) { S.singular = true; S.diagnostic = "fully constrained (no free DOF)"; return ps; }

    // reduced K_ff (free-free block; independent of prescribed values)
    std::vector<Triplet> ftrips;
    for (int c = 0; c < N; ++c)
        for (SpMat::InnerIterator it(S.K, c); it; ++it) {
            const int r = it.row();
            if (S.fmap[r] < 0) continue;
            if (S.fmap[c] >= 0) ftrips.emplace_back(S.fmap[r], S.fmap[c], it.value());
        }
    SpMat Kff(S.nf, S.nf);
    Kff.setFromTriplets(ftrips.begin(), ftrips.end());
    Kff.makeCompressed();

    // factor + mechanism detection (from the factorization, never from connectivity)
    S.ldlt.compute(Kff);
    if (S.ldlt.info() != Eigen::Success) {
        S.singular = true;
        S.diagnostic = "LDLT factorization failed (info != Success): rank-deficient stiffness / mechanism";
        return ps;
    }
    const VecX D = S.ldlt.vectorD();
    real maxAbs = 0;
    for (int i = 0; i < D.size(); ++i) maxAbs = std::max(maxAbs, std::abs(D(i)));
    const real tol = opts.pivotTol * std::max(real(1), maxAbs);
    for (int i = 0; i < D.size(); ++i) {
        if (std::abs(D(i)) < tol) {
            S.singular = true;
            S.diagnostic = "near-zero pivot in LDLT D (mechanism) at reduced dof #" + std::to_string(i);
            return ps;
        }
        if (D(i) < 0) {
            S.singular = true;
            S.diagnostic = "negative pivot in LDLT D (indefinite stiffness / mechanism) at reduced dof #" + std::to_string(i);
            return ps;
        }
    }
    return ps;
}

// ============================================================================
// Phase 2: cheap re-solve. Reuses the factorization; only re-builds the RHS from the
// current model's nodal loads + prescribed values + baked distributed loads.
// ============================================================================
SolveResult solveLoad(const PreparedSystem& prepared, const FrameModel& model) {
    const PreparedSystem::Impl& S = *prepared.impl;
    SolveResult R;
    const int N = S.N;
    R.u.assign((size_t)std::max(0, N), 0.0);
    R.reactions.assign((size_t)std::max(0, N), 0.0);
    if (S.singular) { R.singular = true; R.diagnostic = S.diagnostic; return R; }

    // global load vector: current nodal loads + baked distributed equivalent loads
    VecX F = VecX::Zero(N);
    for (const auto& nl : model.nodalLoads) {
        const int ni = model.nodeIndex(nl.node);
        if (ni < 0) continue;
        for (int d = 0; d < 6; ++d) F(gdof(ni, d)) += nl.comp[d];
    }
    for (const auto& el : S.elems) el->addEquivalentNodalLoads(F);

    // current prescribed (support displacement) values — may differ between solveLoad calls
    std::vector<real> presc((size_t)N, 0.0);
    for (size_t k = 0; k < model.nodes.size(); ++k)
        for (int d = 0; d < 6; ++d)
            if (model.nodes[k].fixed[d]) presc[(size_t)gdof((int)k, d)] = model.nodes[k].prescribed[d];

    // reduced RHS  F_f - K_fc u_c
    VecX Ff = VecX::Zero(S.nf);
    for (int c = 0; c < N; ++c)
        for (SpMat::InnerIterator it(S.K, c); it; ++it) {
            const int r = it.row();
            if (S.fmap[r] < 0) continue;
            if (S.fmap[c] < 0 && presc[(size_t)c] != 0.0) Ff(S.fmap[r]) -= it.value() * presc[(size_t)c];
        }
    for (int g = 0; g < N; ++g) if (S.fmap[g] >= 0) Ff(S.fmap[g]) += F(g);

    // solve (reuse factorization)
    VecX uf = S.ldlt.solve(Ff);
    if (S.ldlt.info() != Eigen::Success || !uf.allFinite()) {
        R.singular = true;
        R.diagnostic = "solve produced non-finite displacements (mechanism)";
        return R;
    }

    VecX u = VecX::Zero(N);
    for (int g = 0; g < N; ++g) u(g) = (S.fmap[g] >= 0) ? uf(S.fmap[g]) : presc[(size_t)g];
    for (int g = 0; g < N; ++g) R.u[(size_t)g] = u(g);

    const VecX Rv = S.K * u - F;
    for (int g = 0; g < N; ++g) R.reactions[(size_t)g] = Rv(g);

    R.memberForces.resize(model.members.size());
    R.shellForces.resize(model.shells.size());
    for (const auto& el : S.elems) el->recover(u, R);
    return R;
}

// One-shot linear static solve (back-compatible thin wrapper).
SolveResult solve(const FrameModel& model, const SolveOptions& opts) {
    PreparedSystem ps = assembleAndFactor(model, opts);
    return solveLoad(ps, model);
}

} // namespace frame
