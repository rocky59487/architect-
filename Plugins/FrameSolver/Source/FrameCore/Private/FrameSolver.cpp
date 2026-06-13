#include "FrameCore/FrameSolver.h"
#include "PreparedSystemImpl.h"
#include "BeamColumnElement.h"
#include "MITC4ShellElement.h"
#include "ElementFactory.h"

#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <memory>
#include <cstring>
#include <cstdint>

namespace frame {

namespace {
// Structural fingerprint: hashes everything solveLoad() must NOT change between an
// assembleAndFactor() and its reuse — node id/positions/support FLAGS, member
// id/connectivity/refVec/releases/matIdx/secIdx/active, shell id/connectivity/thickness/matIdx/active,
// the referenced material VALUES (E/G/nu/rho) and section VALUES (A/Iy/Iz/J/Asy/Asz), the
// baked distributed loads (member UDLs, shell pressures), and the plastic hinges (a hinge
// changes BOTH the condensed stiffness and the baked fixed-end forces). The factorization
// bakes in the stiffness those properties imply, so changing E, Iz, or which material/section
// an element points to would make a reused factorization a SILENT STALE SOLVE — they must be
// fingerprinted (the section/material values were the gap; see PROJECT.txt P1). It deliberately
// EXCLUDES nodal loads and prescribed VALUES, which solveLoad is allowed to vary (the
// interactive / settlement path), and the capacity-side fields fy/Zy/Zz/cap (post-processing
// screens only — they never enter K or the baked loads).
inline uint64_t fpMix(uint64_t h, uint64_t v) {
    return h ^ (v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2));
}
inline uint64_t fpBits(real d) {
    uint64_t u = 0; std::memcpy(&u, &d, sizeof(real)); return u;
}
uint64_t modelFingerprint(const FrameModel& m) {
    uint64_t h = 1469598103934665603ULL;
    h = fpMix(h, m.nodes.size());
    h = fpMix(h, m.members.size());
    h = fpMix(h, m.shells.size());
    h = fpMix(h, m.memberUDLs.size());
    h = fpMix(h, m.shellPressures.size());
    for (const auto& n : m.nodes) {
        h = fpMix(h, static_cast<uint64_t>(n.id));
        h = fpMix(h, fpBits(n.pos.x)); h = fpMix(h, fpBits(n.pos.y)); h = fpMix(h, fpBits(n.pos.z));
        uint64_t fb = 0; for (int d = 0; d < 6; ++d) if (n.fixed[d]) fb |= (1ull << d);
        h = fpMix(h, fb);
    }
    for (const auto& mem : m.members) {
        h = fpMix(h, static_cast<uint64_t>(mem.id));
        h = fpMix(h, static_cast<uint64_t>(mem.i)); h = fpMix(h, static_cast<uint64_t>(mem.j));
        h = fpMix(h, static_cast<uint64_t>(static_cast<int64_t>(mem.matIdx)));
        h = fpMix(h, static_cast<uint64_t>(static_cast<int64_t>(mem.secIdx)));
        h = fpMix(h, fpBits(mem.refVec.x)); h = fpMix(h, fpBits(mem.refVec.y)); h = fpMix(h, fpBits(mem.refVec.z));
        uint64_t rb = 0; for (int d = 0; d < 12; ++d) if (mem.release[d]) rb |= (1ull << d);
        h = fpMix(h, rb);
        h = fpMix(h, mem.active ? 1ull : 0ull);   // toggling active is a structural (remove/restore) change
        h = fpMix(h, mem.tensionOnly ? 1ull : 0ull);   // flipping tension-only changes runTensionOnly semantics
    }
    for (const auto& sh : m.shells) {
        h = fpMix(h, static_cast<uint64_t>(static_cast<int64_t>(sh.id)));
        for (int k = 0; k < 4; ++k) h = fpMix(h, static_cast<uint64_t>(sh.n[k]));
        h = fpMix(h, static_cast<uint64_t>(static_cast<int64_t>(sh.matIdx)));
        h = fpMix(h, fpBits(sh.t));
        h = fpMix(h, sh.active ? 1ull : 0ull);   // toggling active is a structural (remove/restore) change
    }
    // Material / section VALUES the elements reference by index. Changing E / Iz / etc. alters
    // the assembled K, so a reused factorization built on the OLD values is a stale solve. We
    // hash all entries (not just referenced ones) — editing an unused entry is harmless
    // over-conservatism, while missing a used one is a silent correctness bug.
    for (const auto& mat : m.materials) {
        h = fpMix(h, fpBits(mat.E)); h = fpMix(h, fpBits(mat.G));
        h = fpMix(h, fpBits(mat.nu)); h = fpMix(h, fpBits(mat.rho));
    }
    for (const auto& s : m.sections) {
        h = fpMix(h, fpBits(s.A));   h = fpMix(h, fpBits(s.Iy));  h = fpMix(h, fpBits(s.Iz));
        h = fpMix(h, fpBits(s.J));   h = fpMix(h, fpBits(s.Asy)); h = fpMix(h, fpBits(s.Asz));
    }
    for (const auto& u : m.memberUDLs) {
        h = fpMix(h, static_cast<uint64_t>(u.member));
        h = fpMix(h, fpBits(u.w_local.x)); h = fpMix(h, fpBits(u.w_local.y)); h = fpMix(h, fpBits(u.w_local.z));
    }
    for (const auto& sp : m.shellPressures) {
        h = fpMix(h, static_cast<uint64_t>(sp.shell)); h = fpMix(h, fpBits(sp.p));
    }
    for (const auto& ph : m.hinges) {   // a hinge alters K (release) AND the baked Qf (Mp)
        h = fpMix(h, static_cast<uint64_t>(ph.member));
        h = fpMix(h, static_cast<uint64_t>(static_cast<int64_t>(ph.dof)));
        h = fpMix(h, fpBits(ph.Mp));
    }
    return h;
}
}  // namespace

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
        if (model.members[e].active)   // inactive members are excluded from assembly (element removal)
            S.elems.push_back(makeMemberElem((int)e));
    for (size_t s = 0; s < model.shells.size(); ++s)
        if (model.shells[s].active)   // inactive shells are excluded from assembly (element removal)
            S.elems.push_back(makeShellElem((int)s));

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
    const SpMat Kff = reduceFF(S.K, S.fmap, S.nf);

    // factor + mechanism detection (from the factorization, never from connectivity)
    S.ldlt.compute(Kff);
    if (S.ldlt.info() != Eigen::Success) {
        S.singular = true;
        S.diagnostic = "LDLT factorization failed (info != Success): rank-deficient stiffness / mechanism";
        return ps;
    }
    const VecX D = S.ldlt.vectorD();
    real maxAbs = 0, minAbs = 0; bool firstPivot = true;
    for (int i = 0; i < D.size(); ++i) {
        const real a = std::abs(D(i));
        maxAbs = std::max(maxAbs, a);
        if (firstPivot || a < minAbs) { minAbs = a; firstPivot = false; }
    }
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
    S.pivotMargin = (maxAbs > 0) ? minAbs / maxAbs : 0;   // C4: proximity-to-singular (min/max pivot)
    S.fingerprint = modelFingerprint(model);   // baseline for the solveLoad reuse-validity guard
    return ps;
}

real PreparedSystem::pivotMargin() const { return impl ? impl->pivotMargin : real(0); }

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

    // Reuse-validity guard: the factorization + baked distributed loads are only valid while the
    // structural model is unchanged. Reject a model whose fingerprint differs (geometry / topology /
    // support flags / member UDLs / shell pressures changed since assembleAndFactor) instead of
    // silently returning a stale-load solve. Nodal loads and prescribed VALUES may still vary
    // freely — they are excluded from the fingerprint (the interactive / settlement path).
    if (modelFingerprint(model) != S.fingerprint) {
        R.singular = true;
        R.diagnostic = "solveLoad: model changed since assembleAndFactor (geometry/topology/"
                       "support flags/distributed loads). Re-run assembleAndFactor.";
        return R;
    }

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
    // Stamp every row's id, including INACTIVE elements (their forces stay zero but recover()
    // never visits them): leaving the default id 0 in a vector that is parallel to the model
    // would make any consumer keying results by id silently read the wrong element.
    for (size_t e = 0; e < model.members.size(); ++e) R.memberForces[e].member = model.members[e].id;
    for (size_t s = 0; s < model.shells.size(); ++s)  R.shellForces[s].shell   = model.shells[s].id;
    for (const auto& el : S.elems) el->recover(u, R);
    R.pivotMargin = S.pivotMargin;   // C4: report the factorization's criticality margin
    return R;
}

// One-shot linear static solve (back-compatible thin wrapper).
SolveResult solve(const FrameModel& model, const SolveOptions& opts) {
    PreparedSystem ps = assembleAndFactor(model, opts);
    return solveLoad(ps, model);
}

} // namespace frame
