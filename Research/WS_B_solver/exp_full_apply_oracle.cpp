// exp_full_apply_oracle.cpp -- research-only full-global matrix-free oracle.
//
// Free-free Kff*x checks are not enough for adoption: solveLoad also needs the
// constrained-column RHS contribution from prescribed supports and full K*u-F
// reactions. This experiment builds a generic element-by-element full operator
// from PreparedSystem::Impl::elems and validates those paths against the
// assembled S.K / solveLoad oracle.

#include "research_common.h"
#include "FrameTestFixtures.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace research;

namespace {

struct Entry {
    int r = 0;
    int c = 0;
    real v = 0;
};

struct ElementEntries {
    std::vector<Entry> entries;
};

struct FullElementOperator {
    int n = 0;
    std::vector<ElementEntries> blocks;
    size_t entryCount = 0;

    explicit FullElementOperator(const PreparedSystem::Impl& S) : n(S.N) {
        blocks.reserve(S.elems.size());
        for (const auto& el : S.elems) {
            std::vector<Triplet> trips;
            const int ldof = std::max(0, el->localDof());
            trips.reserve(static_cast<size_t>(ldof) * static_cast<size_t>(ldof));
            el->assemble(trips);

            ElementEntries b;
            b.entries.reserve(trips.size());
            for (const Triplet& t : trips) {
                const int r = static_cast<int>(t.row());
                const int c = static_cast<int>(t.col());
                if (r < 0 || c < 0 || r >= n || c >= n) continue;
                b.entries.push_back({ r, c, t.value() });
            }
            if (!b.entries.empty()) {
                entryCount += b.entries.size();
                blocks.push_back(std::move(b));
            }
        }
    }

    void apply(const VecX& x, VecX& y) const {
        y.setZero(n);
        for (const ElementEntries& b : blocks)
            for (const Entry& e : b.entries)
                y(e.r) += e.v * x(e.c);
    }
};

VecX fullLoadVector(const FrameModel& model, const PreparedSystem::Impl& S) {
    VecX F = VecX::Zero(S.N);
    for (const auto& nl : model.nodalLoads) {
        const int ni = model.nodeIndex(nl.node);
        if (ni < 0) continue;
        for (int d = 0; d < 6; ++d) F(gdof(ni, d)) += nl.comp[static_cast<size_t>(d)];
    }
    for (const auto& el : S.elems) el->addEquivalentNodalLoads(F);
    return F;
}

VecX prescribedVector(const FrameModel& model, int N) {
    VecX u = VecX::Zero(N);
    for (size_t k = 0; k < model.nodes.size(); ++k)
        for (int d = 0; d < 6; ++d)
            if (model.nodes[k].fixed[d])
                u(gdof(static_cast<int>(k), d)) = model.nodes[k].prescribed[static_cast<size_t>(d)];
    return u;
}

VecX reduceFree(const VecX& v, const PreparedSystem::Impl& S) {
    VecX out = VecX::Zero(S.nf);
    for (int g = 0; g < S.N; ++g)
        if (S.fmap[static_cast<size_t>(g)] >= 0) out(S.fmap[static_cast<size_t>(g)]) = v(g);
    return out;
}

VecX scatterFull(const VecX& uf, const VecX& prescribed, const PreparedSystem::Impl& S) {
    VecX u = prescribed;
    for (int g = 0; g < S.N; ++g)
        if (S.fmap[static_cast<size_t>(g)] >= 0) u(g) = uf(S.fmap[static_cast<size_t>(g)]);
    return u;
}

VecX solveLoadRhsAssembled(const VecX& F, const VecX& prescribed, const PreparedSystem::Impl& S) {
    VecX Ff = VecX::Zero(S.nf);
    for (int c = 0; c < S.N; ++c) {
        if (S.fmap[static_cast<size_t>(c)] >= 0 || prescribed(c) == real(0)) continue;
        for (SpMat::InnerIterator it(S.K, c); it; ++it) {
            const int r = static_cast<int>(it.row());
            if (S.fmap[static_cast<size_t>(r)] >= 0)
                Ff(S.fmap[static_cast<size_t>(r)]) -= it.value() * prescribed(c);
        }
    }
    for (int g = 0; g < S.N; ++g)
        if (S.fmap[static_cast<size_t>(g)] >= 0) Ff(S.fmap[static_cast<size_t>(g)]) += F(g);
    return Ff;
}

double relL2(const VecX& a, const VecX& b) {
    return static_cast<double>((a - b).norm() / std::max<real>(real(1e-300), b.norm()));
}

double relMaxVec(const VecX& a, const VecX& b) {
    double maxDiff = 0.0;
    double maxRef = 0.0;
    for (int i = 0; i < std::min(a.size(), b.size()); ++i) {
        maxDiff = std::max(maxDiff, std::abs(static_cast<double>(a(i) - b(i))));
        maxRef = std::max(maxRef, std::abs(static_cast<double>(b(i))));
    }
    return maxRef > 0.0 ? maxDiff / maxRef : maxDiff;
}

VecX resultVector(const std::vector<real>& values) {
    VecX v(static_cast<int>(values.size()));
    for (int i = 0; i < v.size(); ++i) v(i) = values[static_cast<size_t>(i)];
    return v;
}

struct CaseOutcome {
    bool ok = false;
    double applyRel = 0;
    double rhsRel = 0;
    double uRel = 0;
    double reactionRel = 0;
    double reactionMaxRel = 0;
    int reactionMaxDof = -1;
    double reactionMaxAbs = 0;
};

CaseOutcome runCase(const char* name, const FrameModel& preparedModel, const FrameModel& solveModel) {
    CaseOutcome out;
    PreparedSystem ps = assembleAndFactor(preparedModel);
    const PreparedSystem::Impl& S = *ps.impl;
    if (S.singular) {
        std::printf("[full_apply] case=%s status=singular diagnostic=\"%s\"\n", name, S.diagnostic.c_str());
        return out;
    }

    const SolveResult engine = solveLoad(ps, solveModel);
    if (engine.singular) {
        std::printf("[full_apply] case=%s status=solveLoad_failed diagnostic=\"%s\"\n", name, engine.diagnostic.c_str());
        return out;
    }

    FullElementOperator A(S);
    VecX probe = VecX::Zero(S.N);
    for (int i = 0; i < S.N; ++i)
        probe(i) = std::sin(0.007 * static_cast<double>(i + 1)) + real(0.25) * std::cos(0.013 * static_cast<double>(i + 5));
    VecX yMf(S.N);
    A.apply(probe, yMf);
    out.applyRel = relL2(yMf, S.K * probe);

    const VecX F = fullLoadVector(solveModel, S);
    const VecX uc = prescribedVector(solveModel, S.N);
    VecX Kuc(S.N);
    A.apply(uc, Kuc);

    const VecX FfMf = reduceFree(F - Kuc, S);
    const VecX FfRef = solveLoadRhsAssembled(F, uc, S);
    out.rhsRel = relL2(FfMf, FfRef);

    const VecX uf = S.ldlt.solve(FfMf);
    const VecX uMf = scatterFull(uf, uc, S);
    const VecX uEngine = resultVector(engine.u);
    out.uRel = relL2(uMf, uEngine);

    VecX Ku(S.N);
    A.apply(uMf, Ku);
    const VecX reactionMf = Ku - F;
    const VecX reactionEngine = resultVector(engine.reactions);
    out.reactionRel = relL2(reactionMf, reactionEngine);
    out.reactionMaxRel = relMaxVec(reactionMf, reactionEngine);
    for (int i = 0; i < std::min(reactionMf.size(), reactionEngine.size()); ++i) {
        const double d = std::abs(static_cast<double>(reactionMf(i) - reactionEngine(i)));
        if (d > out.reactionMaxAbs) {
            out.reactionMaxAbs = d;
            out.reactionMaxDof = i;
        }
    }

    out.ok = out.applyRel <= 1e-10 &&
             out.rhsRel <= 1e-10 &&
             out.uRel <= 1e-9 &&
             (out.reactionRel <= 1e-9 || out.reactionMaxAbs <= 1e-6) &&
             (out.reactionMaxRel <= 1e-8 || out.reactionMaxAbs <= 1e-6);

    std::printf("[full_apply] case=%s status=%s N=%d nf=%d elems=%zu entries=%zu applyRel=%.3e rhsRel=%.3e uRel=%.3e reactionRel=%.3e reactionMaxRel=%.3e reactionMaxDof=%d reactionMaxAbs=%.3e\n",
                name, out.ok ? "ok" : "fail", S.N, S.nf, A.blocks.size(), A.entryCount,
                out.applyRel, out.rhsRel, out.uRel, out.reactionRel, out.reactionMaxRel,
                out.reactionMaxDof, out.reactionMaxAbs);
    return out;
}

FrameModel settlementCase() {
    Section sec = Section::Rectangular(100.0, 100.0);
    Material mat(210000.0, 80769.0, 7850.0);
    FrameModel m;
    fixtures::clampedSettlement(m, 2000.0, 1.0, mat, sec);
    return m;
}

FrameModel udlCase() {
    Section sec = Section::Rectangular(100.0, 100.0);
    Material mat(210000.0, 80769.0, 7850.0);
    FrameModel m;
    fixtures::simplySupportedUDL(m, 5.0, 3000.0, mat, sec);
    return m;
}

FrameModel mbBaseCase() {
    Section sec = Section::Rectangular(100.0, 100.0);
    Material mat(210000.0, 80769.0, 7850.0);
    FrameModel m;
    fixtures::simplySupportedBeamN(m, 8, 4000.0, mat, sec);
    return m;
}

FrameModel towerCase() {
    return makeTower(5, 4, 8);
}

} // namespace

int main() {
    int passed = 0;
    int total = 0;
    double maxApply = 0;
    double maxRhs = 0;
    double maxU = 0;
    double maxReaction = 0;
    double maxReactionAbs = 0;

    auto record = [&](const CaseOutcome& c) {
        ++total;
        if (c.ok) ++passed;
        maxApply = std::max(maxApply, c.applyRel);
        maxRhs = std::max(maxRhs, c.rhsRel);
        maxU = std::max(maxU, c.uRel);
        maxReaction = std::max(maxReaction, c.reactionRel);
        maxReactionAbs = std::max(maxReactionAbs, c.reactionMaxAbs);
    };

    {
        const FrameModel m = settlementCase();
        record(runCase("settlement_prescribed", m, m));
    }
    {
        const FrameModel m = udlCase();
        record(runCase("udl_equiv_load_reaction", m, m));
    }
    {
        const FrameModel prepared = mbBaseCase();
        FrameModel solve = prepared;
        solve.nodalLoads.clear();
        solve.nodes[0].prescribed[Uz] = 1.0;
        record(runCase("muller_breslau_prescribed_after_prepare", prepared, solve));
    }
    {
        const FrameModel prepared = towerCase();
        FrameModel solve = prepared;
        solve.nodes[0].prescribed[Ux] = 1.25;
        solve.nodes[0].prescribed[Rz] = 0.002;
        record(runCase("tower_base_prescribed_with_loads", prepared, solve));
    }

    std::printf("[full_apply_summary] passed=%d total=%d maxApplyRel=%.3e maxRhsRel=%.3e maxURel=%.3e maxReactionRel=%.3e maxReactionAbs=%.3e\n",
                passed, total, maxApply, maxRhs, maxU, maxReaction, maxReactionAbs);
    return passed == total ? 0 : 1;
}
