// exp_mechanism_guard_oracle.cpp -- research-only mechanism/fallback oracle.
//
// Iterative residuals are not mechanism detectors. Any future HP-FEM lane must
// preserve the existing LDLT/element-prepare singular guard and fall back or
// reject before PCG when the prepared system is singular.

#include "research_common.h"
#include "FrameTestFixtures.h"

#include <algorithm>
#include <cstdio>
#include <string>
#include <vector>

using namespace research;

namespace {

struct Outcome {
    bool ok = false;
    bool preparedSingular = false;
    bool solveSingular = false;
    double pivotMargin = 0.0;
};

FrameModel isolatedNodeMechanism(const Material& mat, const Section& sec) {
    FrameModel m;
    fixtures::cantileverTipLoad(m, 1000.0, 1000.0, mat, sec);
    if (!m.members.empty()) m.members[0].active = false;
    return m;
}

Outcome runCase(const char* name, const FrameModel& model, const SolveOptions& opts, bool expectSingular) {
    Outcome out;
    PreparedSystem ps = assembleAndFactor(model, opts);
    const PreparedSystem::Impl& S = *ps.impl;
    SolveResult r = solve(model, opts);
    out.preparedSingular = S.singular;
    out.solveSingular = r.singular;
    out.pivotMargin = S.pivotMargin;

    if (expectSingular) {
        out.ok = out.preparedSingular && out.solveSingular && !r.diagnostic.empty();
    } else {
        out.ok = !out.preparedSingular && !out.solveSingular && out.pivotMargin > 0.0;
    }

    std::printf("[mechanism_guard] case=%s status=%s expectSingular=%d preparedSingular=%d solveSingular=%d pivotMargin=%.3e diagnostic=\"%s\"\n",
                name, out.ok ? "ok" : "fail", expectSingular ? 1 : 0,
                out.preparedSingular ? 1 : 0, out.solveSingular ? 1 : 0,
                out.pivotMargin, r.diagnostic.c_str());
    return out;
}

} // namespace

int main() {
    Section sec = Section::Rectangular(100.0, 100.0);
    Material mat(210000.0, 80769.0, 7850.0);

    int passed = 0;
    int total = 0;
    auto record = [&](const Outcome& o) {
        ++total;
        if (o.ok) ++passed;
    };

    {
        FrameModel m;
        fixtures::mechanism(m, mat, sec);
        record(runCase("underconstrained_member", m, SolveOptions{}, true));
    }
    {
        FrameModel m;
        fixtures::torsionReleaseMechanism(m, mat, sec);
        SolveOptions opts;
        opts.enableReleases = true;
        record(runCase("release_condensation_mechanism", m, opts, true));
    }
    {
        FrameModel m = isolatedNodeMechanism(mat, sec);
        record(runCase("isolated_node_after_removal", m, SolveOptions{}, true));
    }
    {
        FrameModel m;
        fixtures::cantileverTipLoad(m, 1000.0, 1000.0, mat, sec);
        record(runCase("cantilever_control", m, SolveOptions{}, false));
    }

    std::printf("[mechanism_guard_summary] passed=%d total=%d\n", passed, total);
    return passed == total ? 0 : 1;
}
