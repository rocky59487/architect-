//
// HP-FEM seeded-session benchmark (diagnostic, NOT a gate leg).
//
// Measures the per-frame wall-clock of the production HpSession against a reused-LDLT solveLoad on
// a fixed factorized structure, across problem sizes and thread counts. This is the honest evidence
// for the session's performance niche: the in-subspace seeded path (gravity + a known contact set)
// vs the direct back-substitution, plus the out-of-subspace full-PCG frame (the HP worst case).
//
// Timing is hardware-dependent and sensitive to other system load; it reports the MEDIAN over many
// reps. It is not a correctness check (the five-leg gate owns that) — it only quantifies speed.
//
#include "FrameCore/FrameSolver.h"
#include "FrameCore/HpSession.h"
#include "FrameCore/FrameModel.h"
#include "FrameCore/SolveResult.h"

#include <chrono>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <string>

using namespace frame;
using Clock = std::chrono::high_resolution_clock;

static volatile double g_sink = 0;   // defeats dead-code elimination of the timed solve

template <typename Fn>
static double medianMs(int reps, Fn fn) {
    std::vector<double> t;
    t.reserve((size_t)reps);
    for (int i = 0; i < reps; ++i) {
        const auto t0 = Clock::now();
        fn();
        const auto t1 = Clock::now();
        t.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count());
    }
    std::sort(t.begin(), t.end());
    return t[t.size() / 2];
}

// n-segment horizontal cantilever along +X (nf = 6n): a clean, scalable banded system.
static void buildCantilever(FrameModel& m, int n, real L, const Material& mat, const Section& sec) {
    m = FrameModel{};
    m.materials = { mat }; m.sections = { sec };
    Node n0(0, 0, 0, 0); n0.fixAll(); m.nodes = { n0 };
    for (int k = 1; k <= n; ++k) m.nodes.push_back(Node(k, L * real(k) / real(n), 0, 0));
    for (int k = 0; k < n; ++k) m.members.push_back(Member(k, k, k + 1, 0, 0));
}

static std::vector<real> unitUz(const FrameModel& m, int nd, real mag) {
    std::vector<real> lv((size_t)(6 * (int)m.nodes.size()), 0.0);
    lv[(size_t)gdof(nd, Uz)] = mag;
    return lv;
}

int main() {
    std::printf("# HP-FEM seeded-session benchmark | build %s\n", FRAMECORE_BUILD_SHA);
    std::printf("# production HpSession vs reused-LDLT solveLoad, per-frame MEDIAN ms (close other load)\n\n");

    const Material mat(210000.0, 80769.0, 7850.0);
    const Section  sec = Section::Rectangular(100.0, 100.0);
    const real     L   = 3000.0;
    const int      reps = 200;

    std::printf("%-8s %-7s | %-10s %-10s %-10s %-11s | %-7s %-7s %-8s | %-10s\n",
                "scale", "nf", "LDLT ms", "HP-ser", "HP-par8", "HP-disp8", "spd-ser", "spd-par", "spd-disp", "out-sub");
    std::printf("%s\n", "------------------------------------------------------------------------------------------------");

    for (int n : { 64, 256, 1024, 4096 }) {
        FrameModel m; buildCantilever(m, n, L, mat, sec);
        PreparedSystem ps = assembleAndFactor(m);
        if (ps.pivotMargin() <= 0) { std::printf("n=%-6d singular\n", n); continue; }
        const int nf = 6 * n;

        // in-subspace per-frame load: a pure multiple of the node-1 seed
        FrameModel mf = m;
        { NodalLoad nl; nl.node = 1; nl.comp[Uz] = 1234.0; mf.nodalLoads = { nl }; }
        // out-of-subspace per-frame load: an unseeded interior node
        FrameModel mo = m;
        { NodalLoad nl; nl.node = n / 2; nl.comp[Uz] = 1000.0; mo.nodalLoads = { nl }; }

        const double ldltMs = medianMs(reps, [&] {
            SolveResult r = solveLoad(ps, mf);
            g_sink += r.u.size() > 6 ? r.u[6] : 0.0;
        });

        HpSessionOptions so; so.threads = 1;
        HpSession ss(ps, so);
        ss.setLoadBasis({ unitUz(m, 1, 1.0), unitUz(m, n, 1.0) });
        HpSessionStats st;
        const double hpSerMs = medianMs(reps, [&] {
            SolveResult r = ss.solveFrame(mf, &st);
            g_sink += r.u.size() > 6 ? r.u[6] : 0.0;
        });
        const bool proj = st.usedProjection;

        HpSessionOptions po; po.threads = 8;
        HpSession sp(ps, po);
        sp.setLoadBasis({ unitUz(m, 1, 1.0), unitUz(m, n, 1.0) });
        const double hpParMs = medianMs(reps, [&] {
            SolveResult r = sp.solveFrame(mf, &st);
            g_sink += r.u.size() > 6 ? r.u[6] : 0.0;
        });
        const double hpOutMs = medianMs(reps, [&] {
            SolveResult r = sp.solveFrame(mo, &st);
            g_sink += r.u.size() > 6 ? r.u[6] : 0.0;
        });

        // displacements-only fast path (parallel, in-subspace): skips reactions + force recovery
        HpSessionOptions pd = po; pd.displacementsOnly = true;
        HpSession spd(ps, pd);
        spd.setLoadBasis({ unitUz(m, 1, 1.0), unitUz(m, n, 1.0) });
        const double hpDispMs = medianMs(reps, [&] {
            SolveResult r = spd.solveFrame(mf, &st);
            g_sink += r.u.size() > 6 ? r.u[6] : 0.0;
        });

        std::printf("n=%-6d %-7d | %-10.4f %-10.4f %-10.4f %-11.4f | %-7.2f %-7.2f %-8.2f | %-10.4f %s\n",
                    n, nf, ldltMs, hpSerMs, hpParMs, hpDispMs,
                    ldltMs / hpSerMs, ldltMs / hpParMs, ldltMs / hpDispMs, hpOutMs,
                    proj ? "" : "(NOT in-sub!)");
    }

    std::printf("\n# spd-* = LDLT / HP  (>1 => HP faster than a reused LDLT back-substitution).\n");
    std::printf("# HP-disp8 = displacements-only fast path (skips reactions + force recovery, parallel x8).\n");
    std::printf("# out-sub = parallel full-PCG frame for an UNSEEDED load (the HP worst case).\n");
    std::printf("# sink=%g\n", (double)g_sink);
    return 0;
}
