// exp_million_dof.cpp — WS-B scale baseline (scratch, NOT engine).
//
// Runs ONE tower size per process (peak-memory readings stay clean) and prints a
// structured line: build / factor / solve wall times + memory. OOM is caught and
// reported as data (exit 3). A PowerShell driver loops sizes in separate processes.

#define NOMINMAX
#include "research_common.h"
#include <windows.h>
#include <psapi.h>
#pragma comment(lib, "psapi.lib")
#include <cstring>
#include <new>

using namespace research;

namespace {

double peakMiB() {
    PROCESS_MEMORY_COUNTERS_EX pmc{};
    if (GetProcessMemoryInfo(GetCurrentProcess(), reinterpret_cast<PROCESS_MEMORY_COUNTERS*>(&pmc), sizeof(pmc)))
        return static_cast<double>(pmc.PeakWorkingSetSize) / (1024.0 * 1024.0);
    return 0;
}
double privMiB() {
    PROCESS_MEMORY_COUNTERS_EX pmc{};
    if (GetProcessMemoryInfo(GetCurrentProcess(), reinterpret_cast<PROCESS_MEMORY_COUNTERS*>(&pmc), sizeof(pmc)))
        return static_cast<double>(pmc.PrivateUsage) / (1024.0 * 1024.0);
    return 0;
}

} // namespace

int main(int argc, char** argv) {
    int nx = 12, ny = 9, st = 24;
    for (int i = 1; i < argc; ++i) {
        auto next = [&]() -> const char* { return (i + 1 < argc) ? argv[++i] : "0"; };
        if      (!std::strcmp(argv[i], "--nx"))      nx = std::atoi(next());
        else if (!std::strcmp(argv[i], "--ny"))      ny = std::atoi(next());
        else if (!std::strcmp(argv[i], "--stories")) st = std::atoi(next());
    }

    Timer tb;
    FrameModel model = makeTower(nx, ny, st);
    const double tBuild = tb.ms();

    int fixedDof = 0;
    for (const Node& n : model.nodes) for (bool f : n.fixed) if (f) ++fixedDof;
    const int nf = model.dofCount() - fixedDof;
    std::printf("[scale-begin] nx=%d ny=%d st=%d nodes=%zu members=%zu nf=%d buildMs=%.0f memMiB=%.0f\n",
                nx, ny, st, model.nodes.size(), model.members.size(), nf, tBuild, privMiB());
    std::fflush(stdout);

    try {
        Timer tf;
        PreparedSystem ps = assembleAndFactor(model);
        const double tFactor = tf.ms();
        const PreparedSystem::Impl& S = *ps.impl;
        if (S.singular) { std::printf("[scale-fail] singular: %s\n", S.diagnostic.c_str()); return 2; }
        const double memAfterFactor = privMiB();

        Timer ts;
        SolveResult r = solveLoad(ps, model);
        const double tSolve = ts.ms();
        if (r.singular) { std::printf("[scale-fail] solve singular\n"); return 2; }

        double umax = 0;
        for (real v : r.u) umax = std::max(umax, std::abs(static_cast<double>(v)));

        std::printf("[scale] nx=%d ny=%d st=%d nf=%d nnzK=%lld factorMs=%.0f solveMs=%.1f umax=%.6g memFactorMiB=%.0f peakMiB=%.0f\n",
                    nx, ny, st, S.nf, static_cast<long long>(S.K.nonZeros()),
                    tFactor, tSolve, umax, memAfterFactor, peakMiB());
        return 0;
    } catch (const std::bad_alloc&) {
        std::printf("[oom] nx=%d ny=%d st=%d nf=%d peakMiB=%.0f\n", nx, ny, st, nf, peakMiB());
        return 3;
    }
}
