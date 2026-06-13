// exp_threaded_element_apply.cpp -- research-only thread-local element apply benchmark.
//
// Stores each active frame member as one fixed 12x12 dense block. The threaded
// apply partitions blocks across std::thread workers; each worker accumulates
// into its own full y_local vector and the caller reduces those vectors after
// join. This intentionally benchmarks the thread-local reduction path before
// any solver API changes.

#include "research_common.h"

#include <algorithm>
#include <array>
#include <condition_variable>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

using namespace research;

namespace {

volatile double gBenchSink = 0.0;

struct Args {
    std::string preset = "small"; // small | xxl
    int repeat = 30;
    int threads = 0;
};

int defaultThreadCount() {
    const unsigned hc = std::thread::hardware_concurrency();
    return hc > 0 ? static_cast<int>(hc) : 1;
}

Args parseArgs(int argc, char** argv) {
    Args a;
    a.threads = defaultThreadCount();
    for (int i = 1; i < argc; ++i) {
        auto next = [&]() -> const char* { return (i + 1 < argc) ? argv[++i] : ""; };
        if      (!std::strcmp(argv[i], "--preset"))  a.preset = next();
        else if (!std::strcmp(argv[i], "--repeat"))  a.repeat = std::max(1, std::atoi(next()));
        else if (!std::strcmp(argv[i], "--threads")) a.threads = std::max(1, std::atoi(next()));
    }
    if (a.preset != "xxl") a.preset = "small";
    return a;
}

FrameModel makeCase(const std::string& preset) {
    if (preset == "xxl") return makeTower(12, 9, 24);
    return makeTower(5, 4, 8);
}

double checksum(const std::vector<real>& v) {
    double s = 0.0;
    const int step = std::max(1, static_cast<int>(v.size()) / 4096);
    for (int i = 0; i < static_cast<int>(v.size()); i += step)
        s += std::abs(static_cast<double>(v[static_cast<size_t>(i)]));
    return s;
}

double relNormV(const std::vector<real>& a, const std::vector<real>& b) {
    double num = 0.0;
    double den = 0.0;
    const size_t n = std::min(a.size(), b.size());
    for (size_t i = 0; i < n; ++i) {
        const double d = static_cast<double>(a[i] - b[i]);
        num += d * d;
        den += static_cast<double>(b[i]) * static_cast<double>(b[i]);
    }
    return std::sqrt(num / std::max(1e-300, den));
}

std::vector<real> toVector(const VecX& v) {
    std::vector<real> out(static_cast<size_t>(v.size()));
    for (int i = 0; i < v.size(); ++i) out[static_cast<size_t>(i)] = v(i);
    return out;
}

VecX toEigen(const std::vector<real>& v) {
    VecX out(static_cast<int>(v.size()));
    for (int i = 0; i < out.size(); ++i) out(i) = v[static_cast<size_t>(i)];
    return out;
}

struct ElementBlock12 {
    std::array<int, 12> rid{};       // reduced free DOF id, or -1 for fixed DOFs
    std::array<real, 144> k{};       // row-major fixed 12x12 global element block
};

struct ElementApply {
    int nf = 0;
    std::vector<ElementBlock12> blocks;

    bool build(const FrameModel& model, const PreparedSystem::Impl& S, std::string& why) {
        nf = S.nf;
        blocks.clear();
        blocks.reserve(model.members.size());

        SolveOptions opts;
        for (size_t e = 0; e < model.members.size(); ++e) {
            if (!model.members[e].active) continue;

            std::vector<Triplet> trips;
            if (!memberGlobalK(model, static_cast<int>(e), opts, trips, why)) return false;

            int gd[12];
            memberGlobalDofs(model, static_cast<int>(e), gd);

            ElementBlock12 b;
            b.rid.fill(-1);
            b.k.fill(real(0));
            int freeDofs = 0;
            for (int i = 0; i < 12; ++i) {
                const int g = gd[i];
                if (g >= 0 && g < static_cast<int>(S.fmap.size())) {
                    b.rid[static_cast<size_t>(i)] = S.fmap[static_cast<size_t>(g)];
                    if (b.rid[static_cast<size_t>(i)] >= 0) ++freeDofs;
                }
            }
            if (freeDofs == 0) continue;

            auto localIndex = [&](int g) -> int {
                for (int i = 0; i < 12; ++i)
                    if (gd[i] == g) return i;
                return -1;
            };

            for (const Triplet& t : trips) {
                const int lr = localIndex(static_cast<int>(t.row()));
                const int lc = localIndex(static_cast<int>(t.col()));
                if (lr >= 0 && lc >= 0)
                    b.k[static_cast<size_t>(lr * 12 + lc)] += t.value();
            }
            blocks.push_back(std::move(b));
        }
        return true;
    }

    void applySequential(const std::vector<real>& x, std::vector<real>& y) const {
        y.assign(static_cast<size_t>(nf), real(0));
        accumulateRange(0, blocks.size(), x, y);
    }

    int effectiveThreads(int requested) const {
        int t = std::max(1, requested);
        if (!blocks.empty())
            t = std::min(t, static_cast<int>(blocks.size()));
        return t;
    }

    void applyThreadLocal(const std::vector<real>& x, std::vector<real>& y, int requestedThreads,
                          std::vector<std::vector<real>>& localY) const {
        const int nt = effectiveThreads(requestedThreads);
        if (static_cast<int>(localY.size()) != nt)
            localY.assign(static_cast<size_t>(nt), std::vector<real>(static_cast<size_t>(nf), real(0)));
        for (auto& v : localY) {
            if (static_cast<int>(v.size()) != nf) v.resize(static_cast<size_t>(nf));
        }

        const size_t chunk = (blocks.size() + static_cast<size_t>(nt) - 1) / static_cast<size_t>(nt);
        std::vector<std::thread> workers;
        workers.reserve(static_cast<size_t>(nt));
        for (int tid = 0; tid < nt; ++tid) {
            const size_t begin = std::min(blocks.size(), static_cast<size_t>(tid) * chunk);
            const size_t end = std::min(blocks.size(), begin + chunk);
            workers.emplace_back([&, tid, begin, end]() {
                std::vector<real>& yy = localY[static_cast<size_t>(tid)];
                std::fill(yy.begin(), yy.end(), real(0));
                accumulateRange(begin, end, x, yy);
            });
        }
        for (std::thread& w : workers) w.join();

        y.assign(static_cast<size_t>(nf), real(0));
        for (int tid = 0; tid < nt; ++tid) {
            const std::vector<real>& yy = localY[static_cast<size_t>(tid)];
            for (int i = 0; i < nf; ++i) y[static_cast<size_t>(i)] += yy[static_cast<size_t>(i)];
        }
    }

    double benchThreadLocalPooled(const std::vector<real>& x, std::vector<real>& y, int requestedThreads,
                                  int repeat, double& guard) const {
        const int nt = effectiveThreads(requestedThreads);
        std::vector<std::vector<real>> localY(static_cast<size_t>(nt),
                                              std::vector<real>(static_cast<size_t>(nf), real(0)));

        std::mutex mu;
        std::condition_variable cvStart;
        std::condition_variable cvDone;
        int generation = 0;
        int done = 0;
        bool stop = false;

        const size_t chunk = (blocks.size() + static_cast<size_t>(nt) - 1) / static_cast<size_t>(nt);
        std::vector<std::thread> workers;
        workers.reserve(static_cast<size_t>(nt));
        for (int tid = 0; tid < nt; ++tid) {
            const size_t begin = std::min(blocks.size(), static_cast<size_t>(tid) * chunk);
            const size_t end = std::min(blocks.size(), begin + chunk);
            workers.emplace_back([&, tid, begin, end]() {
                int seen = 0;
                for (;;) {
                    std::unique_lock<std::mutex> lock(mu);
                    cvStart.wait(lock, [&]() { return stop || generation != seen; });
                    if (stop) return;
                    seen = generation;
                    lock.unlock();

                    std::vector<real>& yy = localY[static_cast<size_t>(tid)];
                    std::fill(yy.begin(), yy.end(), real(0));
                    accumulateRange(begin, end, x, yy);

                    lock.lock();
                    ++done;
                    if (done == nt) cvDone.notify_one();
                }
            });
        }

        Timer t;
        for (int r = 0; r < repeat; ++r) {
            {
                std::lock_guard<std::mutex> lock(mu);
                done = 0;
                ++generation;
            }
            cvStart.notify_all();
            {
                std::unique_lock<std::mutex> lock(mu);
                cvDone.wait(lock, [&]() { return done == nt; });
            }

            y.assign(static_cast<size_t>(nf), real(0));
            for (int tid = 0; tid < nt; ++tid) {
                const std::vector<real>& yy = localY[static_cast<size_t>(tid)];
                for (int i = 0; i < nf; ++i) y[static_cast<size_t>(i)] += yy[static_cast<size_t>(i)];
            }
            guard += checksum(y);
        }
        const double ms = t.ms() / static_cast<double>(repeat);

        {
            std::lock_guard<std::mutex> lock(mu);
            stop = true;
            ++generation;
        }
        cvStart.notify_all();
        for (std::thread& w : workers) w.join();
        return ms;
    }

private:
    void accumulateRange(size_t begin, size_t end, const std::vector<real>& x, std::vector<real>& y) const {
        for (size_t bi = begin; bi < end; ++bi) {
            const ElementBlock12& b = blocks[bi];
            real xe[12];
            for (int i = 0; i < 12; ++i) {
                const int id = b.rid[static_cast<size_t>(i)];
                xe[i] = id >= 0 ? x[static_cast<size_t>(id)] : real(0);
            }

            for (int r = 0; r < 12; ++r) {
                const int rr = b.rid[static_cast<size_t>(r)];
                if (rr < 0) continue;
                const real* row = &b.k[static_cast<size_t>(r * 12)];
                real acc = real(0);
                for (int c = 0; c < 12; ++c) acc += row[c] * xe[c];
                y[static_cast<size_t>(rr)] += acc;
            }
        }
    }
};

} // namespace

int main(int argc, char** argv) {
    const Args args = parseArgs(argc, argv);
    FrameModel model = makeCase(args.preset);
    assertNodalOnly(model, "exp_threaded_element_apply");

    PreparedSystem ps = assembleAndFactor(model);
    const PreparedSystem::Impl& S = *ps.impl;
    if (S.singular) {
        std::printf("[threaded_apply] status=singular diagnostic=\"%s\"\n", S.diagnostic.c_str());
        return 2;
    }

    std::string why;
    ElementApply op;
    if (!op.build(model, S, why)) {
        std::printf("[threaded_apply] status=build_failed why=\"%s\"\n", why.c_str());
        return 2;
    }

    std::vector<real> x(static_cast<size_t>(S.nf));
    for (int i = 0; i < S.nf; ++i)
        x[static_cast<size_t>(i)] = std::sin(0.00091 * static_cast<double>(i + 1))
                                  + 0.5 * std::cos(0.00037 * static_cast<double>(i + 11));

    const SpMat Kff = research::reduceFF(S.K, S.fmap, S.nf);
    const std::vector<real> yRef = toVector(Kff * toEigen(x));

    std::vector<real> ySeq;
    std::vector<real> ySpawn;
    std::vector<real> yPool;
    std::vector<std::vector<real>> localY;
    op.applySequential(x, ySeq);
    op.applyThreadLocal(x, ySpawn, args.threads, localY);
    double guard = checksum(ySeq) + checksum(ySpawn);
    const double poolMsWarmup = op.benchThreadLocalPooled(x, yPool, args.threads, 1, guard);
    (void)poolMsWarmup;
    const double applyRel = std::max(relNormV(ySeq, yRef),
                                     std::max(relNormV(ySpawn, ySeq), relNormV(yPool, ySeq)));


    Timer tSeq;
    for (int r = 0; r < args.repeat; ++r) {
        op.applySequential(x, ySeq);
        guard += checksum(ySeq);
    }
    const double seqMs = tSeq.ms() / static_cast<double>(args.repeat);

    Timer tSpawn;
    for (int r = 0; r < args.repeat; ++r) {
        op.applyThreadLocal(x, ySpawn, args.threads, localY);
        guard += checksum(ySpawn);
    }
    const double spawnParMs = tSpawn.ms() / static_cast<double>(args.repeat);
    const double parMs = op.benchThreadLocalPooled(x, yPool, args.threads, args.repeat, guard);
    gBenchSink = guard;

    const int threads = op.effectiveThreads(args.threads);
    std::printf("[threaded_apply] preset=%s applyRel=%.3e seqMs=%.6f parMs=%.6f speedup=%.3f spawnParMs=%.6f spawnSpeedup=%.3f threads=%d nf=%d blocks=%zu\n",
                args.preset.c_str(), applyRel, seqMs, parMs, seqMs / std::max(1e-300, parMs),
                spawnParMs, seqMs / std::max(1e-300, spawnParMs), threads, S.nf, op.blocks.size());
    return (applyRel <= 1e-10) ? 0 : 1;
}
