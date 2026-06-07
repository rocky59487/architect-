// In-process FrameCore performance benchmark.
//
// This avoids frame_cli.exe startup and stdout costs, so the reported times are
// closer to the real engine call: build element objects, assemble K/F, factor,
// solve, recover reactions/member forces.
#include "FrameCore/FrameSolver.h"
#include "FrameCore/Section.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

using namespace frame;

namespace {

struct BenchSpec {
    int nx = 12;
    int ny = 9;
    int stories = 24;
    int repeat = 5;
    int warmup = 1;
    bool dry = false;
    std::string name = "XXL-24st-12x9";
};

int nodeId(int i, int j, int k, int nx, int ny) {
    return k * (nx + 1) * (ny + 1) + j * (nx + 1) + i;
}

Section squareSection(real side) {
    Section s = Section::Rectangular(side, side);
    s.J = 0.1406 * side * side * side * side;
    return s;
}

Vec3 refVecFor(const Vec3& pi, const Vec3& pj) {
    Vec3 x = pj - pi;
    const real n = norm(x);
    if (n > 0) x = x * (1.0 / n);
    // Match the Python benchmark: prefer global Z unless nearly parallel.
    return (std::abs(x.z) < 0.92) ? Vec3(0, 0, 1) : Vec3(1, 0, 0);
}

FrameModel makeTower(int nx, int ny, int stories) {
    constexpr real sx = 6000.0, sy = 5000.0, h = 3300.0;
    FrameModel m;
    m.materials.reserve(1);
    m.sections.reserve(3);
    m.nodes.reserve(static_cast<size_t>(nx + 1) * (ny + 1) * (stories + 1));

    m.materials.emplace_back(200000.0, 76923.07692307692, 7850.0);
    m.sections.push_back(squareSection(520.0));
    m.sections.push_back(squareSection(360.0));
    m.sections.push_back(squareSection(220.0));

    for (int k = 0; k <= stories; ++k) {
        const real driftX = 35.0 * k;
        const real driftY = -18.0 * k;
        for (int j = 0; j <= ny; ++j) {
            for (int i = 0; i <= nx; ++i) {
                Node n(nodeId(i, j, k, nx, ny), i * sx + driftX, j * sy + driftY, k * h);
                if (k == 0) n.fixAll();
                m.nodes.push_back(n);
            }
        }
    }

    auto nodePos = [&](int id) -> Vec3 {
        return m.nodes[static_cast<size_t>(id)].pos; // ids are contiguous by construction.
    };
    auto addMember = [&](int i, int j, int sec) {
        Member mem(static_cast<int>(m.members.size()), i, j, 0, sec);   // material index 0
        mem.refVec = refVecFor(nodePos(i), nodePos(j));
        m.members.push_back(mem);
    };

    const size_t columnCount = static_cast<size_t>(nx + 1) * (ny + 1) * stories;
    const size_t beamPerFloor = static_cast<size_t>(nx) * (ny + 1) + static_cast<size_t>(ny) * (nx + 1);
    const size_t bracePerStory = static_cast<size_t>(2 * nx + 2 * ny);
    m.members.reserve(columnCount + beamPerFloor * stories + bracePerStory * stories);

    for (int k = 0; k < stories; ++k)
        for (int j = 0; j <= ny; ++j)
            for (int i = 0; i <= nx; ++i)
                addMember(nodeId(i, j, k, nx, ny), nodeId(i, j, k + 1, nx, ny), 0);

    for (int k = 1; k <= stories; ++k) {
        for (int j = 0; j <= ny; ++j)
            for (int i = 0; i < nx; ++i)
                addMember(nodeId(i, j, k, nx, ny), nodeId(i + 1, j, k, nx, ny), 1);
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i <= nx; ++i)
                addMember(nodeId(i, j, k, nx, ny), nodeId(i, j + 1, k, nx, ny), 1);
    }

    for (int k = 0; k < stories; ++k) {
        for (int i = 0; i < nx; ++i) {
            for (int face = 0; face < 2; ++face) {
                const int jf = face ? ny : 0;
                const int a = ((i + k) % 2 == 0) ? nodeId(i, jf, k, nx, ny)     : nodeId(i + 1, jf, k, nx, ny);
                const int b = ((i + k) % 2 == 0) ? nodeId(i + 1, jf, k + 1, nx, ny) : nodeId(i, jf, k + 1, nx, ny);
                addMember(a, b, 2);
            }
        }
        for (int j = 0; j < ny; ++j) {
            for (int face = 0; face < 2; ++face) {
                const int ifc = face ? nx : 0;
                const int a = ((j + k) % 2 == 0) ? nodeId(ifc, j, k, nx, ny)     : nodeId(ifc, j + 1, k, nx, ny);
                const int b = ((j + k) % 2 == 0) ? nodeId(ifc, j + 1, k + 1, nx, ny) : nodeId(ifc, j, k + 1, nx, ny);
                addMember(a, b, 2);
            }
        }
    }

    m.nodalLoads.reserve(static_cast<size_t>(nx + 1) * (ny + 1) * stories + 2);
    for (int k = 1; k <= stories; ++k) {
        const real sf = static_cast<real>(k) / static_cast<real>(stories);
        for (int j = 0; j <= ny; ++j) {
            for (int i = 0; i <= nx; ++i) {
                const int edge = (i == 0 || i == nx ? 1 : 0) + (j == 0 || j == ny ? 1 : 0);
                const real trib = edge == 0 ? 1.0 : (edge == 1 ? 0.65 : 0.42);
                NodalLoad l;
                l.node = nodeId(i, j, k, nx, ny);
                l.comp[Ux] = 900.0 * sf * (1.0 + 0.08 * (j - ny / 2.0));
                l.comp[Uy] = -350.0 * sf * (1.0 + 0.05 * (i - nx / 2.0));
                l.comp[Uz] = -38000.0 * trib;
                m.nodalLoads.push_back(l);
            }
        }
    }
    NodalLoad t1; t1.node = nodeId(nx, ny, stories, nx, ny); t1.comp[Rz] = 3.0e6; m.nodalLoads.push_back(t1);
    NodalLoad t2; t2.node = nodeId(0, 0, stories, nx, ny);  t2.comp[Ry] = -2.0e6; m.nodalLoads.push_back(t2);
    return m;
}

double checksum(const SolveResult& r) {
    double s = 0.0;
    for (size_t i = 0; i < r.u.size(); i += std::max<size_t>(1, r.u.size() / 4096)) s += std::abs(r.u[i]);
    for (size_t i = 0; i < r.reactions.size(); i += std::max<size_t>(1, r.reactions.size() / 4096)) s += 1e-9 * std::abs(r.reactions[i]);
    for (size_t i = 0; i < r.memberForces.size(); i += std::max<size_t>(1, r.memberForces.size() / 4096)) {
        s += 1e-9 * (std::abs(r.memberForces[i].endI.N) + std::abs(r.memberForces[i].endI.My) + std::abs(r.memberForces[i].endJ.Mz));
    }
    return s;
}

void applyPreset(const std::string& p, BenchSpec& s) {
    if (p == "xxl")      { s = BenchSpec{12, 9, 24, s.repeat, s.warmup, s.dry, "XXL-24st-12x9"}; }
    else if (p == "mega") { s = BenchSpec{18, 14, 36, s.repeat, s.warmup, s.dry, "MEGA-36st-18x14"}; }
    else if (p == "million-dof-est") { s = BenchSpec{60, 50, 54, s.repeat, s.warmup, true, "EST-1MfreeDOF-60x50x54"}; }
}

BenchSpec parseArgs(int argc, char** argv) {
    BenchSpec s;
    for (int i = 1; i < argc; ++i) {
        auto next = [&]() -> const char* { return (i + 1 < argc) ? argv[++i] : ""; };
        if (std::strcmp(argv[i], "--preset") == 0) applyPreset(next(), s);
        else if (std::strcmp(argv[i], "--nx") == 0) s.nx = std::atoi(next());
        else if (std::strcmp(argv[i], "--ny") == 0) s.ny = std::atoi(next());
        else if (std::strcmp(argv[i], "--stories") == 0) s.stories = std::atoi(next());
        else if (std::strcmp(argv[i], "--repeat") == 0) s.repeat = std::max(1, std::atoi(next()));
        else if (std::strcmp(argv[i], "--warmup") == 0) s.warmup = std::max(0, std::atoi(next()));
        else if (std::strcmp(argv[i], "--dry") == 0) s.dry = true;
    }
    return s;
}

} // namespace

int main(int argc, char** argv) {
    BenchSpec spec = parseArgs(argc, argv);
    std::printf("FrameCore precise in-process benchmark\n");
    std::printf("case=%s nx=%d ny=%d stories=%d repeat=%d warmup=%d dry=%d\n",
                spec.name.c_str(), spec.nx, spec.ny, spec.stories, spec.repeat, spec.warmup, spec.dry ? 1 : 0);

    auto tBuild0 = std::chrono::steady_clock::now();
    FrameModel model = makeTower(spec.nx, spec.ny, spec.stories);
    auto tBuild1 = std::chrono::steady_clock::now();

    const int dof = model.dofCount();
    int fixed = 0;
    for (const Node& n : model.nodes) for (bool f : n.fixed) if (f) ++fixed;
    const int freeDof = dof - fixed;
    const size_t tripletContrib = model.members.size() * 144ull;
    const double tripletMiB = static_cast<double>(tripletContrib) * 24.0 / (1024.0 * 1024.0);
    const double buildMs = std::chrono::duration<double, std::milli>(tBuild1 - tBuild0).count();

    std::printf("nodes=%zu members=%zu dof=%d freeDof=%d loads=%zu tripletContrib=%zu (~%.1f MiB raw triplets) buildMs=%.3f\n",
                model.nodes.size(), model.members.size(), dof, freeDof, model.nodalLoads.size(),
                tripletContrib, tripletMiB, buildMs);
    if (spec.dry) {
        std::printf("DRY-RUN only. Estimated factorization memory may be many times raw triplets.\n");
        return 0;
    }

    SolveOptions opt;
    double cs = 0.0;
    for (int i = 0; i < spec.warmup; ++i) {
        const SolveResult r = solve(model, opt);
        if (r.singular) { std::printf("warmup singular: %s\n", r.diagnostic.c_str()); return 2; }
        cs += checksum(r);
    }

    std::vector<double> ms;
    ms.reserve(static_cast<size_t>(spec.repeat));
    for (int i = 0; i < spec.repeat; ++i) {
        const auto t0 = std::chrono::steady_clock::now();
        const SolveResult r = solve(model, opt);
        const auto t1 = std::chrono::steady_clock::now();
        if (r.singular) { std::printf("solve singular: %s\n", r.diagnostic.c_str()); return 2; }
        cs += checksum(r);
        ms.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count());
    }
    std::sort(ms.begin(), ms.end());
    double mean = 0.0;
    for (double v : ms) mean += v;
    mean /= static_cast<double>(ms.size());
    std::printf("solveMs median=%.3f min=%.3f max=%.3f mean=%.3f checksum=%.12g\n",
                ms[ms.size() / 2], ms.front(), ms.back(), mean, cs);
    return 0;
}
