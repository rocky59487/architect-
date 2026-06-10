// exp_size_opt.cpp — WS-H fully-stressed design (FSD) on the classic 10-bar truss
// (scratch, NOT engine). Stress-ratio resizing A <- max(Amin, A * |sigma|/sigma_allow).
//
// Classic benchmark (imperial, converted to the engine's N-mm units):
//   E = 1e7 psi, sigma_allow = +/-25000 psi, P = 100 kip down at nodes 2 and 4,
//   bays 360 in, rho = 0.1 lb/in^3, Amin = 0.1 in^2. Literature FSD/optimum weight
//   ~1593 lb [cross-checked against sources in the research doc, not hardcoded here].
//
// Members are engine beam-columns with moment joints (TrussPin is gated off in the
// engine): slender squares keep bending negligible — the max bending-stress share is
// printed as the honest truss-likeness measure.

#include "research_common.h"
#include <cstring>

using namespace research;

namespace {

constexpr real IN = 25.4;                  // mm per inch
constexpr real LB = 4.4482216152605;       // N per lbf
constexpr real PSI = 0.0068947572931684;   // MPa per psi

struct Bar { int i, j; };

FrameModel buildTruss(const std::vector<real>& A_mm2, const std::vector<Bar>& bars) {
    FrameModel m;
    m.materials.emplace_back(1.0e7 * PSI, 1.0e7 * PSI / 2.6, 0.0);   // G irrelevant-ish, >0
    for (real A : A_mm2) m.sections.push_back(squareSection(std::sqrt(A)));

    const real B = 360.0 * IN;   // bay
    auto N = [&](int id, real xin, real zin, bool support) {
        Node n(id, xin, 0.0, zin);
        if (support) { n.fixed[Ux] = n.fixed[Uy] = n.fixed[Uz] = true; }
        else n.fixed[Uy] = true;   // planar problem
        return n;
    };
    m.nodes = {
        N(1, 2 * B, B, false),   // idx0
        N(2, 2 * B, 0, false),   // idx1
        N(3, B,     B, false),   // idx2
        N(4, B,     0, false),   // idx3
        N(5, 0,     B, true),    // idx4
        N(6, 0,     0, true),    // idx5
    };
    for (size_t k = 0; k < bars.size(); ++k) {
        Member mem(static_cast<int>(k), bars[k].i, bars[k].j, 0, static_cast<int>(k));
        mem.refVec = Vec3(0, 1, 0);   // member axes in XZ plane -> +Y is never parallel
        m.members.push_back(mem);
    }
    const real P = 1.0e5 * LB;
    NodalLoad l2; l2.node = 2; l2.comp[Uz] = -P;
    NodalLoad l4; l4.node = 4; l4.comp[Uz] = -P;
    m.nodalLoads.push_back(l2);
    m.nodalLoads.push_back(l4);
    return m;
}

} // namespace

int main() {
    // classic bar list (node IDs): 1:5-3 2:3-1 3:6-4 4:4-2 5:3-4 6:1-2 7:5-4 8:6-3 9:3-2 10:4-1
    const std::vector<Bar> bars = {
        {5, 3}, {3, 1}, {6, 4}, {4, 2}, {3, 4}, {1, 2}, {5, 4}, {6, 3}, {3, 2}, {4, 1}
    };
    const real Amin = 0.1 * IN * IN;            // 0.1 in^2
    const real sigA = 25000.0 * PSI;            // MPa
    std::vector<real> A(10, 10.0 * IN * IN);    // start: 10 in^2 each

    real maxBendShare = 0;
    for (int it = 1; it <= 40; ++it) {
        FrameModel m = buildTruss(A, bars);
        SolveResult r = solve(m);
        if (r.singular) { std::printf("[fsd] FATAL singular at iter %d\n", it); return 2; }

        real maxDev = 0, weightLb = 0;
        maxBendShare = 0;
        std::vector<real> An(10);
        for (int k = 0; k < 10; ++k) {
            const auto& f = r.memberForces[static_cast<size_t>(k)];
            const real Nax = std::max(std::abs(f.endI.N), std::abs(f.endJ.N));
            const real sig = Nax / A[static_cast<size_t>(k)];
            An[static_cast<size_t>(k)] = std::max(Amin, A[static_cast<size_t>(k)] * sig / sigA);

            const Section& s = m.sections[static_cast<size_t>(k)];
            const real sb = std::abs(f.endI.My) / std::max<real>(1e-300, s.Wy())
                          + std::abs(f.endI.Mz) / std::max<real>(1e-300, s.Wz());
            if (sig > 0) maxBendShare = std::max(maxBendShare, sb / sig);

            const Vec3 pi = m.nodes[static_cast<size_t>(m.nodeIndex(bars[static_cast<size_t>(k)].i))].pos;
            const Vec3 pj = m.nodes[static_cast<size_t>(m.nodeIndex(bars[static_cast<size_t>(k)].j))].pos;
            const real L = norm(pj - pi);
            weightLb += 0.1 * (A[static_cast<size_t>(k)] / (IN * IN)) * (L / IN);
            if (An[static_cast<size_t>(k)] > Amin * 1.001)
                maxDev = std::max(maxDev, std::abs(sig / sigA - 1.0));
        }
        std::printf("[fsd] iter=%d weightLb=%.4f maxDevFromFS=%.3e maxBendShare=%.3e\n",
                    it, static_cast<double>(weightLb), static_cast<double>(maxDev),
                    static_cast<double>(maxBendShare));
        A = An;
        if (maxDev < 1e-10 && it > 3) break;
    }

    std::printf("[fsd-final] areas_in2 =");
    for (int k = 0; k < 10; ++k) std::printf(" %.4f", static_cast<double>(A[static_cast<size_t>(k)] / (IN * IN)));
    std::printf("\n[done]\n");
    return 0;
}
