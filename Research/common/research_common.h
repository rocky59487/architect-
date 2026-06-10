#pragma once
//
// Research-only shared helpers. NOT part of the FrameCore engine.
//
//  * makeTower / squareSection / refVecFor are COPIED from Standalone/frame_perf.cpp
//    (kept byte-equivalent in behaviour) so experiments run on the same benchmark
//    family without touching engine files.
//  * Experiments may include PRIVATE engine headers (PreparedSystemImpl.h, ...) for
//    introspection. That is a research privilege only — production features go through
//    the public API + the adoption ladder (see docs/KARAMBA3D_ROADMAP.md).
//  * Scope guard: helpers assume NODAL loads only (no UDL / shell pressure / self
//    weight), which all experiment models respect. assertNodalOnly() enforces it.
//
#include "FrameCore/FrameSolver.h"
#include "FrameCore/Section.h"
#include "PreparedSystemImpl.h"   // private: SpMat K, fmap, nf, ldlt, elems
#include "BeamColumnElement.h"    // private: per-member element stiffness extraction

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>

namespace research {

using namespace frame;

// ---------------------------------------------------------------- timing
struct Timer {
    std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
    void   reset() { t0 = std::chrono::steady_clock::now(); }
    double ms() const {
        return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    }
};

// ---------------------------------------------------------------- tower benchmark model
// Copied from Standalone/frame_perf.cpp (research copy; engine file untouched).
inline int towerNodeId(int i, int j, int k, int nx, int ny) {
    return k * (nx + 1) * (ny + 1) + j * (nx + 1) + i;
}

inline Section squareSection(real side) {
    Section s = Section::Rectangular(side, side);
    s.J = 0.1406 * side * side * side * side;
    return s;
}

inline Vec3 refVecFor(const Vec3& pi, const Vec3& pj) {
    Vec3 x = pj - pi;
    const real n = norm(x);
    if (n > 0) x = x * (1.0 / n);
    return (std::abs(x.z) < 0.92) ? Vec3(0, 0, 1) : Vec3(1, 0, 0);
}

inline FrameModel makeTower(int nx, int ny, int stories) {
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
                Node n(towerNodeId(i, j, k, nx, ny), i * sx + driftX, j * sy + driftY, k * h);
                if (k == 0) n.fixAll();
                m.nodes.push_back(n);
            }
        }
    }

    auto nodePos = [&](int id) -> Vec3 { return m.nodes[static_cast<size_t>(id)].pos; };
    auto addMember = [&](int i, int j, int sec) {
        Member mem(static_cast<int>(m.members.size()), i, j, 0, sec);
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
                addMember(towerNodeId(i, j, k, nx, ny), towerNodeId(i, j, k + 1, nx, ny), 0);

    for (int k = 1; k <= stories; ++k) {
        for (int j = 0; j <= ny; ++j)
            for (int i = 0; i < nx; ++i)
                addMember(towerNodeId(i, j, k, nx, ny), towerNodeId(i + 1, j, k, nx, ny), 1);
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i <= nx; ++i)
                addMember(towerNodeId(i, j, k, nx, ny), towerNodeId(i, j + 1, k, nx, ny), 1);
    }

    for (int k = 0; k < stories; ++k) {
        for (int i = 0; i < nx; ++i) {
            for (int face = 0; face < 2; ++face) {
                const int jf = face ? ny : 0;
                const int a = ((i + k) % 2 == 0) ? towerNodeId(i, jf, k, nx, ny)         : towerNodeId(i + 1, jf, k, nx, ny);
                const int b = ((i + k) % 2 == 0) ? towerNodeId(i + 1, jf, k + 1, nx, ny) : towerNodeId(i, jf, k + 1, nx, ny);
                addMember(a, b, 2);
            }
        }
        for (int j = 0; j < ny; ++j) {
            for (int face = 0; face < 2; ++face) {
                const int ifc = face ? nx : 0;
                const int a = ((j + k) % 2 == 0) ? towerNodeId(ifc, j, k, nx, ny)         : towerNodeId(ifc, j + 1, k, nx, ny);
                const int b = ((j + k) % 2 == 0) ? towerNodeId(ifc, j + 1, k + 1, nx, ny) : towerNodeId(ifc, j, k + 1, nx, ny);
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
                l.node = towerNodeId(i, j, k, nx, ny);
                l.comp[Ux] = 900.0 * sf * (1.0 + 0.08 * (j - ny / 2.0));
                l.comp[Uy] = -350.0 * sf * (1.0 + 0.05 * (i - nx / 2.0));
                l.comp[Uz] = -38000.0 * trib;
                m.nodalLoads.push_back(l);
            }
        }
    }
    NodalLoad t1; t1.node = towerNodeId(nx, ny, stories, nx, ny); t1.comp[Rz] = 3.0e6;  m.nodalLoads.push_back(t1);
    NodalLoad t2; t2.node = towerNodeId(0, 0, stories, nx, ny);   t2.comp[Ry] = -2.0e6; m.nodalLoads.push_back(t2);
    return m;
}

// ---------------------------------------------------------------- simple cantilever column
// Vertical column along +Z, base fixAll, meshed into nElem beams. Tip loads:
// axialP > 0 pushes DOWN (compression-positive axial force in every element),
// lateralH applies +X at the tip. Used by the buckling / P-Delta experiments.
inline FrameModel makeColumn(int nElem, real L, real side, real axialP, real lateralH) {
    FrameModel m;
    m.materials.emplace_back(200000.0, 76923.07692307692, 7850.0);
    m.sections.push_back(squareSection(side));
    for (int k = 0; k <= nElem; ++k) {
        Node n(k, 0.0, 0.0, L * k / nElem);
        if (k == 0) n.fixAll();
        m.nodes.push_back(n);
    }
    for (int k = 0; k < nElem; ++k) {
        Member mem(k, k, k + 1, 0, 0);
        mem.refVec = Vec3(1, 0, 0);     // axis ~ +Z -> refVec must not be parallel
        m.members.push_back(mem);
    }
    NodalLoad tip; tip.node = nElem;
    tip.comp[Uz] = -axialP;
    tip.comp[Ux] = lateralH;
    m.nodalLoads.push_back(tip);
    return m;
}

// ---------------------------------------------------------------- reduced-system helpers
// Free-free block of a global sparse matrix using the prepared free-DOF map.
inline SpMat reduceFF(const SpMat& A, const std::vector<int>& fmap, int nf, real scale = 1.0) {
    std::vector<Triplet> t;
    t.reserve(static_cast<size_t>(A.nonZeros()));
    for (int c = 0; c < A.outerSize(); ++c)
        for (SpMat::InnerIterator it(A, c); it; ++it) {
            const int r = static_cast<int>(it.row());
            if (fmap[static_cast<size_t>(r)] >= 0 && fmap[static_cast<size_t>(c)] >= 0)
                t.emplace_back(fmap[static_cast<size_t>(r)], fmap[static_cast<size_t>(c)], scale * it.value());
        }
    SpMat R(nf, nf);
    R.setFromTriplets(t.begin(), t.end());
    R.makeCompressed();
    return R;
}

// Global nodal-load vector. Experiments use nodal loads ONLY (no UDL / pressure /
// self-weight), matching what solveLoad would bake for these models.
inline void assertNodalOnly(const FrameModel& m, const char* who) {
    if (!m.memberUDLs.empty() || !m.shellPressures.empty()) {
        std::fprintf(stderr, "[%s] FATAL: experiment helpers support nodal loads only\n", who);
        std::exit(9);
    }
}

inline VecX nodalLoadVector(const FrameModel& m) {
    VecX F = VecX::Zero(m.dofCount());
    for (const auto& L : m.nodalLoads) {
        const int ni = m.nodeIndex(L.node);
        if (ni < 0) continue;
        for (int d = 0; d < 6; ++d) F(gdof(ni, d)) += L.comp[static_cast<size_t>(d)];
    }
    return F;
}

inline VecX reduceVec(const VecX& F, const std::vector<int>& fmap, int nf) {
    VecX f = VecX::Zero(nf);
    for (size_t g = 0; g < fmap.size(); ++g)
        if (fmap[g] >= 0) f(fmap[g]) = F(static_cast<int>(g));
    return f;
}

inline VecX scatterVec(const VecX& uf, const std::vector<int>& fmap, int N) {
    VecX u = VecX::Zero(N);
    for (size_t g = 0; g < fmap.size(); ++g)
        if (fmap[g] >= 0) u(static_cast<int>(g)) = uf(fmap[g]);
    return u;
}

// ---------------------------------------------------------------- single-member stiffness
// Global stiffness triplets of ONE member, exactly as the engine assembles it
// (BeamColumnElement::prepare + assemble on the given model). The model passed in
// must have the member ACTIVE-compatible data (geometry/material/section); the
// member's own `active` flag is irrelevant here because we construct the element
// directly. Returns false (with `why`) if the element is ill-posed.
inline bool memberGlobalK(const FrameModel& m, int memberIdx, const SolveOptions& opts,
                          std::vector<Triplet>& trips, std::string& why) {
    BeamColumnElement el(memberIdx);
    if (!el.prepare(m, opts, why)) return false;
    el.assemble(trips);
    return true;
}

// The 12 global DOF indices of a member (node_i 6, node_j 6).
inline void memberGlobalDofs(const FrameModel& m, int memberIdx, int out[12]) {
    const Member& mem = m.members[static_cast<size_t>(memberIdx)];
    const int ni = m.nodeIndex(mem.i), nj = m.nodeIndex(mem.j);
    for (int d = 0; d < 6; ++d) { out[d] = gdof(ni, d); out[6 + d] = gdof(nj, d); }
}

// ---------------------------------------------------------------- comparisons
inline double relMaxDiff(const VecX& a, const std::vector<real>& b) {
    double maxd = 0, maxref = 0;
    const size_t n = std::min(static_cast<size_t>(a.size()), b.size());
    for (size_t i = 0; i < n; ++i) {
        maxd   = std::max(maxd, std::abs(static_cast<double>(a(static_cast<int>(i)) - b[i])));
        maxref = std::max(maxref, std::abs(static_cast<double>(b[i])));
    }
    return maxref > 0 ? maxd / maxref : maxd;
}

inline double relMaxDiffV(const std::vector<real>& a, const std::vector<real>& b) {
    double maxd = 0, maxref = 0;
    const size_t n = std::min(a.size(), b.size());
    for (size_t i = 0; i < n; ++i) {
        maxd   = std::max(maxd, std::abs(static_cast<double>(a[i] - b[i])));
        maxref = std::max(maxref, std::abs(static_cast<double>(b[i])));
    }
    return maxref > 0 ? maxd / maxref : maxd;
}

} // namespace research
