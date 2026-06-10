// Connected-component analysis + debris mass properties (collapse stage 3b).
// Pure POD math on the model graph — no Eigen, no solver state. The mass geometry is
// closed-form: a member is a uniform slender rod, a shell facet is a uniform thin lamina
// split into two triangles (the triangle second-moment formula is EXACT for polygons, so
// the only modelling approximation is "rod"/"lamina", not the integration).
#include "FrameCore/Connectivity.h"

#include <algorithm>
#include <map>

namespace frame {
namespace {

// Union-find over node INDICES (ids can be sparse; indices are dense).
struct UnionFind {
    std::vector<int> p;
    explicit UnionFind(size_t n) : p(n) { for (size_t i = 0; i < n; ++i) p[i] = (int)i; }
    int find(int a) { while (p[a] != a) { p[a] = p[p[a]]; a = p[a]; } return a; }
    void unite(int a, int b) { a = find(a); b = find(b); if (a != b) p[b] = a; }
};

// Symmetric 3x3 second-moment accumulator S = integral(r r^T dm) about the global origin.
struct Sym3 {
    real xx = 0, yy = 0, zz = 0, xy = 0, xz = 0, yz = 0;
    void addScaledOuter(const Vec3& v, real w) {           // S += w * v v^T
        xx += w * v.x * v.x; yy += w * v.y * v.y; zz += w * v.z * v.z;
        xy += w * v.x * v.y; xz += w * v.x * v.z; yz += w * v.y * v.z;
    }
    void add(const Sym3& o) { xx += o.xx; yy += o.yy; zz += o.zz; xy += o.xy; xz += o.xz; yz += o.yz; }
};

// Running mass aggregate for one cluster: total mass, first moment, S about the origin.
struct MassAgg {
    real mass = 0;
    Vec3 first;          // integral(r dm) = mass * com
    Sym3 S;              // integral(r r^T dm), origin

    // Uniform slender rod between a and b: S_com = (m L^2/12) a_hat a_hat^T, shifted to the
    // origin by the parallel-axis term m c c^T.
    void addRod(const Vec3& a, const Vec3& b, real m) {
        if (m <= 0) return;
        const Vec3 d = b - a;
        const real L = norm(d);
        const Vec3 c = (a + b) * 0.5;
        mass += m;
        first = first + c * m;
        if (L > 0) S.addScaledOuter(Vec3(d.x / L, d.y / L, d.z / L), m * L * L / 12.0);
        S.addScaledOuter(c, m);
    }

    // Uniform thin triangular lamina P1P2P3 with total mass m (exact closed form):
    //   integral(r r^T dm) about the ORIGIN = (m/12) * (sum_i Pi Pi^T + 9 c c^T),  c = centroid.
    void addTriangle(const Vec3& P1, const Vec3& P2, const Vec3& P3, real m) {
        if (m <= 0) return;
        const Vec3 c { (P1.x + P2.x + P3.x) / 3.0, (P1.y + P2.y + P3.y) / 3.0, (P1.z + P2.z + P3.z) / 3.0 };
        mass += m;
        first = first + c * m;
        S.addScaledOuter(P1, m / 12.0);
        S.addScaledOuter(P2, m / 12.0);
        S.addScaledOuter(P3, m / 12.0);
        S.addScaledOuter(c, 9.0 * m / 12.0);
    }
};

real triangleArea(const Vec3& a, const Vec3& b, const Vec3& c) {
    return 0.5 * norm(cross(b - a, c - a));
}

}  // namespace

ConnectivityResult analyzeConnectivity(const FrameModel& model) {
    ConnectivityResult out;
    std::string why;
    if (!model.validate(why)) return out;   // valid stays false
    out.valid = true;

    const size_t nN = model.nodes.size();
    UnionFind uf(nN);
    std::vector<bool> attached(nN, false);

    for (const Member& mem : model.members) {
        if (!mem.active) continue;
        const int a = model.nodeIndex(mem.i), b = model.nodeIndex(mem.j);
        uf.unite(a, b);
        attached[a] = attached[b] = true;
    }
    for (const ShellQuad& sh : model.shells) {
        if (!sh.active) continue;
        int idx[4];
        for (int k = 0; k < 4; ++k) { idx[k] = model.nodeIndex(sh.n[k]); attached[idx[k]] = true; }
        for (int k = 0; k < 4; ++k) uf.unite(idx[k], idx[(k + 1) % 4]);
    }

    // Bare nodes: no active element. Any free DOF there is a zero-stiffness mechanism.
    for (size_t k = 0; k < nN; ++k) {
        if (attached[k]) continue;
        bool fullyFixed = true;
        for (int d = 0; d < DOF_PER_NODE; ++d) fullyFixed = fullyFixed && model.nodes[k].fixed[d];
        if (!fullyFixed) out.looseNodes.push_back(model.nodes[k].id);
    }
    std::sort(out.looseNodes.begin(), out.looseNodes.end());

    // Group attached nodes by component root; a component is grounded iff any node has any
    // fixed DOF. std::map keys are ordered by root index, but the FINAL ordering contract is
    // re-established below by sorting on the smallest contained node id.
    struct Comp { std::vector<int> nodeIdx; bool grounded = false; };
    std::map<int, Comp> comps;
    for (size_t k = 0; k < nN; ++k) {
        if (!attached[k]) continue;
        Comp& c = comps[uf.find((int)k)];
        c.nodeIdx.push_back((int)k);
        for (int d = 0; d < DOF_PER_NODE; ++d) c.grounded = c.grounded || model.nodes[k].fixed[d];
    }

    for (auto& kv : comps) {
        Comp& c = kv.second;
        if (c.grounded) { ++out.groundedComponents; continue; }

        FragmentCluster fc;
        MassAgg agg;
        for (int k : c.nodeIdx) fc.nodes.push_back(model.nodes[(size_t)k].id);
        std::sort(fc.nodes.begin(), fc.nodes.end());

        // Gather the cluster's elements FIRST and sort them by id BEFORE accumulating mass:
        // float addition is order-dependent, so accumulating in storage order would make the
        // output depend on the caller's vector ordering. Id-sorted accumulation keeps the
        // result bit-identical under any shuffle of model.members / model.shells.
        const int root = kv.first;
        std::vector<std::pair<int, size_t>> memIdx, shIdx;   // (id, storage index)
        for (size_t e = 0; e < model.members.size(); ++e) {
            const Member& mem = model.members[e];
            if (mem.active && uf.find(model.nodeIndex(mem.i)) == root) memIdx.push_back({ mem.id, e });
        }
        for (size_t s = 0; s < model.shells.size(); ++s) {
            const ShellQuad& sh = model.shells[s];
            if (sh.active && uf.find(model.nodeIndex(sh.n[0])) == root) shIdx.push_back({ sh.id, s });
        }
        std::sort(memIdx.begin(), memIdx.end());
        std::sort(shIdx.begin(), shIdx.end());

        for (const auto& [id, e] : memIdx) {
            const Member& mem = model.members[e];
            fc.members.push_back(id);
            real rho = 0, A = 0;
            if (mem.matIdx >= 0 && mem.matIdx < (int)model.materials.size()) rho = model.materials[(size_t)mem.matIdx].rho;
            if (mem.secIdx >= 0 && mem.secIdx < (int)model.sections.size())  A   = model.sections[(size_t)mem.secIdx].A;
            const Vec3& pi = model.nodes[(size_t)model.nodeIndex(mem.i)].pos;
            const Vec3& pj = model.nodes[(size_t)model.nodeIndex(mem.j)].pos;
            agg.addRod(pi, pj, rho * 1e-12 * A * norm(pj - pi));   // rho[kg/m^3] -> tonne/mm^3
        }
        for (const auto& [id, s] : shIdx) {
            const ShellQuad& sh = model.shells[s];
            fc.shells.push_back(id);
            real rho = 0;
            if (sh.matIdx >= 0 && sh.matIdx < (int)model.materials.size()) rho = model.materials[(size_t)sh.matIdx].rho;
            const real sigma = rho * 1e-12 * sh.t;                 // surface density (tonne/mm^2)
            const Vec3& P0 = model.nodes[(size_t)model.nodeIndex(sh.n[0])].pos;
            const Vec3& P1 = model.nodes[(size_t)model.nodeIndex(sh.n[1])].pos;
            const Vec3& P2 = model.nodes[(size_t)model.nodeIndex(sh.n[2])].pos;
            const Vec3& P3 = model.nodes[(size_t)model.nodeIndex(sh.n[3])].pos;
            agg.addTriangle(P0, P1, P2, sigma * triangleArea(P0, P1, P2));
            agg.addTriangle(P0, P2, P3, sigma * triangleArea(P0, P2, P3));
        }

        fc.mass = agg.mass;
        if (agg.mass > 0) {
            fc.com = agg.first * (1.0 / agg.mass);
            Sym3 S = agg.S;                                        // shift to the com (parallel axis)
            S.addScaledOuter(fc.com, -agg.mass);
            fc.inertia[0] = S.yy + S.zz;    // Ixx
            fc.inertia[1] = S.xx + S.zz;    // Iyy
            fc.inertia[2] = S.xx + S.yy;    // Izz
            fc.inertia[3] = -S.xy;          // Ixy (tensor matrix entry)
            fc.inertia[4] = -S.xz;          // Ixz
            fc.inertia[5] = -S.yz;          // Iyz
        } else {                                                   // massless (rho==0) chunk
            Vec3 sum;
            for (int k : c.nodeIdx) sum = sum + model.nodes[(size_t)k].pos;
            fc.com = sum * (1.0 / (real)c.nodeIdx.size());
        }
        out.detached.push_back(std::move(fc));
    }

    std::sort(out.detached.begin(), out.detached.end(),
              [](const FragmentCluster& a, const FragmentCluster& b) { return a.nodes.front() < b.nodes.front(); });
    return out;
}

}  // namespace frame
