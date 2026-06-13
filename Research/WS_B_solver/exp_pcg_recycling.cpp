// exp_pcg_recycling.cpp -- research-only multi-RHS recycled PCG experiment.
//
// The experiment targets game/interactive FEM workloads where many linear solves
// share one stiffness matrix and receive correlated load vectors across frames or
// loadcases. Recycling is used only as an initial guess; final accuracy is still
// checked against assembled K_ff and LDLT oracle solutions.

#include "research_common.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

using namespace research;

namespace {

struct Args {
    std::string preset = "small"; // small | xxl | mega
    int rhsCount = 16;
    int basisMax = 6;
    int pcgMaxIter = 5000;
    double pcgTol = 1e-10;
};

Args parseArgs(int argc, char** argv) {
    Args a;
    for (int i = 1; i < argc; ++i) {
        auto next = [&]() -> const char* { return (i + 1 < argc) ? argv[++i] : ""; };
        if      (!std::strcmp(argv[i], "--preset")) a.preset = next();
        else if (!std::strcmp(argv[i], "--rhs")) a.rhsCount = std::max(1, std::atoi(next()));
        else if (!std::strcmp(argv[i], "--basisMax")) a.basisMax = std::max(0, std::atoi(next()));
        else if (!std::strcmp(argv[i], "--pcgMaxIter")) a.pcgMaxIter = std::max(1, std::atoi(next()));
        else if (!std::strcmp(argv[i], "--pcgTol")) a.pcgTol = std::atof(next());
    }
    return a;
}

FrameModel makeCase(const std::string& preset) {
    if (preset == "xxl")  return makeTower(12, 9, 24);
    if (preset == "mega") return makeTower(18, 14, 36);
    return makeTower(5, 4, 8);
}

double relNorm(const VecX& a, const VecX& b) {
    return static_cast<double>((a - b).norm() / std::max<real>(real(1e-300), b.norm()));
}

struct ElementBlock {
    int n = 0;
    std::array<int, 12> rid{};
    std::array<real, 144> k{};
};

struct ElementOperator {
    int nf = 0;
    std::vector<ElementBlock> blocks;
    VecX diag;

    bool build(const FrameModel& model, const PreparedSystem::Impl& S, std::string& why) {
        nf = S.nf;
        diag = VecX::Zero(nf);
        blocks.clear();
        blocks.reserve(model.members.size());

        SolveOptions opts;
        for (size_t e = 0; e < model.members.size(); ++e) {
            if (!model.members[e].active) continue;
            std::vector<Triplet> trips;
            if (!memberGlobalK(model, static_cast<int>(e), opts, trips, why)) return false;

            int gd[12];
            memberGlobalDofs(model, static_cast<int>(e), gd);

            std::vector<int> keep;
            std::vector<int> rid;
            keep.reserve(12);
            rid.reserve(12);
            for (int i = 0; i < 12; ++i) {
                const int f = S.fmap[static_cast<size_t>(gd[i])];
                if (f >= 0) {
                    keep.push_back(i);
                    rid.push_back(f);
                }
            }
            const int n = static_cast<int>(rid.size());
            if (n == 0) continue;

            auto localIndex = [&](int g) -> int {
                for (int i = 0; i < 12; ++i) if (gd[i] == g) return i;
                return -1;
            };

            MatX k = MatX::Zero(n, n);
            for (const Triplet& t : trips) {
                const int lr = localIndex(static_cast<int>(t.row()));
                const int lc = localIndex(static_cast<int>(t.col()));
                if (lr < 0 || lc < 0) continue;
                const auto ir = std::find(keep.begin(), keep.end(), lr);
                const auto ic = std::find(keep.begin(), keep.end(), lc);
                if (ir == keep.end() || ic == keep.end()) continue;
                k(static_cast<int>(ir - keep.begin()), static_cast<int>(ic - keep.begin())) += t.value();
            }

            ElementBlock b;
            b.n = n;
            for (int i = 0; i < n; ++i) {
                b.rid[static_cast<size_t>(i)] = rid[static_cast<size_t>(i)];
                diag(rid[static_cast<size_t>(i)]) += k(i, i);
            }
            for (int r = 0; r < n; ++r)
                for (int c = 0; c < n; ++c)
                    b.k[static_cast<size_t>(r * 12 + c)] = k(r, c);
            blocks.push_back(std::move(b));
        }
        return true;
    }

    void apply(const VecX& x, VecX& y) const {
        y.setZero(nf);
        for (const ElementBlock& b : blocks) {
            real xe[12];
            for (int i = 0; i < b.n; ++i) xe[i] = x(b.rid[static_cast<size_t>(i)]);
            for (int r = 0; r < b.n; ++r) {
                const real* row = &b.k[static_cast<size_t>(r * 12)];
                real acc = 0;
                for (int c = 0; c < b.n; ++c) acc += row[c] * xe[c];
                y(b.rid[static_cast<size_t>(r)]) += acc;
            }
        }
    }

    void jacobi(const VecX& r, VecX& z) const {
        z.resize(nf);
        for (int i = 0; i < nf; ++i)
            z(i) = (diag(i) != real(0)) ? r(i) / diag(i) : r(i);
    }
};

struct PcgStats {
    int iters = 0;
    double relResidual = 0.0;
    double ms = 0.0;
};

PcgStats pcg(const ElementOperator& A, const VecX& b, VecX& x, int maxIter, double tol) {
    PcgStats st;
    VecX Ax(b.size());
    A.apply(x, Ax);
    VecX r = b - Ax;
    const real bnorm = std::max<real>(real(1e-300), b.norm());
    st.relResidual = static_cast<double>(r.norm() / bnorm);
    if (st.relResidual <= tol) return st;

    VecX z;
    A.jacobi(r, z);
    VecX p = z;
    real rzOld = r.dot(z);
    Timer t;
    if (!(rzOld > real(0)) || !std::isfinite(static_cast<double>(rzOld))) {
        st.ms = t.ms();
        return st;
    }

    VecX Ap(b.size());
    for (int it = 0; it < maxIter; ++it) {
        A.apply(p, Ap);
        const real denom = p.dot(Ap);
        if (!(denom > real(0)) || !std::isfinite(static_cast<double>(denom))) break;
        const real alpha = rzOld / denom;
        x += alpha * p;
        r -= alpha * Ap;
        st.iters = it + 1;
        st.relResidual = static_cast<double>(r.norm() / bnorm);
        if (st.relResidual <= tol) break;
        A.jacobi(r, z);
        const real rzNew = r.dot(z);
        if (!(rzNew > real(0)) || !std::isfinite(static_cast<double>(rzNew))) break;
        p = z + (rzNew / rzOld) * p;
        rzOld = rzNew;
    }
    st.ms = t.ms();
    return st;
}

struct RecycleBasis {
    int maxSize = 0;
    std::vector<VecX> v;
    std::vector<VecX> av;
    double projectMs = 0.0;
    double addMs = 0.0;
    int accepted = 0;
    int rejected = 0;

    explicit RecycleBasis(int maxBasis) : maxSize(maxBasis) {}

    VecX initialGuess(const ElementOperator& A, const VecX& b) {
        Timer t;
        VecX x = VecX::Zero(A.nf);
        for (size_t i = 0; i < v.size(); ++i)
            x.noalias() += v[i].dot(b) * v[i];
        projectMs += t.ms();
        return x;
    }

    void add(const ElementOperator& A, const VecX& x) {
        if (maxSize <= 0) return;
        Timer t;

        VecX q = x;
        VecX Aq(A.nf);
        A.apply(q, Aq);
        for (size_t i = 0; i < v.size(); ++i) {
            const real c = v[i].dot(Aq);
            q.noalias() -= c * v[i];
            Aq.noalias() -= c * av[i];
        }

        const real normA2 = q.dot(Aq);
        const double scaleRef = std::max(1e-300, static_cast<double>(x.norm()));
        if (!(normA2 > real(0)) || !std::isfinite(static_cast<double>(normA2)) ||
            static_cast<double>(q.norm()) <= 1e-12 * scaleRef) {
            ++rejected;
            addMs += t.ms();
            return;
        }

        const real inv = real(1) / std::sqrt(normA2);
        q *= inv;
        Aq *= inv;
        if (static_cast<int>(v.size()) >= maxSize) {
            v.erase(v.begin());
            av.erase(av.begin());
        }
        v.push_back(std::move(q));
        av.push_back(std::move(Aq));
        ++accepted;
        addMs += t.ms();
    }
};

std::vector<VecX> makeRhsSequence(const FrameModel& model, const PreparedSystem::Impl& S, int count) {
    const VecX base = reduceVec(nodalLoadVector(model), S.fmap, S.nf);
    std::vector<VecX> rhs;
    rhs.reserve(static_cast<size_t>(count));

    real minX = std::numeric_limits<real>::max(), maxX = -std::numeric_limits<real>::max();
    real minY = std::numeric_limits<real>::max(), maxY = -std::numeric_limits<real>::max();
    real minZ = std::numeric_limits<real>::max(), maxZ = -std::numeric_limits<real>::max();
    for (const Node& n : model.nodes) {
        minX = std::min(minX, n.pos.x); maxX = std::max(maxX, n.pos.x);
        minY = std::min(minY, n.pos.y); maxY = std::max(maxY, n.pos.y);
        minZ = std::min(minZ, n.pos.z); maxZ = std::max(maxZ, n.pos.z);
    }
    const real spanX = std::max<real>(real(1), maxX - minX);
    const real spanY = std::max<real>(real(1), maxY - minY);
    const real spanZ = std::max<real>(real(1), maxZ - minZ);

    VecX windX = VecX::Zero(S.nf);
    VecX windY = VecX::Zero(S.nf);
    VecX liveZ = VecX::Zero(S.nf);
    VecX torsion = VecX::Zero(S.nf);

    for (size_t ni = 0; ni < model.nodes.size(); ++ni) {
        const Node& node = model.nodes[ni];
        const real xn = (node.pos.x - minX) / spanX - real(0.5);
        const real yn = (node.pos.y - minY) / spanY - real(0.5);
        const real zn = (node.pos.z - minZ) / spanZ;
        for (int d = 0; d < 6; ++d) {
            const int f = S.fmap[static_cast<size_t>(gdof(static_cast<int>(ni), d))];
            if (f < 0) continue;
            if (d == Ux) windX(f) = real(1600.0) * zn * (real(1) + real(0.15) * yn);
            if (d == Uy) windY(f) = real(1200.0) * zn * (real(1) - real(0.10) * xn);
            if (d == Uz) liveZ(f) = real(-4500.0) * (real(0.25) + zn) * (real(1) + real(0.10) * xn);
            if (d == Rz) torsion(f) = real(2.5e5) * zn * (xn - yn);
        }
    }

    for (int k = 0; k < count; ++k) {
        const double t = static_cast<double>(k);
        VecX b = (1.0 + 0.035 * std::sin(0.37 * t)) * base;
        b.noalias() += (0.65 + 0.35 * std::sin(0.23 * t + 0.4)) * windX;
        b.noalias() += (0.45 + 0.25 * std::cos(0.31 * t)) * windY;
        b.noalias() += (0.20 + 0.15 * std::sin(0.17 * t + 1.1)) * liveZ;
        b.noalias() += (0.25 * std::sin(0.41 * t)) * torsion;
        rhs.push_back(std::move(b));
    }
    return rhs;
}

} // namespace

int main(int argc, char** argv) {
    const Args args = parseArgs(argc, argv);
    FrameModel model = makeCase(args.preset);
    assertNodalOnly(model, "exp_pcg_recycling");

    Timer tPrep;
    PreparedSystem ps = assembleAndFactor(model);
    const double prepMs = tPrep.ms();
    const PreparedSystem::Impl& S = *ps.impl;
    if (S.singular) {
        std::printf("[recycle] status=singular diagnostic=\"%s\"\n", S.diagnostic.c_str());
        return 2;
    }

    const SpMat Kff = research::reduceFF(S.K, S.fmap, S.nf);
    std::string why;
    ElementOperator A;
    Timer tBuild;
    if (!A.build(model, S, why)) {
        std::printf("[recycle] status=build_failed why=\"%s\"\n", why.c_str());
        return 2;
    }
    const double opBuildMs = tBuild.ms();

    VecX xProbe = VecX::Zero(S.nf);
    for (int i = 0; i < xProbe.size(); ++i) xProbe(i) = std::sin(0.003 * static_cast<double>(i + 1));
    VecX yMf(S.nf);
    A.apply(xProbe, yMf);
    const double applyRel = relNorm(yMf, Kff * xProbe);

    const std::vector<VecX> rhs = makeRhsSequence(model, S, args.rhsCount);
    RecycleBasis basis(args.basisMax);

    int baseIters = 0;
    int recIters = 0;
    int recItersSkip1 = 0;
    double baseMs = 0.0;
    double recMs = 0.0;
    double maxRelXBase = 0.0;
    double maxRelXRec = 0.0;
    double maxTrueBase = 0.0;
    double maxTrueRec = 0.0;
    int validSkip = 0;

    for (int i = 0; i < static_cast<int>(rhs.size()); ++i) {
        const VecX& b = rhs[static_cast<size_t>(i)];
        const VecX xRef = S.ldlt.solve(b);

        VecX xBase = VecX::Zero(S.nf);
        const PcgStats bs = pcg(A, b, xBase, args.pcgMaxIter, args.pcgTol);
        baseIters += bs.iters;
        baseMs += bs.ms;
        maxRelXBase = std::max(maxRelXBase, relNorm(xBase, xRef));
        maxTrueBase = std::max(maxTrueBase, relNorm(Kff * xBase, b));

        VecX xRec = basis.initialGuess(A, b);
        const PcgStats rs = pcg(A, b, xRec, args.pcgMaxIter, args.pcgTol);
        recIters += rs.iters;
        recMs += rs.ms;
        if (i > 0) {
            recItersSkip1 += rs.iters;
            ++validSkip;
        }
        maxRelXRec = std::max(maxRelXRec, relNorm(xRec, xRef));
        maxTrueRec = std::max(maxTrueRec, relNorm(Kff * xRec, b));
        basis.add(A, xRec);
    }

    const double n = static_cast<double>(rhs.size());
    const double skipN = static_cast<double>(std::max(1, validSkip));
    const double baseItersAvg = static_cast<double>(baseIters) / n;
    const double recItersAvg = static_cast<double>(recIters) / n;
    const double recItersAvgSkip1 = static_cast<double>(recItersSkip1) / skipN;
    const double baseMsAvg = baseMs / n;
    const double recMsAvg = recMs / n;
    const double projectMsAvg = basis.projectMs / n;
    const bool ok = applyRel <= 1e-10 && maxRelXBase <= 1e-6 && maxRelXRec <= 1e-6 &&
                    maxTrueBase <= 1e-8 && maxTrueRec <= 1e-8;

    std::printf("[recycle] preset=%s rhs=%d basisMax=%d basisAccepted=%d basisRejected=%d nf=%d blocks=%zu prepMs=%.3f opBuildMs=%.3f applyRel=%.3e baseItersAvg=%.2f recItersAvg=%.2f recItersSkip1Avg=%.2f speedupIters=%.3f baseMsAvg=%.3f recMsAvg=%.3f projectMsAvg=%.6f basisAddMs=%.3f speedupMs=%.3f maxRelXBase=%.3e maxRelXRec=%.3e maxTrueBase=%.3e maxTrueRec=%.3e\n",
                args.preset.c_str(), static_cast<int>(rhs.size()), args.basisMax,
                basis.accepted, basis.rejected, S.nf, A.blocks.size(), prepMs, opBuildMs,
                applyRel, baseItersAvg, recItersAvg, recItersAvgSkip1,
                baseItersAvg / std::max(1e-300, recItersAvg),
                baseMsAvg, recMsAvg, projectMsAvg, basis.addMs,
                baseMsAvg / std::max(1e-300, recMsAvg),
                maxRelXBase, maxRelXRec, maxTrueBase, maxTrueRec);
    return ok ? 0 : 1;
}
