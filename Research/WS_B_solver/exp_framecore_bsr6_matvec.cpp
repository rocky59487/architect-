// exp_framecore_bsr6_matvec.cpp -- FrameCore K_ff to 6x6 BSR matvec prototype.
//
// This differs from exp_bsr6_matvec.cpp, which is a synthetic kernel benchmark.
// Here the sparse matrix is the engine-produced reduced stiffness matrix, so the
// result answers whether FrameCore's real sparsity pattern benefits from BSR.

#include "research_common.h"

#include <algorithm>
#include <array>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>

using namespace research;

namespace {

struct Args {
    std::string preset = "small"; // small | xxl | mega
    int repeat = 50;
};

Args parseArgs(int argc, char** argv) {
    Args a;
    for (int i = 1; i < argc; ++i) {
        auto next = [&]() -> const char* { return (i + 1 < argc) ? argv[++i] : ""; };
        if      (!std::strcmp(argv[i], "--preset")) a.preset = next();
        else if (!std::strcmp(argv[i], "--repeat")) a.repeat = std::max(1, std::atoi(next()));
    }
    return a;
}

FrameModel makeCase(const std::string& preset) {
    if (preset == "xxl")  return makeTower(12, 9, 24);
    if (preset == "mega") return makeTower(18, 14, 36);
    return makeTower(5, 4, 8);
}

double checksum(const VecX& v) {
    double s = 0.0;
    const int step = std::max(1, static_cast<int>(v.size()) / 4096);
    for (int i = 0; i < v.size(); i += step) s += std::abs(static_cast<double>(v(i)));
    return s;
}

double relNorm(const VecX& a, const VecX& b) {
    return static_cast<double>((a - b).norm() / std::max<real>(real(1e-300), b.norm()));
}

struct Block6 {
    int col = 0;
    std::array<real, 36> a{};
};

struct Bsr6 {
    int nBlocks = 0;
    std::vector<int> rowPtr;
    std::vector<Block6> blocks;

    bool build(const SpMat& K) {
        if (K.rows() != K.cols() || K.rows() % 6 != 0) return false;
        nBlocks = static_cast<int>(K.rows()) / 6;
        std::vector<std::unordered_map<int, std::array<real, 36>>> rows(static_cast<size_t>(nBlocks));
        for (int c = 0; c < K.outerSize(); ++c) {
            for (SpMat::InnerIterator it(K, c); it; ++it) {
                const int r = static_cast<int>(it.row());
                const int br = r / 6, bc = c / 6;
                const int lr = r - 6 * br, lc = c - 6 * bc;
                rows[static_cast<size_t>(br)][bc][static_cast<size_t>(lr * 6 + lc)] += it.value();
            }
        }
        rowPtr.assign(static_cast<size_t>(nBlocks + 1), 0);
        for (int br = 0; br < nBlocks; ++br) {
            rowPtr[static_cast<size_t>(br + 1)] =
                rowPtr[static_cast<size_t>(br)] + static_cast<int>(rows[static_cast<size_t>(br)].size());
        }
        blocks.resize(static_cast<size_t>(rowPtr.back()));
        for (int br = 0; br < nBlocks; ++br) {
            std::vector<int> cols;
            cols.reserve(rows[static_cast<size_t>(br)].size());
            for (const auto& kv : rows[static_cast<size_t>(br)]) cols.push_back(kv.first);
            std::sort(cols.begin(), cols.end());
            int out = rowPtr[static_cast<size_t>(br)];
            for (int bc : cols) {
                blocks[static_cast<size_t>(out)].col = bc;
                blocks[static_cast<size_t>(out)].a = rows[static_cast<size_t>(br)][bc];
                ++out;
            }
        }
        return true;
    }

    void apply(const VecX& x, VecX& y) const {
        y.setZero(nBlocks * 6);
        for (int br = 0; br < nBlocks; ++br) {
            real acc[6] = {0, 0, 0, 0, 0, 0};
            for (int bi = rowPtr[static_cast<size_t>(br)]; bi < rowPtr[static_cast<size_t>(br + 1)]; ++bi) {
                const Block6& b = blocks[static_cast<size_t>(bi)];
                const int xc = b.col * 6;
                for (int r = 0; r < 6; ++r) {
                    const real* row = &b.a[static_cast<size_t>(r * 6)];
                    acc[r] += row[0] * x(xc + 0) + row[1] * x(xc + 1) + row[2] * x(xc + 2)
                            + row[3] * x(xc + 3) + row[4] * x(xc + 4) + row[5] * x(xc + 5);
                }
            }
            for (int r = 0; r < 6; ++r) y(br * 6 + r) = acc[r];
        }
    }
};

} // namespace

int main(int argc, char** argv) {
    const Args args = parseArgs(argc, argv);
    FrameModel model = makeCase(args.preset);
    assertNodalOnly(model, "exp_framecore_bsr6_matvec");

    Timer tPrep;
    PreparedSystem ps = assembleAndFactor(model);
    const double prepMs = tPrep.ms();
    const PreparedSystem::Impl& S = *ps.impl;
    if (S.singular) {
        std::printf("[framecore_bsr6] status=singular diagnostic=\"%s\"\n", S.diagnostic.c_str());
        return 2;
    }
    const SpMat Kff = research::reduceFF(S.K, S.fmap, S.nf);
    Bsr6 bsr;
    Timer tBuild;
    if (!bsr.build(Kff)) {
        std::printf("[framecore_bsr6] status=unsupported nf=%d\n", S.nf);
        return 2;
    }
    const double buildMs = tBuild.ms();

    VecX x = VecX::Zero(S.nf);
    for (int i = 0; i < x.size(); ++i) x(i) = std::cos(0.0007 * static_cast<double>(i + 3));
    VecX ySparse = Kff * x;
    VecX yBsr(S.nf);
    bsr.apply(x, yBsr);
    const double applyRel = relNorm(yBsr, ySparse);

    Timer tSparse;
    double csSparse = 0.0;
    for (int r = 0; r < args.repeat; ++r) {
        ySparse.noalias() = Kff * x;
        csSparse += checksum(ySparse);
    }
    const double sparseMs = tSparse.ms() / args.repeat;

    Timer tBsr;
    double csBsr = 0.0;
    for (int r = 0; r < args.repeat; ++r) {
        bsr.apply(x, yBsr);
        csBsr += checksum(yBsr);
    }
    const double bsrMs = tBsr.ms() / args.repeat;

    std::printf("[framecore_bsr6] preset=%s nodes=%zu members=%zu nf=%d scalarNnz=%lld blockRows=%d blockNnz=%zu prepMs=%.3f buildMs=%.3f applyRel=%.3e sparseApplyMs=%.6f bsrApplyMs=%.6f speedup=%.3f checksumRel=%.3e\n",
                args.preset.c_str(), model.nodes.size(), model.members.size(), S.nf,
                static_cast<long long>(Kff.nonZeros()), bsr.nBlocks, bsr.blocks.size(),
                prepMs, buildMs, applyRel, sparseMs, bsrMs, sparseMs / std::max(1e-300, bsrMs),
                std::abs(csBsr - csSparse) / std::max(1e-300, std::abs(csSparse)));
    return (applyRel <= 1e-12) ? 0 : 1;
}
