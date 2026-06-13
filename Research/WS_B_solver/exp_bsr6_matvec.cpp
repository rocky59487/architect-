// exp_bsr6_matvec.cpp -- standalone 6x6 BSR-like matvec microbenchmark.
//
// Research-only prototype. This file does not include or link the FrameCore engine.
// It builds one deterministic 6-DOF-per-node sparse matrix, expands it to scalar
// CSR, and compares scalar CSR matvec against a 6x6 block-row matvec.

#include <Eigen/Core>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr int kBlockSize = 6;
constexpr int kBlockEntries = kBlockSize * kBlockSize;

volatile double g_sink = 0.0;

struct Options {
    int blockRows = 16384;
    int nx = 0;
    int ny = 0;
    int nz = 0;
    int reps = 30;
    int warmup = 3;
    bool help = false;
};

struct Grid {
    int nx = 0;
    int ny = 0;
    int nz = 0;
};

struct Bsr6Matrix {
    int blockRows = 0;
    std::vector<int> rowPtr;
    std::vector<int> colIdx;
    std::vector<std::array<double, kBlockEntries>> values;
};

struct CsrMatrix {
    int rows = 0;
    std::vector<int> rowPtr;
    std::vector<int> colIdx;
    std::vector<double> values;
};

struct Timing {
    double totalMs = 0.0;
    double guard = 0.0;
};

long long scalarNnzEquivalent(const Bsr6Matrix& bsr) {
    return static_cast<long long>(bsr.values.size()) * kBlockEntries;
}

double nowMs() {
    using Clock = std::chrono::steady_clock;
    return std::chrono::duration<double, std::milli>(Clock::now().time_since_epoch()).count();
}

int parseInt(const char* text, const char* name) {
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (!text[0] || *end != '\0' || value <= 0 || value > std::numeric_limits<int>::max()) {
        throw std::runtime_error(std::string("invalid positive integer for ") + name + ": " + text);
    }
    return static_cast<int>(value);
}

void printUsage(const char* exe) {
    std::fprintf(stderr,
        "usage: %s [--blocks N] [--nx N --ny N --nz N] [--reps N] [--warmup N]\n"
        "\n"
        "Defaults: --blocks 16384 --reps 30 --warmup 3\n"
        "If --nx/--ny/--nz are all supplied, block rows become nx*ny*nz.\n",
        exe);
}

Options parseArgs(int argc, char** argv) {
    Options opt;
    for (int i = 1; i < argc; ++i) {
        auto needValue = [&](const char* flag) -> const char* {
            if (i + 1 >= argc) throw std::runtime_error(std::string("missing value after ") + flag);
            return argv[++i];
        };

        if (std::strcmp(argv[i], "--help") == 0 || std::strcmp(argv[i], "-h") == 0) {
            opt.help = true;
        } else if (std::strcmp(argv[i], "--blocks") == 0) {
            opt.blockRows = parseInt(needValue(argv[i]), "--blocks");
        } else if (std::strcmp(argv[i], "--nx") == 0) {
            opt.nx = parseInt(needValue(argv[i]), "--nx");
        } else if (std::strcmp(argv[i], "--ny") == 0) {
            opt.ny = parseInt(needValue(argv[i]), "--ny");
        } else if (std::strcmp(argv[i], "--nz") == 0) {
            opt.nz = parseInt(needValue(argv[i]), "--nz");
        } else if (std::strcmp(argv[i], "--reps") == 0) {
            opt.reps = parseInt(needValue(argv[i]), "--reps");
        } else if (std::strcmp(argv[i], "--warmup") == 0) {
            opt.warmup = parseInt(needValue(argv[i]), "--warmup");
        } else {
            throw std::runtime_error(std::string("unknown argument: ") + argv[i]);
        }
    }

    const bool anyGrid = opt.nx || opt.ny || opt.nz;
    const bool allGrid = opt.nx && opt.ny && opt.nz;
    if (anyGrid && !allGrid) {
        throw std::runtime_error("--nx, --ny, and --nz must be supplied together");
    }
    if (allGrid) {
        const long long product = 1LL * opt.nx * opt.ny * opt.nz;
        if (product > std::numeric_limits<int>::max()) {
            throw std::runtime_error("grid is too large for int-index prototype");
        }
        opt.blockRows = static_cast<int>(product);
    }
    return opt;
}

Grid deriveGrid(const Options& opt) {
    if (opt.nx && opt.ny && opt.nz) return Grid{opt.nx, opt.ny, opt.nz};

    const int n = opt.blockRows;
    int nx = static_cast<int>(std::ceil(std::cbrt(static_cast<double>(n))));
    nx = std::max(1, nx);
    int ny = nx;
    int nz = (n + nx * ny - 1) / (nx * ny);
    nz = std::max(1, nz);
    return Grid{nx, ny, nz};
}

double signedHash01(int row, int col, int r, int c) {
    unsigned int x = 2166136261u;
    auto mix = [&](unsigned int v) {
        x ^= v + 0x9e3779b9u + (x << 6) + (x >> 2);
        x *= 16777619u;
    };
    mix(static_cast<unsigned int>(row));
    mix(static_cast<unsigned int>(col));
    mix(static_cast<unsigned int>(r));
    mix(static_cast<unsigned int>(c));
    const double u = static_cast<double>(x & 0xffffu) / 65535.0;
    return 2.0 * u - 1.0;
}

std::array<double, kBlockEntries> makeBlock(int row, int col) {
    std::array<double, kBlockEntries> a{};
    const bool diagonalBlock = (row == col);
    const double distance = static_cast<double>(std::abs(row - col) + 1);
    for (int r = 0; r < kBlockSize; ++r) {
        for (int c = 0; c < kBlockSize; ++c) {
            const double h = signedHash01(row, col, r, c);
            double v = 0.0;
            if (diagonalBlock) {
                v = (r == c) ? (7.0 + 0.05 * (r + 1)) : (0.015 * h);
            } else {
                const double base = -0.08 / std::sqrt(distance);
                v = base * ((r == c) ? (1.0 + 0.01 * h) : (0.04 * h));
            }
            a[static_cast<size_t>(r * kBlockSize + c)] = v;
        }
    }
    return a;
}

void pushIfValid(std::vector<int>& cols, int col, int blockRows) {
    if (0 <= col && col < blockRows) cols.push_back(col);
}

Bsr6Matrix buildBsr6(const Options& opt, const Grid& grid) {
    Bsr6Matrix bsr;
    bsr.blockRows = opt.blockRows;
    bsr.rowPtr.reserve(static_cast<size_t>(opt.blockRows) + 1);
    bsr.colIdx.reserve(static_cast<size_t>(opt.blockRows) * 7);
    bsr.values.reserve(static_cast<size_t>(opt.blockRows) * 7);
    bsr.rowPtr.push_back(0);

    const int plane = grid.nx * grid.ny;
    for (int row = 0; row < opt.blockRows; ++row) {
        const int x = row % grid.nx;
        const int y = (row / grid.nx) % grid.ny;
        const int z = row / plane;

        std::vector<int> cols;
        cols.reserve(7);
        cols.push_back(row);
        if (x > 0) pushIfValid(cols, row - 1, opt.blockRows);
        if (x + 1 < grid.nx) pushIfValid(cols, row + 1, opt.blockRows);
        if (y > 0) pushIfValid(cols, row - grid.nx, opt.blockRows);
        if (y + 1 < grid.ny) pushIfValid(cols, row + grid.nx, opt.blockRows);
        if (z > 0) pushIfValid(cols, row - plane, opt.blockRows);
        if (z + 1 < grid.nz) pushIfValid(cols, row + plane, opt.blockRows);

        std::sort(cols.begin(), cols.end());
        cols.erase(std::unique(cols.begin(), cols.end()), cols.end());

        for (int col : cols) {
            bsr.colIdx.push_back(col);
            bsr.values.push_back(makeBlock(row, col));
        }
        bsr.rowPtr.push_back(static_cast<int>(bsr.colIdx.size()));
    }
    return bsr;
}

CsrMatrix expandToScalarCsr(const Bsr6Matrix& bsr) {
    CsrMatrix csr;
    csr.rows = bsr.blockRows * kBlockSize;
    csr.rowPtr.resize(static_cast<size_t>(csr.rows) + 1, 0);
    csr.colIdx.reserve(bsr.values.size() * kBlockEntries);
    csr.values.reserve(bsr.values.size() * kBlockEntries);

    int scalarRow = 0;
    for (int br = 0; br < bsr.blockRows; ++br) {
        for (int lr = 0; lr < kBlockSize; ++lr) {
            csr.rowPtr[static_cast<size_t>(scalarRow)] = static_cast<int>(csr.values.size());
            for (int bi = bsr.rowPtr[static_cast<size_t>(br)];
                 bi < bsr.rowPtr[static_cast<size_t>(br + 1)]; ++bi) {
                const int bc = bsr.colIdx[static_cast<size_t>(bi)];
                const auto& block = bsr.values[static_cast<size_t>(bi)];
                for (int lc = 0; lc < kBlockSize; ++lc) {
                    csr.colIdx.push_back(bc * kBlockSize + lc);
                    csr.values.push_back(block[static_cast<size_t>(lr * kBlockSize + lc)]);
                }
            }
            ++scalarRow;
        }
    }
    csr.rowPtr[static_cast<size_t>(csr.rows)] = static_cast<int>(csr.values.size());
    return csr;
}

std::vector<double> makeInput(int n) {
    std::vector<double> x(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        const double a = std::sin(0.013 * static_cast<double>(i + 1));
        const double b = std::cos(0.017 * static_cast<double>((i % 97) + 1));
        x[static_cast<size_t>(i)] = 0.75 * a + 0.25 * b;
    }
    return x;
}

void matvecCsr(const CsrMatrix& csr, const std::vector<double>& x, std::vector<double>& y) {
    for (int row = 0; row < csr.rows; ++row) {
        double sum = 0.0;
        const int begin = csr.rowPtr[static_cast<size_t>(row)];
        const int end = csr.rowPtr[static_cast<size_t>(row + 1)];
        for (int k = begin; k < end; ++k) {
            sum += csr.values[static_cast<size_t>(k)] * x[static_cast<size_t>(csr.colIdx[static_cast<size_t>(k)])];
        }
        y[static_cast<size_t>(row)] = sum;
    }
}

void matvecBsr6(const Bsr6Matrix& bsr, const std::vector<double>& x, std::vector<double>& y) {
    for (int br = 0; br < bsr.blockRows; ++br) {
        double y0 = 0.0, y1 = 0.0, y2 = 0.0, y3 = 0.0, y4 = 0.0, y5 = 0.0;
        const int begin = bsr.rowPtr[static_cast<size_t>(br)];
        const int end = bsr.rowPtr[static_cast<size_t>(br + 1)];
        for (int bi = begin; bi < end; ++bi) {
            const auto& block = bsr.values[static_cast<size_t>(bi)];
            const double* a = block.data();
            const double* xv = &x[static_cast<size_t>(bsr.colIdx[static_cast<size_t>(bi)] * kBlockSize)];
            const double x0 = xv[0], x1 = xv[1], x2 = xv[2], x3 = xv[3], x4 = xv[4], x5 = xv[5];

            y0 += a[ 0] * x0 + a[ 1] * x1 + a[ 2] * x2 + a[ 3] * x3 + a[ 4] * x4 + a[ 5] * x5;
            y1 += a[ 6] * x0 + a[ 7] * x1 + a[ 8] * x2 + a[ 9] * x3 + a[10] * x4 + a[11] * x5;
            y2 += a[12] * x0 + a[13] * x1 + a[14] * x2 + a[15] * x3 + a[16] * x4 + a[17] * x5;
            y3 += a[18] * x0 + a[19] * x1 + a[20] * x2 + a[21] * x3 + a[22] * x4 + a[23] * x5;
            y4 += a[24] * x0 + a[25] * x1 + a[26] * x2 + a[27] * x3 + a[28] * x4 + a[29] * x5;
            y5 += a[30] * x0 + a[31] * x1 + a[32] * x2 + a[33] * x3 + a[34] * x4 + a[35] * x5;
        }
        double* yy = &y[static_cast<size_t>(br * kBlockSize)];
        yy[0] = y0; yy[1] = y1; yy[2] = y2; yy[3] = y3; yy[4] = y4; yy[5] = y5;
    }
}

double checksum(const std::vector<double>& v) {
    long double sum = 0.0;
    for (size_t i = 0; i < v.size(); ++i) {
        const long double weight = 1.0L + static_cast<long double>(i % 1009) * 0.0009765625L;
        sum += static_cast<long double>(v[i]) * weight;
    }
    return static_cast<double>(sum);
}

double maxAbsDiff(const std::vector<double>& a, const std::vector<double>& b) {
    double m = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        m = std::max(m, std::abs(a[i] - b[i]));
    }
    return m;
}

double relL2Diff(const std::vector<double>& a, const std::vector<double>& b) {
    Eigen::Map<const Eigen::VectorXd> av(a.data(), static_cast<Eigen::Index>(a.size()));
    Eigen::Map<const Eigen::VectorXd> bv(b.data(), static_cast<Eigen::Index>(b.size()));
    return (av - bv).norm() / std::max(1.0, av.norm());
}

template <typename Fn>
Timing runTimed(const char* name, Fn&& fn, std::vector<double>& y, int warmup, int reps) {
    for (int i = 0; i < warmup; ++i) fn();

    Timing t;
    const double t0 = nowMs();
    for (int i = 0; i < reps; ++i) {
        fn();
        const size_t guardIndex = static_cast<size_t>((1LL * i * 9973) % static_cast<long long>(y.size()));
        t.guard += y[guardIndex];
    }
    t.totalMs = nowMs() - t0;
    g_sink += t.guard;
    (void)name;
    return t;
}

void printResult(const char* kernel, const Timing& t, int reps, long long logicalNnz, double csum) {
    const double avgMs = t.totalMs / static_cast<double>(reps);
    const double seconds = t.totalMs / 1000.0;
    const double gnnzPerSec = seconds > 0.0 ? (static_cast<double>(logicalNnz) * reps) / seconds / 1.0e9 : 0.0;
    std::printf(
        "{\"event\":\"result\",\"kernel\":\"%s\",\"total_ms\":%.6f,\"avg_ms\":%.6f,"
        "\"logical_nnz\":%lld,\"gnnz_per_s\":%.6f,\"checksum\":%.17g,\"guard\":%.17g}\n",
        kernel, t.totalMs, avgMs, logicalNnz, gnnzPerSec, csum, t.guard);
}

} // namespace

int main(int argc, char** argv) {
    try {
        const Options opt = parseArgs(argc, argv);
        if (opt.help) {
            printUsage(argv[0]);
            return 0;
        }

        const Grid grid = deriveGrid(opt);
        const Bsr6Matrix bsr = buildBsr6(opt, grid);
        const CsrMatrix csr = expandToScalarCsr(bsr);
        const std::vector<double> x = makeInput(csr.rows);
        std::vector<double> yCsr(static_cast<size_t>(csr.rows), 0.0);
        std::vector<double> yBsr(static_cast<size_t>(csr.rows), 0.0);

        matvecCsr(csr, x, yCsr);
        matvecBsr6(bsr, x, yBsr);

        const double rel = relL2Diff(yCsr, yBsr);
        const double maxAbs = maxAbsDiff(yCsr, yBsr);
        const double csrChecksum = checksum(yCsr);
        const double bsrChecksum = checksum(yBsr);
        const long long scalarNnz = static_cast<long long>(csr.values.size());
        const long long blockNnz = static_cast<long long>(bsr.values.size());

        std::printf(
            "{\"event\":\"config\",\"block_size\":%d,\"block_rows\":%d,\"scalar_rows\":%d,"
            "\"grid_nx\":%d,\"grid_ny\":%d,\"grid_nz\":%d,\"block_nnz\":%lld,"
            "\"scalar_nnz\":%lld,\"avg_block_nnz_per_row\":%.6f,\"reps\":%d,\"warmup\":%d}\n",
            kBlockSize, bsr.blockRows, csr.rows, grid.nx, grid.ny, grid.nz, blockNnz, scalarNnz,
            static_cast<double>(blockNnz) / static_cast<double>(bsr.blockRows), opt.reps, opt.warmup);
        std::printf(
            "{\"event\":\"check\",\"rel_l2\":%.17g,\"max_abs\":%.17g,"
            "\"csr_checksum\":%.17g,\"bsr_checksum\":%.17g}\n",
            rel, maxAbs, csrChecksum, bsrChecksum);

        Timing csrTiming = runTimed("scalar_csr", [&]() { matvecCsr(csr, x, yCsr); }, yCsr, opt.warmup, opt.reps);
        Timing bsrTiming = runTimed("bsr6", [&]() { matvecBsr6(bsr, x, yBsr); }, yBsr, opt.warmup, opt.reps);

        const double finalRel = relL2Diff(yCsr, yBsr);
        const double finalMaxAbs = maxAbsDiff(yCsr, yBsr);
        const double finalCsrChecksum = checksum(yCsr);
        const double finalBsrChecksum = checksum(yBsr);

        printResult("scalar_csr", csrTiming, opt.reps, scalarNnz, finalCsrChecksum);
        printResult("bsr6", bsrTiming, opt.reps, scalarNnzEquivalent(bsr), finalBsrChecksum);

        const double speedup = (bsrTiming.totalMs > 0.0) ? (csrTiming.totalMs / bsrTiming.totalMs) : 0.0;
        std::printf(
            "{\"event\":\"compare\",\"speedup_bsr6_vs_scalar_csr\":%.6f,"
            "\"final_rel_l2\":%.17g,\"final_max_abs\":%.17g,\"sink\":%.17g}\n",
            speedup, finalRel, finalMaxAbs, g_sink);
        return finalRel < 1e-12 ? 0 : 1;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "{\"event\":\"error\",\"message\":\"%s\"}\n", e.what());
        printUsage(argv[0]);
        return 2;
    }
}
