#include "FrameCore/SizeOpt.h"
#include "FrameCore/FrameSolver.h"
#include "FrameCore/ElasticAllowable.h"
#include "MemberQuery.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <set>
#include <string>

namespace frame {
namespace {

// Similarity (shape-preserving) rescale of a section to a new area. Linear dimensions scale by
// lambda = sqrt(Anew/Aold), so A ~ lambda^2, I/J ~ lambda^4, extreme-fibre c ~ lambda, shear area
// ~ lambda^2, plastic modulus Z ~ lambda^3. Then W = I/c ~ lambda^3. Rescaling from the ORIGINAL
// section (not the previous iterate) avoids compounding round-off. For a square this reproduces
// Section::Rectangular(sqrt(Anew), sqrt(Anew)) bit-for-bit.
Section scaledFrom(const Section& orig, real Anew) {
    Section s = orig;
    const real Aold = orig.A;
    if (Aold <= 0) return s;        // defensive: leave untouched
    if (Anew <= 0) {
        s.A = s.Iy = s.Iz = s.J = s.cy = s.cz = s.Asy = s.Asz = s.Zy = s.Zz = 0;
        return s;
    }
    const real r   = Anew / Aold;
    const real lam = std::sqrt(r);
    s.A   = Anew;
    s.Iy  = orig.Iy * r * r;
    s.Iz  = orig.Iz * r * r;
    s.J   = orig.J  * r * r;
    s.cy  = orig.cy * lam;
    s.cz  = orig.cz * lam;
    s.Asy = orig.Asy * r;
    s.Asz = orig.Asz * r;
    s.Zy  = orig.Zy * (r * lam);
    s.Zz  = orig.Zz * (r * lam);
    // shape preserved
    return s;
}

// Round UP to the smallest table area >= A (conservative); clamp to the largest if A exceeds it.
// `table` is assumed ascending and non-empty (the caller only calls this when it is set).
real snapUp(const std::vector<real>& table, real A) {
    for (real t : table) if (t >= A) return t;
    return table.back();
}

// FNV-1a hash of the area vector's raw bit patterns -- a finite discrete state space (snapped to a
// table) revisiting a state means a flip-flop, which the guard turns into finite termination.
uint64_t hashAreas(const std::vector<real>& a) {
    uint64_t h = 1469598103934665603ull;
    for (real x : a) {
        uint64_t bits = 0;
        std::memcpy(&bits, &x, sizeof(bits) < sizeof(x) ? sizeof(bits) : sizeof(x));
        for (int b = 0; b < 8; ++b) { h ^= (bits >> (b * 8)) & 0xffu; h *= 1099511628211ull; }
    }
    return h;
}

}  // namespace

SizeOptResult runSizeOptimization(const FrameModel& model, const SizeOptOptions& opts,
                                  const std::vector<int>& sizableMembers) {
    SizeOptResult R;
    std::string why;
    if (!model.validate(why)) { R.singular = true; return R; }

    FrameModel work = model;                 // the caller's model is never mutated
    const size_t nMem = work.members.size();

    // Resolve the sizable set (member indices). Empty request = every screenable active member.
    std::vector<int> sized;
    if (sizableMembers.empty()) {
        for (size_t e = 0; e < nMem; ++e) if (screenable(work, e)) sized.push_back((int)e);
    } else {
        for (int e : sizableMembers)
            if (e >= 0 && e < (int)nMem && screenable(work, (size_t)e)) sized.push_back(e);
    }

    // Give each sizable member a PRIVATE section copy so resizing one never affects another that
    // happened to share a section index. Track the original section + current area per sizable member.
    std::vector<Section> origSec(sized.size());
    std::vector<real>    area(sized.size());
    std::vector<int>     sizePos(nMem, -1);
    for (size_t k = 0; k < sized.size(); ++k) {
        const int e = sized[k];
        const Section orig = work.sections[(size_t)work.members[(size_t)e].secIdx];
        origSec[k] = orig;
        area[k]    = orig.A;
        work.sections.push_back(orig);
        work.members[(size_t)e].secIdx = (int)work.sections.size() - 1;
        sizePos[(size_t)e] = (int)k;
    }

    const ElasticAllowable screen;
    const int nCases = opts.cases.empty() ? 1 : (int)opts.cases.size();

    auto sizedVolume = [&]() -> real {
        real vol = 0;
        for (size_t k = 0; k < sized.size(); ++k)
            vol += area[k] * memberLen(work, sized[k]);
        return vol;
    };

    // Nothing to resize: report the model as-is (converged, no work).
    if (sized.empty()) {
        R.converged = true;
    } else {
        std::set<uint64_t> seen;
        for (int it = 1; it <= opts.maxIter; ++it) {
            R.iterations = it;

            // 1) per-sizable-member worst D/C across both ends and all load cases (envelope)
            std::vector<real> dc(sized.size(), 0);
            bool sawSingular = false;
            for (int ci = 0; ci < nCases; ++ci) {
                if (!opts.cases.empty()) {
                    work.nodalLoads = opts.cases[(size_t)ci].nodalLoads;
                    work.memberUDLs = opts.cases[(size_t)ci].memberUDLs;
                }
                const SolveResult r = solve(work, opts.solve);
                if (r.singular) { sawSingular = true; break; }
                const size_t nf = std::min(nMem, r.memberForces.size());
                for (size_t e = 0; e < nf; ++e) {
                    if (!screenable(work, e)) continue;
                    const int p = sizePos[e];
                    if (p < 0) continue;
                    dc[(size_t)p] = std::max(dc[(size_t)p], memberDC(screen, work, r, e));
                }
            }
            if (sawSingular) { R.singular = true; break; }

            bool sawNonFinite = false;
            real worstDC = 0;
            for (real d : dc) {
                if (!std::isfinite(d)) {
                    sawNonFinite = true;
                    worstDC = std::numeric_limits<real>::infinity();
                } else if (!sawNonFinite) {
                    worstDC = std::max(worstDC, d);
                }
            }
            if (sawNonFinite) {
                R.invalidDemand = true;
                R.dcHistory.push_back(worstDC);
                R.weightHistory.push_back(sizedVolume());
                break;
            }

            // 2) resize + bookkeeping
            real maxDevFS = 0, vol = 0;
            int changed = 0;
            for (size_t k = 0; k < sized.size(); ++k) {
                const real factor = dc[k];
                real Anew = std::max(opts.Amin, area[k] * factor);
                if (!opts.sectionTable.empty()) Anew = snapUp(opts.sectionTable, Anew);
                if (Anew != area[k]) ++changed;
                area[k] = Anew;
                Member& mem = work.members[(size_t)sized[k]];
                if (Anew <= 0) {
                    mem.active = false;
                } else {
                    mem.active = true;
                    work.sections[(size_t)mem.secIdx] = scaledFrom(origSec[k], Anew);
                }
                if (Anew > opts.Amin * (real(1) + real(1e-9)))
                    maxDevFS = std::max(maxDevFS, std::fabs(factor - real(1)));
                vol += Anew * memberLen(work, sized[k]);
            }
            R.dcHistory.push_back(worstDC);
            R.weightHistory.push_back(vol);

            // 3) convergence
            if (opts.sectionTable.empty()) {
                if (it >= 3 && maxDevFS < opts.dcTol) { R.converged = true; break; }
            } else {
                if (changed == 0) { R.converged = true; break; }
                if (!seen.insert(hashAreas(area)).second) { R.cycled = true; R.converged = true; break; }
            }
        }
    }

    // 4) fill per-member results from the converged working model
    R.finalAreas.assign(nMem, 0);
    R.finalSections.assign(nMem, Section{});
    R.finalDC.assign(nMem, 0);
    for (size_t e = 0; e < nMem; ++e) {
        const int p = sizePos[e];
        if (p >= 0) {
            R.finalSections[e] = scaledFrom(origSec[(size_t)p], area[(size_t)p]);
            R.finalAreas[e] = area[(size_t)p];
            continue;
        }
        const Member& m = work.members[e];
        if (m.secIdx >= 0 && m.secIdx < (int)work.sections.size()) {
            R.finalSections[e] = work.sections[(size_t)m.secIdx];
            R.finalAreas[e]    = work.sections[(size_t)m.secIdx].A;
        }
    }

    // finalDC: worst D/C per member over all load cases on the CONVERGED sections (so the audit can
    // assert the sized structure is fully-stressed / safe under every case).
    if (!R.singular) {
        std::vector<real> fdc(nMem, 0);
        bool ok = true;
        for (int ci = 0; ci < nCases && ok; ++ci) {
            if (!opts.cases.empty()) {
                work.nodalLoads = opts.cases[(size_t)ci].nodalLoads;
                work.memberUDLs = opts.cases[(size_t)ci].memberUDLs;
            }
            const SolveResult r = solve(work, opts.solve);
            if (r.singular) { ok = false; break; }
            const size_t nf = std::min(nMem, r.memberForces.size());
            for (size_t e = 0; e < nf; ++e)
                if (screenable(work, e)) fdc[e] = std::max(fdc[e], memberDC(screen, work, r, e));
        }
        if (ok) R.finalDC = fdc;
    }
    return R;
}

}  // namespace frame
