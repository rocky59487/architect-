#include "FrameCore/ElasticAllowable.h"
#include <algorithm>
#include <cmath>
#include <limits>

namespace frame {

// Elastic / allowable-stress combined-stress screen (spec PFSFv2-to-UE5 §1.4).
// NOT RC ultimate strength.
DemandResult ElasticAllowable::checkSection(const MemberEndForces& f, const Section& s, const Capacity& c) const {
    const real A  = std::max(s.A,   real(1e-12));
    const real Wy = std::max(s.Wy(), real(1e-12));
    const real Wz = std::max(s.Wz(), real(1e-12));
    const real Jt = std::max(s.J,   real(1e-12));

    const real sN    = f.N / A;                                      // compression-positive
    // Biaxial bending stress. A round section has no worst corner, so the resultant
    // moment sqrt(My^2+Mz^2)/W is exact; a rectangle uses the conservative corner sum
    // |My|/Wy + |Mz|/Wz (report PFSFv2-to-UE5 §0c).
    const real sM    = (s.shape == Section::Shape::Circular)
                       ? std::sqrt(f.My * f.My + f.Mz * f.Mz) / Wz   // Wy == Wz for a circle
                       : std::abs(f.My) / Wy + std::abs(f.Mz) / Wz;
    const real sComp = std::max(sM + sN, real(0));
    const real sTens = std::max(sM - sN, real(0));
    // Transverse shear stress. V/A is the cross-section AVERAGE; the peak (at the neutral
    // axis) is k*V/A with k = 1.5 for a rectangle and 4/3 for a circle. We screen on the
    // peak so a shear-controlled member is not under-checked.
    const real shearPeak = (s.shape == Section::Shape::Circular) ? real(4.0 / 3.0) : real(1.5);
    const real tau   = shearPeak * std::sqrt(f.Vy * f.Vy + f.Vz * f.Vz) / A;
    // Torsional shear ~ T*c/J  (c = extreme-fibre distance). Units N*mm*mm/mm^4 = MPa,
    // comparable to the shear capacity. For a CIRCLE this is the exact max shear T*r/J at
    // the surface (c = r); for a RECTANGLE we use the diagonal corner hypot(cy,cz) as a
    // conservative heuristic (true St-Venant rectangular torsion needs warping, deferred).
    const real cTor  = (s.shape == Section::Shape::Circular) ? s.cy : std::hypot(s.cy, s.cz);
    const real sTor  = std::abs(f.T) * cTor / Jt;

    auto ratio = [](real demand, real cap) -> real {
        if (cap > 0) return demand / cap;
        return demand > 0 ? std::numeric_limits<real>::infinity() : real(0);
    };
    const real r[5] = { ratio(sComp, c.comp), ratio(sTens, c.tens),
                        ratio(tau, c.shear),  ratio(sM, c.bend), ratio(sTor, c.tors) };
    int k = 0;
    for (int i = 1; i < 5; ++i) if (r[i] > r[k]) k = i;
    static const FailMode M[5] = { FailMode::Crush, FailMode::Tension,
                                   FailMode::Shear, FailMode::Bending, FailMode::Torsion };

    DemandResult d;
    d.risk  = r[k];
    d.mode  = r[k] > 0 ? M[k] : FailMode::None;
    d.sComp = sComp; d.sTens = sTens; d.tau = tau; d.sTor = sTor;
    return d;
}

// Stage 3d: shell surface von Mises screen. Membrane stress from the centre resultants
// (element-constant approximation), bending stress from the centre AND per-corner moments,
// both faces; the worst sample governs. Mirrors checkSection's ratio() semantics, including
// "zero capacity under demand = infinite D/C".
ShellDemandResult checkShellSurface(const ShellElementForces& f, real t, const Capacity& c) {
    ShellDemandResult out;
    if (!(t > 0)) return out;   // validate() rejects t <= 0; defensive zero here

    auto vonMises = [](real sx, real sy, real txy) {
        return std::sqrt(std::max(real(0), sx * sx - sx * sy + sy * sy + 3.0 * txy * txy));
    };
    auto ratio = [](real demand, real cap) -> real {
        if (cap > 0) return demand / cap;
        return demand > 0 ? std::numeric_limits<real>::infinity() : real(0);
    };

    const real bend = 6.0 / (t * t);              // sigma_bend = 6*M/t^2 per unit moment
    const real mx = f.Nxx / t, my = f.Nyy / t, mxy = f.Nxy / t;   // membrane (centre)

    for (int kc = -1; kc < 4; ++kc) {             // -1 = centre, 0..3 = corners
        const real Mx  = (kc < 0) ? f.Mxx : f.MxxC[kc];
        const real My  = (kc < 0) ? f.Myy : f.MyyC[kc];
        const real Mxy = (kc < 0) ? f.Mxy : f.MxyC[kc];
        for (int face = 0; face < 2; ++face) {    // 0 = top (+bending), 1 = bottom (-bending)
            const real s = (face == 0) ? real(1) : real(-1);
            const real r = ratio(vonMises(mx + s * bend * Mx, my + s * bend * My, mxy + s * bend * Mxy), c.vm);
            if (r > out.risk) { out.risk = r; out.corner = kc; out.top = (face == 0); }
        }
    }
    return out;
}

// Shell counterpart of worstUtilization (same skip rules: inactive / out-of-range matIdx).
ShellDemandSummary worstShellUtilization(const FrameModel& model, const SolveResult& r) {
    ShellDemandSummary out;
    const size_t nS = std::min(model.shells.size(), r.shellForces.size());
    bool any = false;
    real maxDC = 0;
    for (size_t s = 0; s < nS; ++s) {
        const ShellQuad& sh = model.shells[s];
        if (!sh.active) continue;
        if (sh.matIdx < 0 || sh.matIdx >= (int)model.materials.size()) continue;
        const ShellDemandResult d = checkShellSurface(r.shellForces[s], sh.t, model.materials[(size_t)sh.matIdx].cap);
        if (!any || d.risk > maxDC) {
            maxDC = d.risk;
            out.governingShell = sh.id;
        }
        any = true;
    }
    out.valid = any;
    out.maxDC = any ? maxDC : real(0);
    return out;
}

// C3: worst Demand/Capacity over all ACTIVE members (both ends) + the elastic safety factor.
// Inactive members (element removal) are skipped; members with an out-of-range material/section
// index are skipped (validate() should have caught those, but be defensive). Pure post-process.
DemandSummary worstUtilization(const FrameModel& model, const SolveResult& r) {
    DemandSummary out;
    const ElasticAllowable screen;
    const size_t nM = std::min(model.members.size(), r.memberForces.size());
    bool any = false;
    real maxDC = 0;
    for (size_t e = 0; e < nM; ++e) {
        const Member& mem = model.members[e];
        if (!mem.active) continue;
        if (mem.matIdx < 0 || mem.matIdx >= (int)model.materials.size()) continue;
        if (mem.secIdx < 0 || mem.secIdx >= (int)model.sections.size())  continue;
        const Section&  s = model.sections[mem.secIdx];
        const Capacity& c = model.materials[mem.matIdx].cap;
        const MemberForcePair& mf = r.memberForces[e];
        const DemandResult di = screen.checkSection(mf.endI, s, c);
        const DemandResult dj = screen.checkSection(mf.endJ, s, c);
        const DemandResult& d = (di.risk >= dj.risk) ? di : dj;   // worse of the two ends
        if (!any || d.risk > maxDC) {
            maxDC = d.risk;
            out.governingMember = mem.id;
            out.mode = d.mode;
        }
        any = true;
    }
    out.valid = any;
    out.maxDC = any ? maxDC : real(0);
    out.safetyFactor = any ? ((maxDC > 0) ? real(1) / maxDC : std::numeric_limits<real>::infinity())
                           : real(0);
    return out;
}

} // namespace frame
