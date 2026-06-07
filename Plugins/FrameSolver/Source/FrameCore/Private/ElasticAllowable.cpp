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

} // namespace frame
