#pragma once
#include "FrameCore/ISectionStrength.h"
#include "FrameCore/FrameModel.h"
#include "FrameCore/SolveResult.h"

namespace frame {

// Elastic / allowable-stress combined-stress screen (spec PFSFv2-to-UE5 §1.4).
//   sigma = N/A +/- M/W ,  tau = |V|/A ,  demand/capacity per mode, argmax = dominant.
// This is a real-time SCREENING layer, NOT RC ultimate strength (no P-M interaction,
// no concrete nonlinearity). Biaxial bending uses the conservative |My|/Wy + |Mz|/Wz
// worst-corner sum (valid for rectangular; circular sections need resultant moment).
struct ElasticAllowable final : ISectionStrength {
    FRAMECORE_API DemandResult checkSection(const MemberEndForces& f, const Section& s, const Capacity& c) const override;
};

// Structural utilization aggregate (C3 / collapse foundation). The worst Demand/Capacity over all
// ACTIVE members (both ends), and the elastic SAFETY FACTOR = 1/maxDC: the linear load multiplier
// that brings the most-utilized member to its first allowable-stress limit ("how far from first
// failure"). This is the same ElasticAllowable SCREEN aggregated per structure — NOT an RC
// ultimate-strength margin. Kept as a free post-process (the solver stays capacity-free), matching
// combine()/envelope(): a pure function of (model, result).
struct DemandSummary {
    real     maxDC          = 0;            // worst Demand/Capacity (0 if no screenable active member)
    real     safetyFactor   = 0;            // 1/maxDC; +infinity when maxDC == 0; 0 when invalid
    MemberId governingMember = 0;           // id of the member carrying maxDC
    FailMode mode           = FailMode::None;
    bool     valid          = false;        // false when no active member could be screened
};
FRAMECORE_API DemandSummary worstUtilization(const FrameModel& model, const SolveResult& r);

// Shell SURFACE-STRESS von Mises screen (stage 3d; elastic, screening-grade):
//   sigma_x = Nxx/t +/- 6*Mxx/t^2 (top/bottom faces), same for yy / xy,
//   sigma_vM = sqrt(sx^2 - sx*sy + sy^2 + 3*txy^2),  D/C = sigma_vM / cap.vm.
// Evaluated at the element CENTRE and at the four CORNERS (membrane taken at the centre --
// the element-constant approximation -- bending from the per-corner recovery), on both
// faces; risk = the max over all samples. Honest boundary: transverse shear Qx/Qy is NOT
// screened (thin-plate regime), no plate buckling, no plate ultimate strength.
struct ShellDemandResult {
    real risk   = 0;
    int  corner = -1;     // governing sample: -1 = centre, 0..3 = corner (ShellQuad::n order)
    bool top    = true;   // governing face (+bending side)
};
FRAMECORE_API ShellDemandResult checkShellSurface(const ShellElementForces& f, real t, const Capacity& c);

// Shell counterpart of worstUtilization: worst surface von Mises D/C over all ACTIVE shells.
// Same free-post-process contract: a pure function of (model, result).
struct ShellDemandSummary {
    real maxDC          = 0;       // worst D/C (0 if no screenable active shell)
    int  governingShell = 0;       // id of the shell carrying maxDC
    bool valid          = false;   // false when no active shell could be screened
};
FRAMECORE_API ShellDemandSummary worstShellUtilization(const FrameModel& model, const SolveResult& r);

} // namespace frame
