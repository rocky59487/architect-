#pragma once
#include "FrameCore/FrameTypes.h"
#include "FrameCore/Section.h"
#include "FrameCore/Material.h"

namespace frame {

// Governing screening mode. ShellVonMises is the facet surface-stress screen (stage 3d);
// the others are member section modes. Appended only -- existing switches stay valid.
enum class FailMode { None, Crush, Tension, Shear, Bending, Torsion, ShellVonMises };

// Member end forces in LOCAL coordinates. N is COMPRESSION-POSITIVE.
struct MemberEndForces {
    real N = 0, Vy = 0, Vz = 0, T = 0, My = 0, Mz = 0;
};

struct DemandResult {
    real     risk = 0;
    FailMode mode = FailMode::None;
    real     sComp = 0, sTens = 0, tau = 0, sTor = 0;   // intermediate stresses (diagnostics)
};

// Seam between the elastic screening layer (now: ElasticAllowable) and the
// future fiber-section nonlinear layer (RC ultimate). Do not conflate them.
struct ISectionStrength {
    virtual ~ISectionStrength() = default;
    // NOTE: named checkSection (not check) — `check` is a UnrealEngine assertion macro.
    virtual DemandResult checkSection(const MemberEndForces& f, const Section& s, const Capacity& c) const = 0;
};

} // namespace frame
