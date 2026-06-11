#pragma once
#include "FrameCore/FrameTypes.h"

namespace frame {

// POD options threaded into solve(). MUST stay engine-agnostic: only real / bool /
// int / plain enums -- never UE or Eigen types -- so FrameCore stays pure.
struct SolveOptions {
    real pivotTol       = 1e-12;   // mechanism-detection pivot tolerance (relative to max|D|)
    bool enableReleases = false;   // honor Member.release[12] via per-element static condensation
    bool useTimoshenko  = false;   // include shear flexibility (needs Section.Asy/Asz > 0);
                                   // false keeps the Euler-Bernoulli element bit-for-bit.

    // S8-8a: opt-in QM6 incompatible-mode membrane (Wilson Q6 1973 + Taylor 1976 correction).
    // false keeps the original bilinear Q4 shell membrane BIT-FOR-BIT, so the OpenSees
    // ShellMITC4 plate gate (~1e-10) is preserved. true adds two bubble modes (1-xi^2,
    // 1-eta^2) on the in-plane u,v and condenses them out per element, defeating in-plane
    // membrane locking. The drilling (Rz / Hughes-Brezzi) block is untouched either way.
    bool useIncompatibleMembrane = false;

    // S8-8b: opt-in DKQ discrete-Kirchhoff THIN-plate bending (Batoz & Tahar 1982), replacing
    // the MITC4 assumed-shear bending block. false keeps the MITC4 Reissner-Mindlin plate
    // BIT-FOR-BIT (OpenSees ShellMITC4 gate). true = no transverse-shear DOF -> thin plate
    // only (t/L < ~1/20); mid/thick plates MUST keep this false (MITC4). Membrane + drilling
    // are shared -- this flag ONLY swaps the bending block. Recovered Qx=Qy=0 (Kirchhoff).
    bool useDKQPlate = false;
};

} // namespace frame
