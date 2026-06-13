#pragma once
#include "FrameCore/FrameTypes.h"

namespace frame {

// S10 -- axially-reduced plastic moment for a plastic hinge forming under combined axial
// force N and bending. Returns
//     Mp_eff(N) = Mp * max(0, 1 - (N/Ny)^2),     Mp = fy*Z,  Ny = fy*A.
//
// EXACTNESS / SCOPE (honest):
//  * RECTANGULAR solid section -> EXACT. The plastic neutral axis shifts so a central band of
//    depth 2c carries the axial force (c = N/(2 b fy)); the surviving moment is
//    Mp - N^2/(4 b fy) which is algebraically identical to Mp(1-(N/Ny)^2). This is the classic
//    textbook plastic N-M envelope (EC3 EN1993-1-1 sec 6.2.9 for a solid rectangle).
//  * CIRCULAR solid section -> CONSERVATIVE. Its true plastic envelope is FULLER than this
//    parabola, so a hinge predicted here forms no later than reality (safe for collapse).
//  * AISC H1.1 is a more conservative BILINEAR design check, not this exact envelope -- do not
//    conflate; this is the cross-section plastic capacity, not a code acceptance ratio.
//  * UNIAXIAL: the reduction uses |N| only. Tension and compression yield at the same |N|=Ny
//    for a compact section, and the square makes the sign of N irrelevant. It does NOT model
//    My-Mz BIAXIAL moment coupling (each bending axis is reduced by the same axial term).
//  * It is a section-capacity ENVELOPE, NOT true elastoplasticity (no unloading/reversal, no
//    N-M tangent coupling). The collapse driver freezes Mp_eff at hinge formation.
//
// Degenerate guards: Mp<=0 -> 0 (not hinge-capable); Ny<=0 -> Mp unreduced (no axial-yield
// reference). |N|>=Ny -> 0 (the section is fully plastified in axial; no moment capacity left).
inline real reducedPlasticMoment(real Mp, real N, real Ny) {
    if (!(Mp > real(0))) return real(0);
    if (!(Ny > real(0))) return Mp;            // no axial-yield reference -> no reduction
    const real n = N / Ny;
    const real f = real(1) - n * n;
    return f > real(0) ? Mp * f : real(0);     // beyond the squash load: no moment capacity
}

}  // namespace frame
