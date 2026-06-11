#pragma once
#include "FrameCore/FrameModel.h"
#include "FrameCore/SolveResult.h"
#include "FrameCore/SolveOptions.h"

namespace frame {

// POD options for runPDelta (S3). POD/std only — no Eigen, no UE — so this stays on the pure
// public boundary. Two solution paths select via refactorPath:
//   * FROZEN (default, refactorPath=false): pseudo-load iteration  u <- K_e^-1 (F - Kg u)
//     that REUSES the existing K_e factorization (the assembleAndFactor LDLT). Zero re-factor,
//     so it is the natural client of the factorize-once architecture (one expensive factor,
//     many cheap forward/back-substitutions). Converges while P < P_cr at rate ~ P/P_cr.
//   * REFERENCE (refactorPath=true): one fresh LDLT of the tangent K_T = K_e + Kg (Wilson &
//     Habibullah 1987 form). A single solve; no iteration. Used as the independent cross-check
//     of the frozen path (they share one fixed point) and when the frozen iteration is slow
//     (P near P_cr) so the iteration count would blow past maxIter.
struct PDeltaOptions {
    int  maxIter     = 200;     // frozen-path iteration cap (P near P_cr needs more; raise it)
    real tolU        = 1e-12;   // frozen-path convergence on ||du|| / ||u||
    bool accelerate  = true;    // PROTECTED geometric (Aitken-style) extrapolation; see note below
    bool refactorPath = false;  // false = frozen reuse path (default); true = K_T refactor reference
    SolveOptions solve;         // threaded into assembleAndFactor / the underlying linear solve
};

// POD result. finalState is the FULL recovered second-order state (displacements scaled by the
// P-Delta amplification + member end forces evaluated at those displacements = the second-order
// internal forces). converged / diverged are mutually exclusive on a healthy run; both false
// means the cap (maxIter) was hit without a verdict (report lastIncrement to the caller).
struct PDeltaResult {
    bool        converged    = false;
    bool        diverged     = false;   // P past critical: frozen increments grow (sliding-window
                                        // detector) / reference K_T not positive-definite.
    int         iterations   = 0;       // frozen path: pseudo-load steps; reference path: 1; P=0: 0
    real        lastIncrement = 0;      // last ||du|| / ||u|| (frozen path)
    SolveResult finalState;             // recovered second-order result (singular flag forwarded)
};

// Linearized (Theory-II) P-Delta second-order analysis. Builds the geometric stiffness Kg from
// the FIRST-ORDER member axial forces (N FROZEN at the linear solve — the standard Theory-II
// convention; NOT updated as the structure deflects, which would be the large-displacement
// corotational analysis reserved for S9), then either iterates the pseudo-load with the reused
// K_e factorization (frozen, default) or factors K_T = K_e + Kg once (reference). Both converge
// to the same second-order displacement. For a tip-loaded cantilever column under axial P this
// reproduces the beam-column amplification  delta = H (tan(kL) - kL)/(P k),  k = sqrt(P/EI), and
// recovers the linear solution bit-for-bit at P = 0.
//
// Honest scope: small-sway linearization, axial force frozen at first order, beam-column
// geometric stiffness only (shells do not contribute Kg — the existing buckling limitation).
// Convergence (frozen path) slows as P -> P_cr (iteration count ~ log(tol)/log(P/P_cr)); past
// P_cr the analysis is unstable and is reported as diverged, not silently wrong.
FRAMECORE_API PDeltaResult runPDelta(const FrameModel& model, const PDeltaOptions& opts = {});

}  // namespace frame
