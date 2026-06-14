#pragma once
#include "FrameCore/FrameSolver.h"      // PreparedSystem (opaque PIMPL)
#include "FrameCore/FrameModel.h"
#include "FrameCore/SolveResult.h"
#include "FrameCore/HpSolveOptions.h"

namespace frame {

// Opt-in HP-FEM iterative solve lane (research-validated; see Research/WS_B_solver and
// HPFEM_RESEARCH_NOTES.md). Computes the SAME SolveResult as solveLoad() — identical RHS
// assembly, reactions (K*u - F) and element force recovery — but replaces the LDLT
// forward/back-substitution with a matrix-free preconditioned conjugate gradient over the
// prepared elements.
//
// One-shot (stateless): serial element apply + Jacobi preconditioner, frame-only. It falls back to
// the PreparedSystem's LDLT (the oracle) whenever the model is singular, contains a non-beam (shell)
// element, or the PCG does not converge — so the numeric result never deviates from the direct
// solve. With opts.enabled == false it goes straight to LDLT, making solveLoadHP(prepared, model) a
// drop-in equal to solveLoad(prepared, model).
//
// PERFORMANCE: being stateless it has NO seeded subspace, so for an arbitrary RHS it is at best on
// par with — and usually an order of magnitude or more SLOWER than — a reused LDLT back-substitution.
// It exists for completeness / single-shot use; the repeated-solve speedup lives in HpSession (a
// seeded session). solve() / solveLoad() remain the default and the oracle.
//
// This is a parallel entry point: solve(), solveLoad() and assembleAndFactor() are NOT
// changed, and the public boundary stays Eigen-free (POD SolveResult).
FRAMECORE_API SolveResult solveLoadHP(const PreparedSystem& prepared,
                                      const FrameModel& model,
                                      const HpSolveOptions& opts = {});

}  // namespace frame
