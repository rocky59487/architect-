#pragma once
#include "FrameCore/SolveResult.h"
#include <vector>

namespace frame {

// Linear superposition of per-load-case results:  out = sum_i factors[i] * cases[i].
//
// This is the engine's load-COMBINATION primitive (e.g. 1.2*Dead + 1.6*Live). It is a
// pure post-process and is valid ONLY for LINEAR-elastic analysis (the engine's regime) —
// superposition does NOT hold once geometric (P-Delta) or material nonlinearity is added.
// All result fields are combined component-wise: displacements, reactions, member end
// forces and shell stress resultants. `singular` is true if any combined case was singular.
//
// `cases` and `factors` are paired by index; extra entries in either are ignored. All cases
// must come from the SAME model (same DOF/member/shell counts).
FRAMECORE_API SolveResult combine(const std::vector<SolveResult>& cases,
                                  const std::vector<real>& factors);

}  // namespace frame
