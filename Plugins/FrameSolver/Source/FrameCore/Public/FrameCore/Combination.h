#pragma once
#include "FrameCore/SolveResult.h"
#include <vector>
#include <string>

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

// Component-wise ENVELOPE (max/min) of result quantities across many cases/combinations —
// the standard way to capture the most-unfavourable live-load PATTERN (e.g. chessboard /
// skip loading): solve each pattern, then take the worst value of each quantity. Pure
// post-process. All results must come from the same model.
//
// CAUTION (component-wise): each component is enveloped INDEPENDENTLY, so a min/max of one
// component need not come from the same load case as another. A RESULTANT magnitude formed
// from envelope components (e.g. sqrt(My_min^2 + Mz_min^2)) therefore does NOT correspond to
// any single real load case and can be non-conservative/meaningless under biaxial action.
// Envelope per-component; do not combine enveloped components into a resultant.
struct ResultEnvelope {
    std::vector<real>            uMax, uMin;            // per global DOF
    std::vector<real>            reactMax, reactMin;    // per global DOF
    std::vector<MemberEndForces> endIMax, endIMin;      // per member, component-wise
    std::vector<MemberEndForces> endJMax, endJMin;
    std::vector<ShellElementForces> shellMax, shellMin; // per shell, component-wise

    // Set true when the cases are NOT envelopable together (results from different models —
    // mismatched DOF / member / shell counts). Mirrors SolveResult/combine(): rather than let
    // the min()-bounded loop silently fold partial/garbage data, the envelope is flagged invalid
    // and `diagnostic` explains why. Empty input is a no-op (singular stays false).
    bool        singular = false;
    std::string diagnostic;
};
FRAMECORE_API ResultEnvelope envelope(const std::vector<SolveResult>& cases);

}  // namespace frame
