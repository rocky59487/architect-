#pragma once
#include "FrameCore/FrameTypes.h"
#include <vector>
#include <string>

namespace frame {

struct BucklingResult {
    bool              singular = false;
    std::string       diagnostic;
    real              criticalFactor = 0;   // lambda_cr: P_cr = lambda_cr * applied reference load
    std::vector<real> mode;                 // 6N buckling mode shape (constrained DOFs = 0)
};

}  // namespace frame
