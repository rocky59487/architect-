#pragma once
#include "FrameCore/FrameTypes.h"
#include <vector>
#include <string>

namespace frame {

struct ModalOptions {
    int numModes = 3;   // number of lowest modes to extract
};

struct ModeShape {
    real              omega  = 0;   // circular frequency (rad/s)
    real              freqHz = 0;   // = omega / (2*pi)
    std::vector<real> shape;        // 6N mass-normalized eigenvector (constrained DOFs = 0)
};

struct ModalResult {
    bool                   singular = false;
    std::string            diagnostic;
    std::vector<ModeShape> modes;   // ascending frequency
};

}  // namespace frame
