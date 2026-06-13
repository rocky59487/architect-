#pragma once
//
// Element factory: build the concrete IElement for a member / shell index. Shared by
// FrameSolver.cpp (assembleAndFactor) and Reanalysis.cpp so element construction lives in one
// place — changing how an element is constructed (e.g. a new ctor argument) is then a one-line
// edit. Header-only inline; pulls in the concrete element headers, so include it only from TUs
// that already build elements.
//
#include "IElement.h"
#include "BeamColumnElement.h"
#include "MITC4ShellElement.h"
#include <memory>

namespace frame {

inline std::unique_ptr<IElement> makeMemberElem(int e) { return std::make_unique<BeamColumnElement>(e); }
inline std::unique_ptr<IElement> makeShellElem(int s)  { return std::make_unique<MITC4ShellElement>(s); }

}  // namespace frame
