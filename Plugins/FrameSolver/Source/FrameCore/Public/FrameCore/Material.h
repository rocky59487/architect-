#pragma once
#include "FrameCore/FrameTypes.h"
#include <algorithm>

namespace frame {

// ALLOWABLE / ELASTIC capacities (MPa). This is the screening layer — it is
// NOT RC ultimate strength (P-M interaction + concrete nonlinearity). The real
// RC strength belongs to a future fiber-section precision layer; do not conflate.
struct Capacity {
    real comp  = 0;   // allowable compressive stress
    real tens  = 0;   // allowable tensile stress
    real shear = 0;   // allowable shear stress
    real bend  = 0;   // elastic proxy = min(comp, tens)  (NOT a real Mn)
    real tors  = 0;   // elastic proxy = shear

    static Capacity make(real comp_, real tens_, real shear_) {
        Capacity c;
        c.comp  = comp_;
        c.tens  = tens_;
        c.shear = shear_;
        c.bend  = std::min(comp_, tens_);
        c.tors  = shear_;
        return c;
    }
};

struct Material {
    real     E   = 0;     // Young's modulus (MPa)
    real     G   = 0;     // shear modulus (MPa)
    real     nu  = 0;     // Poisson ratio — used by the shell plane-stress/bending
                          // constitutive (Dm, Db). Beams use E,G directly and ignore
                          // nu, so this defaults to 0 and leaves existing fixtures
                          // (which never set it) byte-for-byte unchanged.
    real     rho = 0;     // density (kg/m^3). Used by self-weight (addSelfWeight) and the
                          // mass matrix; both apply the unit bridge rho[tonne/mm^3] =
                          // rho[kg/m^3]*1e-12 to stay in the engine's N-mm-tonne-s system.
    Capacity cap;

    Material() = default;
    Material(real E_, real G_, real rho_ = 0) : E(E_), G(G_), rho(rho_) {}
};

} // namespace frame
