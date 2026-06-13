#include "FrameCore/Section.h"
#include <algorithm>
#include <cmath>

namespace frame {

// b = width (along local z), d = depth (along local y).
Section Section::Rectangular(real b, real d) {
    Section s;
    s.A  = b * d;
    s.Iz = b * d * d * d / 12.0;   // about local z (depth d in local y) -> governs Mz
    s.Iy = d * b * b * b / 12.0;   // about local y (width b in local z) -> governs My
    s.cz = d / 2.0;                // -> Wz = Iz/cz = b*d^2/6  (textbook section modulus)
    s.cy = b / 2.0;                // -> Wy = Iy/cy = d*b^2/6
    // St-Venant torsion constant for a solid rectangle (a = long side, t = short side).
    const real a = std::max(b, d);
    const real t = std::min(b, d);
    const real r = t / a;
    s.J = a * t * t * t * (1.0 / 3.0 - 0.21 * r * (1.0 - (r * r * r * r) / 12.0));
    // Timoshenko shear coefficient for a rectangle: k = 5/6 -> A_s = 5/6 * A.
    s.Asy = s.Asz = (5.0 / 6.0) * s.A;
    // Plastic section moduli (full-section plastic moment Mp = fy * Z).
    s.Zz = b * d * d / 4.0;
    s.Zy = d * b * b / 4.0;
    s.shape = Shape::Rectangular;
    return s;
}

// Solid circular section of radius r.
Section Section::Circular(real r) {
    Section s;
    const real r2 = r * r;
    s.A  = kPi * r2;
    s.Iy = s.Iz = kPi * r2 * r2 / 4.0;   // I = pi r^4 / 4 (about any centroidal axis)
    s.J  = kPi * r2 * r2 / 2.0;          // J = pi r^4 / 2 (polar = 2I)
    s.cy = s.cz = r;                    // extreme fibre at radius r -> Wy = Wz = pi r^3 / 4
    // Timoshenko shear coefficient for a solid circle ~ 0.9 -> A_s = 0.9 * A.
    s.Asy = s.Asz = 0.9 * s.A;
    // Plastic section modulus of a solid circle: Z = 4 r^3 / 3 (both axes).
    s.Zy = s.Zz = 4.0 * r * r2 / 3.0;
    s.shape = Shape::Circular;
    return s;
}

} // namespace frame
