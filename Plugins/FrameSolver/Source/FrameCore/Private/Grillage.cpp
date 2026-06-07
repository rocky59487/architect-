#include "FrameCore/Grillage.h"
#include <cmath>
#include <string>

namespace frame { namespace grillage {

// A slab strip of width w, thickness t, as a beam cross-section, sized so the grillage
// reproduces the ISOTROPIC plate's deflection. Writing the grillage as an equivalent
// orthotropic plate, matching Dx = Dy = 2H = D = E t^3/[12(1-nu^2)] requires
//   bending  i = t^3 / [12 (1 - nu^2)]        -> I = w * i        (so E*i = D)
//   torsion  c = t^3 / [6  (1 - nu )]         -> J = w * c        (so G*c = D, with
//                                                G = E/[2(1+nu)], giving 2*Dxy = D and,
//                                                with no Poisson cross-term, 2H = D)
// The naive Hambly recipe (i=t^3/12, c=t^3/6) instead gives Dx=Dy=(1-nu^2)D and
// 2H=(1-nu)D, i.e. ~26% too flexible for nu=0.3 — these nu-inflations remove that.
// Deflection then tracks plate theory to discretization error; transverse bending
// moments remain over-estimated (the grillage trades the Poisson cross-moment for
// extra twisting moment) — the known grillage-analogy caveat. Iy==Iz so the result is
// independent of the refVec local-axis mapping.
static Section stripSection(real w, real t, real nu) {
    const real i = t * t * t / (12.0 * (1.0 - nu * nu));   // bending inertia per unit width
    const real c = t * t * t / (6.0  * (1.0 - nu));        // torsion constant per unit width
    Section s;
    s.A   = w * t;
    s.Iy  = w * i;
    s.Iz  = w * i;
    s.J   = w * c;
    s.cy  = t / 2.0;
    s.cz  = t / 2.0;
    s.Asy = 0.0;                    // Euler-Bernoulli (thin-plate / Kirchhoff)
    s.Asz = 0.0;
    s.shape = Section::Shape::Rectangular;
    return s;
}

bool buildGrillage(FrameModel& m, const PlateSpec& spec, const Material& mat, std::string& why) {
    m = FrameModel{};
    if (spec.nx < 1 || spec.ny < 1)            { why = "grillage needs nx>=1, ny>=1"; return false; }
    if (spec.a <= 0 || spec.b <= 0 || spec.t <= 0) { why = "grillage needs a,b,t > 0"; return false; }

    const int  nx = spec.nx, ny = spec.ny;
    const real dx = spec.a / nx, dy = spec.b / ny, t = spec.t;

    // Poisson ratio recovered from the isotropic material (G = E/[2(1+nu)]); the strip
    // section inflations above depend on it. Guard to a physical range.
    const real nu = (mat.G > 0) ? (mat.E / (2.0 * mat.G) - 1.0) : 0.0;
    if (!(nu > -0.999 && nu < 0.5)) { why = "grillage: implausible Poisson ratio from E,G (need G=E/2(1+nu))"; return false; }

    // one material (index 0); four strip sections (full / half tributary width in each
    // direction) at indices 0..3. Members reference them by index (no pointer capture).
    m.materials.push_back(mat);                            // material index 0
    m.sections.push_back(stripSection(dy,        t, nu));  // 0 longitudinal interior (width dy)
    m.sections.push_back(stripSection(dy / 2.0,  t, nu));  // 1 longitudinal edge
    m.sections.push_back(stripSection(dx,        t, nu));  // 2 transverse interior (width dx)
    m.sections.push_back(stripSection(dx / 2.0,  t, nu));  // 3 transverse edge
    constexpr int kLongFull = 0, kLongHalf = 1, kTranFull = 2, kTranHalf = 3;

    // grid nodes: lock the in-plane DOFs everywhere (grillage carries only out-of-plane
    // action); a simply-supported boundary fixes Uz on the four edges.
    m.nodes.reserve(static_cast<size_t>(nx + 1) * (ny + 1));
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nx; ++i) {
            Node n(static_cast<NodeId>(gridNode(i, j, nx)), i * dx, j * dy, 0.0);
            n.fixed[Ux] = n.fixed[Uy] = n.fixed[Rz] = true;
            const bool onEdge = (i == 0 || i == nx || j == 0 || j == ny);
            if (onEdge && spec.simplySupported) n.fixed[Uz] = true;
            m.nodes.push_back(n);
        }

    // longitudinal (x) + transverse (y) beams; default refVec=+Z is valid (none vertical).
    MemberId mid = 0;
    for (int j = 0; j <= ny; ++j) {
        const int sIdx = (j == 0 || j == ny) ? kLongHalf : kLongFull;
        for (int i = 0; i < nx; ++i)
            m.members.push_back(Member(mid++, gridNode(i, j, nx), gridNode(i + 1, j, nx), 0, sIdx));
    }
    for (int i = 0; i <= nx; ++i) {
        const int sIdx = (i == 0 || i == nx) ? kTranHalf : kTranFull;
        for (int j = 0; j < ny; ++j)
            m.members.push_back(Member(mid++, gridNode(i, j, nx), gridNode(i, j + 1, nx), 0, sIdx));
    }

    // uniform pressure q -> consistent nodal loads (q * tributary area), downward (-z).
    // Tributary areas partition a*b, so the total applied load is exactly q*a*b.
    for (int j = 0; j <= ny; ++j) {
        const real wy = (j == 0 || j == ny) ? dy / 2.0 : dy;
        for (int i = 0; i <= nx; ++i) {
            const real wx = (i == 0 || i == nx) ? dx / 2.0 : dx;
            NodalLoad nl;
            nl.node = static_cast<NodeId>(gridNode(i, j, nx));
            nl.comp[Uz] = -spec.q * wx * wy;
            m.nodalLoads.push_back(nl);
        }
    }

    why.clear();
    return true;
}

}} // namespace frame::grillage
