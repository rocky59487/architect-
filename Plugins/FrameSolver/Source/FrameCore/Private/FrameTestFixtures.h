#pragma once
//
// Pure-`frame` fixture builders (NO Eigen, NO test framework) shared by the
// standalone harness (Standalone/main.cpp) AND the UE automation tests
// (Private/Tests/*.cpp), so both exercise the identical solver path.
//
// Each builder FILLS a caller-owned FrameModel& (never returns by value) so the
// Material*/Section* pointers stored in Members remain valid for the model's
// lifetime. Geometry/units: N, mm, MPa.
//
#include "FrameCore/FrameModel.h"

namespace frame { namespace fixtures {

// Resets the model and stores ONE material + ONE section, both at index 0 (matIdx = secIdx = 0).
// Members/shells reference them by index, so there is no pointer capture / reserve-before-push.
inline void prepMatSec(FrameModel& m, const Material& mat, const Section& sec) {
    m = FrameModel{};                 // reset
    m.materials.push_back(mat);       // index 0
    m.sections.push_back(sec);        // index 0
}

// F1: horizontal cantilever along +X. node0 encastre at origin, node1 free at
// (L,0,0). Tip point load P in global +Z at node1.  delta = P L^3 / (3 E Iz).
inline void cantileverTipLoad(FrameModel& m, real P, real L, const Material& mat, const Section& sec) {
    prepMatSec(m, mat, sec);
    Node n0(0, 0, 0, 0); n0.fixAll();
    Node n1(1, L, 0, 0);
    m.nodes   = { n0, n1 };
    m.members = { Member(0, 0, 1, 0, 0) };
    NodalLoad nl; nl.node = 1; nl.comp[Uz] = P;
    m.nodalLoads = { nl };
}

// F2: simply-supported beam along +X, split into 3 nodes so midspan is a node.
// node0 = pin (Ux,Uy,Uz,Rx), node2 = roller (Uy,Uz). UDL w downward (local -y =
// global -Z) on both half-members.  midspan delta = 5 w L^4 / (384 E Iz).
inline void simplySupportedUDL(FrameModel& m, real w, real L, const Material& mat, const Section& sec) {
    prepMatSec(m, mat, sec);
    Node n0(0, 0,     0, 0); n0.fixed[Ux] = n0.fixed[Uy] = n0.fixed[Uz] = n0.fixed[Rx] = true;
    Node n1(1, L / 2, 0, 0);
    Node n2(2, L,     0, 0); n2.fixed[Uy] = n2.fixed[Uz] = true;
    m.nodes   = { n0, n1, n2 };
    m.members = { Member(0, 0, 1, 0, 0), Member(1, 1, 2, 0, 0) };
    MemberUDL u0; u0.member = 0; u0.w_local = { 0, -w, 0 };
    MemberUDL u1; u1.member = 1; u1.w_local = { 0, -w, 0 };
    m.memberUDLs = { u0, u1 };
}

// F4 (extra): vertical column along +Z (forces the refVec degeneracy fallback).
// node0 encastre, node1 at (0,0,h). Gravity load P downward at node1.
// axial shortening = P h / (E A); axial force must read N > 0 (compression).
inline void axialColumn(FrameModel& m, real P, real h, const Material& mat, const Section& sec) {
    prepMatSec(m, mat, sec);
    Node n0(0, 0, 0, 0); n0.fixAll();
    Node n1(1, 0, 0, h);
    m.nodes   = { n0, n1 };
    m.members = { Member(0, 0, 1, 0, 0) };
    NodalLoad nl; nl.node = 1; nl.comp[Uz] = -P;
    m.nodalLoads = { nl };
}

// F3: under-constrained single member. node0 translations fixed but rotations
// free (ball joint), node1 fully free -> rigid-body rotation -> K singular.
inline void mechanism(FrameModel& m, const Material& mat, const Section& sec) {
    prepMatSec(m, mat, sec);
    Node n0(0, 0, 0, 0); n0.pinTranslations();
    Node n1(1, 1000.0, 0, 0);
    m.nodes   = { n0, n1 };
    m.members = { Member(0, 0, 1, 0, 0) };
    NodalLoad nl; nl.node = 1; nl.comp[Uz] = 10.0;
    m.nodalLoads = { nl };
}

// Propped cantilever via a MEMBER RELEASE: a clamped-clamped beam (3 nodes, midspan
// free) whose bending moment is released at node2's end -> effectively fixed at node0,
// pinned at node2. Under UDL w: |M(node0)| = wL^2/8, M(released end) = 0,
// R_z(node0) = 5wL/8, R_z(node2) = 3wL/8. Requires solve with enableReleases=true.
inline void proppedCantileverRelease(FrameModel& m, real w, real L, const Material& mat, const Section& sec) {
    prepMatSec(m, mat, sec);
    Node n0(0, 0,     0, 0); n0.fixAll();
    Node n1(1, L / 2, 0, 0);                 // midspan free DOF
    Node n2(2, L,     0, 0); n2.fixAll();
    m.nodes = { n0, n1, n2 };
    Member mm0(0, 0, 1, 0, 0);
    Member mm1(1, 1, 2, 0, 0);
    mm1.release = makeRelease(ReleasePreset::HingeJ);   // moment hinge at node2
    m.members = { mm0, mm1 };
    MemberUDL u0; u0.member = 0; u0.w_local = { 0, -w, 0 };
    MemberUDL u1; u1.member = 1; u1.w_local = { 0, -w, 0 };
    m.memberUDLs = { u0, u1 };
}

// A release that leaves a FREE MECHANISM: both torsional ends released with no other
// torsional restraint -> the released sub-block kcc is singular -> solve must flag
// singular (the condensation guard). Requires enableReleases=true.
inline void torsionReleaseMechanism(FrameModel& m, const Material& mat, const Section& sec) {
    prepMatSec(m, mat, sec);
    Node n0(0, 0,      0, 0); n0.fixAll();
    Node n1(1, 1000.0, 0, 0); n1.fixed[Uy] = n1.fixed[Uz] = true;
    m.nodes = { n0, n1 };
    Member mem(0, 0, 1, 0, 0);
    mem.release = { { false, false, false, true,  false, false,    // Rx_i released
                      false, false, false, true,  false, false } };// Rx_j released
    m.members = { mem };
    NodalLoad nl; nl.node = 1; nl.comp[Rx] = 100.0;
    m.nodalLoads = { nl };
}

// #4 curved member: a quarter-circle cantilever in the global X-Y plane, discretized
// into nSeg straight segments. Fixed (encastre) at (R,0,0), free at (0,R,0); in-plane
// tip load P in global +Y. Bending-only Castigliano for a thin arc gives the closed
// form  delta_Y(tip) = (pi/4) P R^3 / (E I).  The faceted polygon approximates the
// arc and the tip deflection converges to that value as nSeg grows. Use a section with
// Iy == Iz (square/circular) so "E I" is unambiguous (in-plane bending uses Iy here).
inline void circularArchCantilever(FrameModel& m, real R, int nSeg, real P,
                                   const Material& mat, const Section& sec) {
    prepMatSec(m, mat, sec);
    const real kPi = 3.14159265358979323846;   // NOTE: 'PI' is a UE macro -> use a local name
    m.nodes.clear();
    m.members.clear();
    for (int k = 0; k <= nSeg; ++k) {
        const real th = (kPi / 2.0) * (real(k) / real(nSeg));
        Node n(k, R * std::cos(th), R * std::sin(th), 0);
        if (k == 0) n.fixAll();
        m.nodes.push_back(n);
    }
    for (int k = 0; k < nSeg; ++k)
        m.members.push_back(Member(k, k, k + 1, 0, 0));
    NodalLoad nl; nl.node = nSeg; nl.comp[Uy] = P;
    m.nodalLoads = { nl };
}

// Bare geometry (NO loads) for self-weight tests — the caller applies addSelfWeight().
// cantileverBare: horizontal cantilever along +X, encastre at origin.
inline void cantileverBare(FrameModel& m, real L, const Material& mat, const Section& sec) {
    prepMatSec(m, mat, sec);
    Node n0(0, 0, 0, 0); n0.fixAll();
    Node n1(1, L, 0, 0);
    m.nodes   = { n0, n1 };
    m.members = { Member(0, 0, 1, 0, 0) };
}

// simplySupportedBare: pin at x=0 (Ux,Uy,Uz,Rx), roller at x=L (Uy,Uz), midspan node.
inline void simplySupportedBare(FrameModel& m, real L, const Material& mat, const Section& sec) {
    prepMatSec(m, mat, sec);
    Node n0(0, 0,     0, 0); n0.fixed[Ux] = n0.fixed[Uy] = n0.fixed[Uz] = n0.fixed[Rx] = true;
    Node n1(1, L / 2, 0, 0);
    Node n2(2, L,     0, 0); n2.fixed[Uy] = n2.fixed[Uz] = true;
    m.nodes   = { n0, n1, n2 };
    m.members = { Member(0, 0, 1, 0, 0), Member(1, 1, 2, 0, 0) };
}

// Simply-supported beam discretized into n equal segments (n+1 nodes along +X). Pin at
// node 0 (Ux,Uy,Uz,Rx), roller at node n (Uy,Uz). For influence-line / moving-load tests.
inline void simplySupportedBeamN(FrameModel& m, int n, real L, const Material& mat, const Section& sec) {
    prepMatSec(m, mat, sec);
    m.nodes.clear(); m.members.clear();
    for (int i = 0; i <= n; ++i) {
        Node nd(i, L * i / n, 0, 0);
        if (i == 0) { nd.fixed[Ux] = nd.fixed[Uy] = nd.fixed[Uz] = nd.fixed[Rx] = true; }
        if (i == n) { nd.fixed[Uy] = nd.fixed[Uz] = true; }
        m.nodes.push_back(nd);
    }
    for (int i = 0; i < n; ++i) m.members.push_back(Member(i, i, i + 1, 0, 0));
}

// Cantilever beam discretized into n equal segments (encastre at node 0). For modal /
// buckling tests. omega_1 = 1.875^2 * sqrt(EI/(rho*A*L^4)).
inline void cantileverBeamN(FrameModel& m, int n, real L, const Material& mat, const Section& sec) {
    prepMatSec(m, mat, sec);
    m.nodes.clear(); m.members.clear();
    for (int i = 0; i <= n; ++i) { Node nd(i, L * i / n, 0, 0); if (i == 0) nd.fixAll(); m.nodes.push_back(nd); }
    for (int i = 0; i < n; ++i) m.members.push_back(Member(i, i, i + 1, 0, 0));
}

// Two-span continuous beam (A - midL - B - midR - C), equal spans L, NO loads (the caller
// applies live-load PATTERNS via member UDLs). Members 0,1 = left span; 2,3 = right span.
// Supports: A pin (Ux,Uy,Uz,Rx), interior B and end C restrain Uy,Uz (mirror of the F2
// simply-supported convention). Node B (index 2) is the interior support.
inline void twoSpanContinuous(FrameModel& m, real L, const Material& mat, const Section& sec) {
    prepMatSec(m, mat, sec);
    Node a (0, 0.0,         0, 0); a.fixed[Ux] = a.fixed[Uy] = a.fixed[Uz] = a.fixed[Rx] = true;
    Node ml(1, L / 2.0,     0, 0);
    Node b (2, L,           0, 0); b.fixed[Uy] = b.fixed[Uz] = true;
    Node mr(3, 3.0 * L / 2, 0, 0);
    Node c (4, 2.0 * L,     0, 0); c.fixed[Uy] = c.fixed[Uz] = true;
    m.nodes   = { a, ml, b, mr, c };
    m.members = { Member(0, 0, 1, 0, 0), Member(1, 1, 2, 0, 0),
                  Member(2, 2, 3, 0, 0), Member(3, 3, 4, 0, 0) };
}

// Fixed-fixed beam with a prescribed SETTLEMENT delta (downward) imposed at the far end.
// 3 nodes (free midspan node) so there ARE free DOFs to solve; the beam element is cubic-
// exact for a load-free span, so the end forces are exact regardless of the mesh.
// Support-settlement oracle: end moment 6*E*I*delta/L^2, reaction 12*E*I*delta/L^3.
inline void clampedSettlement(FrameModel& m, real L, real delta, const Material& mat, const Section& sec) {
    prepMatSec(m, mat, sec);
    Node n0(0, 0,     0, 0); n0.fixAll();
    Node n1(1, L / 2, 0, 0);                                   // free midspan node
    Node n2(2, L,     0, 0); n2.fixAll(); n2.prescribed[Uz] = -delta;
    m.nodes   = { n0, n1, n2 };
    m.members = { Member(0, 0, 1, 0, 0), Member(1, 1, 2, 0, 0) };
}

// P-Delta cantilever COLUMN along +Z (vertical -> exercises the refVec degeneracy fallback).
// base (node0) encastre, tip (node nElem) free; meshed into nElem equal beams. The tip carries an
// AXIAL compression P (global -Z) plus a LATERAL load H (global +X), so the whole column is in
// uniform compression N = P. Second-order tip sway (beam-column theory):
//   delta = H (tan(kL) - kL)/(P k),  k = sqrt(P/EI),  limit H L^3/(3 E I) as P -> 0.
// Euler critical load  P_cr = pi^2 E I / (4 L^2).  Use a SQUARE section so EI is unambiguous
// (Iy == Iz). Tip sway is read on Ux. For modal/buckling/P-Delta tests.
inline void pdeltaColumn(FrameModel& m, int nElem, real L, real P, real H,
                         const Material& mat, const Section& sec) {
    prepMatSec(m, mat, sec);
    m.nodes.clear(); m.members.clear();
    for (int k = 0; k <= nElem; ++k) {
        Node nd(k, 0, 0, L * real(k) / real(nElem));
        if (k == 0) nd.fixAll();
        m.nodes.push_back(nd);
    }
    for (int k = 0; k < nElem; ++k) m.members.push_back(Member(k, k, k + 1, 0, 0));
    NodalLoad nl; nl.node = nElem; nl.comp[Uz] = -P; nl.comp[Ux] = H;
    m.nodalLoads = { nl };
}

// S9 elastica: PLANAR cantilever along +X in the global XY plane, n equal beam elements. Encastre at
// node 0; tip transverse point load P (+Y) at node n. Out-of-plane DOFs (Uz,Rx,Ry) restrained at every
// node so the problem is purely planar (large-deflection bending about global Z via Iz). The end-loaded
// elastica (Bisshopp-Drucker / Mattiasson 1981) gives tip delta_v/L (transverse) and delta_h/L
// (axial shortening) vs alpha = P L^2 / (E Iz). Used by runCorotational (S9).
inline void cantileverPlanarTipShearN(FrameModel& m, int n, real L, real P,
                                      const Material& mat, const Section& sec) {
    prepMatSec(m, mat, sec);
    m.nodes.clear(); m.members.clear();
    for (int k = 0; k <= n; ++k) {
        Node nd(k, L * real(k) / real(n), 0, 0);
        nd.fixed[Uz] = nd.fixed[Rx] = nd.fixed[Ry] = true;   // planar: out-of-plane restrained
        if (k == 0) nd.fixAll();
        m.nodes.push_back(nd);
    }
    for (int k = 0; k < n; ++k) m.members.push_back(Member(k, k, k + 1, 0, 0));
    NodalLoad nl; nl.node = n; nl.comp[Uy] = P;
    m.nodalLoads = { nl };
}

// X-braced portal in the global X-Z plane (out-of-plane Uy pinned at the free top nodes), for
// tension-only tests. Stocky columns/beam (section index 0) form a stable moment frame; the two
// SLENDER diagonals (section index 1) are flagged tension-only. node0/1 = base (encastre),
// node2/3 = top (Uy pinned). Lateral H (+X) and downward V at the top. Members: 0=col L, 1=col R,
// 2=beam, 3=brace A (0->3), 4=brace B (1->2). Under pure lateral H one diagonal compresses and, as
// a tension-only member, drops out -> the converged state equals the model with that brace omitted.
inline void xBracedPortal(FrameModel& m, real H, real V,
                          const Material& mat, const Section& stocky, const Section& brace) {
    m = FrameModel{};
    m.materials.push_back(mat);
    m.sections.push_back(stocky);   // index 0: columns + beam
    m.sections.push_back(brace);    // index 1: slender diagonals
    Node n0(0, 0,    0, 0);    n0.fixAll();
    Node n1(1, 6000, 0, 0);    n1.fixAll();
    Node n2(2, 0,    0, 3000); n2.fixed[Uy] = true;
    Node n3(3, 6000, 0, 3000); n3.fixed[Uy] = true;
    m.nodes = { n0, n1, n2, n3 };
    auto add = [&](int id, int i, int j, int sec, bool tonly) {
        Member mm(id, i, j, 0, sec);
        mm.refVec = Vec3(0, 1, 0);     // +Y is out of plane for every member here (no degeneracy)
        mm.tensionOnly = tonly;
        m.members.push_back(mm);
    };
    add(0, 0, 2, 0, false);   // column L
    add(1, 1, 3, 0, false);   // column R
    add(2, 2, 3, 0, false);   // beam
    add(3, 0, 3, 1, true);    // brace A (tension-only diagonal)
    add(4, 1, 2, 1, true);    // brace B (tension-only diagonal)
    NodalLoad l2; l2.node = 2; l2.comp[Ux] = H; l2.comp[Uz] = -V;
    NodalLoad l3; l3.node = 3; l3.comp[Uz] = -V;
    m.nodalLoads = { l2, l3 };
}

// Classic 10-bar cantilever truss (Schmit/Berke sizing benchmark) in the global X-Z plane. Two
// bays of `bay`, height `bay`; nodes 5,6 are supported at x=0; downward (global -Z) point loads P
// at nodes 2 and 4. ALL ten bars SHARE section index 0 (a uniform start area A0) -- runSizeOptimization
// expands a private section per member. Planar: Uy is pinned everywhere, supports pin Ux,Uy,Uz.
// Bars (standard numbering, by node id): 1:5-3 2:3-1 3:6-4 4:4-2 5:3-4 6:1-2 7:5-4 8:6-3 9:3-2 10:4-1.
// FSD optimum (Haftka & Gurdal): weight ~ 1593.2 lb, areas in^2
//   [7.94, 0.10, 8.06, 3.94, 0.10, 0.10, 5.74, 5.57, 5.57, 0.10] (bars 2,5,6,10 at the A_min bound).
inline void tenBarTruss(FrameModel& m, real bay, real A0, real P, const Material& mat) {
    m = FrameModel{};
    m.materials.push_back(mat);                                                  // index 0
    m.sections.push_back(Section::Rectangular(std::sqrt(A0), std::sqrt(A0)));    // index 0: square, area A0
    auto MakeNode = [](int id, real x, real z, bool support) {
        Node n(id, x, 0, z);
        if (support) { n.fixed[Ux] = n.fixed[Uy] = n.fixed[Uz] = true; }
        else         { n.fixed[Uy] = true; }                                     // planar X-Z problem
        return n;
    };
    m.nodes = {
        MakeNode(1, 2 * bay, bay, false), MakeNode(2, 2 * bay, 0, false),
        MakeNode(3, bay,     bay, false), MakeNode(4, bay,     0, false),
        MakeNode(5, 0,       bay, true),  MakeNode(6, 0,       0, true),
    };
    const int bars[10][2] = { {5,3},{3,1},{6,4},{4,2},{3,4},{1,2},{5,4},{6,3},{3,2},{4,1} };
    for (int k = 0; k < 10; ++k) {
        Member mem(k, bars[k][0], bars[k][1], 0, 0);
        mem.refVec = Vec3(0, 1, 0);            // member local axes lie in the X-Z plane (+Y out of plane)
        m.members.push_back(mem);
    }
    NodalLoad l2; l2.node = 2; l2.comp[Uz] = -P;
    NodalLoad l4; l4.node = 4; l4.comp[Uz] = -P;
    m.nodalLoads = { l2, l4 };
}

// ---------------------------------------------------------------------------
// Shell (MITC4) fixtures. Geometry in the global X-Y plane (facet normal +Z), so
// at milestone 2 (plate bending only) the in-plane DOFs (Ux,Uy,Rz) are restrained
// at every node and only the bending action (Uz,Rx,Ry) is exercised.
// node(i,j) = j*(nx+1) + i, i in [0,nx], j in [0,ny].
// ---------------------------------------------------------------------------

// Simply-supported (soft: w=0 on the boundary) square plate under uniform pressure q
// (downward). Kirchhoff thin-plate theory: w_c = 0.00406 q a^4 / D, D = E t^3/[12(1-nu^2)]
// for nu = 0.3. The MITC4 mesh converges to this from a faceted bilinear field.
inline void squarePlateShell(FrameModel& m, real a, real t, int n, real q,
                             const Material& mat) {
    m = FrameModel{};
    m.materials.reserve(1);
    m.materials.push_back(mat);
    const real h = a / n;
    auto gid = [n](int i, int j) { return j * (n + 1) + i; };
    for (int j = 0; j <= n; ++j)
        for (int i = 0; i <= n; ++i) {
            Node nd(gid(i, j), i * h, j * h, 0);
            nd.fixed[Ux] = nd.fixed[Uy] = nd.fixed[Rz] = true;   // in-plane restrained (milestone 2)
            const bool edge = (i == 0 || i == n || j == 0 || j == n);
            if (edge) nd.fixed[Uz] = true;                        // simple support: w = 0
            m.nodes.push_back(nd);
        }
    int sid = 0;
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) {
            ShellQuad sh(sid, gid(i, j), gid(i + 1, j), gid(i + 1, j + 1), gid(i, j + 1), 0, t);
            m.shells.push_back(sh);
            ShellPressure sp; sp.shell = sid; sp.p = -q;          // downward (local -z)
            m.shellPressures.push_back(sp);
            ++sid;
        }
}

// Plate bending PATCH TEST on a 2x2 mesh. `skew` shears the whole grid into a uniform
// PARALLELOGRAM mesh (x += skew*y), so every element is a congruent parallelogram
// (constant Jacobian) — the regime in which MITC4 is GUARANTEED to reproduce a
// constant-curvature state exactly. The exact Kirchhoff cylindrical-bending field
//   w = 0.5*c*x^2,  Ry = bx = -c*x,  Rx = 0
// is prescribed on all 8 boundary nodes (using each node's actual x); the interior
// node is free. A correct element reproduces a CONSTANT moment field everywhere:
//   Mxx = -D*c,  Myy = -nu*D*c,  Mxy = 0,  Q = 0    (D = E t^3/[12(1-nu^2)]).
// (General, non-parallelogram quads instead show an O(h) patch residual that vanishes
//  under refinement — see the square-plate convergence and the milestone-4 benchmarks.)
inline void platePatchCylindrical(FrameModel& m, real a, real t, real skew, real c,
                                  const Material& mat) {
    m = FrameModel{};
    m.materials.reserve(1);
    m.materials.push_back(mat);
    const real h = a / 2;
    auto isBoundary = [](int k) { return k != 4; };   // only centre (k=4) is interior
    for (int j = 0; j < 3; ++j)
        for (int i = 0; i < 3; ++i) {
            const int k = j * 3 + i;
            const real y = j * h;
            const real x = i * h + skew * y;            // parallelogram shear
            Node nd(k, x, y, 0);
            nd.fixed[Ux] = nd.fixed[Uy] = nd.fixed[Rz] = true;     // in-plane restrained
            if (isBoundary(k)) {
                nd.fixed[Uz] = nd.fixed[Rx] = nd.fixed[Ry] = true; // prescribe the exact field
                nd.prescribed[Uz] = 0.5 * c * x * x;               // w  = 0.5 c x^2
                nd.prescribed[Rx] = 0.0;                           // Rx = -by = 0
                nd.prescribed[Ry] = -c * x;                        // Ry =  bx = -c x
            }
            m.nodes.push_back(nd);
        }
    auto gid = [](int i, int j) { return j * 3 + i; };
    int sid = 0;
    for (int j = 0; j < 2; ++j)
        for (int i = 0; i < 2; ++i)
            m.shells.push_back(ShellQuad(sid++, gid(i, j), gid(i + 1, j),
                                         gid(i + 1, j + 1), gid(i, j + 1), 0, t));
}

// Membrane PATCH TEST on a 2x2 parallelogram mesh (x += skew*y). A constant-strain
// in-plane field u = gx*x, v = 0 (eps_x = gx, others 0; true rotation omega = 0) is
// prescribed on the 8 boundary nodes; the interior node's in-plane DOFs are free; the
// out-of-plane DOFs (Uz,Rx,Ry) are restrained. A correct membrane reproduces constant
// stress resultants: Nxx = t*E/(1-nu^2)*gx, Nyy = nu*Nxx, Nxy = 0.
inline void membranePatch(FrameModel& m, real a, real t, real skew, real gx,
                          const Material& mat) {
    m = FrameModel{};
    m.materials.reserve(1);
    m.materials.push_back(mat);
    const real h = a / 2;
    auto isBoundary = [](int k) { return k != 4; };
    for (int j = 0; j < 3; ++j)
        for (int i = 0; i < 3; ++i) {
            const int k = j * 3 + i;
            const real y = j * h;
            const real x = i * h + skew * y;
            Node nd(k, x, y, 0);
            nd.fixed[Uz] = nd.fixed[Rx] = nd.fixed[Ry] = true;     // out-of-plane restrained
            if (isBoundary(k)) {
                nd.fixed[Ux] = nd.fixed[Uy] = nd.fixed[Rz] = true; // prescribe constant-strain field
                nd.prescribed[Ux] = gx * x;
                nd.prescribed[Uy] = 0.0;
                nd.prescribed[Rz] = 0.0;                            // omega = 0
            }
            m.nodes.push_back(nd);
        }
    auto gid = [](int i, int j) { return j * 3 + i; };
    int sid = 0;
    for (int j = 0; j < 2; ++j)
        for (int i = 0; i < 2; ++i)
            m.shells.push_back(ShellQuad(sid++, gid(i, j), gid(i + 1, j),
                                         gid(i + 1, j + 1), gid(i, j + 1), 0, t));
}

// Slender in-plane cantilever MEMBRANE (S8-8a: QM6-vs-Q4 membrane-locking pin). A thin
// rectangle L x H in the xy plane (L >> H), clamped at x=0, end-loaded by a transverse
// in-plane shear P (+y) lumped uniformly on the x=L edge. Out-of-plane Uz,Rx,Ry are
// restrained everywhere (pure membrane); the clamped edge fixes Ux,Uy,Rz. The bilinear Q4
// membrane LOCKS in coarse in-plane bending (parasitic shear -> far too stiff); the QM6
// bubbles release it. Reference: Euler-Bernoulli tip P L^3 / (3 E I), I = t H^3 / 12.
inline void slenderMembraneCantilever(FrameModel& m, real L, real H, real t,
                                      int nx, int ny, real P, const Material& mat) {
    m = FrameModel{};
    m.materials.reserve(1);
    m.materials.push_back(mat);
    const real hx = L / nx, hy = H / ny;
    auto gid = [nx](int i, int j) { return j * (nx + 1) + i; };
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nx; ++i) {
            Node nd(gid(i, j), i * hx, j * hy, 0);
            nd.fixed[Uz] = nd.fixed[Rx] = nd.fixed[Ry] = true;          // pure membrane
            if (i == 0) { nd.fixed[Ux] = nd.fixed[Uy] = nd.fixed[Rz] = true; }  // clamped end
            m.nodes.push_back(nd);
        }
    int sid = 0;
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i)
            m.shells.push_back(ShellQuad(sid++, gid(i, j), gid(i + 1, j),
                                         gid(i + 1, j + 1), gid(i, j + 1), 0, t));
    const real per = P / (ny + 1);                                      // lump shear uniformly
    for (int j = 0; j <= ny; ++j) {
        NodalLoad nl; nl.node = gid(nx, j); nl.comp[Uy] = per;
        m.nodalLoads.push_back(nl);
    }
}

// Cook's skew membrane (standard QM6 / distorted-membrane benchmark). Trapezoid with corners
// A(0,0) D(0,44) on the clamped left edge, B(48,44) C(48,60) on the loaded right edge; total
// in-plane shear P (+y) lumped uniformly on the right edge. Out-of-plane restrained (membrane).
// The tip (corner C = node(n,n)) vertical displacement converges to ~23.96 (E=1, nu=1/3, t=1,
// P=1). Coarse Q4 underestimates badly (membrane + skew locking); QM6 is far closer.
inline void cooksMembrane(FrameModel& m, int n, real P, real t, const Material& mat) {
    m = FrameModel{};
    m.materials.reserve(1);
    m.materials.push_back(mat);
    auto gid = [n](int i, int j) { return j * (n + 1) + i; };
    for (int j = 0; j <= n; ++j)
        for (int i = 0; i <= n; ++i) {
            const real xi = real(i) / n, eta = real(j) / n;
            const real x = 48.0 * xi;
            const real y = (1.0 - xi) * (44.0 * eta) + xi * (44.0 + 16.0 * eta);
            Node nd(gid(i, j), x, y, 0);
            nd.fixed[Uz] = nd.fixed[Rx] = nd.fixed[Ry] = true;          // pure membrane
            if (i == 0) { nd.fixed[Ux] = nd.fixed[Uy] = nd.fixed[Rz] = true; }  // clamped edge
            m.nodes.push_back(nd);
        }
    int sid = 0;
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            m.shells.push_back(ShellQuad(sid++, gid(i, j), gid(i + 1, j),
                                         gid(i + 1, j + 1), gid(i, j + 1), 0, t));
    // Consistent lumping of a uniform edge traction (total P): half weight at the two end nodes.
    for (int j = 0; j <= n; ++j) {
        const real wj = (j == 0 || j == n) ? 0.5 : 1.0;
        NodalLoad nl; nl.node = gid(n, j); nl.comp[Uy] = wj * P / n;
        m.nodalLoads.push_back(nl);
    }
}

// Fully CLAMPED square plate (boundary nodes encastre, ALL 6 DOF) under uniform
// pressure q. The interior nodes leave every DOF free — including the in-plane
// translations (held by the membrane) and the drilling Rz (held only by the drilling
// stiffness). So this is the gate that the drilling treatment removes the in-plane
// rotational zero-energy mode: a flat clamped shell must solve NON-SINGULAR.
inline void clampedPlateShell(FrameModel& m, real a, real t, int n, real q,
                              const Material& mat) {
    m = FrameModel{};
    m.materials.reserve(1);
    m.materials.push_back(mat);
    const real h = a / n;
    auto gid = [n](int i, int j) { return j * (n + 1) + i; };
    for (int j = 0; j <= n; ++j)
        for (int i = 0; i <= n; ++i) {
            Node nd(gid(i, j), i * h, j * h, 0);
            const bool edge = (i == 0 || i == n || j == 0 || j == n);
            if (edge) nd.fixAll();          // clamped boundary; interior all-free
            m.nodes.push_back(nd);
        }
    int sid = 0;
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) {
            m.shells.push_back(ShellQuad(sid, gid(i, j), gid(i + 1, j),
                                         gid(i + 1, j + 1), gid(i, j + 1), 0, t));
            ShellPressure sp; sp.shell = sid; sp.p = -q;
            m.shellPressures.push_back(sp);
            ++sid;
        }
}

// Rigidly rotate a whole model's node coordinates by R (rows = R[0],R[1],R[2]).
// Used by the shell rotation-invariance check: ShellPressure follows the facet normal
// automatically and an encastre boundary is frame-invariant, so the solution must
// simply rotate with the model (displacement magnitudes preserved).
inline void rotateModelRigid(FrameModel& m, const real R[3][3]) {
    for (auto& nd : m.nodes) {
        const real x = nd.pos.x, y = nd.pos.y, z = nd.pos.z;
        nd.pos.x = R[0][0] * x + R[0][1] * y + R[0][2] * z;
        nd.pos.y = R[1][0] * x + R[1][1] * y + R[1][2] * z;
        nd.pos.z = R[2][0] * x + R[2][1] * y + R[2][2] * z;
    }
}

// ---------------------------------------------------------------------------
// Curved-shell benchmarks (MacNeal-Harder). These exercise the flat-facet
// approximation of a curved surface + the membrane/bending coupling end-to-end.
// ---------------------------------------------------------------------------

// Scordelis-Lo roof (a classic membrane-dominated shell benchmark). A cylindrical
// barrel-vault segment under self-weight, supported by rigid diaphragms at the curved
// ends, free along the straight edges. Modelled as a QUARTER with two symmetry planes:
//   crown   phi=0      (x=0 plane):     u_x=0, Ry=0, Rz=0
//   midspan y=L/2      (y=const plane): u_y=0, Rx=0, Rz=0
//   diaphragm y=0:                       u_x=0, u_z=0      (free axial + rotations)
//   free edge phi=phi0:                  unrestrained
// Reference: vertical deflection at the free-edge midspan = 0.3024 (downward).
// Self-weight g (force/area) is lumped to nodes as global -Z nodal loads.
inline void scordelisLoRoof(FrameModel& m, real R, real L, real phi0, real t, real g,
                            int nf, int ny, const Material& mat) {
    m = FrameModel{};
    m.materials.reserve(1);
    m.materials.push_back(mat);
    const int NN = (nf + 1) * (ny + 1);
    auto gid = [nf](int i, int j) { return j * (nf + 1) + i; };
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nf; ++i) {
            const real phi = phi0 * (real)i / nf;
            const real y   = (L * 0.5) * (real)j / ny;
            Node nd(gid(i, j), R * std::sin(phi), y, R * std::cos(phi));
            if (j == 0)  { nd.fixed[Ux] = nd.fixed[Uz] = true; }                 // diaphragm
            if (j == ny) { nd.fixed[Uy] = nd.fixed[Rx] = nd.fixed[Rz] = true; }  // midspan symmetry
            if (i == 0)  { nd.fixed[Ux] = nd.fixed[Ry] = nd.fixed[Rz] = true; }  // crown symmetry
            m.nodes.push_back(nd);
        }
    std::vector<real> fz((size_t)NN, 0.0);
    int sid = 0;
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nf; ++i) {
            const int a = gid(i, j), b = gid(i + 1, j), c = gid(i + 1, j + 1), d = gid(i, j + 1);
            m.shells.push_back(ShellQuad(sid++, a, b, c, d, 0, t));
            const Vec3 pa = m.nodes[a].pos, pb = m.nodes[b].pos, pc = m.nodes[c].pos, pd = m.nodes[d].pos;
            const real area = 0.5 * norm(cross(pc - pa, pd - pb));
            const real fn = -g * area * 0.25;                  // self-weight per node (global -Z)
            fz[a] += fn; fz[b] += fn; fz[c] += fn; fz[d] += fn;
        }
    for (int k = 0; k < NN; ++k)
        if (fz[k] != 0.0) { NodalLoad nl; nl.node = k; nl.comp[Uz] = fz[k]; m.nodalLoads.push_back(nl); }
}

// Pinched cylinder with rigid end diaphragms (a very demanding inextensional-bending
// benchmark). A short cylinder is pinched by two opposite radial point loads P at
// mid-length. Modelled as a 1/8 OCTANT (axis = z); the loaded generator is theta=0.
//   theta=0   (x-z plane, normal y):  Uy=0, Rx=0, Rz=0
//   theta=90  (y-z plane, normal x):  Ux=0, Ry=0, Rz=0
//   midspan z=L/2 (normal z):         Uz=0, Rx=0, Ry=0
//   diaphragm z=0:                     Ux=0, Uy=0      (rigid diaphragm; axial free)
// The load node (theta=0, z=L/2) sits on two symmetry planes, so it carries P/4 inward
// (-x). Reference radial deflection under the load = 1.8248e-5. MITC4 flat facets
// converge to it from BELOW (this test is notoriously slow-converging).
inline void pinchedCylinder(FrameModel& m, real R, real L, real t, real P,
                            int nth, int nz, const Material& mat) {
    m = FrameModel{};
    m.materials.reserve(1);
    m.materials.push_back(mat);
    const real kPi = 3.14159265358979323846;
    auto gid = [nth](int i, int j) { return j * (nth + 1) + i; };
    for (int j = 0; j <= nz; ++j)
        for (int i = 0; i <= nth; ++i) {
            const real th = (0.5 * kPi) * (real)i / nth;       // theta in [0, 90deg]
            const real z  = (L * 0.5) * (real)j / nz;          // z in [0, L/2]
            Node nd(gid(i, j), R * std::cos(th), R * std::sin(th), z);
            if (i == 0)   { nd.fixed[Uy] = nd.fixed[Rx] = nd.fixed[Rz] = true; }  // theta=0 symmetry
            if (i == nth) { nd.fixed[Ux] = nd.fixed[Ry] = nd.fixed[Rz] = true; }  // theta=90 symmetry
            if (j == nz)  { nd.fixed[Uz] = nd.fixed[Rx] = nd.fixed[Ry] = true; }  // midspan symmetry
            if (j == 0)   { nd.fixed[Ux] = nd.fixed[Uy] = true; }                 // diaphragm
            m.nodes.push_back(nd);
        }
    int sid = 0;
    for (int j = 0; j < nz; ++j)
        for (int i = 0; i < nth; ++i)
            m.shells.push_back(ShellQuad(sid++, gid(i, j), gid(i + 1, j),
                                         gid(i + 1, j + 1), gid(i, j + 1), 0, t));
    NodalLoad nl; nl.node = gid(0, nz); nl.comp[Ux] = -0.25 * P;   // P/4 inward at load node
    m.nodalLoads.push_back(nl);
}

}} // namespace frame::fixtures
