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

inline void prepMatSec(FrameModel& m, const Material& mat, const Section& sec,
                       const Material*& pm, const Section*& ps) {
    m = FrameModel{};                 // reset
    m.materials.reserve(1);
    m.sections.reserve(1);
    m.materials.push_back(mat);
    m.sections.push_back(sec);
    pm = &m.materials.back();
    ps = &m.sections.back();
}

// F1: horizontal cantilever along +X. node0 encastre at origin, node1 free at
// (L,0,0). Tip point load P in global +Z at node1.  delta = P L^3 / (3 E Iz).
inline void cantileverTipLoad(FrameModel& m, real P, real L, const Material& mat, const Section& sec) {
    const Material* pm; const Section* ps; prepMatSec(m, mat, sec, pm, ps);
    Node n0(0, 0, 0, 0); n0.fixAll();
    Node n1(1, L, 0, 0);
    m.nodes   = { n0, n1 };
    m.members = { Member(0, 0, 1, pm, ps) };
    NodalLoad nl; nl.node = 1; nl.comp[Uz] = P;
    m.nodalLoads = { nl };
}

// F2: simply-supported beam along +X, split into 3 nodes so midspan is a node.
// node0 = pin (Ux,Uy,Uz,Rx), node2 = roller (Uy,Uz). UDL w downward (local -y =
// global -Z) on both half-members.  midspan delta = 5 w L^4 / (384 E Iz).
inline void simplySupportedUDL(FrameModel& m, real w, real L, const Material& mat, const Section& sec) {
    const Material* pm; const Section* ps; prepMatSec(m, mat, sec, pm, ps);
    Node n0(0, 0,     0, 0); n0.fixed[Ux] = n0.fixed[Uy] = n0.fixed[Uz] = n0.fixed[Rx] = true;
    Node n1(1, L / 2, 0, 0);
    Node n2(2, L,     0, 0); n2.fixed[Uy] = n2.fixed[Uz] = true;
    m.nodes   = { n0, n1, n2 };
    m.members = { Member(0, 0, 1, pm, ps), Member(1, 1, 2, pm, ps) };
    MemberUDL u0; u0.member = 0; u0.w_local = { 0, -w, 0 };
    MemberUDL u1; u1.member = 1; u1.w_local = { 0, -w, 0 };
    m.memberUDLs = { u0, u1 };
}

// F4 (extra): vertical column along +Z (forces the refVec degeneracy fallback).
// node0 encastre, node1 at (0,0,h). Gravity load P downward at node1.
// axial shortening = P h / (E A); axial force must read N > 0 (compression).
inline void axialColumn(FrameModel& m, real P, real h, const Material& mat, const Section& sec) {
    const Material* pm; const Section* ps; prepMatSec(m, mat, sec, pm, ps);
    Node n0(0, 0, 0, 0); n0.fixAll();
    Node n1(1, 0, 0, h);
    m.nodes   = { n0, n1 };
    m.members = { Member(0, 0, 1, pm, ps) };
    NodalLoad nl; nl.node = 1; nl.comp[Uz] = -P;
    m.nodalLoads = { nl };
}

// F3: under-constrained single member. node0 translations fixed but rotations
// free (ball joint), node1 fully free -> rigid-body rotation -> K singular.
inline void mechanism(FrameModel& m, const Material& mat, const Section& sec) {
    const Material* pm; const Section* ps; prepMatSec(m, mat, sec, pm, ps);
    Node n0(0, 0, 0, 0); n0.pinTranslations();
    Node n1(1, 1000.0, 0, 0);
    m.nodes   = { n0, n1 };
    m.members = { Member(0, 0, 1, pm, ps) };
    NodalLoad nl; nl.node = 1; nl.comp[Uz] = 10.0;
    m.nodalLoads = { nl };
}

// Propped cantilever via a MEMBER RELEASE: a clamped-clamped beam (3 nodes, midspan
// free) whose bending moment is released at node2's end -> effectively fixed at node0,
// pinned at node2. Under UDL w: |M(node0)| = wL^2/8, M(released end) = 0,
// R_z(node0) = 5wL/8, R_z(node2) = 3wL/8. Requires solve with enableReleases=true.
inline void proppedCantileverRelease(FrameModel& m, real w, real L, const Material& mat, const Section& sec) {
    const Material* pm; const Section* ps; prepMatSec(m, mat, sec, pm, ps);
    Node n0(0, 0,     0, 0); n0.fixAll();
    Node n1(1, L / 2, 0, 0);                 // midspan free DOF
    Node n2(2, L,     0, 0); n2.fixAll();
    m.nodes = { n0, n1, n2 };
    Member mm0(0, 0, 1, pm, ps);
    Member mm1(1, 1, 2, pm, ps);
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
    const Material* pm; const Section* ps; prepMatSec(m, mat, sec, pm, ps);
    Node n0(0, 0,      0, 0); n0.fixAll();
    Node n1(1, 1000.0, 0, 0); n1.fixed[Uy] = n1.fixed[Uz] = true;
    m.nodes = { n0, n1 };
    Member mem(0, 0, 1, pm, ps);
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
    const Material* pm; const Section* ps; prepMatSec(m, mat, sec, pm, ps);
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
        m.members.push_back(Member(k, k, k + 1, pm, ps));
    NodalLoad nl; nl.node = nSeg; nl.comp[Uy] = P;
    m.nodalLoads = { nl };
}

}} // namespace frame::fixtures
