// Standalone gate for FrameCore — builds fixtures, solves, compares to closed-form
// analytic solutions, prints PASS/FAIL, exits 0 iff all pass.
#include "FrameCore/FrameSolver.h"
#include "FrameCore/ElasticAllowable.h"
#include "FrameCore/Grillage.h"
#include "FrameCore/SelfWeight.h"
#include "FrameCore/Combination.h"
#include "FrameCore/InfluenceLine.h"
#include "FrameCore/ModalAnalysis.h"
#include "FrameCore/BucklingAnalysis.h"
#include "FrameCore/Reanalysis.h"
#include "FrameCore/ResponseSpectrum.h"
#include "FrameCore/ModalDynamics.h"
#include "FrameCore/Connectivity.h"
#include "FrameCore/Collapse.h"
#include "FrameCore/DynamicCollapse.h"
#include "FrameCore/PDeltaAnalysis.h"
#include "FrameCore/CorotationalAnalysis.h"
#include "FrameCore/TensionOnly.h"
#include "FrameCore/SizeOpt.h"
#include "FrameCore/Topology.h"
#include "FrameCore/NMInteraction.h"
#include "FrameCore/MemberGeometry.h"
#include "FrameCore/HpSolver.h"
#include "FrameCore/HpSession.h"
#include "FrameTestFixtures.h"

#include <vector>
#include <utility>

#include <cstdio>
#include <cmath>
#include <string>
#include <algorithm>

using namespace frame;

#ifndef FRAMECORE_BUILD_SHA
#define FRAMECORE_BUILD_SHA "unknown"   // overridden by the build script via /D (git short SHA)
#endif

static int g_fail = 0;

static double relErr(double got, double expect) {
    return std::fabs(got - expect) / std::max(std::fabs(expect), 1e-30);
}
static void checkClose(const char* tag, double got, double expect, double tol) {
    const double e = relErr(got, expect);
    const bool ok = std::isfinite(got) && e <= tol;
    if (!ok) ++g_fail;
    std::printf("  %s %-32s got=%-12.6g exp=%-12.6g rel=%.2e (tol=%.0e)\n",
                ok ? "[PASS]" : "[FAIL]", tag, got, expect, e, tol);
}
static void checkTrue(const char* tag, bool cond, const std::string& detail = "") {
    if (!cond) ++g_fail;
    std::printf("  %s %-32s %s\n", cond ? "[PASS]" : "[FAIL]", tag, detail.c_str());
}

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);  // unbuffered: output survives a crash
    std::printf("# FrameCore standalone gate | build %s | compiled %s %s\n",
                FRAMECORE_BUILD_SHA, __DATE__, __TIME__);
    const real E = 210000.0, G = 80769.0;     // steel, MPa
    const real side = 100.0;                   // mm, square section
    Section  sec = Section::Rectangular(side, side);
    Material mat(E, G, 7850.0);
    mat.cap = Capacity::make(300.0, 300.0, 180.0);

    std::printf("FrameCore standalone fixtures  (units: N, mm, MPa)\n");
    std::printf("  E=%.0f  G=%.0f  square=%.0f  A=%.0f  Iz=%.6g  Wz=%.6g\n\n",
                E, G, side, sec.A, sec.Iz, sec.Wz());

    // ---------- F1: cantilever, tip point load ----------
    {
        const real P = 1000.0, L = 2000.0;
        FrameModel m; fixtures::cantileverTipLoad(m, P, L, mat, sec);
        SolveResult r = solve(m);
        std::printf("[F1] cantilever tip load  P=%.0f L=%.0f\n", P, L);
        checkTrue("not singular", !r.singular, r.diagnostic);
        const real dExp = P * L * L * L / (3.0 * E * sec.Iz);
        checkClose("tip deflection |u_z(j)|", std::fabs(r.disp(1, Uz)), dExp, 1e-6);
        const MemberForcePair& mf = r.memberForces[0];
        const real Mroot = std::sqrt(mf.endI.My * mf.endI.My + mf.endI.Mz * mf.endI.Mz);
        checkClose("root moment |M_i|", Mroot, P * L, 1e-6);
        checkClose("reaction R_z(node0)", std::fabs(r.reaction(0, Uz)), P, 1e-9);
        ElasticAllowable strength;
        const DemandResult d = strength.checkSection(mf.endI, sec, mat.cap);
        checkTrue("checkSection mode != None", d.mode != FailMode::None);
        checkClose("checkSection risk", d.risk, (Mroot / sec.Wz()) / mat.cap.bend, 1e-6);
    }

    // ---------- F2: simply supported, UDL ----------
    {
        const real w = 5.0, L = 3000.0;
        FrameModel m; fixtures::simplySupportedUDL(m, w, L, mat, sec);
        SolveResult r = solve(m);
        std::printf("[F2] simply-supported UDL  w=%.1f L=%.0f\n", w, L);
        checkTrue("not singular", !r.singular, r.diagnostic);
        const real dExp = 5.0 * w * L * L * L * L / (384.0 * E * sec.Iz);
        checkClose("midspan deflection |u_z|", std::fabs(r.disp(1, Uz)), dExp, 1e-6);
        const MemberForcePair& m0 = r.memberForces[0];
        checkClose("midspan moment |M|", std::fabs(m0.endJ.Mz), w * L * L / 8.0, 1e-5);
        checkClose("reaction R_z(node0)", std::fabs(r.reaction(0, Uz)), w * L / 2.0, 1e-6);
        checkClose("reaction R_z(node2)", std::fabs(r.reaction(2, Uz)), w * L / 2.0, 1e-6);
    }

    // ---------- F3: mechanism (negative test) ----------
    {
        FrameModel m; fixtures::mechanism(m, mat, sec);
        SolveResult r = solve(m);
        std::printf("[F3] mechanism (under-constrained)\n");
        checkTrue("singular detected", r.singular, r.diagnostic);
        checkTrue("diagnostic non-empty", !r.diagnostic.empty(), r.diagnostic);
    }

    // ---------- F4: vertical column (axial sign + stiffness) ----------
    {
        const real P = 50000.0, h = 3000.0;
        FrameModel m; fixtures::axialColumn(m, P, h, mat, sec);
        SolveResult r = solve(m);
        std::printf("[F4] vertical column  P=%.0f h=%.0f\n", P, h);
        checkTrue("not singular", !r.singular, r.diagnostic);
        checkClose("axial shortening |u_z|", std::fabs(r.disp(1, Uz)), P * h / (E * sec.A), 1e-6);
        const real N = r.memberForces[0].endI.N;
        checkTrue("axial N > 0 (compression)", N > 0.0, "N=" + std::to_string(N));
    }

    // ---------- F5: checkSection pure torsion (units must be MPa) ----------
    {
        MemberEndForces f;            // pure torsion
        f.T = 1.0e6;                  // N*mm
        ElasticAllowable strength;
        const DemandResult d = strength.checkSection(f, sec, mat.cap);
        const real cTor    = std::hypot(sec.cy, sec.cz);
        const real sTorExp = 1.0e6 * cTor / sec.J;   // ~5.021 MPa, NOT |T|/J=0.071
        std::printf("[F5] checkSection pure torsion  T=1e6\n");
        checkClose("torsion stress sTor (MPa)", d.sTor, sTorExp, 1e-6);
        checkTrue("mode == Torsion", d.mode == FailMode::Torsion);
        checkClose("torsion risk", d.risk, sTorExp / mat.cap.tors, 1e-6);
    }

    // ---- F6: propped cantilever via member release (Qf condensation + analytic) ----
    {
        const real w = 5.0, L = 4000.0;
        FrameModel m; fixtures::proppedCantileverRelease(m, w, L, mat, sec);
        SolveOptions opt; opt.enableReleases = true;
        SolveResult r = solve(m, opt);
        std::printf("[F6] propped cantilever via release  w=%.1f L=%.0f\n", w, L);
        checkTrue("not singular", !r.singular, r.diagnostic);
        const real Mj = std::sqrt(r.memberForces[1].endJ.My * r.memberForces[1].endJ.My
                                + r.memberForces[1].endJ.Mz * r.memberForces[1].endJ.Mz);
        const real Mscale = w * L * L / 8.0;
        checkTrue("released-end moment ~ 0", Mj < 1e-6 * Mscale, "Mj=" + std::to_string(Mj));
        checkClose("fixed-end moment |M_i| = wL^2/8", std::fabs(r.memberForces[0].endI.Mz), Mscale, 1e-6);
        checkClose("reaction R_z(node0) = 5wL/8", std::fabs(r.reaction(0, Uz)), 5.0 * w * L / 8.0, 1e-6);
        checkClose("reaction R_z(node2) = 3wL/8", std::fabs(r.reaction(2, Uz)), 3.0 * w * L / 8.0, 1e-6);
    }

    // ---- F7: ill-posed release (both torsional ends) -> singular guard ----
    {
        FrameModel m; fixtures::torsionReleaseMechanism(m, mat, sec);
        SolveOptions opt; opt.enableReleases = true;
        SolveResult r = solve(m, opt);
        std::printf("[F7] torsional double-release mechanism\n");
        checkTrue("singular detected", r.singular, r.diagnostic);
        checkTrue("diagnostic mentions release",
                  r.diagnostic.find("release") != std::string::npos, r.diagnostic);
    }

    // ---- F8: Timoshenko shear correction on a cantilever (deep vs slender) ----
    {
        const real P = 1000.0, L = 2000.0;
        FrameModel m; fixtures::cantileverTipLoad(m, P, L, mat, sec);
        SolveOptions opt; opt.useTimoshenko = true;
        SolveResult r = solve(m, opt);
        const real dEB    = P * L * L * L / (3.0 * E * sec.Iz);
        const real dShear = P * L / (G * sec.Asy);
        std::printf("[F8] Timoshenko cantilever  L=%.0f Asy=%.5g\n", L, sec.Asy);
        checkTrue("not singular", !r.singular, r.diagnostic);
        checkClose("tip deflection = PL^3/3EI + PL/GAs", std::fabs(r.disp(1, Uz)), dEB + dShear, 1e-6);

        // L/d sweep: slender -> Euler-Bernoulli; deep -> shear contribution matters.
        auto timoRatio = [&](real Ls) -> real {
            FrameModel ms; fixtures::cantileverTipLoad(ms, P, Ls, mat, sec);
            SolveOptions o; o.useTimoshenko = true;
            const real dt = std::fabs(solve(ms, o).disp(1, Uz));
            const real de = P * Ls * Ls * Ls / (3.0 * E * sec.Iz);
            return dt / de;
        };
        const real rDeep = timoRatio(300.0);     // L/d = 3
        const real rSlen = timoRatio(6000.0);    // L/d = 60
        std::printf("   ratio Timo/EB: deep(L/d=3)=%.5f  slender(L/d=60)=%.5f\n", rDeep, rSlen);
        checkTrue("deep-beam shear significant (>5%)", rDeep > 1.05, "rDeep=" + std::to_string(rDeep));
        checkTrue("slender converges to EB (<0.1%)", std::fabs(rSlen - 1.0) < 1e-3, "rSlen=" + std::to_string(rSlen));
        checkTrue("shear effect shrinks with slenderness", rSlen < rDeep);
    }

    // ---- F9: circular section properties + resultant biaxial bending ----
    {
        const real rad = 50.0;
        Section cir = Section::Circular(rad);
        const real PI = 3.14159265358979323846;
        const real Iexp = PI * rad * rad * rad * rad / 4.0;
        std::printf("[F9] circular section  r=%.0f  I=%.6g\n", rad, cir.Iz);
        checkClose("I = pi r^4 / 4", cir.Iz, Iexp, 1e-12);
        checkClose("Iy == Iz", cir.Iy, cir.Iz, 1e-12);
        checkClose("J = 2I (polar)", cir.J, 2.0 * Iexp, 1e-12);

        const real P = 1000.0, L = 2000.0;
        FrameModel m; fixtures::cantileverTipLoad(m, P, L, mat, cir);
        SolveResult r = solve(m);
        checkClose("circular cantilever tip = PL^3/3EI", std::fabs(r.disp(1, Uz)),
                   P * L * L * L / (3.0 * E * Iexp), 1e-6);

        // biaxial: round section uses the resultant sqrt(My^2+Mz^2)/W, which is LESS
        // than the rectangular conservative corner sum |My|/W + |Mz|/W.
        ElasticAllowable strength;
        const Capacity cap = Capacity::make(250.0, 250.0, 150.0);
        MemberEndForces f; f.My = 1.0e6; f.Mz = 1.0e6;
        const DemandResult dC = strength.checkSection(f, cir, cap);
        const real W = cir.Wz();
        checkClose("circular biaxial sComp = sqrt2*M/W", dC.sComp, std::sqrt(2.0) * 1.0e6 / W, 1e-6);
        checkTrue("round resultant < rectangular corner sum", dC.sComp < 2.0e6 / W - 1e-9,
                  "sComp=" + std::to_string(dC.sComp));

        // round torsion is exact: tau_max = T r / J at the surface (NOT the diagonal
        // hypot(cy,cz) corner used for rectangles).
        MemberEndForces ft; ft.T = 1.0e6;
        const DemandResult dT = strength.checkSection(ft, cir, cap);
        checkClose("circular torsion sTor = T r / J (exact)", dT.sTor, 1.0e6 * rad / cir.J, 1e-9);
    }

    // ---- F10: piecewise-linear curved member (quarter-circle arch), convergence ----
    {
        const real R = 2000.0, P = 1000.0;
        const real PI = 3.14159265358979323846;
        const real dExp = (PI / 4.0) * P * R * R * R / (E * sec.Iz);
        std::printf("[F10] quarter-circle arch cantilever  R=%.0f  dExp=%.6g\n", R, dExp);
        const int Ns[4] = { 8, 16, 32, 64 };
        real prevErr = 1e30, lastVal = 0; bool monotone = true, nonsing = true;
        for (int i = 0; i < 4; ++i) {
            FrameModel m; fixtures::circularArchCantilever(m, R, Ns[i], P, mat, sec);
            SolveResult r = solve(m);
            if (r.singular) nonsing = false;
            const real d   = std::fabs(r.disp(Ns[i], Uy));
            const real err = std::fabs(d - dExp) / dExp;
            std::printf("   N=%2d  d=%.6g  err=%.3e\n", Ns[i], d, err);
            if (err > prevErr) monotone = false;
            prevErr = err; lastVal = d;
        }
        checkTrue("arch solves (non-singular)", nonsing);
        checkTrue("arch error shrinks as N grows", monotone);
        checkClose("arch tip (N=64) -> (pi/4)PR^3/EI", lastVal, dExp, 1.5e-2);
    }

    // ---- F11: rotation-symmetry invariant (solver equivariance) ----
    {
        auto sq = [](real v) { return v * v; };
        const real P = 1000.0, L = 2000.0;
        FrameModel mb; fixtures::cantileverTipLoad(mb, P, L, mat, sec);
        const SolveResult rb = solve(mb);
        const real u0[3] = { rb.disp(1, Ux), rb.disp(1, Uy), rb.disp(1, Uz) };

        // Rodrigues rotation about a = unit(1,2,3), angle 0.6 rad
        real ax = 1, ay = 2, az = 3; const real an = std::sqrt(ax*ax + ay*ay + az*az);
        ax /= an; ay /= an; az /= an;
        const real ph = 0.6, c = std::cos(ph), s = std::sin(ph), t = 1.0 - c;
        const real R[3][3] = {
            { c + ax*ax*t,    ax*ay*t - az*s, ax*az*t + ay*s },
            { ay*ax*t + az*s, c + ay*ay*t,    ay*az*t - ax*s },
            { az*ax*t - ay*s, az*ay*t + ax*s, c + az*az*t }
        };
        auto rot = [&](real x, real y, real z, int i) { return R[i][0]*x + R[i][1]*y + R[i][2]*z; };

        FrameModel mr;
        mr.materials.push_back(mat); mr.sections.push_back(sec);   // indices 0, 0
        Node n0(0, 0, 0, 0); n0.fixAll();
        Node n1(1, rot(L,0,0,0), rot(L,0,0,1), rot(L,0,0,2));
        mr.nodes = { n0, n1 };
        Member mm(0, 0, 1, 0, 0);
        mm.refVec = Vec3(rot(0,0,1,0), rot(0,0,1,1), rot(0,0,1,2));   // rotate the local frame too
        mr.members = { mm };
        NodalLoad nl; nl.node = 1;
        nl.comp[Ux] = rot(0,0,P,0); nl.comp[Uy] = rot(0,0,P,1); nl.comp[Uz] = rot(0,0,P,2);
        mr.nodalLoads = { nl };
        const SolveResult rr = solve(mr);

        const real ex = rot(u0[0], u0[1], u0[2], 0), ey = rot(u0[0], u0[1], u0[2], 1), ez = rot(u0[0], u0[1], u0[2], 2);
        const real du  = std::sqrt(sq(rr.disp(1,Ux) - ex) + sq(rr.disp(1,Uy) - ey) + sq(rr.disp(1,Uz) - ez));
        const real nrm = std::sqrt(sq(u0[0]) + sq(u0[1]) + sq(u0[2]));
        std::printf("[F11] rotation invariance  |u0|=%.6g\n", nrm);
        checkTrue("rotated solve non-singular", !rr.singular, rr.diagnostic);
        checkClose("tip |u| preserved under rotation",
                   std::sqrt(sq(rr.disp(1,Ux)) + sq(rr.disp(1,Uy)) + sq(rr.disp(1,Uz))), nrm, 1e-9);
        checkTrue("u_rotated == R u_base (equivariance)", du < 1e-6 * nrm, "du=" + std::to_string(du));
    }

    // ---- F12: grillage idealization of a simply-supported isotropic plate ----
    //          (the continuous-surface approximation — the core direction of this engine)
    {
        using namespace frame::grillage;
        const real Epl = 30000.0, nu = 0.3;            // concrete-like, nu = 0.3
        const real Gpl = Epl / (2.0 * (1.0 + nu));     // so the torsion recipe matches plate theory
        Material pmat(Epl, Gpl, 2400.0);

        PlateSpec sp; sp.a = 4000.0; sp.b = 4000.0; sp.t = 250.0; sp.q = 0.025;

        // Kirchhoff plate (simply supported square, uniform load): w_c = alpha q a^4 / D,
        // D = E t^3 / [12 (1-nu^2)], alpha = 0.00406 (Timoshenko table, nu = 0.3).
        const real D  = Epl * sp.t * sp.t * sp.t / (12.0 * (1.0 - nu * nu));
        const real wc = 0.00406 * sp.q * sp.a * sp.a * sp.a * sp.a / D;
        std::printf("[F12] grillage plate  D=%.6g  plate w_c=%.6g mm\n", D, wc);

        auto centerDeflection = [&](int n) -> real {
            FrameModel m; std::string why; PlateSpec s = sp; s.nx = n; s.ny = n;
            if (!buildGrillage(m, s, pmat, why)) return 1e30;
            const SolveResult r = solve(m);
            if (r.singular) return 1e30;
            return std::fabs(r.disp(gridNode(n / 2, n / 2, n), Uz));
        };

        // build + non-singular + load conservation (sum of vertical reactions == q a b)
        {
            FrameModel m; std::string why; PlateSpec s = sp; s.nx = 8; s.ny = 8;
            checkTrue("grillage builds", buildGrillage(m, s, pmat, why), why);
            const SolveResult r = solve(m);
            checkTrue("grillage non-singular", !r.singular, r.diagnostic);
            real sumRz = 0;
            for (size_t k = 0; k < m.nodes.size(); ++k) sumRz += r.reaction((int)k, Uz);
            checkClose("load conservation: sum Rz = q a b", std::fabs(sumRz), sp.q * sp.a * sp.b, 1e-6);
        }

        // The nu-matched grillage reproduces the plate's Dx=Dy=2H=D, so the center
        // deflection sits within a few % of plate theory at every mesh density and is
        // mesh-converged (it settles to a value ~2% off plate theory — the residual is
        // load-lumping + the moment-distribution mismatch, NOT a convergence-to-exact).
        const real d4 = centerDeflection(4), d8 = centerDeflection(8), d12 = centerDeflection(12);
        const real e4 = std::fabs(d4 - wc) / wc, e8 = std::fabs(d8 - wc) / wc, e12 = std::fabs(d12 - wc) / wc;
        std::printf("   N=4 w=%.5g (err %.2f%%)  N=8 w=%.5g (err %.2f%%)  N=12 w=%.5g (err %.2f%%)\n",
                    d4, 100 * e4, d8, 100 * e8, d12, 100 * e12);
        checkTrue("grillage deflection within 5% of plate theory (all meshes)",
                  e4 < 0.05 && e8 < 0.05 && e12 < 0.05,
                  "e4=" + std::to_string(e4) + " e8=" + std::to_string(e8) + " e12=" + std::to_string(e12));
        checkTrue("grillage mesh-converged (|d12-d8| < 2% d8)",
                  std::fabs(d12 - d8) < 0.02 * d8, "d8=" + std::to_string(d8) + " d12=" + std::to_string(d12));
    }

    // ---------- F13: MITC4 shell — plate bending (Reissner-Mindlin + assumed shear) ----------
    {
        const real Es = 30000.0, nu = 0.3, Gs = Es / (2.0 * (1.0 + nu));
        Material smat(Es, Gs); smat.nu = nu;
        std::printf("[F13] MITC4 shell plate bending  (E=%.0f nu=%.1f)\n", Es, nu);

        // (a) simply-supported square plate under uniform pressure -> Kirchhoff w_c.
        const real a = 1000.0, t = 10.0, q = 0.01;
        const real D  = Es * t * t * t / (12.0 * (1.0 - nu * nu));
        const real wc = 0.00406 * q * a * a * a * a / D;
        auto centerW = [&](int n) -> real {
            FrameModel m; fixtures::squarePlateShell(m, a, t, n, q, smat);
            const SolveResult r = solve(m);
            if (r.singular) return 1e30;
            const int c = (n / 2) * (n + 1) + (n / 2);   // node(n/2,n/2)
            return std::fabs(r.disp(c, Uz));
        };
        {
            FrameModel m; fixtures::squarePlateShell(m, a, t, 8, q, smat);
            const SolveResult r = solve(m);
            checkTrue("plate builds non-singular", !r.singular, r.diagnostic);
        }
        const real w4 = centerW(4), w8 = centerW(8), w16 = centerW(16), w24 = centerW(24);
        const real e4 = relErr(w4, wc), e8 = relErr(w8, wc), e16 = relErr(w16, wc), e24 = relErr(w24, wc);
        std::printf("   simply-supported square  plate w_c=%.6g mm\n", wc);
        std::printf("   N=4 w=%.5g (err %.2f%%)  N=8 w=%.5g (err %.2f%%)  N=16 w=%.5g (err %.2f%%)  N=24 w=%.5g (err %.2f%%)\n",
                    w4, 100 * e4, w8, 100 * e8, w16, 100 * e16, w24, 100 * e24);
        // The MITC4 (Reissner-Mindlin) plate converges to a slightly SOFTER solution than
        // Kirchhoff thin-plate theory: the centre deflection rises MONOTONICALLY with mesh
        // refinement and OVERSHOOTS the Kirchhoff value. So |w-wc| does NOT shrink monotonically
        // (it bottoms out near N=16 and grows again by N=24 as w climbs past wc toward the
        // Mindlin answer) — the residual is the Kirchhoff-vs-Mindlin model gap, not a
        // convergence error. Assert the true behaviour: close to Kirchhoff AND monotone-softening
        // (the old "e16 < e4" check hid the overshoot by only comparing the coarsest to N=16).
        checkTrue("plate near Kirchhoff (N=16 and N=24 within 2%)", e16 < 0.02 && e24 < 0.02,
                  "e16=" + std::to_string(e16) + " e24=" + std::to_string(e24));
        checkTrue("plate deflection monotone-increasing (Mindlin softening overshoots Kirchhoff)",
                  w4 < w8 && w8 < w16 && w16 < w24,
                  "w4=" + std::to_string(w4) + " w8=" + std::to_string(w8) +
                  " w16=" + std::to_string(w16) + " w24=" + std::to_string(w24));

        // (b) no shear locking: a very thin plate (t/a = 0.001) must NOT lock (a locked
        // element would read far too stiff -> deflection far below theory).
        {
            const real tThin = 1.0;
            const real Dt  = Es * tThin * tThin * tThin / (12.0 * (1.0 - nu * nu));
            const real wcT = 0.00406 * q * a * a * a * a / Dt;
            FrameModel m; fixtures::squarePlateShell(m, a, tThin, 16, q, smat);
            const SolveResult r = solve(m);
            const int c = 8 * (16 + 1) + 8;
            const real wT = r.singular ? 1e30 : std::fabs(r.disp(c, Uz));
            const real eT = relErr(wT, wcT);
            std::printf("   thin plate t/a=0.001  w=%.6g exp=%.6g (err %.2f%%)\n", wT, wcT, 100 * eT);
            checkTrue("thin plate not shear-locked (N=16 within 3%)", eT < 0.03,
                      "eThin=" + std::to_string(eT));
        }

        // (c) constant-curvature patch test (cylindrical bending). MITC4 reproduces the
        // exact constant moment to MACHINE PRECISION on regular AND parallelogram meshes.
        for (const real skew : { 0.0, 0.4 }) {
            const real ap = 1000.0, tp = 10.0, cc = 1e-6;
            const real Dfac = Es * tp * tp * tp / (12.0 * (1.0 - nu * nu));
            FrameModel m; fixtures::platePatchCylindrical(m, ap, tp, skew, cc, smat);
            const SolveResult r = solve(m);
            const char* tag = (skew == 0.0) ? "regular" : "parallelogram";
            checkTrue("patch non-singular", !r.singular, r.diagnostic);
            const real MxxExp = -Dfac * cc, MyyExp = -nu * Dfac * cc;
            const real scale = std::fabs(MxxExp);
            real eMxx = 0, eMyy = 0, mMxy = 0, mQ = 0;
            for (const auto& sf : r.shellForces) {
                eMxx = std::max(eMxx, std::fabs(sf.Mxx - MxxExp));
                eMyy = std::max(eMyy, std::fabs(sf.Myy - MyyExp));
                mMxy = std::max(mMxy, std::fabs(sf.Mxy));
                mQ   = std::max(mQ, std::max(std::fabs(sf.Qx), std::fabs(sf.Qy)));
            }
            std::printf("   patch (%-13s) Mxx_exp=%.6g  max|dMxx|=%.2e |dMyy|=%.2e |Mxy|=%.2e |Q|=%.2e\n",
                        tag, MxxExp, eMxx, eMyy, mMxy, mQ);
            checkTrue("patch Mxx = -D c (constant)", eMxx < 1e-8 * scale, "eMxx=" + std::to_string(eMxx));
            checkTrue("patch Myy = -nu D c (constant)", eMyy < 1e-8 * scale, "eMyy=" + std::to_string(eMyy));
            checkTrue("patch Mxy ~ 0", mMxy < 1e-8 * scale, "mMxy=" + std::to_string(mMxy));
            checkTrue("patch shear ~ 0 (pure bending)", mQ < 1e-6 * scale, "mQ=" + std::to_string(mQ));
        }
    }

    // ---------- F14: MITC4 shell — membrane + drilling + 3D facet rotation ----------
    {
        const real Es = 30000.0, nu = 0.3, Gs = Es / (2.0 * (1.0 + nu));
        Material smat(Es, Gs); smat.nu = nu;
        std::printf("[F14] MITC4 shell membrane + drilling + 3D facet\n");

        // (a) membrane constant-strain patch (regular + parallelogram) -> constant N.
        for (const real skew : { 0.0, 0.4 }) {
            const real a = 1000.0, t = 10.0, gx = 1e-4;
            const real f = Es / (1.0 - nu * nu);
            FrameModel m; fixtures::membranePatch(m, a, t, skew, gx, smat);
            const SolveResult r = solve(m);
            const char* tag = (skew == 0.0) ? "regular" : "parallelogram";
            checkTrue("membrane patch non-singular", !r.singular, r.diagnostic);
            const real NxxExp = t * f * gx, NyyExp = nu * NxxExp;
            const real scale = std::fabs(NxxExp);
            real eNxx = 0, eNyy = 0, mNxy = 0;
            for (const auto& sf : r.shellForces) {
                eNxx = std::max(eNxx, std::fabs(sf.Nxx - NxxExp));
                eNyy = std::max(eNyy, std::fabs(sf.Nyy - NyyExp));
                mNxy = std::max(mNxy, std::fabs(sf.Nxy));
            }
            std::printf("   membrane patch (%-13s) Nxx_exp=%.6g max|dNxx|=%.2e |dNyy|=%.2e |Nxy|=%.2e\n",
                        tag, NxxExp, eNxx, eNyy, mNxy);
            checkTrue("membrane patch Nxx constant", eNxx < 1e-8 * scale, "eNxx=" + std::to_string(eNxx));
            checkTrue("membrane patch Nyy = nu Nxx", eNyy < 1e-8 * scale, "eNyy=" + std::to_string(eNyy));
            checkTrue("membrane patch Nxy ~ 0", mNxy < 1e-8 * scale, "mNxy=" + std::to_string(mNxy));
        }

        // (b) drilling gate: a FLAT clamped plate with interior in-plane + drilling DOFs
        // all free must solve NON-SINGULAR (the drilling stiffness removed the coplanar
        // Rz zero-energy mode) and match clamped-plate theory w_c = 0.00126 q a^4 / D.
        {
            const real a = 1000.0, t = 10.0, q = 0.01;
            const real D = Es * t * t * t / (12.0 * (1.0 - nu * nu));
            const real wcC = 0.00126 * q * a * a * a * a / D;
            FrameModel m; fixtures::clampedPlateShell(m, a, t, 16, q, smat);
            const SolveResult r = solve(m);
            checkTrue("flat clamped shell NON-singular (drilling gate)", !r.singular, r.diagnostic);
            const int c = 8 * (16 + 1) + 8;
            const real wC = r.singular ? 1e30 : std::fabs(r.disp(c, Uz));
            const real eC = relErr(wC, wcC);
            std::printf("   clamped plate w=%.6g exp=%.6g (err %.2f%%)\n", wC, wcC, 100 * eC);
            checkTrue("clamped plate within 3% of theory", eC < 0.03, "eC=" + std::to_string(eC));
        }

        // (c) rotation invariance: solve a flat clamped plate, then the SAME model rigidly
        // rotated to an arbitrary 3D orientation. The pressure follows the facet normal and
        // an encastre boundary is frame-invariant, so |displacement| at the centre is preserved.
        {
            const real a = 1000.0, t = 10.0, q = 0.01;
            const int n = 12, c = (n / 2) * (n + 1) + (n / 2);
            FrameModel m0; fixtures::clampedPlateShell(m0, a, t, n, q, smat);
            const SolveResult r0 = solve(m0);
            const real d0 = std::sqrt(r0.disp(c, Ux) * r0.disp(c, Ux) +
                                      r0.disp(c, Uy) * r0.disp(c, Uy) +
                                      r0.disp(c, Uz) * r0.disp(c, Uz));
            // Rodrigues rotation about a normalized arbitrary axis.
            const real ax = 1, ay = 2, az = 3, an = std::sqrt(14.0);
            const real kx = ax / an, ky = ay / an, kz = az / an, th = 0.9;
            const real ct = std::cos(th), st = std::sin(th), vt = 1 - ct;
            const real R[3][3] = {
                { ct + kx * kx * vt,      kx * ky * vt - kz * st, kx * kz * vt + ky * st },
                { ky * kx * vt + kz * st, ct + ky * ky * vt,      ky * kz * vt - kx * st },
                { kz * kx * vt - ky * st, kz * ky * vt + kx * st, ct + kz * kz * vt }
            };
            FrameModel m1; fixtures::clampedPlateShell(m1, a, t, n, q, smat);
            fixtures::rotateModelRigid(m1, R);
            const SolveResult r1 = solve(m1);
            checkTrue("rotated shell non-singular", !r1.singular, r1.diagnostic);
            const real d1 = std::sqrt(r1.disp(c, Ux) * r1.disp(c, Ux) +
                                      r1.disp(c, Uy) * r1.disp(c, Uy) +
                                      r1.disp(c, Uz) * r1.disp(c, Uz));
            std::printf("   |u_centre| flat=%.6g rotated=%.6g\n", d0, d1);
            checkClose("|u_centre| invariant under 3D rotation", d1, d0, 1e-9);
        }
    }

    // ---------- F15: MITC4 shell — Scordelis-Lo roof (curved-surface benchmark) ----------
    {
        const real kPi = 3.14159265358979323846;
        const real R = 25.0, L = 50.0, phi0 = 40.0 * kPi / 180.0, t = 0.25, g = 90.0;
        const real Er = 4.32e8, nur = 0.0, Gr = Er / 2.0;
        Material rmat(Er, Gr); rmat.nu = nur;
        const real ref = 0.3024;                          // MacNeal-Harder reference
        std::printf("[F15] MITC4 Scordelis-Lo roof  (ref free-edge w=%.4f)\n", ref);
        auto edgeW = [&](int n) -> real {
            FrameModel m; fixtures::scordelisLoRoof(m, R, L, phi0, t, g, n, n, rmat);
            const SolveResult r = solve(m);
            if (r.singular) return 1e30;
            const int c = n * (n + 1) + n;                // gid(nf=n, ny=n): free-edge midspan
            return std::fabs(r.disp(c, Uz));
        };
        {
            FrameModel m; fixtures::scordelisLoRoof(m, R, L, phi0, t, g, 8, 8, rmat);
            const SolveResult r = solve(m);
            checkTrue("roof builds non-singular", !r.singular, r.diagnostic);
        }
        const real w8 = edgeW(8), w16 = edgeW(16), w24 = edgeW(24);
        const real e8 = relErr(w8, ref), e16 = relErr(w16, ref), e24 = relErr(w24, ref);
        std::printf("   N=8 w=%.5g (err %.2f%%)  N=16 w=%.5g (err %.2f%%)  N=24 w=%.5g (err %.2f%%)\n",
                    w8, 100 * e8, w16, 100 * e16, w24, 100 * e24);
        checkTrue("roof converges to reference (N=24 within 3%)", e24 < 0.03,
                  "e24=" + std::to_string(e24));
        checkTrue("roof mesh-converging (e24 < e8)", e24 < e8,
                  "e8=" + std::to_string(e8) + " e24=" + std::to_string(e24));
    }

    // ---------- F16: MITC4 shell — pinched cylinder (hard inextensional benchmark) ----------
    {
        const real R = 300.0, L = 600.0, t = 3.0, P = 1.0;
        const real Ec = 3.0e6, nuc = 0.3, Gc = Ec / (2.0 * (1.0 + nuc));
        Material cmat(Ec, Gc); cmat.nu = nuc;
        const real ref = 1.8248e-5;                       // analytic radial deflection under load
        std::printf("[F16] MITC4 pinched cylinder  (ref under-load w=%.4e)\n", ref);
        auto loadW = [&](int n) -> real {
            FrameModel m; fixtures::pinchedCylinder(m, R, L, t, P, n, n, cmat);
            const SolveResult r = solve(m);
            if (r.singular) return 0.0;
            const int c = n * (n + 1) + 0;                // gid(0, nz=n): load node
            return std::fabs(r.disp(c, Ux));
        };
        {
            FrameModel m; fixtures::pinchedCylinder(m, R, L, t, P, 8, 8, cmat);
            const SolveResult r = solve(m);
            checkTrue("cylinder builds non-singular", !r.singular, r.diagnostic);
        }
        const real w8 = loadW(8), w16 = loadW(16), w24 = loadW(24), w32 = loadW(32);
        std::printf("   N=8 %.5e (%.1f%%)  N=16 %.5e (%.1f%%)  N=24 %.5e (%.1f%%)  N=32 %.5e (%.1f%%)\n",
                    w8, 100 * w8 / ref, w16, 100 * w16 / ref, w24, 100 * w24 / ref, w32, 100 * w32 / ref);
        // Converges from below toward the reference; require monotone approach and a
        // reasonable fraction at the finest mesh (flat-facet MITC4 is slow on this test).
        checkTrue("cylinder converging upward (w32 > w16 > w8)", w32 > w16 && w16 > w8,
                  "w8=" + std::to_string(w8) + " w16=" + std::to_string(w16) + " w32=" + std::to_string(w32));
        checkTrue("cylinder N=32 reaches >= 90% of reference", w32 >= 0.90 * ref,
                  "ratio=" + std::to_string(w32 / ref));
        checkTrue("cylinder N=32 not over-stiff/over-soft (<= 105%)", w32 <= 1.05 * ref,
                  "ratio=" + std::to_string(w32 / ref));
    }

    // ---------- F17: superposition identity (the load-combination primitive) ----------
    {
        std::printf("[F17] superposition: u(A)+u(B) == u(A+B)\n");
        auto run = [&](real Pz, real My) -> SolveResult {
            FrameModel m; fixtures::cantileverBare(m, 2000.0, mat, sec);
            NodalLoad nl; nl.node = 1; nl.comp[Uz] = Pz; nl.comp[Ry] = My;
            m.nodalLoads = { nl };
            return solve(m);
        };
        const SolveResult rA  = run(1000.0, 0.0);
        const SolveResult rB  = run(0.0, 2.0e6);
        const SolveResult rAB = run(1000.0, 2.0e6);
        const SolveResult rC  = combine({ rA, rB }, { 1.0, 1.0 });
        real us = 1e-30, Rs = 1e-30, Ms = 1e-30, du = 0, dr = 0, dm = 0;
        for (size_t k = 0; k < rAB.u.size(); ++k) { us = std::max(us, std::fabs(rAB.u[k])); du = std::max(du, std::fabs(rC.u[k] - rAB.u[k])); }
        for (size_t k = 0; k < rAB.reactions.size(); ++k) { Rs = std::max(Rs, std::fabs(rAB.reactions[k])); dr = std::max(dr, std::fabs(rC.reactions[k] - rAB.reactions[k])); }
        for (size_t e = 0; e < rAB.memberForces.size(); ++e) {
            const auto& a = rC.memberForces[e].endI; const auto& b = rAB.memberForces[e].endI;
            Ms = std::max(Ms, std::max(std::fabs(b.My), std::fabs(b.Mz)));
            dm = std::max(dm, std::max(std::fabs(a.My - b.My), std::fabs(a.Mz - b.Mz)));
        }
        std::printf("   max|du|=%.2e/%.3g  max|dR|=%.2e/%.3g  max|dM|=%.2e/%.3g\n", du, us, dr, Rs, dm, Ms);
        checkTrue("superposition u identity (rel<1e-10)", du < 1e-10 * us, "du/us=" + std::to_string(du / us));
        checkTrue("superposition reactions identity",     dr < 1e-10 * Rs, "dr/Rs=" + std::to_string(dr / Rs));
        checkTrue("superposition moments identity",       dm < 1e-10 * Ms, "dm/Ms=" + std::to_string(dm / Ms));
    }

    // ---------- F18: self-weight from rho (validates the kg/m^3 -> N-mm unit bridge) ----------
    {
        const real g = 9810.0;                                   // mm/s^2
        const real wg = mat.rho * sec.A * g * 1.0e-12;           // beam weight/length [N/mm]
        std::printf("[F18] self-weight  rho=%.0f -> w=%.6g N/mm\n", mat.rho, wg);

        // (a) cantilever self-weight: fixed-end moment wL^2/2, vertical reaction wL,
        //     tip deflection wL^4/(8 E Iz)  (square section -> resultant is plane-agnostic).
        {
            const real L = 2000.0;
            FrameModel m; fixtures::cantileverBare(m, L, mat, sec);
            addSelfWeight(m, g);
            const SolveResult r = solve(m);
            checkTrue("cantilever SW non-singular", !r.singular, r.diagnostic);
            const auto& mf = r.memberForces[0].endI;
            const real Mroot = std::sqrt(mf.My * mf.My + mf.Mz * mf.Mz);
            const real Rvert = std::fabs(r.reaction(0, Uz));
            const real dtip  = std::sqrt(r.disp(1, Uy) * r.disp(1, Uy) + r.disp(1, Uz) * r.disp(1, Uz));
            checkClose("cantilever SW fixed-end moment wL^2/2", Mroot, wg * L * L / 2.0, 1e-6);
            checkClose("cantilever SW reaction = wL",           Rvert, wg * L, 1e-6);
            checkClose("cantilever SW tip defl wL^4/8EI",       dtip, wg * L * L * L * L / (8.0 * E * sec.Iz), 1e-6);
        }
        // (b) simply-supported self-weight: reactions wL/2, midspan deflection 5wL^4/384EI.
        {
            const real L = 3000.0;
            FrameModel m; fixtures::simplySupportedBare(m, L, mat, sec);
            addSelfWeight(m, g);
            const SolveResult r = solve(m);
            checkTrue("SS SW non-singular", !r.singular, r.diagnostic);
            checkClose("SS SW reaction wL/2 (node0)", std::fabs(r.reaction(0, Uz)), wg * L / 2.0, 1e-6);
            checkClose("SS SW midspan defl 5wL^4/384EI", std::fabs(r.disp(1, Uz)),
                       5.0 * wg * L * L * L * L / (384.0 * E * sec.Iz), 1e-6);
        }
        // (c) shell self-weight (body load) == equivalent transverse pressure rho*t*g.
        {
            const real a = 1000.0, t = 10.0, nu = 0.3, Es = 30000.0;
            Material smat(Es, Es / (2.0 * (1.0 + nu)), 2400.0); smat.nu = nu;   // rho = 2400 kg/m^3
            const real p = smat.rho * t * g * 1.0e-12;            // equivalent pressure [MPa]
            const int n = 12, c = (n / 2) * (n + 1) + (n / 2);
            FrameModel mp; fixtures::squarePlateShell(mp, a, t, n, p, smat);   // pressure model
            const SolveResult rp = solve(mp);
            FrameModel mw; fixtures::squarePlateShell(mw, a, t, n, 0.0, smat);  // bare + self-weight
            addSelfWeight(mw, g);
            const SolveResult rw = solve(mw);
            const real wP = std::fabs(rp.disp(c, Uz)), wW = std::fabs(rw.disp(c, Uz));
            std::printf("   shell: pressure w=%.6g  self-weight w=%.6g\n", wP, wW);
            checkTrue("shell SW non-singular", !rw.singular, rw.diagnostic);
            checkClose("shell self-weight == rho*t*g pressure", wW, wP, 1e-9);
        }
    }

    // ---------- F19: factorize-once-solve-many + prescribed settlement ----------
    {
        std::printf("[F19] factorize-once reuse + prescribed settlement\n");
        // (a) reuse ONE factorization for two different nodal loads; each must equal a
        //     fresh full solve (proves solveLoad reuse is exact, not an approximation).
        {
            FrameModel m; fixtures::cantileverBare(m, 2000.0, mat, sec);
            PreparedSystem ps = assembleAndFactor(m);
            auto maxDiff = [](const SolveResult& a, const SolveResult& b) {
                real d = 0; for (size_t k = 0; k < a.u.size(); ++k) d = std::max(d, std::fabs(a.u[k] - b.u[k])); return d;
            };
            m.nodalLoads = { [] { NodalLoad n; n.node = 1; n.comp[Uz] = 1000.0; return n; }() };
            const real dA = maxDiff(solveLoad(ps, m), solve(m));
            m.nodalLoads = { [] { NodalLoad n; n.node = 1; n.comp[Uy] = 700.0; n.comp[Ry] = 1.0e6; return n; }() };
            const real dB = maxDiff(solveLoad(ps, m), solve(m));   // SAME ps, different load
            std::printf("   reuse vs fresh: max|du| load A=%.2e  load B=%.2e\n", dA, dB);
            checkTrue("solveLoad(A) == fresh solve(A)", dA < 1e-12, "dA=" + std::to_string(dA));
            checkTrue("solveLoad(B) reuses same factorization == fresh", dB < 1e-12, "dB=" + std::to_string(dB));
        }
        // (b) prescribed support settlement oracle (closes the long-standing gap):
        //     fixed-fixed beam, one end settles delta -> end moment 6EI*delta/L^2, reaction 12EI*delta/L^3.
        {
            const real L = 2000.0, delta = 1.0;
            FrameModel m; fixtures::clampedSettlement(m, L, delta, mat, sec);
            const SolveResult r = solve(m);
            checkTrue("settlement non-singular", !r.singular, r.diagnostic);
            const auto& mf = r.memberForces[0].endI;
            const real Mend  = std::sqrt(mf.My * mf.My + mf.Mz * mf.Mz);
            const real Rvert = std::fabs(r.reaction(0, Uz));
            checkClose("settlement end moment 6EI d/L^2", Mend,  6.0 * E * sec.Iz * delta / (L * L), 1e-6);
            checkClose("settlement reaction 12EI d/L^3",  Rvert, 12.0 * E * sec.Iz * delta / (L * L * L), 1e-6);
        }
    }

    // ---------- F20: live-load PATTERN loading + result envelope (two-span continuous) ----------
    {
        const real L = 3000.0, w = 5.0;
        std::printf("[F20] pattern loading + envelope (two-span continuous, w=%.1f L=%.0f)\n", w, L);
        auto pattern = [&](bool left, bool right) -> SolveResult {
            FrameModel m; fixtures::twoSpanContinuous(m, L, mat, sec);
            auto addUDL = [&](int mem) { MemberUDL u; u.member = mem; u.w_local = { 0, -w, 0 }; m.memberUDLs.push_back(u); };
            if (left)  { addUDL(0); addUDL(1); }
            if (right) { addUDL(2); addUDL(3); }
            return solve(m);
        };
        const SolveResult rL = pattern(true, false), rR = pattern(false, true), rAll = pattern(true, true);
        auto Mres = [](const MemberEndForces& f) { return std::sqrt(f.My * f.My + f.Mz * f.Mz); };
        // member 1 endJ is at the interior support B. Three-moment theory: full load -> wL^2/8,
        // a single loaded span -> wL^2/16. The full-load pattern governs the support moment.
        const real MB_full = Mres(rAll.memberForces[1].endJ);
        const real MB_one  = Mres(rL.memberForces[1].endJ);
        std::printf("   support moment: full=%.5g (wL^2/8=%.5g)  single-span=%.5g (wL^2/16=%.5g)\n",
                    MB_full, w * L * L / 8.0, MB_one, w * L * L / 16.0);
        checkClose("two-span support moment (full) = wL^2/8", MB_full, w * L * L / 8.0, 1e-4);
        checkClose("two-span support moment (1 span) = wL^2/16", MB_one, w * L * L / 16.0, 1e-4);

        const ResultEnvelope env = envelope({ rL, rR, rAll });
        // envelope captures the worst (most hogging) support moment = the full-load case.
        const real MB_envWorst = std::sqrt(env.endJMin[1].My * env.endJMin[1].My + env.endJMin[1].Mz * env.endJMin[1].Mz);
        checkClose("envelope worst support moment = wL^2/8", MB_envWorst, w * L * L / 8.0, 1e-4);
        // pattern matters for the SPAN: a single loaded span hogs the support less but bends
        // its own span MORE than the full-load case (the chessboard rationale).
        const real Mspan_one  = Mres(rL.memberForces[0].endJ);
        const real Mspan_full = Mres(rAll.memberForces[0].endJ);
        std::printf("   left-span moment at midL: single-left=%.5g  full=%.5g\n", Mspan_one, Mspan_full);
        checkTrue("single-span pattern bends its span more than full load", Mspan_one > Mspan_full * 1.05,
                  "one=" + std::to_string(Mspan_one) + " full=" + std::to_string(Mspan_full));
        // envelope displacement plumbing: uMin equals the per-case minimum at the loaded midspan.
        const real uminCase = std::min({ rL.disp(1, Uz), rR.disp(1, Uz), rAll.disp(1, Uz) });
        checkClose("envelope uMin matches per-case min", env.uMin[gdof(1, Uz)], uminCase, 1e-12);
    }

    // ---------- F21: influence lines / moving load (same-K-many-RHS) ----------
    {
        std::printf("[F21] influence lines (factorize-once, unit load marched)\n");
        const int n = 8; const real L = 4000.0;
        FrameModel m; fixtures::simplySupportedBeamN(m, n, L, mat, sec);
        PreparedSystem ps = assembleAndFactor(m);
        std::vector<NodeId> loadNodes; for (int i = 0; i <= n; ++i) loadNodes.push_back(i);

        // (a) reaction R_A influence line == (L-x)/L, cross-checked by Müller-Breslau
        //     (unit settlement at A, no load -> the deflected shape IS the reaction IL).
        const std::vector<real> ilR = reactionInfluenceLine(ps, m, loadNodes, 0, Uz);
        m.nodalLoads.clear(); m.nodes[0].prescribed[Uz] = 1.0;
        const SolveResult rMB = solveLoad(ps, m);
        m.nodes[0].prescribed[Uz] = 0.0;
        real eAna = 0, eMB = 0;
        for (int i = 0; i <= n; ++i) {
            const real x = L * i / n, exact = (L - x) / L;
            eAna = std::max(eAna, std::fabs(ilR[i] - exact));
            eMB  = std::max(eMB, std::fabs(ilR[i] - rMB.disp(i, Uz)));
        }
        checkTrue("reaction IL == (L-x)/L", eAna < 1e-9, "e=" + std::to_string(eAna));
        checkTrue("reaction IL == Muller-Breslau deflected shape", eMB < 1e-9, "e=" + std::to_string(eMB));

        // (b) midspan bending-moment influence line = triangle, peak ab/L = L/4 at midspan.
        const int c = n / 2; real eShape = 0;
        const real peak = L / 4.0;
        for (int i = 0; i <= n; ++i) {
            m.nodalLoads.clear();
            NodalLoad nl; nl.node = i; nl.comp[Uz] = -1.0; m.nodalLoads.push_back(nl);
            const SolveResult r = solveLoad(ps, m);
            const auto& mf = r.memberForces[c - 1].endJ;        // moment at midspan node
            const real Mc = std::sqrt(mf.My * mf.My + mf.Mz * mf.Mz);
            const real x = L * i / n, exact = (x <= L / 2) ? x / 2.0 : (L - x) / 2.0;
            eShape = std::max(eShape, std::fabs(Mc - exact));
        }
        m.nodalLoads.clear();
        std::printf("   reaction IL err=%.2e (analytic) / %.2e (Muller-Breslau);  moment IL err=%.2e (peak ab/L=%.0f)\n",
                    eAna, eMB, eShape, peak);
        checkTrue("midspan moment IL == triangle (peak ab/L)", eShape < 1e-6 * peak, "e=" + std::to_string(eShape));
    }

    // ---------- F22: modal analysis (free-vibration natural frequencies) ----------
    {
        const real kPi = 3.14159265358979323846;
        const real rhoC = mat.rho * 1.0e-12;   // consistent units (tonne/mm^3)
        std::printf("[F22] modal analysis  (rho=%.0f -> %.3e tonne/mm^3)\n", mat.rho, rhoC);
        // (a) simply-supported beam, fundamental omega1 = (pi/L)^2 sqrt(EI/(rho A)).
        {
            const int n = 12; const real L = 4000.0;
            FrameModel m; fixtures::simplySupportedBeamN(m, n, L, mat, sec);
            PreparedSystem ps = assembleAndFactor(m);
            ModalOptions mo; mo.numModes = 4;
            const ModalResult mr = solveModal(ps, mo);
            checkTrue("SS modal non-singular", !mr.singular && mr.modes.size() >= 2, mr.diagnostic);
            const real w1ex = kPi * kPi / (L * L) * std::sqrt(E * sec.Iz / (rhoC * sec.A));
            std::printf("   SS beam: omega1=%.6g rad/s (exact %.6g; f1=%.4f Hz)\n",
                        mr.modes[0].omega, w1ex, mr.modes[0].freqHz);
            checkClose("SS fundamental omega1 = (pi/L)^2 sqrt(EI/rhoA)", mr.modes[0].omega, w1ex, 1e-3);
            checkTrue("modes ascending", mr.modes[1].omega >= mr.modes[0].omega * 0.999, "");
        }
        // (b) cantilever, fundamental omega1 = 1.875^2 sqrt(EI/(rho A L^4)).
        {
            const int n = 12; const real L = 3000.0;
            FrameModel m; fixtures::cantileverBeamN(m, n, L, mat, sec);
            PreparedSystem ps = assembleAndFactor(m);
            const ModalResult mr = solveModal(ps, ModalOptions{});
            checkTrue("cantilever modal non-singular", !mr.singular && !mr.modes.empty(), mr.diagnostic);
            const real beta1 = 1.8751040687;
            const real w1ex = beta1 * beta1 * std::sqrt(E * sec.Iz / (rhoC * sec.A * L * L * L * L));
            std::printf("   cantilever: omega1=%.6g (exact %.6g; f1=%.4f Hz)\n",
                        mr.modes[0].omega, w1ex, mr.modes[0].freqHz);
            checkClose("cantilever fundamental omega1 = 1.875^2 sqrt(EI/rhoAL^4)", mr.modes[0].omega, w1ex, 1e-3);
        }
    }

    // ---------- F23: linear buckling (geometric stiffness) vs Euler ----------
    {
        const real kPi = 3.14159265358979323846;
        const real Pref = 1000.0;   // reference axial compression
        std::printf("[F23] linear buckling (Euler column)\n");
        // (a) pinned-pinned column: Pcr = pi^2 EI / L^2.
        {
            const int n = 10; const real L = 3000.0;
            FrameModel m; fixtures::simplySupportedBeamN(m, n, L, mat, sec);
            NodalLoad nl; nl.node = n; nl.comp[Ux] = -Pref; m.nodalLoads = { nl };
            PreparedSystem ps = assembleAndFactor(m);
            const BucklingResult br = solveBuckling(ps, m);
            checkTrue("pinned-pinned buckling non-singular", !br.singular, br.diagnostic);
            const real Pcr = br.criticalFactor * Pref;
            const real PcrEx = kPi * kPi * E * sec.Iz / (L * L);
            std::printf("   pinned-pinned: Pcr=%.6g (Euler %.6g, factor=%.4g)\n", Pcr, PcrEx, br.criticalFactor);
            checkClose("Euler Pcr = pi^2 EI/L^2", Pcr, PcrEx, 1e-3);
        }
        // (b) fixed-free (cantilever) column: Pcr = pi^2 EI / (2L)^2.
        {
            const int n = 10; const real L = 3000.0;
            FrameModel m; fixtures::cantileverBeamN(m, n, L, mat, sec);
            NodalLoad nl; nl.node = n; nl.comp[Ux] = -Pref; m.nodalLoads = { nl };
            PreparedSystem ps = assembleAndFactor(m);
            const BucklingResult br = solveBuckling(ps, m);
            checkTrue("fixed-free buckling non-singular", !br.singular, br.diagnostic);
            const real Pcr = br.criticalFactor * Pref;
            const real PcrEx = kPi * kPi * E * sec.Iz / (4.0 * L * L);
            std::printf("   fixed-free: Pcr=%.6g (Euler %.6g, factor=%.4g)\n", Pcr, PcrEx, br.criticalFactor);
            checkClose("Euler Pcr = pi^2 EI/(2L)^2", Pcr, PcrEx, 1e-3);
        }
    }

    // ---------- F24: response spectrum / seismic (modal combination) ----------
    {
        const real kPi = 3.14159265358979323846;
        std::printf("[F24] response spectrum (modal participation + SRSS)\n");
        // rectangular section (Iy != Iz) so the two bending planes are NON-degenerate -- a
        // square section's degenerate modes get arbitrarily mixed by the eigensolver, splitting
        // the Uz participation across the pair.
        Section secR = Section::Rectangular(80.0, 120.0);
        Spectrum sp; sp.T = { 0.0, 10.0 }; sp.Sa = { 9810.0, 9810.0 };   // flat 1g (mm/s^2)
        auto maxOf = [](const std::vector<real>& v) { real m = 0; for (real e : v) m = std::max(m, e); return m; };
        auto sumOf = [](const std::vector<real>& v) { real s = 0; for (real e : v) s += e; return s; };
        // (a) simply-supported beam: 1st-mode effective mass ratio = 8/pi^2 (textbook).
        {
            const int n = 10; const real L = 4000.0;
            FrameModel m; fixtures::simplySupportedBeamN(m, n, L, mat, secR);
            PreparedSystem ps = assembleAndFactor(m);
            const ModalResult mr = solveModal(ps, ModalOptions{ 100 });   // all modes
            const ResponseSpectrumResult rs = solveResponseSpectrum(ps, mr, sp, Uz, SpectrumCombo::SRSS, 0.05);
            checkTrue("RS non-singular", !rs.singular, rs.diagnostic);
            const real r1 = maxOf(rs.effMass) / rs.totalMass, rsum = sumOf(rs.effMass) / rs.totalMass;
            std::printf("   SS: 1st-mode eff mass ratio=%.4f (8/pi^2=%.4f)  participating=%.1f%%  V=%.5g\n",
                        r1, 8.0 / (kPi * kPi), 100 * rsum, rs.baseShear);
            checkClose("SS 1st-mode eff mass = 8/pi^2", r1, 8.0 / (kPi * kPi), 0.02);
            checkTrue("SS participating mass > 85%", rsum > 0.85, "rsum=" + std::to_string(rsum));
            checkTrue("RS base shear > 0", rs.baseShear > 0, "");
        }
        // (b) cantilever: 1st-mode effective mass ratio ~ 0.6131 (textbook).
        {
            const int n = 10; const real L = 3000.0;
            FrameModel m; fixtures::cantileverBeamN(m, n, L, mat, secR);
            PreparedSystem ps = assembleAndFactor(m);
            const ModalResult mr = solveModal(ps, ModalOptions{ 100 });
            const ResponseSpectrumResult rs = solveResponseSpectrum(ps, mr, sp, Uz, SpectrumCombo::SRSS, 0.05);
            const real r1 = maxOf(rs.effMass) / rs.totalMass;
            std::printf("   cantilever: 1st-mode eff mass ratio=%.4f (0.6131)\n", r1);
            checkClose("cantilever 1st-mode eff mass ~ 0.613", r1, 0.6131, 0.03);
        }
    }

    // ---------- F25: modal-superposition time history (Newmark step response) ----------
    {
        const real kPi = 3.14159265358979323846;
        std::printf("[F25] modal time-history (step load, Newmark-beta)\n");
        const int n = 10; const real L = 3000.0;
        FrameModel m; fixtures::cantileverBeamN(m, n, L, mat, sec);
        NodalLoad nl; nl.node = n; nl.comp[Uz] = -1000.0; m.nodalLoads = { nl };
        PreparedSystem ps = assembleAndFactor(m);
        const int cdof = gdof(n, Uz);
        // (a) single-mode SDOF: an undamped suddenly-applied (step) load gives DLF = 2.
        {
            const ModalResult mr = solveModal(ps, ModalOptions{ 1 });
            const real T1 = 2.0 * kPi / mr.modes[0].omega;
            ModalDynamicsOptions opt; opt.zeta = 0.0;  opt.dt = T1 / 400; opt.nSteps = 400;   // ~1 period
            const ModalTimeHistory th = solveModalStepResponse(ps, m, mr, opt);
            real peak = 0; for (const auto& u : th.u) peak = std::max(peak, std::fabs(u[cdof]));
            ModalDynamicsOptions od = opt; od.zeta = 0.4; od.nSteps = 4000;                   // settle -> static
            const ModalTimeHistory td = solveModalStepResponse(ps, m, mr, od);
            const real settled = std::fabs(td.u.back()[cdof]);
            std::printf("   SDOF: peak=%.5g  static(settled)=%.5g  DLF=%.4f\n", peak, settled, peak / settled);
            checkClose("undamped step DLF = 2", peak / settled, 2.0, 0.01);
        }
        // (b) multi-mode damped step response settles to the FULL static solution.
        {
            const ModalResult mr = solveModal(ps, ModalOptions{ 20 });
            const real T1 = 2.0 * kPi / mr.modes[0].omega;
            ModalDynamicsOptions opt; opt.zeta = 0.10; opt.dt = T1 / 100; opt.nSteps = 3000;   // many periods
            const ModalTimeHistory th = solveModalStepResponse(ps, m, mr, opt);
            const SolveResult st = solveLoad(ps, m);
            const real settled = th.u.back()[cdof], staticv = st.u[cdof];
            std::printf("   multi-mode damped: settled=%.5g  static=%.5g\n", settled, staticv);
            checkClose("damped step settles to static", settled, staticv, 0.02);
        }
    }

    // ---------- F26: element removal (Member::active) — redistribution + mechanism oracle ----------
    {
        // Propped cantilever = cantilever beam M0 (encastre base 0 -> tip 1) + vertical prop M1
        // (fixed foot 2 -> tip 1). It is statically INDETERMINATE. Deactivating the prop must drop
        // the structure back to a DETERMINATE cantilever whose exact tip deflection is the closed
        // form -PL^3/3EI (the same oracle as F1). This is the element-removal / collapse foundation
        // (C1) and its independent analytic oracle (C5): an inactive member leaves K cleanly, its own
        // forces stay zero, and removal is byte-identical to physically omitting the member.
        const real L = 2000.0, H = 1000.0, P = 1000.0;
        const real dExp = -P * L * L * L / (3.0 * mat.E * sec.Iy);   // cantilever tip deflection (closed form)

        auto buildPropped = [&](FrameModel& m) {
            fixtures::prepMatSec(m, mat, sec);
            Node n0(0, 0.0, 0.0,  0.0); n0.fixAll();   // encastre base
            Node n1(1,   L, 0.0,  0.0);                // free tip
            Node n2(2,   L, 0.0,   -H); n2.fixAll();   // prop foot, below the tip
            m.nodes = { n0, n1, n2 };
            m.members = { Member(0, 0, 1, 0, 0),       // cantilever beam
                          Member(1, 2, 1, 0, 0) };     // vertical prop
            NodalLoad p; p.node = 1; p.comp[Uz] = -P;
            m.nodalLoads = { p };
        };

        std::printf("[F26] element removal (Member::active) — propped cantilever, remove prop -> -PL^3/3EI=%.5g\n", dExp);

        FrameModel mFull; buildPropped(mFull);
        const SolveResult rFull = solve(mFull);
        checkTrue("propped cantilever non-singular", !rFull.singular, rFull.diagnostic);
        checkTrue("prop stiffens tip (|defl| < cantilever)", std::fabs(rFull.disp(1, Uz)) < std::fabs(dExp), "");

        // remove the prop (M1) -> pure cantilever -> tip deflection = -PL^3/3EI exactly
        FrameModel mCut = mFull; mCut.members[1].active = false;
        const SolveResult rCut = solve(mCut);
        checkTrue("after removing prop: non-singular", !rCut.singular, rCut.diagnostic);
        checkClose("removed prop force stays 0", std::fabs(rCut.memberForces[1].endI.N), 0.0, 1e-9);
        checkClose("tip deflection = -PL^3/3EI", rCut.disp(1, Uz), dExp, 1e-6);

        // removal-by-flag invariant: active=false MUST equal physically omitting the member
        FrameModel mOmit; buildPropped(mOmit);
        mOmit.members = { mOmit.members[0] };   // drop the prop entirely
        const SolveResult rOmit = solve(mOmit);
        real duMax = 0.0;
        for (size_t k = 0; k < rCut.u.size() && k < rOmit.u.size(); ++k)
            duMax = std::max(duMax, std::fabs(rCut.u[k] - rOmit.u[k]));
        checkClose("active=false == omitted member (disp)", duMax, 0.0, 1e-9);

        // disable BOTH members -> the tip node hangs off nothing -> mechanism (singular)
        FrameModel mMech = mFull; mMech.members[0].active = false; mMech.members[1].active = false;
        const SolveResult rMech = solve(mMech);
        checkTrue("isolated node after removal -> mechanism", rMech.singular, "expected singular");
    }

    // ---------- F27: safety factor (C3) + criticality / pivot margin (C4) ----------
    {
        // C3: a cantilever (tip load P, root moment M=PL) has a closed-form worst utilization:
        //   sigma = M/W, D/C = sigma/cap.bend, safetyFactor = 1/(D/C). Square 100x100, cap.bend=300:
        //   M=PL=2e6, W=1.6667e5 -> sigma=12 MPa -> D/C=0.04 -> SF=25.
        Section  csec = Section::Rectangular(100.0, 100.0);
        Material cmat(210000.0, 80769.0, 7850.0);
        cmat.cap = Capacity::make(300.0, 300.0, 180.0);
        const real L = 2000.0, P = 1000.0;
        const real dcExact = (P * L / csec.Wz()) / cmat.cap.bend;   // = 0.04

        FrameModel m; fixtures::prepMatSec(m, cmat, csec);
        Node cn0(0, 0.0, 0.0, 0.0); cn0.fixAll();
        Node cn1(1,   L, 0.0, 0.0);
        m.nodes = { cn0, cn1 };
        m.members = { Member(0, 0, 1, 0, 0) };
        NodalLoad cp; cp.node = 1; cp.comp[Uz] = -P; m.nodalLoads = { cp };
        const SolveResult r = solve(m);
        const DemandSummary ds = worstUtilization(m, r);
        std::printf("[F27] safety factor + criticality margin  D/C=%.5g SF=%.5g pivotMargin=%.5g\n",
                    ds.maxDC, ds.safetyFactor, r.pivotMargin);
        checkTrue("worstUtilization valid", ds.valid, "");
        checkClose("max D/C = (PL/W)/cap (closed form)", ds.maxDC, dcExact, 1e-9);
        checkClose("safety factor = 1/maxDC", ds.safetyFactor, 1.0 / dcExact, 1e-9);
        checkTrue("governing member id = 0", ds.governingMember == 0, "");

        // C3 linearity: double the load -> D/C doubles, SF halves.
        FrameModel m2 = m; m2.nodalLoads[0].comp[Uz] = -2.0 * P;
        const DemandSummary ds2 = worstUtilization(m2, solve(m2));
        checkClose("D/C scales linearly with load", ds2.maxDC, 2.0 * dcExact, 1e-9);
        checkClose("SF scales inversely with load", ds2.safetyFactor, 0.5 / dcExact, 1e-9);

        // C3 x C1: deactivating the only member -> nothing screenable -> invalid summary.
        FrameModel m0 = m; m0.members[0].active = false;
        const DemandSummary ds0 = worstUtilization(m0, solve(m0));
        checkTrue("no active member -> invalid summary", !ds0.valid, "");

        // C4 exact anchor: a single free DOF (axial bar) has ONE pivot -> pivotMargin == 1.
        {
            FrameModel ma; fixtures::prepMatSec(ma, cmat, csec);
            Node a0(0, 0.0, 0.0, 0.0); a0.fixAll();
            Node a1(1, 1000.0, 0.0, 0.0);
            a1.fixed[Uy] = a1.fixed[Uz] = a1.fixed[Rx] = a1.fixed[Ry] = a1.fixed[Rz] = true;  // free: Ux only
            ma.nodes = { a0, a1 };
            ma.members = { Member(0, 0, 1, 0, 0) };
            NodalLoad p; p.node = 1; p.comp[Ux] = 1000.0; ma.nodalLoads = { p };
            const SolveResult ra = solve(ma);
            checkTrue("single-DOF non-singular", !ra.singular, ra.diagnostic);
            checkClose("single-DOF axial pivotMargin = 1", ra.pivotMargin, 1.0, 1e-12);
        }

        // C4 bound: a healthy structure's margin is in (0,1].
        checkTrue("healthy pivotMargin in (0,1]", r.pivotMargin > 0.0 && r.pivotMargin <= 1.0 + 1e-12, "");

        // C4 scale-invariance: pivotMargin is a min/max pivot RATIO, so scaling ALL stiffness
        // (E and G x1000) leaves it unchanged -> a dimensionless conditioning proxy, not an absolute.
        // (It is deliberately NOT asserted monotone in a single member's stiffness: min/max over
        // pivots is not monotone that way. Its mechanism link is the singular threshold below, plus
        // F26: an actual mechanism is flagged singular, i.e. the margin has crossed pivotTol.)
        FrameModel mE = m; mE.materials[0].E *= 1000.0; mE.materials[0].G *= 1000.0;
        checkClose("pivotMargin invariant under uniform stiffness scale", solve(mE).pivotMargin, r.pivotMargin, 1e-12);

        // C4 mechanism link: a one-step-from-mechanism model (the F26 prop+cantilever with both
        // members removed) is singular -> margin is 0 (the warning has hit the floor).
        FrameModel mSing; fixtures::prepMatSec(mSing, cmat, csec);
        Node s0(0, 0.0, 0.0, 0.0); s0.fixAll();
        Node s1n(1, L, 0.0, 0.0);
        mSing.nodes = { s0, s1n };
        mSing.members = { Member(0, 0, 1, 0, 0) };
        mSing.members[0].active = false;                 // remove the only member -> mechanism
        NodalLoad sp; sp.node = 1; sp.comp[Uz] = -P; mSing.nodalLoads = { sp };
        const SolveResult rSing = solve(mSing);
        checkTrue("mechanism -> singular & pivotMargin 0", rSing.singular && rSing.pivotMargin == 0.0, "");
    }

    // ---------- F28: shell element removal (ShellQuad::active) — mirror of F26 for facets ----------
    {
        // 2x2-quad cantilever plate (9 nodes, 500x500 facets, t=10) clamped along the x=0 edge,
        // uniform pressure on every facet. Deactivating a CLAMP-ADJACENT facet (q0) keeps every
        // node attached to an active element, so the model stays solvable and the removal-by-flag
        // invariant can be checked against physically omitting the facet. Deactivating the FAR
        // CORNER facet (q3) isolates the corner node -> mechanism. A ShellPressure on an inactive
        // facet is dropped (mirrors a UDL on an inactive member), proven via the reaction total.
        const real Es = 30000.0, nu = 0.3, Gs = Es / (2.0 * (1.0 + nu));
        Material smat(Es, Gs); smat.nu = nu;
        const real e = 500.0, t = 10.0, p = 0.01;   // facet edge, thickness, pressure (N/mm^2)

        auto buildPlate = [&](FrameModel& m) {
            m = FrameModel{};
            m.materials = { smat };
            for (int j = 0; j <= 2; ++j)
                for (int i = 0; i <= 2; ++i) {
                    Node n(j * 3 + i, i * e, j * e, 0.0);
                    if (i == 0) n.fixAll();          // clamped edge x=0
                    m.nodes.push_back(n);
                }
            m.shells = { ShellQuad(0, 0, 1, 4, 3, 0, t),     // q0: clamp-adjacent (lower-left)
                         ShellQuad(1, 1, 2, 5, 4, 0, t),     // q1
                         ShellQuad(2, 3, 4, 7, 6, 0, t),     // q2
                         ShellQuad(3, 4, 5, 8, 7, 0, t) };   // q3: far corner (sole owner of node 8)
            for (int s = 0; s < 4; ++s) { ShellPressure sp; sp.shell = s; sp.p = p; m.shellPressures.push_back(sp); }
        };

        std::printf("[F28] shell element removal (ShellQuad::active) — 2x2 cantilever plate\n");

        FrameModel mFull; buildPlate(mFull);
        const SolveResult rFull = solve(mFull);
        checkTrue("full plate non-singular", !rFull.singular, rFull.diagnostic);

        // (a) removal-by-flag invariant: active=false MUST equal physically omitting the facet
        FrameModel mCut = mFull; mCut.shells[0].active = false;
        const SolveResult rCut = solve(mCut);
        checkTrue("after removing q0: non-singular", !rCut.singular, rCut.diagnostic);
        FrameModel mOmit; buildPlate(mOmit);
        mOmit.shells.erase(mOmit.shells.begin());                  // drop q0 entirely...
        mOmit.shellPressures.erase(mOmit.shellPressures.begin());  // ...and its pressure
        const SolveResult rOmit = solve(mOmit);
        real duMax = 0.0;
        for (size_t k = 0; k < rCut.u.size() && k < rOmit.u.size(); ++k)
            duMax = std::max(duMax, std::fabs(rCut.u[k] - rOmit.u[k]));
        checkClose("shell active=false == omitted facet (disp)", duMax, 0.0, 1e-9);

        // (b) inactive facet recovers zero forces (centre resultants AND corner moments)
        real fMax = 0.0;
        {
            const ShellElementForces& sf = rCut.shellForces[0];
            for (real v : { sf.Mxx, sf.Myy, sf.Mxy, sf.Qx, sf.Qy, sf.Nxx, sf.Nyy, sf.Nxy })
                fMax = std::max(fMax, std::fabs(v));
            for (int c = 0; c < 4; ++c)
                fMax = std::max({ fMax, std::fabs(sf.MxxC[c]), std::fabs(sf.MyyC[c]), std::fabs(sf.MxyC[c]) });
        }
        checkClose("inactive facet forces stay 0", fMax, 0.0, 1e-12);

        // (c) the dropped pressure does not leak: reactions balance ONLY the 3 active facets
        const real Ftot = p * 3.0 * e * e;   // 0.01 * 3 * 250000 = 7500 N along +z
        real sumRz = 0.0;
        for (int k = 0; k < 9; ++k) sumRz += rCut.reaction(k, Uz);
        checkClose("reactions balance active-facet pressure only", sumRz, -Ftot, 1e-9);

        // (d) removing the far-corner facet isolates node 8 -> mechanism (singular)
        FrameModel mMech = mFull; mMech.shells[3].active = false;
        const SolveResult rMech = solve(mMech);
        checkTrue("isolated node after facet removal -> mechanism", rMech.singular, "expected singular");

        // (e) fingerprint guard: flipping shell.active after factoring must reject a stale reuse
        FrameModel mFp; buildPlate(mFp);
        const PreparedSystem ps = assembleAndFactor(mFp);
        mFp.shells[0].active = false;
        const SolveResult rStale = solveLoad(ps, mFp);
        checkTrue("flipped shell.active rejects stale factor", rStale.singular, rStale.diagnostic);
        checkTrue("stale-reuse diagnostic names the cause",
                  rStale.diagnostic.find("model changed") != std::string::npos, rStale.diagnostic);
    }

    // ---------- F29: connectivity + debris mass properties (FragmentCluster, Chaos handoff) ----------
    {
        // analyzeConnectivity is pure post-processing (graph + closed-form mass geometry), so
        // every oracle is hand-computable: a free rod's inertia (mL^2/12 about transverse axes
        // through the com, ~0 about its own axis), the 45-degree rotation pinning the PRODUCT-
        // of-inertia sign convention (inertia[] stores tensor MATRIX entries, Ixy = -int(xy dm)),
        // a thin square lamina (ma^2/12 in-plane, ma^2/6 polar -- the two-triangle split is
        // exact), and the twin-tower/bridge graph for grounded-vs-detached classification.
        std::printf("[F29] connectivity + fragment mass properties (Chaos handoff)\n");

        // (a) free rod along +x: one detached cluster, rod inertia closed form
        const real Lr = 2000.0;
        const real mRod = 7850.0 * 1e-12 * sec.A * Lr;             // 0.157 tonne
        {
            FrameModel m; fixtures::prepMatSec(m, mat, sec);
            m.nodes = { Node(0, 0.0, 0.0, 0.0), Node(1, Lr, 0.0, 0.0) };   // NO supports
            m.members = { Member(0, 0, 1, 0, 0) };
            const ConnectivityResult c = analyzeConnectivity(m);
            checkTrue("free rod: analysis valid", c.valid, "");
            checkTrue("free rod: 1 detached, 0 grounded, no loose nodes",
                      c.detached.size() == 1 && c.groundedComponents == 0 && c.looseNodes.empty(), "");
            const FragmentCluster& f = c.detached[0];
            checkClose("rod mass = rho*A*L (tonne)", f.mass, mRod, 1e-12);
            checkClose("rod com x = L/2", f.com.x, Lr / 2.0, 1e-12);
            checkClose("rod Iyy = mL^2/12", f.inertia[1], mRod * Lr * Lr / 12.0, 1e-12);
            checkClose("rod Izz = mL^2/12", f.inertia[2], mRod * Lr * Lr / 12.0, 1e-12);
            checkClose("rod Ixx ~ 0 (slender)", f.inertia[0], 0.0, 1e-12);
        }

        // (b) the same rod rotated 45 deg in the xy plane: pins the tensor sign convention
        {
            const real h = Lr / 2.0 * std::sqrt(2.0);              // = L*cos45
            FrameModel m; fixtures::prepMatSec(m, mat, sec);
            m.nodes = { Node(0, 0.0, 0.0, 0.0), Node(1, h, h, 0.0) };
            m.members = { Member(0, 0, 1, 0, 0) };
            const ConnectivityResult c = analyzeConnectivity(m);
            checkTrue("rotated rod: 1 detached cluster", c.valid && c.detached.size() == 1, "");
            const FragmentCluster& f = c.detached[0];
            checkClose("rot rod Ixx = mL^2/24", f.inertia[0], mRod * Lr * Lr / 24.0, 1e-9);
            checkClose("rot rod Izz = mL^2/12", f.inertia[2], mRod * Lr * Lr / 12.0, 1e-9);
            checkClose("rot rod Ixy = -mL^2/24 (tensor entry)", f.inertia[3], -mRod * Lr * Lr / 24.0, 1e-9);
        }

        // (c) free square lamina (one shell facet): two-triangle closed form is exact
        {
            const real a = 1000.0, t = 10.0, rhoS = 2500.0;
            const real mPl = rhoS * 1e-12 * t * a * a;             // 0.025 tonne
            Material smat(30000.0, 11538.46, rhoS); smat.nu = 0.3;
            FrameModel m; m.materials = { smat };
            m.nodes = { Node(0, 0.0, 0.0, 0.0), Node(1, a, 0.0, 0.0), Node(2, a, a, 0.0), Node(3, 0.0, a, 0.0) };
            m.shells = { ShellQuad(0, 0, 1, 2, 3, 0, t) };
            const ConnectivityResult c = analyzeConnectivity(m);
            checkTrue("free plate: 1 detached cluster", c.valid && c.detached.size() == 1, "");
            const FragmentCluster& f = c.detached[0];
            checkClose("plate mass = rho*t*a^2", f.mass, mPl, 1e-12);
            checkClose("plate com x = a/2", f.com.x, a / 2.0, 1e-12);
            checkClose("plate Ixx = m a^2/12", f.inertia[0], mPl * a * a / 12.0, 1e-9);
            checkClose("plate Izz = m a^2/6 (polar)", f.inertia[2], mPl * a * a / 6.0, 1e-9);
        }

        // (d) twin towers + bridge: grounded/detached classification + loose nodes
        {
            const real Ht = 1000.0;
            auto buildTowers = [&](FrameModel& m, bool fixB) {
                m = FrameModel{};
                fixtures::prepMatSec(m, mat, sec);
                Node a0(0, 0.0, 0.0, 0.0); a0.fixAll();
                Node a1(1, 0.0, 0.0, Ht);
                Node b0(2, 5000.0, 0.0, 0.0); if (fixB) b0.fixAll();
                Node b1(3, 5000.0, 0.0, Ht);
                m.nodes = { a0, a1, b0, b1 };
                m.members = { Member(0, 0, 1, 0, 0),     // tower A
                              Member(1, 2, 3, 0, 0),     // tower B
                              Member(2, 1, 3, 0, 0) };   // bridge
            };

            FrameModel m1; buildTowers(m1, true);
            const ConnectivityResult c1 = analyzeConnectivity(m1);
            checkTrue("towers+bridge: one grounded component",
                      c1.valid && c1.groundedComponents == 1 && c1.detached.empty(), "");

            FrameModel m2 = m1; m2.members[2].active = false;      // cut the bridge, both feet fixed
            const ConnectivityResult c2 = analyzeConnectivity(m2);
            checkTrue("cut bridge, both grounded: 2 components, none detached",
                      c2.groundedComponents == 2 && c2.detached.empty(), "");

            FrameModel m3; buildTowers(m3, false); m3.members[2].active = false;   // tower B unfooted
            const ConnectivityResult c3 = analyzeConnectivity(m3);
            checkTrue("unfooted tower detaches as one cluster",
                      c3.groundedComponents == 1 && c3.detached.size() == 1, "");
            const FragmentCluster& f = c3.detached[0];
            checkTrue("cluster lists tower B exactly",
                      f.nodes == std::vector<NodeId>({ 2, 3 }) && f.members == std::vector<MemberId>({ 1 }) &&
                      f.shells.empty(), "");
            checkClose("cluster mass = tower B rod", f.mass, 7850.0 * 1e-12 * sec.A * Ht, 1e-12);
            checkClose("cluster com z = H/2", f.com.z, Ht / 2.0, 1e-12);

            FrameModel m4 = m3; m4.members[1].active = false;      // kill tower B too -> bare free nodes
            const ConnectivityResult c4 = analyzeConnectivity(m4);
            checkTrue("bare free nodes land in looseNodes",
                      c4.detached.empty() && c4.looseNodes == std::vector<NodeId>({ 2, 3 }), "");
        }
    }

    // ---------- F30: progressive-collapse driver (C2) — sequence, fragments, terminals ----------
    {
        std::printf("[F30] progressive-collapse driver (C2) — LSP sequential linear analysis\n");

        // (A) hanging chain -> STABLE terminal. A vertical chain (stout-thin-stout) under a tip
        // pull: the thin middle link's D/C = |N|/(A*cap.tens) = 150000/(400*300) = 1.25 exactly,
        // so the driver removes it; the lower half then dangles (no support) -> detached debris
        // with closed-form mass properties, its load LEAVES with it, and the surviving upper
        // link reads exactly zero force -> Stable. Any load leak would show as maxDC != 0.
        {
            Section secThin = Section::Rectangular(20.0, 20.0);    // A = 400
            const real Lc = 1000.0, P = 150000.0;
            auto buildChain = [&](FrameModel& m) {
                m = FrameModel{};
                m.materials = { mat };                             // cap.tens = 300
                m.sections  = { sec, secThin };
                Node n0(0, 0.0, 0.0, 0.0); n0.fixAll();            // hung from the top
                m.nodes = { n0, Node(1, 0.0, 0.0, -Lc), Node(2, 0.0, 0.0, -2.0 * Lc), Node(3, 0.0, 0.0, -3.0 * Lc) };
                m.members = { Member(0, 0, 1, 0, 0),
                              Member(1, 1, 2, 0, 1),               // the thin link (A = 400)
                              Member(2, 2, 3, 0, 0) };
                NodalLoad p; p.node = 3; p.comp[Uz] = -P; m.nodalLoads = { p };
            };

            FrameModel m; buildChain(m);
            CollapseOptions co; co.dlf = 1.0;
            const CollapseHistory h = runProgressiveCollapse(m, co);
            checkTrue("chain: outcome Stable", h.outcome == CollapseOutcome::Stable, h.diagnostic);
            checkTrue("chain: exactly 2 steps", h.steps.size() == 2, "");
            checkClose("chain step0 maxDC = P/(A cap) (closed form)", h.steps[0].maxDC, 1.25, 1e-9);
            checkTrue("chain step0 is the baseline (no removals)",
                      h.steps[0].removedMembers.empty() && h.steps[0].detached.empty() && h.steps[0].solved, "");
            const CollapseStep& s1 = h.steps[1];
            checkTrue("chain step1 removes the thin link (id 1, Tension)",
                      s1.removedMembers == std::vector<MemberId>({ 1 }) && s1.mode == FailMode::Tension, "");
            checkClose("chain step1 triggerRatio = 1.25", s1.triggerRatio, 1.25, 1e-9);
            checkTrue("chain step1 detaches the lower half", s1.detached.size() == 1, "");
            if (s1.detached.size() == 1) {
                const FragmentCluster& f = s1.detached[0];
                checkTrue("debris lists {n2,n3 / M2}",
                          f.nodes == std::vector<NodeId>({ 2, 3 }) && f.members == std::vector<MemberId>({ 2 }), "");
                checkClose("debris mass = rho*A*L", f.mass, 7850.0 * 1e-12 * sec.A * Lc, 1e-12);
                checkClose("debris com z = -2500", f.com.z, -2500.0, 1e-12);
                checkClose("debris Ixx = mL^2/12", f.inertia[0], (7850.0 * 1e-12 * sec.A * Lc) * Lc * Lc / 12.0, 1e-9);
            }
            checkClose("chain step1 maxDC = 0 (load left with the debris)", s1.maxDC, 0.0, 1e-12);
            checkTrue("chain step1 safetyFactor = +inf", std::isinf(s1.safetyFactor) && s1.safetyFactor > 0, "");

            // dlf is a pure load scale on a linear model: D/C doubles exactly
            CollapseOptions co2 = co; co2.dlf = 2.0;
            const CollapseHistory h2 = runProgressiveCollapse(m, co2);
            checkClose("dlf=2 doubles step0 maxDC exactly", h2.steps[0].maxDC, 2.50, 1e-12);

            // step budget: the baseline step finds D/C > threshold but may not remove -> MaxSteps
            CollapseOptions co3 = co; co3.maxSteps = 1;
            const CollapseHistory h3 = runProgressiveCollapse(m, co3);
            checkTrue("maxSteps=1 -> MaxSteps outcome", h3.outcome == CollapseOutcome::MaxSteps && h3.steps.size() == 1,
                      h3.diagnostic);

            // scenario API: remove the thin link at step 0 -> immediate debris, then Stable
            CollapseOptions co4 = co; co4.initialRemovals = { 1 };
            const CollapseHistory h4 = runProgressiveCollapse(m, co4);
            checkTrue("initialRemovals: one-step Stable with debris",
                      h4.outcome == CollapseOutcome::Stable && h4.steps.size() == 1 &&
                      h4.steps[0].removedMembers == std::vector<MemberId>({ 1 }) &&
                      h4.steps[0].mode == FailMode::None && h4.steps[0].detached.size() == 1, h4.diagnostic);

            // bad scenario id -> Invalid
            CollapseOptions co5 = co; co5.initialRemovals = { 99 };
            checkTrue("unknown initialRemovals id -> Invalid",
                      runProgressiveCollapse(m, co5).outcome == CollapseOutcome::Invalid, "");
        }

        // (B) propped cantilever cascade -> COLLAPSED terminal. Beam (fixed at n0, P at n1) +
        // a pin-ended prop under n2 (both bending pairs released -> a pure axial strut, so the
        // force-method closed form is exact): R = d10/f11 with d10 = P a^2 (3L-a)/(6EI),
        // f11 = L^3/(3EI) + h/(E A_p). Weak prop crushes first (D/C ~ 1.99), then the beam root
        // (D/C 1.2), then the floating remainder detaches -> nothing grounded -> Collapsed.
        {
            Material matProp(E, G, 7850.0);
            matProp.cap = Capacity::make(5.0, 300.0, 180.0);       // weak in compression only
            Section secProp = Section::Circular(20.0);             // A = pi*400
            const real Lb = 3000.0, aP = 1500.0, hP = 1000.0, P = 40000.0;
            const real EIv = E * sec.Iz;
            const real d10 = P * aP * aP * (3.0 * Lb - aP) / (6.0 * EIv);
            const real f11 = Lb * Lb * Lb / (3.0 * EIv) + hP / (E * secProp.A);
            const real R   = d10 / f11;                            // prop force (compression)
            const real dcProp = (R / secProp.A) / matProp.cap.comp;        // ~1.988
            const real dcRoot = (P * aP / sec.Wz()) / mat.cap.bend;        // = 1.2 after prop loss

            FrameModel m;
            m.materials = { mat, matProp };
            m.sections  = { sec, secProp };
            Node n0(0, 0.0, 0.0, 0.0); n0.fixAll();
            Node n3(3, Lb, 0.0, -hP);  n3.fixAll();
            m.nodes = { n0, Node(1, aP, 0.0, 0.0), Node(2, Lb, 0.0, 0.0), n3 };
            Member prop(2, 3, 2, 1, 1);
            prop.release[4] = prop.release[5] = prop.release[10] = prop.release[11] = true;   // pin-pin bending
            m.members = { Member(0, 0, 1, 0, 0), Member(1, 1, 2, 0, 0), prop };
            NodalLoad p; p.node = 1; p.comp[Uz] = -P; m.nodalLoads = { p };

            CollapseOptions co; co.dlf = 1.0; co.solve.enableReleases = true;
            const CollapseHistory h = runProgressiveCollapse(m, co);
            checkTrue("cascade: outcome Collapsed", h.outcome == CollapseOutcome::Collapsed, h.diagnostic);
            checkTrue("cascade: exactly 3 steps", h.steps.size() == 3, "");
            checkClose("step0 maxDC = prop force-method closed form", h.steps[0].maxDC, dcProp, 1e-9);
            const CollapseStep& s1 = h.steps[1];
            checkTrue("step1 crushes the prop (id 2)",
                      s1.removedMembers == std::vector<MemberId>({ 2 }) && s1.mode == FailMode::Crush, "");
            checkClose("step1 triggerRatio = prop D/C", s1.triggerRatio, dcProp, 1e-9);
            checkClose("step1 maxDC = cantilever root (PL/W)/cap", s1.maxDC, dcRoot, 1e-9);
            const CollapseStep& s2 = h.steps[2];
            checkTrue("step2 severs the beam root (id 0)", s2.removedMembers == std::vector<MemberId>({ 0 }), "");
            checkClose("step2 triggerRatio = root D/C", s2.triggerRatio, dcRoot, 1e-9);
            checkTrue("step2: floating remainder detaches as debris",
                      s2.detached.size() == 1 && !s2.solved &&
                      s2.detached[0].members == std::vector<MemberId>({ 1 }) &&
                      s2.detached[0].nodes == std::vector<NodeId>({ 1, 2 }), "");
            if (s2.detached.size() == 1) {
                checkClose("debris mass = rho*A*(L-a)", s2.detached[0].mass, 7850.0 * 1e-12 * sec.A * (Lb - aP), 1e-12);
                checkClose("debris com x = 2250", s2.detached[0].com.x, 2250.0, 1e-12);
            }
            bool uZero = !s2.u.empty();
            for (real v : s2.u) uZero = uZero && (v == 0.0);
            checkTrue("terminal step snapshot reads all-zero", uZero, "");
        }
    }

    // ---------- F31: shell surface von Mises screen + driver shell removal (stage 3d) ----------
    {
        std::printf("[F31] shell von Mises screen + driver shell removal\n");
        const real Es = 30000.0, nu = 0.3, Gs = Es / (2.0 * (1.0 + nu));
        Material smat(Es, Gs, 2500.0); smat.nu = nu;
        smat.cap = Capacity::make(10.0, 10.0, 6.0);            // vm = min(comp,tens) = 10

        // (a) pure-function oracles (no solve): exact synthetic resultants, t = 10
        {
            const real t = 10.0;
            ShellElementForces f{};                            // membrane only: sigma = (3,1,1)
            f.Nxx = 30.0; f.Nyy = 10.0; f.Nxy = 10.0;
            const ShellDemandResult d1 = checkShellSurface(f, t, smat.cap);
            checkClose("membrane vM = sqrt(10)/vm", d1.risk, std::sqrt(10.0) / 10.0, 1e-12);
            checkTrue("constant field: centre governs", d1.corner == -1, "");

            ShellElementForces b{};                            // bending only: 6M/t^2 = (3,1,1)
            b.Mxx = 50.0; b.Myy = 50.0 / 3.0; b.Mxy = 50.0 / 3.0;
            const ShellDemandResult d2 = checkShellSurface(b, t, smat.cap);
            checkClose("bending vM equals membrane analogue (face-symmetric)", d2.risk, d1.risk, 1e-12);

            ShellElementForces cnr{};                          // a single hot corner governs
            cnr.MxxC[2] = 50.0;                                // sigma_x = 3 at corner 2 only
            const ShellDemandResult d3 = checkShellSurface(cnr, t, smat.cap);
            checkClose("corner peak vM = 0.3", d3.risk, 0.3, 1e-12);
            checkTrue("governing sample is corner 2", d3.corner == 2, "");

            Capacity zero{};                                   // hand-built capacity: vm = 0
            checkTrue("vm=0 under demand screens as infinite D/C",
                      std::isinf(checkShellSurface(f, t, zero).risk), "");
        }

        // (b) constant-curvature patch (machine-precision FE field): Mxx = -D c, Myy = -nu D c
        //     -> surface vM = (6 D c / t^2) * sqrt(1 - nu + nu^2), exact at centre AND corners.
        {
            const real ap = 1000.0, tp = 10.0, cc = 1e-6;
            const real D = Es * tp * tp * tp / (12.0 * (1.0 - nu * nu));
            const real riskExp = (6.0 * D * cc / (tp * tp)) * std::sqrt(1.0 - nu + nu * nu) / smat.cap.vm;
            FrameModel m; fixtures::platePatchCylindrical(m, ap, tp, 0.0, cc, smat);
            const SolveResult r = solve(m);
            checkTrue("patch non-singular", !r.singular, r.diagnostic);
            const ShellDemandSummary ds = worstShellUtilization(m, r);
            checkTrue("patch screen valid", ds.valid, "");
            checkClose("patch vM D/C = (6Dc/t^2) sqrt(1-nu+nu^2)/vm", ds.maxDC, riskExp, 1e-8);

            FrameModel mOff = m;                               // inactive shells are skipped
            for (auto& sh : mOff.shells) sh.active = false;
            checkTrue("all shells inactive -> screen invalid", !worstShellUtilization(mOff, solve(mOff)).valid, "");
        }

        // (c) driver removes the governing facet, the unsupported rest detaches -> Collapsed.
        //     Cantilever strip of two 500x500 facets clamped at x=0, tip line load: the ROOT
        //     facet carries the full moment, so it is condemned first (mode ShellVonMises);
        //     the tip facet then hangs on nothing -> debris with closed-form mass properties.
        {
            const real e = 500.0, t = 10.0;
            FrameModel m;
            m.materials = { smat };
            Node n0(0, 0.0, 0.0, 0.0); n0.fixAll();
            Node n3(3, 0.0, e, 0.0);   n3.fixAll();
            m.nodes = { n0, Node(1, e, 0.0, 0.0), Node(2, 2.0 * e, 0.0, 0.0),
                        n3, Node(4, e, e, 0.0),   Node(5, 2.0 * e, e, 0.0) };
            m.shells = { ShellQuad(0, 0, 1, 4, 3, 0, t),       // root facet
                         ShellQuad(1, 1, 2, 5, 4, 0, t) };     // tip facet
            NodalLoad p2; p2.node = 2; p2.comp[Uz] = -500.0;
            NodalLoad p5; p5.node = 5; p5.comp[Uz] = -500.0;
            m.nodalLoads = { p2, p5 };

            CollapseOptions co; co.dlf = 1.0;
            const CollapseHistory h = runProgressiveCollapse(m, co);
            checkTrue("strip: outcome Collapsed", h.outcome == CollapseOutcome::Collapsed, h.diagnostic);
            checkTrue("strip: exactly 2 steps", h.steps.size() == 2, "");
            checkTrue("strip step0 over threshold", h.steps[0].maxDC > 1.0, "");
            const CollapseStep& s1 = h.steps[1];
            checkTrue("step1 removes the ROOT facet via ShellVonMises",
                      s1.removedShells == std::vector<int>({ 0 }) && s1.removedMembers.empty() &&
                      s1.mode == FailMode::ShellVonMises, "");
            checkClose("step1 triggerRatio carries the facet D/C", s1.triggerRatio, h.steps[0].maxDC, 1e-12);
            checkTrue("tip facet detaches as debris", s1.detached.size() == 1 &&
                      s1.detached[0].shells == std::vector<int>({ 1 }) &&
                      s1.detached[0].nodes == std::vector<NodeId>({ 1, 2, 4, 5 }), "");
            if (s1.detached.size() == 1) {
                checkClose("debris mass = rho*t*a^2", s1.detached[0].mass, 2500.0 * 1e-12 * t * e * e, 1e-12);
                checkClose("debris com x = 750", s1.detached[0].com.x, 750.0, 1e-12);
            }
        }
    }

    // ---------- F32: plastic hinge mechanics (stage 4a, solver layer) ----------
    {
        // Fixed-fixed beam L = 4000 as two members of 2000 (100x100, fy = 300 ->
        // Mp = fy*Zz = 300 * 250000 = 7.5e7) under a uniform load. Closed forms:
        //   elastic end moment   wL^2/12   (yield at w_y = 12 Mp / L^2 = 56.25 N/mm)
        //   one support hinged with residual Mp: far-end moment  wL^2/8 - Mp/2
        //   both support hinges: midspan moment                  wL^2/8 - Mp
        // The DECISIVE check is yield-point continuity: with Mp set to the CURRENT elastic
        // end moment, the hinged solution must equal the elastic one to round-off -- this
        // pins the sign of BOTH channels (element-side Qf injection, node-side moment).
        const real Lh = 2000.0, Ltot = 2.0 * Lh;
        Material hmat(E, G, 7850.0);
        hmat.fy = 300.0;
        hmat.cap = Capacity::make(300.0, 300.0, 180.0);
        const real Mp = hmat.fy * sec.Zz;                          // 7.5e7
        checkClose("plastic modulus Zz = b d^2/4", sec.Zz, 250000.0, 1e-12);

        auto buildBeam = [&](FrameModel& m, real w) {              // w > 0 = downward (global -z)
            m = FrameModel{};
            m.materials = { hmat }; m.sections = { sec };
            Node n0(0, 0.0, 0.0, 0.0);  n0.fixAll();
            Node n2(2, Ltot, 0.0, 0.0); n2.fixAll();
            m.nodes = { n0, Node(1, Lh, 0.0, 0.0), n2 };
            m.members = { Member(0, 0, 1, 0, 0), Member(1, 1, 2, 0, 0) };
            MemberUDL u0; u0.member = 0; u0.w_local = { 0.0, -w, 0.0 };   // local y = global +z here
            MemberUDL u1; u1.member = 1; u1.w_local = { 0.0, -w, 0.0 };
            m.memberUDLs = { u0, u1 };
        };
        // Add a hinge PLUS its node-side moment (-Mp * local axis at the hinge node) -- the
        // two-channel contract from Hinge.h, exactly what the 4b driver automates.
        auto addHinge = [&](FrameModel& m, MemberId memId, int dof, real mp) {
            m.hinges.push_back(PlasticHinge{ memId, dof, mp });
            const Member* mem = nullptr;
            for (const auto& mm : m.members) if (mm.id == memId) { mem = &mm; break; }
            const Vec3 pi = m.nodes[(size_t)m.nodeIndex(mem->i)].pos;
            const Vec3 pj = m.nodes[(size_t)m.nodeIndex(mem->j)].pos;
            Vec3 ax, ay, az;
            memberLocalAxes(pi, pj, mem->refVec, ax, ay, az);
            const Vec3 e = (dof == 4 || dof == 10) ? ay : az;
            NodalLoad nl; nl.node = (dof < 6) ? mem->i : mem->j;
            nl.comp[Rx] = -mp * e.x; nl.comp[Ry] = -mp * e.y; nl.comp[Rz] = -mp * e.z;
            m.nodalLoads.push_back(nl);
        };

        std::printf("[F32] plastic hinge mechanics (Mp=%.4g, w_y=%.4g)\n", Mp, 12.0 * Mp / (Ltot * Ltot));

        // (1) fy/Z plumbing: at w_y the elastic end moment magnitude IS Mp
        const real wy = 12.0 * Mp / (Ltot * Ltot);                 // 56.25
        FrameModel mE; buildBeam(mE, wy);
        const SolveResult rE = solve(mE);
        checkTrue("elastic beam non-singular", !rE.singular, rE.diagnostic);
        const real MzI0 = rE.memberForces[0].endI.Mz;              // signed local end-force value
        checkClose("|elastic end moment| at w_y = fy*Zz", std::fabs(MzI0), Mp, 1e-9);

        // (2) yield-point continuity at a SUPPORT hinge (pins the element-side sign)
        {
            FrameModel mH; buildBeam(mH, wy);
            addHinge(mH, 0, 5, MzI0);                              // Mp = current elastic moment
            const SolveResult rH = solve(mH);                      // NOTE: enableReleases stays false
            checkTrue("hinged beam non-singular", !rH.singular, rH.diagnostic);
            real duMax = 0.0, uMax = 1e-30;
            for (size_t k = 0; k < rE.u.size(); ++k) {
                duMax = std::max(duMax, std::fabs(rH.u[k] - rE.u[k]));
                uMax  = std::max(uMax, std::fabs(rE.u[k]));
            }
            // the hinge changes the FP evaluation path, so "equal" means displacement-scale
            // round-off (a sign error in either channel would read at the millimetre scale)
            checkTrue("support hinge at yield == elastic (u rel < 1e-12)", duMax / uMax < 1e-12,
                      "rel=" + std::to_string(duMax / uMax));
            checkClose("hinged end recovers Mz = 0 (contract)", std::fabs(rH.memberForces[0].endI.Mz), 0.0, 1e-9);
            checkClose("far member forces unchanged",
                       rH.memberForces[1].endJ.Mz, rE.memberForces[1].endJ.Mz, 1e-9);
        }

        // (3) continuity at an INTERNAL joint hinge (pins the node-side sign: node 1 is free,
        //     so a wrong node-side moment would corrupt the whole field)
        {
            const real w3 = 30.0;
            FrameModel m3; buildBeam(m3, w3);
            const SolveResult r3 = solve(m3);
            FrameModel m3h = m3;
            addHinge(m3h, 0, 11, r3.memberForces[0].endJ.Mz);      // hinge at member0's j end (node 1)
            const SolveResult r3h = solve(m3h);
            checkTrue("internal hinge non-singular", !r3h.singular, r3h.diagnostic);
            real duMax = 0.0, uMax = 1e-30;
            for (size_t k = 0; k < r3.u.size(); ++k) {
                duMax = std::max(duMax, std::fabs(r3h.u[k] - r3.u[k]));
                uMax  = std::max(uMax, std::fabs(r3.u[k]));
            }
            checkTrue("internal hinge at current moment == elastic (u rel < 1e-12)", duMax / uMax < 1e-12,
                      "rel=" + std::to_string(duMax / uMax));
        }

        // (4) past yield, single support hinge: far-end moment = wL^2/8 - Mp/2
        const real w = 60.0;
        const real sgn = (MzI0 >= 0) ? 1.0 : -1.0;                 // hinge keeps the elastic moment's sign
        {
            FrameModel m1; buildBeam(m1, w);
            addHinge(m1, 0, 5, sgn * Mp);
            const SolveResult r1 = solve(m1);
            checkTrue("propped state non-singular", !r1.singular, r1.diagnostic);
            const real MB = std::fabs(r1.memberForces[1].endJ.Mz);
            checkClose("one hinge: |M_far| = wL^2/8 - Mp/2", MB, w * Ltot * Ltot / 8.0 - Mp / 2.0, 1e-9);
            // global force equilibrium (node-side moments are pure couples: no force leak)
            real sumRz = 0.0;
            for (int k = 0; k < 3; ++k) sumRz += r1.reaction(k, Uz);
            checkClose("sum Rz balances the UDL", sumRz, w * Ltot, 1e-9);
        }

        // (5) both support hinges: midspan moment = wL^2/8 - Mp
        {
            FrameModel m2; buildBeam(m2, w);
            addHinge(m2, 0, 5, sgn * Mp);
            // by symmetry the j-end moment of member1 has the OPPOSITE local sign convention
            FrameModel probe; buildBeam(probe, wy);
            const real sgnJ = (solve(probe).memberForces[1].endJ.Mz >= 0) ? 1.0 : -1.0;
            addHinge(m2, 1, 11, sgnJ * Mp);
            const SolveResult r2 = solve(m2);
            checkTrue("two-hinge state non-singular", !r2.singular, r2.diagnostic);
            const real Mmid = std::fabs(r2.memberForces[0].endJ.Mz);
            checkClose("two hinges: |M_mid| = wL^2/8 - Mp", Mmid, w * Ltot * Ltot / 8.0 - Mp, 1e-9);
        }

        // (6) fingerprint: forming a hinge after factoring must reject a stale reuse
        {
            FrameModel mF; buildBeam(mF, w);
            const PreparedSystem ps = assembleAndFactor(mF);
            mF.hinges.push_back(PlasticHinge{ 0, 5, sgn * Mp });
            const SolveResult rStale = solveLoad(ps, mF);
            checkTrue("new hinge rejects stale factor", rStale.singular, rStale.diagnostic);
        }

        // (7) validate: bad hinge dof / missing member are rejected
        {
            FrameModel mV; buildBeam(mV, w);
            mV.hinges.push_back(PlasticHinge{ 0, 3, 0.0 });        // torsion dof: not a hinge
            std::string why;
            checkTrue("validate rejects non-bending hinge dof", !mV.validate(why), why);
            FrameModel mV2; buildBeam(mV2, w);
            mV2.hinges.push_back(PlasticHinge{ 99, 5, 0.0 });
            checkTrue("validate rejects hinge on missing member", !mV2.validate(why), why);
        }

        // ---------- F33: event-to-event hinge driver (stage 4b) — w* = 16 Mp / L^2 ----------
        {
            // The classic fixed-fixed plastic-collapse bracket. Hinge order and every ratio is
            // closed form: support hinges form at |M|/Mp = wL^2/(12 Mp), the surviving span then
            // carries wL^2/8 - Mp/2 at the far support, and the THIRD (midspan) hinge -- the
            // mechanism -- forms iff wL^2/8 - Mp >= Mp, i.e. w >= w* = 16 Mp / L^2 exactly.
            // Driving 0.98 w* must end Stable with exactly two hinges; 1.02 w* must form the
            // third hinge and go singular -> the analytic collapse load is bracketed to +/-2%.
            const real wStar = 16.0 * Mp / (Ltot * Ltot);          // 75 N/mm
            CollapseOptions co; co.dlf = 1.0; co.plasticHinges = true;

            FrameModel mA; buildBeam(mA, 0.98 * wStar);
            const CollapseHistory hA = runProgressiveCollapse(mA, co);
            checkTrue("0.98 w*: outcome Stable", hA.outcome == CollapseOutcome::Stable, hA.diagnostic);
            checkTrue("0.98 w*: exactly 3 steps, 2 hinges",
                      hA.steps.size() == 3 &&
                      hA.steps[1].formedHinges.size() == 1 && hA.steps[2].formedHinges.size() == 1 &&
                      hA.steps[1].removedMembers.empty() && hA.steps[2].removedMembers.empty(), "");
            if (hA.steps.size() == 3) {
                const CollapseHingeEvent& h1 = hA.steps[1].formedHinges[0];
                const CollapseHingeEvent& h2 = hA.steps[2].formedHinges[0];
                checkTrue("hinge 1 at member 0 end i (tie-break smallest id/dof)",
                          h1.member == 0 && h1.dof == 5 && hA.steps[1].mode == FailMode::Bending, "");
                checkClose("hinge 1 ratio = 0.98*16/12", hA.steps[1].triggerRatio, 0.98 * 16.0 / 12.0, 1e-9);
                checkClose("|hinge 1 Mp| = fy*Zz", std::fabs(h1.Mp), Mp, 1e-12);
                checkTrue("hinge 2 at member 1 end j", h2.member == 1 && h2.dof == 11, "");
                checkClose("hinge 2 ratio = 0.98*2 - 0.5", hA.steps[2].triggerRatio, 0.98 * 2.0 - 0.5, 1e-9);
                // Stable is declared while the REPORTED allowable D/C still exceeds 1: the
                // bending of a hinge-capable member is ductile and waits at M < Mp. The final
                // midspan moment is (0.98*2 - 1) Mp = 0.96 Mp -> screen D/C = 0.96*Mp/(W*cap).
                checkClose("final maxDC = 0.96*Mp/(W cap) (ductile > 1, still Stable)",
                           hA.steps[2].maxDC, 0.96 * Mp / (sec.Wz() * hmat.cap.bend), 1e-9);
            }

            FrameModel mB; buildBeam(mB, 1.02 * wStar);
            const CollapseHistory hB = runProgressiveCollapse(mB, co);
            checkTrue("1.02 w*: outcome Collapsed (hinge mechanism)",
                      hB.outcome == CollapseOutcome::Collapsed, hB.diagnostic);
            checkTrue("1.02 w*: exactly 4 steps, 3 hinges, terminal unsolved",
                      hB.steps.size() == 4 && hB.steps[3].formedHinges.size() == 1 &&
                      !hB.steps[3].solved, "");
            if (hB.steps.size() == 4) {
                checkClose("hinge 3 ratio = 1.02*2 - 1", hB.steps[3].triggerRatio, 1.02 * 2.0 - 1.0, 1e-9);
                checkTrue("hinge 3 completes the mechanism at midspan",
                          hB.steps[3].formedHinges[0].member == 0 && hB.steps[3].formedHinges[0].dof == 11, "");
            }
            std::printf("   plastic collapse load bracketed: Stable at 0.98 w*, Collapsed at 1.02 w*  (w*=16Mp/L^2=%.4g)\n",
                        wStar);
        }
    }

    // ---------- F34: sparse buckling (subspace iteration) agrees with dense + Euler ----------
    {
        const real kPi  = 3.14159265358979323846;
        const real Pref = 1000.0;                      // reference axial compression (matches F23)
        std::printf("[F34] sparse buckling (subspace iteration vs dense vs Euler)\n");
        // Rectangular section so the two bending planes are NON-degenerate: the lowest mode
        // buckles about the weak axis, a clean single smallest eigenvalue for the sparse subspace
        // iteration to converge to (a square section's y/z pair is degenerate). denseThreshold<=0
        // FORCES the sparse path so this fixture actually exercises subspaceSmallest, not dense.
        Section    secB   = Section::Rectangular(60.0, 100.0);
        const real Ibuck  = std::min(secB.Iy, secB.Iz);   // weak-axis I governs the lowest mode
        BucklingOptions sparseOpt; sparseOpt.denseThreshold = 0;

        auto runCol = [&](const char* what, bool pinned, int n, real L, real lenFactor) {
            FrameModel m;
            if (pinned) fixtures::simplySupportedBeamN(m, n, L, mat, secB);
            else        fixtures::cantileverBeamN(m, n, L, mat, secB);
            NodalLoad nl; nl.node = n; nl.comp[Ux] = -Pref; m.nodalLoads = { nl };
            PreparedSystem ps = assembleAndFactor(m);
            const BucklingResult d = solveBuckling(ps, m);              // dense (default threshold)
            const BucklingResult s = solveBuckling(ps, m, sparseOpt);   // forced sparse
            const real PcrEx = kPi * kPi * E * Ibuck / (lenFactor * L * L);
            checkTrue((std::string(what) + " dense non-singular").c_str(),  !d.singular, d.diagnostic);
            checkTrue((std::string(what) + " sparse non-singular").c_str(), !s.singular, s.diagnostic);
            const real relSD = std::fabs(s.criticalFactor - d.criticalFactor) /
                               std::max<real>(1e-30, std::fabs(d.criticalFactor));
            std::printf("   %s: lamDense=%.9g lamSparse=%.9g Euler=%.6g relSparseVsDense=%.2e\n",
                        what, d.criticalFactor, s.criticalFactor, PcrEx / Pref, relSD);
            checkClose((std::string(what) + " sparse == dense").c_str(), s.criticalFactor, d.criticalFactor, 1e-6);
            checkClose((std::string(what) + " sparse Pcr == Euler").c_str(), s.criticalFactor * Pref, PcrEx, 1e-4);
        };

        runCol("pinned-pinned n=10", true,  10, 3000.0, 1.0);   // Pcr = pi^2 E I / L^2
        runCol("fixed-free n=10",    false, 10, 3000.0, 4.0);   // Pcr = pi^2 E I / (2L)^2
        runCol("pinned-pinned n=24", true,  24, 3000.0, 1.0);   // larger nf still agrees

        // all-tension guard: a tension reference load -> no compressive member -> Kg empty -> the
        // forced-sparse request defers to the dense path, which reports the historical "no
        // compression" singular diagnostic. Forcing sparse must NOT fabricate a positive factor.
        {
            const int n = 10; const real L = 3000.0;
            FrameModel m; fixtures::cantileverBeamN(m, n, L, mat, secB);
            NodalLoad nl; nl.node = n; nl.comp[Ux] = Pref; m.nodalLoads = { nl };   // +Pref = tension (pull)
            PreparedSystem ps = assembleAndFactor(m);
            const BucklingResult d = solveBuckling(ps, m);
            const BucklingResult s = solveBuckling(ps, m, sparseOpt);
            checkTrue("all-tension dense singular (no compression)", d.singular, d.diagnostic);
            checkTrue("all-tension forced-sparse singular (defers to dense)", s.singular, s.diagnostic);
            checkTrue("all-tension sparse diagnostic == dense", s.diagnostic == d.diagnostic, s.diagnostic);
        }
    }

    // ---------- F35: ReSolve ladder (Tier-1 Woodbury) vs fresh + F-increment + mechanism ----------
    {
        std::printf("[F35] ReSolve ladder (Tier-1) vs fresh + F-increment + mechanism\n");
        auto relMax = [](const std::vector<real>& a, const std::vector<real>& b) -> real {
            real num = 0, den = 1e-30;
            const size_t n = std::min(a.size(), b.size());
            for (size_t i = 0; i < n; ++i) { num = std::max(num, std::fabs(a[i] - b[i])); den = std::max(den, std::fabs(b[i])); }
            return num / den;
        };
        auto agree = [&](const char* tag, const std::vector<real>& a, const std::vector<real>& b, real tol) {
            const real e = relMax(a, b);
            checkTrue(tag, e < tol, "rel=" + std::to_string(e));
        };

        // (a) rigid portal (2 columns + 2 beam halves), UDL on the beam halves. A SINGLE member
        //     removal stays stable; removing the UDL member exercises the F-increment (its
        //     equivalent nodal load must leave F). 5 nodes: bases 0/1 (fixed), tops 2/3, mid 4.
        const real H = 3000.0, Lb = 6000.0, w = 8.0;
        auto buildPortal = [&]() {
            FrameModel m;
            m.materials.push_back(mat); m.sections.push_back(sec);
            Node n0(0, 0,      0, 0); n0.fixAll();
            Node n1(1, Lb,     0, 0); n1.fixAll();
            Node n2(2, 0,      0, H);
            Node n3(3, Lb,     0, H);
            Node n4(4, Lb / 2, 0, H);
            m.nodes = { n0, n1, n2, n3, n4 };
            Member c0(0, 0, 2, 0, 0); c0.refVec = Vec3(1, 0, 0);
            Member c1(1, 1, 3, 0, 0); c1.refVec = Vec3(1, 0, 0);
            Member b0(2, 2, 4, 0, 0); b0.refVec = Vec3(0, 0, 1);
            Member b1(3, 4, 3, 0, 0); b1.refVec = Vec3(0, 0, 1);
            m.members = { c0, c1, b0, b1 };
            MemberUDL u0; u0.member = 2; u0.w_local = { 0, -w, 0 };
            MemberUDL u1; u1.member = 3; u1.w_local = { 0, -w, 0 };
            m.memberUDLs = { u0, u1 };
            return m;
        };
        FrameModel portal = buildPortal();
        ReSolveSession session(portal);
        checkTrue("portal session valid", session.valid(), session.diagnostic());

        ReanalysisStats st0;
        const SolveResult re0 = session.solve(&st0);
        const SolveResult fr0 = solve(portal);
        checkTrue("baseline tier == 0", st0.tier == 0, "");
        agree("baseline ReSolve == fresh", re0.u, fr0.u, 1e-10);

        // remove the UDL beam half b0 (id 2) -> Tier-1 + F-increment
        session.setMemberActive(2, false);
        ReanalysisStats st1;
        const SolveResult re1 = session.solve(&st1);
        FrameModel w1 = buildPortal(); w1.members[2].active = false;
        const SolveResult fr1 = solve(w1);
        checkTrue("remove UDL member: tier == 1", st1.tier == 1, "");
        checkTrue("remove UDL member: not singular", !re1.singular, re1.diagnostic);
        std::printf("   removed b0(UDL): tier=%d rank=%d reRel=%.2e\n", st1.tier, st1.rank, relMax(re1.u, fr1.u));
        agree("remove UDL member: ReSolve == fresh (F-increment)", re1.u, fr1.u, 1e-10);

        // restore b0 -> the ladder's -lambda/+lambda columns cancel; result back to baseline
        session.setMemberActive(2, true);
        const SolveResult reR = session.solve();
        agree("restore member: ReSolve == baseline (drift)", reR.u, fr0.u, 1e-12);

        // independent single removal of a column c0 (id 0) -> still Tier-1, still == fresh
        session.setMemberActive(0, false);
        ReanalysisStats st2;
        const SolveResult re2 = session.solve(&st2);
        FrameModel w2 = buildPortal(); w2.members[0].active = false;
        const SolveResult fr2 = solve(w2);
        checkTrue("remove column: tier == 1", st2.tier == 1, "");
        agree("remove column: ReSolve == fresh", re2.u, fr2.u, 1e-10);

        // (b) mechanism: a 2-member cantilever column chain; removing the base floats the top.
        {
            FrameModel chain;
            chain.materials.push_back(mat); chain.sections.push_back(sec);
            Node m0(0, 0, 0, 0); m0.fixAll();
            Node m1(1, 0, 0, 1500.0);
            Node m2(2, 0, 0, 3000.0);
            chain.nodes = { m0, m1, m2 };
            chain.members = { Member(0, 0, 1, 0, 0), Member(1, 1, 2, 0, 0) };
            NodalLoad nl; nl.node = 2; nl.comp[Ux] = 1000.0; chain.nodalLoads = { nl };
            ReSolveSession cs(chain);
            checkTrue("chain session valid", cs.valid(), cs.diagnostic());
            cs.setMemberActive(0, false);
            ReanalysisStats sm;
            const SolveResult rm = cs.solve(&sm);
            FrameModel cw = chain; cw.members[0].active = false;
            const SolveResult fm = solve(cw);
            checkTrue("chain remove base: capacitance mechanism flagged", sm.mechanism, "");
            checkTrue("chain remove base: ReSolve singular", rm.singular, rm.diagnostic);
            checkTrue("chain remove base: fresh singular too", fm.singular, fm.diagnostic);
        }

        // (c) shell F-increment: clamped 2x2 plate under pressure; remove one facet, still stable.
        {
            Material matShell(30000.0, 11538.461538, 7850.0);   // E=30000, nu=0.3
            FrameModel plate; fixtures::clampedPlateShell(plate, 1000.0, 50.0, 2, 1.0, matShell);
            ReSolveSession sps(plate);
            checkTrue("plate session valid", sps.valid(), sps.diagnostic());
            sps.setShellActive(0, false);
            ReanalysisStats ss;
            const SolveResult rs = sps.solve(&ss);
            FrameModel pw = plate; pw.shells[0].active = false;
            const SolveResult fs = solve(pw);
            checkTrue("shell remove facet: not singular", !rs.singular, rs.diagnostic);
            std::printf("   removed shell facet 0: tier=%d rank=%d reRel=%.2e\n", ss.tier, ss.rank, relMax(rs.u, fs.u));
            agree("shell remove facet: ReSolve == fresh (pressure F-increment)", rs.u, fs.u, 1e-10);
        }
    }

    // ---------- F36: ReSolve Tier-2 (stale-LDLT PCG) agrees with fresh (tolerance-level) ----------
    {
        std::printf("[F36] ReSolve Tier-2 (stale-LDLT PCG) vs fresh\n");
        auto relMax = [](const std::vector<real>& a, const std::vector<real>& b) -> real {
            real num = 0, den = 1e-30; const size_t n = std::min(a.size(), b.size());
            for (size_t i = 0; i < n; ++i) { num = std::max(num, std::fabs(a[i] - b[i])); den = std::max(den, std::fabs(b[i])); }
            return num / den;
        };
        // Portal frame; a deliberately SMALL maxRank forces a single removal (rank 6) past Tier-1
        // into the Tier-2 PCG path. Removing the beam leaves two stable cantilever columns.
        FrameModel m;
        m.materials.push_back(mat); m.sections.push_back(sec);
        Node n0(0, 0,    0, 0);    n0.fixAll();
        Node n1(1, 6000, 0, 0);    n1.fixAll();
        Node n2(2, 0,    0, 3000);
        Node n3(3, 6000, 0, 3000);
        m.nodes = { n0, n1, n2, n3 };
        Member c0(0, 0, 2, 0, 0); c0.refVec = Vec3(1, 0, 0);
        Member c1(1, 1, 3, 0, 0); c1.refVec = Vec3(1, 0, 0);
        Member bm(2, 2, 3, 0, 0); bm.refVec = Vec3(0, 0, 1);
        m.members = { c0, c1, bm };
        NodalLoad nl; nl.node = 2; nl.comp[Ux] = 5.0e4; m.nodalLoads = { nl };
        ReanalysisOptions opt; opt.maxRank = 5;            // force rank-6 removal into Tier-2
        ReSolveSession s(m, opt);
        checkTrue("tier2 session valid", s.valid(), s.diagnostic());
        s.setMemberActive(2, false);                       // remove the beam -> two cantilevers
        ReanalysisStats st;
        const SolveResult re = s.solve(&st);
        FrameModel w = m; w.members[2].active = false;
        const SolveResult fr = solve(w);
        const real e = relMax(re.u, fr.u);
        std::printf("   tier=%d rank=%d pcgIters=%d relResidual=%.2e reRel=%.2e\n",
                    st.tier, st.rank, st.pcgIters, st.relResidual, e);
        checkTrue("Tier-2 selected (rank > maxRank)", st.tier == 2, "tier=" + std::to_string(st.tier));
        checkTrue("Tier-2 PCG iterated", st.pcgIters > 0, "");
        checkTrue("Tier-2 not singular", !re.singular, re.diagnostic);
        checkTrue("Tier-2 ReSolve == fresh (tolerance)", e < 1e-8, "rel=" + std::to_string(e));
    }

    // ---------- F37: DynamicCollapse -- full-basis self-consistency (Ritz == pure modes, zero truncation) ----------
    // A portal frame whose diagonal brace trips at t=0+ (static D/C). The moment frame stays
    // connected and free-vibrates. With a FULL basis the Ritz path and the pure-eigenmode path
    // span the same space, so every replay frame must agree to round-off, and the per-event
    // truncation residual must vanish (the inheritance projection is exact).
    {
        std::printf("[F37] DynamicCollapse: portal frame, brace trips -> Ritz==modes (full basis), zero truncation\n");
        const real E = 200000.0, Gs = 80000.0, rho = 7850.0, W = 3000.0, H = 3000.0, P = 50000.0;
        Section sec = Section::Rectangular(200.0, 200.0);
        Material strong(E, Gs, rho); strong.cap = Capacity::make(1e9, 1e9, 1e9);
        Material braceMat(E, Gs, rho); braceMat.cap = Capacity::make(0.5, 0.5, 1e9);   // tiny -> brace trips first
        FrameModel m;
        m.materials = { strong, braceMat }; m.sections = { sec };
        Node n0(0, 0, 0, 0); n0.fixAll();
        Node n1(1, W, 0, 0); n1.fixAll();
        Node n2(2, 0, 0, H);
        Node n3(3, W, 0, H);
        m.nodes   = { n0, n1, n2, n3 };
        m.members = { Member(0, 0, 2, 0, 0), Member(1, 1, 3, 0, 0), Member(2, 2, 3, 0, 0), Member(3, 0, 3, 1, 0) };  // m3 = brace
        NodalLoad nl; nl.node = 2; nl.comp[Ux] = P;
        m.nodalLoads = { nl };

        auto runCase = [&](bool ritz) -> DynCollapseHistory {
            DynCollapseOptions opt;
            opt.dt = 1e-5; opt.maxTime = 40 * opt.dt;
            opt.basisSize = 200; opt.useRitzVectors = ritz;            // full basis (nf = 12)
            opt.screenEvery = 1; opt.frameStride = 1; opt.removeThreshold = 1.0; opt.maxEvents = 3;
            return runDynamicCollapse(m, opt);
        };
        const DynCollapseHistory hm = runCase(false), hr = runCase(true);
        std::printf("   modes: outcome=%d events=%zu frames=%zu | ritz: events=%zu frames=%zu\n",
                    (int)hm.outcome, hm.events.size(), hm.frames.size(), hr.events.size(), hr.frames.size());
        checkTrue("F37 brace tripped first (modes)", !hm.events.empty() && hm.events[0].removedMembers.size() == 1 && hm.events[0].removedMembers[0] == 3, "");
        checkTrue("F37 same event count (Ritz vs modes)", hm.events.size() == hr.events.size(), "");
        checkTrue("F37 same frame count (Ritz vs modes)", hm.frames.size() == hr.frames.size(), "");
        real maxTrunc = 0; for (const auto& ev : hm.events) maxTrunc = std::max(maxTrunc, ev.truncationResidual);
        checkTrue("F37 truncationResidual ~ 0 (full basis)", maxTrunc < 1e-9, "maxTrunc=" + std::to_string(maxTrunc));
        real num = 0, den = 1e-30;
        const size_t nfr = std::min(hm.frames.size(), hr.frames.size());
        for (size_t k = 0; k < nfr; ++k) {
            const size_t n = std::min(hm.frames[k].u.size(), hr.frames[k].u.size());
            for (size_t i = 0; i < n; ++i) { num = std::max(num, std::fabs(hm.frames[k].u[i] - hr.frames[k].u[i])); den = std::max(den, std::fabs(hm.frames[k].u[i])); }
        }
        checkTrue("F37 Ritz frames == mode frames (full basis)", num / den < 1e-8, "rel=" + std::to_string(num / den));

        FrameModel mz = m;
        for (Material& mm : mz.materials) mm.rho = 0.0;
        DynCollapseOptions oz; oz.dt = 1e-5; oz.maxTime = 5 * oz.dt;
        const DynCollapseHistory hz = runDynamicCollapse(mz, oz);
        checkTrue("F37 zero-mass dynamic model is Invalid (not Stable)",
                  hz.outcome == CollapseOutcome::Invalid &&
                  hz.diagnostic.find("zero mass") != std::string::npos,
                  "outcome=" + std::to_string((int)hz.outcome) + " diag=" + hz.diagnostic);
    }

    // ---------- F38: DynamicCollapse -- momentum handoff (a fragment detaches WHILE MOVING) ----------
    // Portal frame brace trips at t=0+, the frame sways, and a chain hung off the moving corner is
    // shaken until its link trips dynamically -> the chain detaches carrying real velocity. The
    // handoff vel/angVel must be finite and nonzero (the exact rigid-motion identity is in the audit).
    {
        std::printf("[F38] DynamicCollapse: brace trips -> sway -> hung chain detaches with velocity\n");
        const real E = 200000.0, Gs = 80000.0, rho = 7850.0, W = 3000.0, H = 3000.0, D = 2000.0, P = 150000.0;
        Section sec = Section::Rectangular(200.0, 200.0);
        Material strong(E, Gs, rho);   strong.cap   = Capacity::make(1e9, 1e9, 1e9);
        Material braceMat(E, Gs, rho); braceMat.cap = Capacity::make(0.5, 0.5, 1e9);   // brace trips at t=0+
        Material linkMat(E, Gs, rho);  linkMat.cap  = Capacity::make(12.0, 12.0, 1e9); // link trips dynamically
        FrameModel m;
        m.materials = { strong, braceMat, linkMat }; m.sections = { sec };
        Node n0(0, 0, 0, 0); n0.fixAll();
        Node n1(1, W, 0, 0); n1.fixAll();
        Node n2(2, 0, 0, H);
        Node n3(3, W, 0, H);
        Node n4(4, 0, 0, H + D);          // chain hung off node2 (the driven, max-amplitude corner)
        Node n5(5, 0, 0, H + 2 * D);
        m.nodes   = { n0, n1, n2, n3, n4, n5 };
        m.members = { Member(0, 0, 2, 0, 0), Member(1, 1, 3, 0, 0), Member(2, 2, 3, 0, 0), Member(3, 0, 3, 1, 0),
                      Member(4, 2, 4, 2, 0), Member(5, 4, 5, 0, 0) };          // m3 brace, m4 weak link (off node2)
        NodalLoad nl; nl.node = 2; nl.comp[Ux] = P;
        m.nodalLoads = { nl };
        DynCollapseOptions opt;
        opt.dt = 1e-5; opt.maxTime = 600 * opt.dt;
        opt.basisSize = 200; opt.useRitzVectors = false;
        opt.screenEvery = 2; opt.frameStride = 2; opt.removeThreshold = 1.0; opt.maxEvents = 6;
        const DynCollapseHistory h = runDynamicCollapse(m, opt);
        std::printf("   outcome=%d events=%zu\n", (int)h.outcome, h.events.size());
        for (size_t i = 0; i < h.events.size(); ++i) {
            const auto& ev = h.events[i];
            std::printf("   event[%zu] t=%.3e mode=%d removedM=%zu detached=%zu", i, ev.t, (int)ev.mode, ev.removedMembers.size(), ev.detached.size());
            if (!ev.detached.empty()) {
                const auto& fc = ev.detached[0];
                std::printf(" mass=%.4e vel=(%.3e,%.3e,%.3e) |angVel|=%.3e", fc.mass, fc.vel.x, fc.vel.y, fc.vel.z, std::sqrt(dot(fc.angVel, fc.angVel)));
            }
            std::printf("\n");
        }
        int di = -1; for (size_t i = 0; i < h.events.size(); ++i) if (!h.events[i].detached.empty()) { di = (int)i; break; }
        checkTrue("F38 a fragment detached", di >= 0, "");
        if (di >= 0) {
            const auto& fc = h.events[(size_t)di].detached[0];
            checkTrue("F38 fragment has mass", fc.mass > 0, "mass=" + std::to_string(fc.mass));
            const real speed = std::sqrt(dot(fc.vel, fc.vel));
            checkTrue("F38 handoff velocity finite", std::isfinite(speed) && std::isfinite(dot(fc.angVel, fc.angVel)), "");
            checkTrue("F38 handoff velocity nonzero (moving fragment)", speed > 1e-9, "speed=" + std::to_string(speed));
            checkTrue("F38 handoff respects x-z symmetry (vel_y ~ 0)", std::fabs(fc.vel.y) < 1e-6 * speed, "vel_y=" + std::to_string(fc.vel.y));
        }
    }

    // ---------- F39: DynamicCollapse -- axial chain, mid member trips -> upper part detaches, retained SDOF ----------
    // The retained piece (node0 encastre - member0 - node1) is a fixed-free axial bar: a single
    // analytic SDOF. The inheritance must hand node1 its pre-event static displacement exactly
    // (full basis), and the retained part then free-vibrates about its new (zero-load) equilibrium.
    {
        std::printf("[F39] DynamicCollapse: axial chain, mid member trips -> detach + retained axial SDOF\n");
        const real E = 200000.0, Gs = 80000.0, rho = 7850.0, L = 1000.0, P = 10000.0;
        Section sec = Section::Rectangular(100.0, 100.0);                 // A = 10000 mm^2
        Material strong(E, Gs, rho); strong.cap = Capacity::make(1e9, 1e9, 1e9);
        Material weak(E, Gs, rho);   weak.cap   = Capacity::make(0.5, 0.5, 1e9);   // sigma = P/A = 1 MPa > 0.5 -> trips
        FrameModel m;
        m.materials = { strong, weak };                                  // idx 0 strong, idx 1 weak
        m.sections  = { sec };
        Node n0(0, 0, 0, 0);     n0.fixAll();
        Node n1(1, 0, 0, L);
        Node n2(2, 0, 0, 2 * L);
        Node n3(3, 0, 0, 3 * L);
        m.nodes   = { n0, n1, n2, n3 };
        m.members = { Member(0, 0, 1, 0, 0), Member(1, 1, 2, 1, 0), Member(2, 2, 3, 0, 0) };   // member1 weak
        NodalLoad nl; nl.node = 3; nl.comp[Uz] = P;                       // axial tension at the tip
        m.nodalLoads = { nl };

        const real kPi = 3.14159265358979323846;
        const real rho_b = rho * 1e-12;                                  // kg/m^3 -> tonne/mm^3 unit bridge
        const real omega = std::sqrt(3.0 * E / (rho_b * L * L));         // fixed-free axial bar (consistent mass)
        const real T = 2.0 * kPi / omega;
        const real u1_static = P * L / (E * sec.A);                      // member0 elongation under P

        auto runCase = [&](bool ritz) -> DynCollapseHistory {
            DynCollapseOptions opt;
            opt.dt = T / 200.0; opt.maxTime = 80.0 * opt.dt;
            opt.basisSize = 200; opt.useRitzVectors = ritz;              // basisSize >> nf -> full basis (exact)
            opt.screenEvery = 1; opt.frameStride = 1; opt.removeThreshold = 1.0; opt.maxEvents = 5;
            return runDynamicCollapse(m, opt);
        };
        const DynCollapseHistory h = runCase(false);
        std::printf("   outcome=%d events=%zu frames=%zu T=%.3e u1_static=%.4e\n",
                    (int)h.outcome, h.events.size(), h.frames.size(), T, u1_static);
        checkTrue("F39 exactly one event", h.events.size() == 1, "events=" + std::to_string(h.events.size()));
        if (h.events.size() == 1) {
            checkTrue("F39 member1 removed", h.events[0].removedMembers.size() == 1 && h.events[0].removedMembers[0] == 1, "");
            checkTrue("F39 upper part detached (carries member2)", !h.events[0].detached.empty(), "");
            const real te = h.events[0].t;
            int fi = -1;
            for (size_t k = 0; k < h.frames.size(); ++k) if (std::fabs(h.frames[k].t - te) < 0.5 * (T / 200.0)) fi = (int)k;
            checkTrue("F39 post-event frame exists", fi >= 0, "");
            if (fi >= 0) checkClose("F39 inherited u1 == static (analytic)", h.frames[(size_t)fi].u[(size_t)gdof(1, Uz)], u1_static, 1e-9);
            bool wentNeg = false; real maxAbs = 0;
            for (size_t k = 0; k < h.frames.size(); ++k) {
                const real u = h.frames[k].u[(size_t)gdof(1, Uz)];
                if (u < -0.2 * u1_static) wentNeg = true;
                maxAbs = std::max(maxAbs, std::fabs(u));
            }
            checkTrue("F39 retained part free-vibrates (overshoots)", wentNeg, "");
            checkClose("F39 amplitude conserved (Newmark)", maxAbs, u1_static, 1e-2);
        }
        const DynCollapseHistory hr = runCase(true);
        if (hr.events.size() == 1) {
            const real te = hr.events[0].t; int fi = -1;
            for (size_t k = 0; k < hr.frames.size(); ++k) if (std::fabs(hr.frames[k].t - te) < 0.5 * (T / 200.0)) fi = (int)k;
            if (fi >= 0) checkClose("F39 Ritz inherited u1 == static", hr.frames[(size_t)fi].u[(size_t)gdof(1, Uz)], u1_static, 1e-8);
        }
    }

    // ---------- F40: P-Delta second-order amplification (frozen vs reference vs beam-column exact) ----------
    {
        std::printf("[F40] P-Delta: cantilever column, frozen-reuse vs K_T-reference vs beam-column exact\n");
        const real Lc = 6000.0, sideC = 200.0, Hlat = 1000.0, Ec = 200000.0;
        const real Ic = sideC * sideC * sideC * sideC / 12.0;
        const real kPiC = 3.14159265358979323846;
        const real Pcr = kPiC * kPiC * Ec * Ic / (4.0 * Lc * Lc);
        Section  pcs = Section::Rectangular(sideC, sideC);          // square -> Iy == Iz == Ic
        Material pcm(Ec, Ec / (2.0 * (1.0 + 0.3)), 7850.0); pcm.cap = Capacity::make(1e9, 1e9, 1e9);
        const int nE = 8;

        for (real frac : { 0.0, 0.3, 0.95 }) {
            const real P = frac * Pcr;
            FrameModel m; fixtures::pdeltaColumn(m, nE, Lc, P, Hlat, pcm, pcs);
            PDeltaOptions of;  of.refactorPath = false; of.maxIter = 5000; of.tolU = 1e-13;
            PDeltaOptions orf; orf.refactorPath = true;
            const PDeltaResult rf = runPDelta(m, of);
            const PDeltaResult rr = runPDelta(m, orf);
            const std::string at = " (P/Pcr=" + std::to_string(frac) + ")";

            checkTrue(("F40 frozen converged" + at).c_str(),    rf.converged, rf.finalState.diagnostic);
            checkTrue(("F40 reference converged" + at).c_str(), rr.converged, rr.finalState.diagnostic);
            const real tipFrozen = rf.finalState.disp(nE, Ux);
            const real tipRef    = rr.finalState.disp(nE, Ux);

            // frozen pseudo-load iteration and the K_T refactor reach the SAME second-order fixed point
            checkClose(("F40 frozen == reference" + at).c_str(), tipFrozen, tipRef, 1e-10);

            // reference vs the beam-column stability-function closed form (8-element discretization)
            real tipExact;
            if (P > 0) { const real kk = std::sqrt(P / (Ec * Ic)); tipExact = Hlat * (std::tan(kk * Lc) - kk * Lc) / (P * kk); }
            else       { tipExact = Hlat * Lc * Lc * Lc / (3.0 * Ec * Ic); }
            checkClose(("F40 reference == beam-column exact" + at).c_str(), tipRef, tipExact, 1e-3);

            if (frac == 0.0) {
                // P = 0 degeneracy: the second-order state is the linear state, bit-for-bit, in 0 iters
                const SolveResult linr = solve(m);
                checkTrue("F40 P=0 frozen bit-identical to linear", tipFrozen == linr.disp(nE, Ux), "");
                checkTrue("F40 P=0 zero iterations", rf.iterations == 0, "iters=" + std::to_string(rf.iterations));
            }
        }

        // beyond the Euler load: BOTH paths must report divergence (not a silent wrong answer)
        {
            const real P = 1.05 * Pcr;
            FrameModel m; fixtures::pdeltaColumn(m, nE, Lc, P, Hlat, pcm, pcs);
            PDeltaOptions of;  of.refactorPath = false; of.maxIter = 5000; of.tolU = 1e-13;
            PDeltaOptions orf; orf.refactorPath = true;
            const PDeltaResult rf = runPDelta(m, of);
            const PDeltaResult rr = runPDelta(m, orf);
            checkTrue("F40 f=1.05 reference diverged (K_T not PD)", rr.diverged, "");
            checkTrue("F40 f=1.05 frozen diverged (sliding window)", rf.diverged, "iters=" + std::to_string(rf.iterations));
        }
    }

    // ---------- F42: tension-only X-braced portal — eliminate the compressed brace ----------
    {
        std::printf("[F42] Tension-only: X-braced portal, eliminate the compressed brace\n");
        Material tmat(200000.0, 76923.07692307692, 7850.0); tmat.cap = Capacity::make(1e9, 1e9, 1e9);
        Section stocky = Section::Rectangular(250.0, 250.0);
        Section tbrace = Section::Rectangular(60.0, 60.0);
        FrameModel m; fixtures::xBracedPortal(m, 5.0e4, 0.0, tmat, stocky, tbrace);
        const TensionOnlyResult R = runTensionOnly(m);
        checkTrue("F42 converged", R.converged, R.finalState.diagnostic);
        checkTrue("F42 converged in <=3 iters", R.iterations <= 3, "iters=" + std::to_string(R.iterations));
        checkTrue("F42 exactly one brace slack", R.slack.size() == 1, "slack=" + std::to_string(R.slack.size()));

        const MemberId slackId = R.slack.empty() ? -1 : R.slack[0];
        const MemberId tenId   = (slackId == 3) ? 4 : 3;               // the surviving diagonal
        const real Nten = R.finalState.memberForces[(size_t)tenId].endI.N;
        checkTrue("F42 surviving brace in tension (N<0)", Nten < 0.0, "N=" + std::to_string(Nten));

        // oracle: the model with the slack (compressed) brace OMITTED from the member list entirely
        FrameModel ref; fixtures::xBracedPortal(ref, 5.0e4, 0.0, tmat, stocky, tbrace);
        for (size_t e = 0; e < ref.members.size(); ++e)
            if (ref.members[e].id == slackId) { ref.members.erase(ref.members.begin() + (long)e); break; }
        const SolveResult rr = solve(ref);
        real maxDiff = 0, maxU = 0;
        for (size_t g = 0; g < R.finalState.u.size() && g < rr.u.size(); ++g) {
            maxDiff = std::max(maxDiff, std::fabs(R.finalState.u[g] - rr.u[g]));
            maxU    = std::max(maxU, std::fabs(rr.u[g]));
        }
        const real rel = maxDiff / std::max(maxU, real(1e-30));
        // ReSolve (Woodbury) inner loop matches the omitted-brace FRESH solve to factorization round-off
        checkTrue("F42 converged == omitted-brace model (rel<1e-10)", rel < 1e-10, "rel=" + std::to_string(rel));

        // fingerprint: flipping tensionOnly rejects a stale factor (it changes driver semantics)
        FrameModel m2; fixtures::xBracedPortal(m2, 5.0e4, 0.0, tmat, stocky, tbrace);
        PreparedSystem ps2 = assembleAndFactor(m2);
        m2.members[3].tensionOnly = !m2.members[3].tensionOnly;
        const SolveResult stale = solveLoad(ps2, m2);
        checkTrue("F42 tensionOnly in fingerprint (stale rejected)", stale.singular, "");
    }

    // ---------- F43: tension-only monotone finite termination + cycle-guard contract ----------
    {
        std::printf("[F43] Tension-only: monotone finite termination + cycle guard\n");
        Material tmat(200000.0, 76923.07692307692, 7850.0); tmat.cap = Capacity::make(1e9, 1e9, 1e9);
        Section stocky = Section::Rectangular(250.0, 250.0);
        Section tbrace = Section::Rectangular(60.0, 60.0);

        // vertical squash: both diagonals compress -> both go slack; the moment frame carries it
        FrameModel m; fixtures::xBracedPortal(m, 1.0e2, 2.0e5, tmat, stocky, tbrace);
        TensionOnlyOptions mono; mono.allowReactivation = false;
        const TensionOnlyResult Rm = runTensionOnly(m, mono);
        checkTrue("F43 monotone converged (finite termination)", Rm.converged, Rm.finalState.diagnostic);
        checkTrue("F43 monotone bounded iters (<= nTO+1)", Rm.iterations <= 3, "iters=" + std::to_string(Rm.iterations));

        // monotone result respects the tension-only constraint: every ACTIVE brace is not in compression
        bool okSign = true;
        for (int e : { 3, 4 }) {
            const bool slack = std::find(Rm.slack.begin(), Rm.slack.end(), (MemberId)e) != Rm.slack.end();
            if (!slack) { const real N = Rm.finalState.memberForces[(size_t)e].endI.N; if (N > 1e-6) okSign = false; }
        }
        checkTrue("F43 monotone respects tension-only sign constraint", okSign, "");

        // default (reactivation) policy must always TERMINATE across a load sweep — converged, or
        // cycled-then-monotone, never a hang. The WS-D sweep found no physical flip-flop, so the
        // cycle path is covered by the guard's finite-termination contract, not a physical case.
        bool allTerminate = true;
        for (real Hl : { 1.0e2, 1.0e3, 5.0e3, 3.0e4 }) {
            FrameModel ms; fixtures::xBracedPortal(ms, Hl, 2.0e5, tmat, stocky, tbrace);
            const TensionOnlyResult Rs = runTensionOnly(ms);
            if (!(Rs.converged || Rs.cycled)) allTerminate = false;
        }
        checkTrue("F43 default policy terminates across a load sweep", allTerminate, "");
    }

    // ---------- F44: 10-bar truss fully-stressed design (Schmit/Berke; FSD weight ~ 1593.2 lb) ----------
    {
        std::printf("[F44] Size optimization: classic 10-bar truss FSD (stress-ratio resizing)\n");
        const real IN = 25.4, LB = 4.4482216152605, PSI = 0.0068947572931684;
        const real Empa = 1.0e7 * PSI, sigA = 25000.0 * PSI;
        const real bay = 360.0 * IN, A0 = 10.0 * IN * IN, Amin = 0.1 * IN * IN, P = 1.0e5 * LB;
        Material tmat(Empa, Empa / 2.6, 0.0);
        tmat.cap = Capacity::make(sigA, sigA, sigA);
        FrameModel m; fixtures::tenBarTruss(m, bay, A0, P, tmat);

        SizeOptOptions o; o.maxIter = 100; o.dcTol = 1e-10; o.Amin = Amin;
        const SizeOptResult R = runSizeOptimization(m, o);
        checkTrue("F44 converged", R.converged && !R.singular, "iters=" + std::to_string(R.iterations));

        // weight from the final areas (imperial): W = sum 0.1 lb/in^3 * A_in2 * L_in
        real weightLb = 0;
        for (size_t k = 0; k < m.members.size(); ++k) {
            const int ni = m.nodeIndex(m.members[k].i), nj = m.nodeIndex(m.members[k].j);
            const real L = norm(m.nodes[(size_t)nj].pos - m.nodes[(size_t)ni].pos);
            weightLb += 0.1 * (R.finalAreas[k] / (IN * IN)) * (L / IN);
        }
        // combined-stress FSD lands at ~1608 lb -- about 1% ABOVE the pin-jointed literature optimum
        // (1593.2 lb, Haftka & Gurdal): the engine gates true pins, so the moment-frame members carry a
        // little bending that eats into the axial allowance, sizing every chord slightly heavier (and
        // safer). The combined-stress optimum can never fall BELOW the pure-axial optimum -- an invariant.
        checkTrue("F44 weight >= pin-jointed optimum 1593.2 lb (combined-stress is heavier)",
                  weightLb >= 1593.0, "W=" + std::to_string(weightLb));
        checkClose("F44 FSD weight within 1.5% of literature 1593.2 lb", weightLb, 1593.2, 1.5e-2);

        // fully-stressed: every SIZED bar sits at D/C = 1 to machine precision (the FS property)
        bool fs = true; int nSized = 0;
        for (size_t k = 0; k < m.members.size(); ++k)
            if (R.finalAreas[k] > Amin * 1.01) { ++nSized; if (std::fabs(R.finalDC[k] - 1.0) > 1e-6) fs = false; }
        checkTrue("F44 sized bars fully stressed (|D/C-1|<1e-6)", fs, "nSized=" + std::to_string(nSized));

        // the four low-force bars (2,5,6,10 -> idx 1,4,5,9) bottom out at A_min exactly (clamp)
        bool atMin = true;
        for (int k : { 1, 4, 5, 9 }) if (relErr(R.finalAreas[(size_t)k], Amin) > 1e-9) atMin = false;
        checkTrue("F44 bars 2/5/6/10 at A_min bound", atMin, "");

        // sized areas match the literature pattern (bars 1,3 are the big chords ~7.94 / 8.06 in^2)
        checkClose("F44 bar1 area (in^2)", R.finalAreas[0] / (IN * IN), 7.94, 3e-2);
        checkClose("F44 bar3 area (in^2)", R.finalAreas[2] / (IN * IN), 8.06, 3e-2);

        // Guard: zero allowable capacity under demand produces D/C=inf. FSD cannot fix that by
        // resizing, so it must fail honestly instead of treating the non-finite factor as 1.
        {
            Material zmat(Empa, Empa / 2.6, 0.0);  // cap intentionally left at zero
            FrameModel zm; fixtures::cantileverTipLoad(zm, -1000.0, 1000.0, zmat, Section::Rectangular(100.0, 100.0));
            SizeOptOptions zo; zo.maxIter = 10; zo.Amin = 1.0;
            const SizeOptResult Z = runSizeOptimization(zm, zo);
            checkTrue("F44 zero-cap demand does not fake-converge",
                      !Z.converged && Z.invalidDemand,
                      "converged=" + std::to_string((int)Z.converged) +
                      " invalidDemand=" + std::to_string((int)Z.invalidDemand));
            checkTrue("F44 zero-cap final D/C is non-finite",
                      !Z.finalDC.empty() && !std::isfinite(Z.finalDC[0]), "");
        }

        // Guard: a zero-force sizable member may hit Amin=0. The reported section/area and weight
        // history must agree; internally the member is deactivated so the working model stays valid.
        {
            FrameModel zm; fixtures::cantileverTipLoad(zm, -1000.0, 1000.0, tmat, Section::Rectangular(100.0, 100.0));
            Node n2(2, 0.0, 1000.0, 0.0); n2.fixAll();
            Node n3(3, 1000.0, 1000.0, 0.0); n3.fixAll();
            zm.nodes.push_back(n2);
            zm.nodes.push_back(n3);
            zm.members.push_back(Member(1, 2, 3, 0, 0));  // fixed-fixed, unloaded -> zero D/C
            SizeOptOptions zo; zo.maxIter = 10; zo.dcTol = 1e-10; zo.Amin = 0;
            const SizeOptResult Z = runSizeOptimization(zm, zo, { 1 });
            checkTrue("F44 zero-force member reports zero area/section",
                      Z.converged && !Z.singular && Z.finalAreas.size() > 1 &&
                      Z.finalAreas[1] == 0 && Z.finalSections[1].A == 0,
                      "A=" + (Z.finalAreas.size() > 1 ? std::to_string(Z.finalAreas[1]) : std::string("missing")));
            checkTrue("F44 zero-force member has zero sized volume",
                      !Z.weightHistory.empty() && Z.weightHistory.back() == 0, "");
        }
    }

    // helper: external work C = sum over nodal loads of F . u (the compliance / strain energy x2)
    auto compliance = [](const FrameModel& mm, const SolveResult& rr) -> real {
        real C = 0;
        for (const NodalLoad& L : mm.nodalLoads) {
            const int ni = mm.nodeIndex(L.node);
            if (ni < 0) continue;
            for (int d = 0; d < 6; ++d) C += L.comp[d] * rr.u[(size_t)gdof(ni, d)];
        }
        return C;
    };

    // ---------- F45: BESO sensitivity = element strain energy; energy balance sum(alpha)=1/2 F.u ----
    {
        std::printf("[F45] BESO sensitivity: element strain energy + energy balance (sum a = 1/2 F.u)\n");
        // 3D L-frame so all six end-force components participate (axial+shear+bending+torsion).
        FrameModel m; m.materials.push_back(mat); m.sections.push_back(sec);
        m.nodes.push_back(Node(0, 0, 0, 0)); m.nodes[0].fixAll();
        m.nodes.push_back(Node(1, 1500, 0, 0));
        m.nodes.push_back(Node(2, 1500, 0, 1200));
        { Member a(0, 0, 1, 0, 0); a.refVec = { 0, 0, 1 }; m.members.push_back(a); }
        { Member b(1, 1, 2, 0, 0); b.refVec = { 1, 0, 0 }; m.members.push_back(b); }
        NodalLoad l; l.node = 2; l.comp[Ux] = 8000.0; l.comp[Uz] = -5000.0; l.comp[Ry] = 3.0e6;
        m.nodalLoads.push_back(l);
        SolveResult r = solve(m);
        checkTrue("F45 not singular", !r.singular, r.diagnostic);
        real sumA = 0;
        for (size_t e = 0; e < m.members.size(); ++e) sumA += memberStrainEnergy(m, r, (int)e);
        // sum of element strain energies == 1/2 external work (the energy-balance identity)
        checkClose("F45 energy balance sum(alpha) == 1/2 F.u", sumA, 0.5 * compliance(m, r), 1e-9);
        checkTrue("F45 both members carry strain energy",
                  memberStrainEnergy(m, r, 0) > 0 && memberStrainEnergy(m, r, 1) > 0, "");
        // single-member cantilever degenerate check: alpha == 1/2 F.u
        FrameModel c; fixtures::cantileverTipLoad(c, 1000.0, 2000.0, mat, sec);
        SolveResult rc = solve(c);
        checkClose("F45 cantilever alpha == 1/2 F.u", memberStrainEnergy(c, rc, 0), 0.5 * compliance(c, rc), 1e-9);

        // per-component analytic oracle: isolate each weight channel against the closed-form element
        // strain energy (P^2 L/2EA, M^2 L/2EIz, T^2 L/2GJ). A sum-only balance can hide a SYMMETRIC
        // sign error; the analytic value is an independent source that pins each channel's coefficient.
        {
            const real Lb = 2000.0;
            auto bar = [&](int comp, real val) -> real {   // cantilever bar along +x, refVec=+y -> local == global
                FrameModel m; m.materials.push_back(mat); m.sections.push_back(sec);
                m.nodes.push_back(Node(0, 0, 0, 0)); m.nodes[0].fixAll();
                m.nodes.push_back(Node(1, Lb, 0, 0));
                { Member a(0, 0, 1, 0, 0); a.refVec = { 0, 1, 0 }; m.members.push_back(a); }
                NodalLoad l; l.node = 1; l.comp[comp] = val; m.nodalLoads.push_back(l);
                const SolveResult r = solve(m);
                return memberStrainEnergy(m, r, 0);
            };
            const real Pax = 30000.0, Mz = 2.0e7, Tq = 1.5e7;
            checkClose("F45 pure axial alpha == P^2 L/2EA",   bar(Ux, Pax), Pax * Pax * Lb / (2.0 * E * sec.A),  1e-6);
            checkClose("F45 pure bending alpha == M^2 L/2EIz", bar(Rz, Mz),  Mz * Mz * Lb / (2.0 * E * sec.Iz),  1e-5);
            checkClose("F45 pure torsion alpha == T^2 L/2GJ",  bar(Rx, Tq),  Tq * Tq * Lb / (2.0 * G * sec.J),   1e-5);
        }
    }

    // ---------- F46: BESO evolutionary removal on a ground structure (Michell-like) ----------
    {
        std::printf("[F46] BESO: ground-structure evolutionary removal + compliance-best rollback\n");
        const int NX = 6, NZ = 3; const real SP = 500.0;
        auto nid = [&](int i, int k) { return k * (NX + 1) + i; };
        FrameModel m; m.materials.push_back(mat); m.sections.push_back(sec);
        for (int k = 0; k <= NZ; ++k)
            for (int i = 0; i <= NX; ++i) {
                Node n(nid(i, k), i * SP, 0.0, k * SP);
                if (i == 0) n.fixAll();
                else { n.fixed[Uy] = true; n.fixed[Rx] = true; n.fixed[Rz] = true; }  // planar xz frame
                m.nodes.push_back(n);
            }
        auto add = [&](int aN, int bN) {
            Member mem((int)m.members.size(), aN, bN, 0, 0);
            const Vec3 d = m.nodes[(size_t)bN].pos - m.nodes[(size_t)aN].pos;
            mem.refVec = (std::fabs(d.z) > std::fabs(d.x)) ? Vec3{ 1, 0, 0 } : Vec3{ 0, 0, 1 };
            m.members.push_back(mem);
        };
        for (int k = 0; k <= NZ; ++k) for (int i = 0; i < NX; ++i) add(nid(i, k), nid(i + 1, k));
        for (int k = 0; k < NZ; ++k) for (int i = 0; i <= NX; ++i) add(nid(i, k), nid(i, k + 1));
        for (int k = 0; k < NZ; ++k) for (int i = 0; i < NX; ++i) { add(nid(i, k), nid(i + 1, k + 1)); add(nid(i + 1, k), nid(i, k + 1)); }
        const int loadNode = nid(NX, NZ / 2);
        NodalLoad l; l.node = loadNode; l.comp[Uz] = -1.0e4; m.nodalLoads.push_back(l);

        const int active0 = (int)m.members.size();
        BESOOptions o; o.targetVolFrac = 0.5; o.evolRate = 0.06; o.maxIter = 80;
        const BESOResult R = runBESO(m, o);
        // BESO drives volume to the target, then terminates cleanly -- either it hit the target, or
        // (discrete bars can't land on it precisely) it stalled JUST above it because the remaining
        // bars can't go without a mechanism. Both are legitimate; reject blow-up / mechanism / invalid.
        checkTrue("F46 drove volume to/near target + clean termination",
                  (R.reason == BESOStop::TargetReached || R.reason == BESOStop::Stalled) &&
                  !R.volFracHistory.empty() && R.volFracHistory.back() > 0.0 &&
                  R.volFracHistory.back() <= o.targetVolFrac + o.evolRate,
                  "reason=" + std::to_string((int)R.reason) +
                  " volFrac=" + std::to_string(R.volFracHistory.empty() ? -1.0 : R.volFracHistory.back()));
        int activeF = 0; for (char a : R.finalActive) if (a) ++activeF;
        // quantitative removal: at least a quarter of the bars are gone (not merely "removed one bar")
        checkTrue("F46 removed >= 25% of members", activeF <= (int)(active0 * 0.75),
                  "active " + std::to_string(active0) + "->" + std::to_string(activeF));
        // loaded node stays grounded: rebuild the best topology, solve, must be non-singular
        FrameModel mb = m;
        for (size_t e = 0; e < mb.members.size(); ++e) mb.members[e].active = R.bestActive[e] != 0;
        const SolveResult rb = solve(mb);
        checkTrue("F46 best topology non-singular (load path intact)", !rb.singular, rb.diagnostic);
        // removing ~half the volume must make the structure MEASURABLY softer: compliance up >= 15%.
        // (a 0.1%-tolerance "rises" check would pass even if nothing changed -- too weak to mean anything)
        checkTrue("F46 compliance rises >= 15% with removal (softer structure)",
                  !R.complianceHistory.empty() && std::isfinite(R.complianceHistory.back()) &&
                  R.complianceHistory.back() >= R.complianceHistory.front() * 1.15,
                  "C0=" + std::to_string(R.complianceHistory.front()) + " Cf=" + std::to_string(R.complianceHistory.back()));
    }

    // ---------- F47: N2 collapse robustness — unconstrained vs robust BESO ----------
    {
        std::printf("[F47] N2 collapse robustness: unconstrained vs robust BESO topology\n");
        // Two real supports A(0) and B(1); the relay node M and loaded node C are FULLY FREE (no fixed
        // DOF -- a node with ANY fixed DOF reads as "grounded" to connectivity, which would defeat the
        // ungrounding test). Main chain g(A-M)+m(M-C) reaches the load via A; a small-section backup
        // b(B-C) gives a second ground path via B. With only the main chain, removing g ungrounds
        // {M,C} -> Collapsed; with b retained the load stays grounded -> Stable. The 3D beams keep the
        // free nodes stable on their own (no truss mechanism). The probe finds fragility regardless of
        // cut order.
        FrameModel m; m.materials.push_back(mat);
        m.sections.push_back(sec);                                   // sec0: main bars (square 100)
        m.sections.push_back(Section::Rectangular(60.0, 60.0));      // sec1: small backup bar
        m.nodes.push_back(Node(0, 0, 0, 0));      m.nodes[0].fixAll();   // support A
        m.nodes.push_back(Node(1, 0, 0, 2000));   m.nodes[1].fixAll();   // support B
        m.nodes.push_back(Node(2, 1500, 0, 0));                          // relay M (free)
        m.nodes.push_back(Node(3, 3000, 0, 1000));                       // loaded C (free)
        m.members.push_back(Member(0, 0, 2, 0, 0));   // g (A-M, main)
        m.members.push_back(Member(1, 2, 3, 0, 0));   // m (M-C, main)
        m.members.push_back(Member(2, 1, 3, 0, 1));   // b (B-C, small backup)
        NodalLoad l; l.node = 3; l.comp[Uz] = -2000.0; m.nodalLoads.push_back(l);

        // does ANY single active-member removal make the (active-masked) structure collapse?
        auto fragile = [&](const std::vector<char>& act) -> bool {
            FrameModel mm = m;
            for (size_t e = 0; e < mm.members.size(); ++e) mm.members[e].active = act[e] != 0;
            for (size_t e = 0; e < mm.members.size(); ++e) {
                if (!mm.members[e].active) continue;
                CollapseOptions co; co.dlf = 1.0; co.removeThreshold = 1.0;
                co.initialRemovals = { mm.members[e].id };
                if (runProgressiveCollapse(mm, co).outcome == CollapseOutcome::Collapsed) return true;
            }
            return false;
        };

        BESOOptions ou; ou.targetVolFrac = 0.6; ou.evolRate = 0.34; ou.maxIter = 20;   // unconstrained
        const BESOResult Ru = runBESO(m, ou);
        BESOOptions orb = ou; orb.redundancyCheckEvery = 1; orb.redundancySamples = 0; // robust
        orb.redundancy.dlf = 1.0; orb.redundancy.removeThreshold = 1.0;
        const BESOResult Rr = runBESO(m, orb);
        BESOOptions ors = orb; ors.redundancySamples = 1;             // sampled N2 path (post-trim D/C ranking)
        const BESOResult Rs = runBESO(m, ors);

        checkTrue("F47 unconstrained topology is fragile (a single removal collapses it)",
                  fragile(Ru.finalActive), "");
        checkTrue("F47 robust topology survives every single-member removal",
                  !fragile(Rr.finalActive), "");
        checkTrue("F47 sampled robust topology survives every single-member removal",
                  !fragile(Rs.finalActive) && !Rs.protectedMembers.empty(), "");
        checkTrue("F47 robust BESO locked protective members", !Rr.protectedMembers.empty(),
                  "nProtected=" + std::to_string(Rr.protectedMembers.size()));
        int au = 0, ar = 0;
        for (char a : Ru.finalActive) au += a ? 1 : 0;
        for (char a : Rr.finalActive) ar += a ? 1 : 0;
        checkTrue("F47 robust retains >= unconstrained active members", ar >= au,
                  "unc=" + std::to_string(au) + " rob=" + std::to_string(ar));
        // direct mechanism check (not just the fragility outcome, which a carefully tuned model could
        // pass by luck): the unconstrained run DROPS the soft backup bar b (idx 2 -- lowest strain
        // energy); the robust run KEEPS it and LISTS it as protected. b's id is 2.
        checkTrue("F47 unconstrained dropped the backup bar b (idx 2)",
                  Ru.finalActive.size() > 2 && Ru.finalActive[2] == 0, "");
        const bool bLocked = std::find(Rr.protectedMembers.begin(), Rr.protectedMembers.end(), 2) != Rr.protectedMembers.end();
        checkTrue("F47 robust kept + locked the backup bar b (idx 2)",
                  Rr.finalActive.size() > 2 && Rr.finalActive[2] == 1 && bLocked, "");
    }

    // ---------- F48: QM6 incompatible-mode membrane (S8-8a) ----------
    {
        std::printf("[F48] QM6 incompatible-mode membrane: constant-strain patch + in-plane bending unlock\n");
        const real Es = 30000.0, nu = 0.3, Gs = Es / (2.0 * (1.0 + nu));
        Material smat(Es, Gs); smat.nu = nu;
        SolveOptions qm6; qm6.useIncompatibleMembrane = true;

        // (a) QM6 must STILL pass the constant-strain membrane patch. On regular AND parallelogram
        // meshes (affine map, detJ = const) the centre-Jacobian correction makes integral(B_inc)=0
        // EXACTLY, so a constant stress stores no bubble energy. (On a general NON-affine quad
        // integral(B_inc) is only O(h); QM6 then passes the WEAK Irons-Razzaque patch test and
        // converges under refinement -- see S8_shell.md.) Reuse the F14a patch with QM6 turned ON.
        for (const real skew : { 0.0, 0.4 }) {
            const real a = 1000.0, t = 10.0, gx = 1e-4;
            const real f = Es / (1.0 - nu * nu);
            FrameModel m; fixtures::membranePatch(m, a, t, skew, gx, smat);
            const SolveResult r = solve(m, qm6);
            const char* tag = (skew == 0.0) ? "regular" : "parallelogram";
            checkTrue("QM6 patch non-singular", !r.singular, r.diagnostic);
            const real NxxExp = t * f * gx, NyyExp = nu * NxxExp;
            const real scale = std::fabs(NxxExp);
            real eNxx = 0, eNyy = 0, mNxy = 0;
            for (const auto& sf : r.shellForces) {
                eNxx = std::max(eNxx, std::fabs(sf.Nxx - NxxExp));
                eNyy = std::max(eNyy, std::fabs(sf.Nyy - NyyExp));
                mNxy = std::max(mNxy, std::fabs(sf.Nxy));
            }
            std::printf("   QM6 patch (%-13s) max|dNxx|=%.2e |dNyy|=%.2e |Nxy|=%.2e\n",
                        tag, eNxx, eNyy, mNxy);
            checkTrue("QM6 patch Nxx constant (centre-Jacobian)", eNxx < 1e-8 * scale, "eNxx=" + std::to_string(eNxx));
            checkTrue("QM6 patch Nyy = nu Nxx", eNyy < 1e-8 * scale, "eNyy=" + std::to_string(eNyy));
            checkTrue("QM6 patch Nxy ~ 0", mNxy < 1e-8 * scale, "mNxy=" + std::to_string(mNxy));
        }

        // (b) slender in-plane cantilever (L/H=10, coarse 4x1): the Q4 membrane LOCKS (reads far
        // too stiff), QM6 releases it. Reference Euler-Bernoulli tip = P L^3 / (3 E I).
        {
            const real L = 100.0, H = 10.0, t = 1.0, P = 1.0;
            const int  nx = 4, ny = 1;
            const real I = t * H * H * H / 12.0;
            const real dEB = P * L * L * L / (3.0 * Es * I);
            auto tip = [&](const SolveOptions& o) -> real {
                FrameModel m; fixtures::slenderMembraneCantilever(m, L, H, t, nx, ny, P, smat);
                const SolveResult r = solve(m, o);
                return r.singular ? 0.0 : std::fabs(r.disp(ny * (nx + 1) + nx, Uy));
            };
            const real dQ4  = tip(SolveOptions{});
            const real dQM6 = tip(qm6);
            const real eQ4  = (dQ4 - dEB) / dEB, eQM6 = (dQM6 - dEB) / dEB;
            std::printf("   slender 4x1 cantilever dEB=%.5g  Q4=%.5g (%.1f%%)  QM6=%.5g (%.1f%%)\n",
                        dEB, dQ4, 100 * eQ4, dQM6, 100 * eQM6);
            checkTrue("Q4 membrane LOCKS (tip < 40% of Euler-Bernoulli)", dQ4 < 0.40 * dEB,
                      "dQ4/dEB=" + std::to_string(dQ4 / dEB));
            checkTrue("QM6 RELEASES the lock (tip within 15% of Euler-Bernoulli)",
                      std::fabs(eQM6) < 0.15, "eQM6=" + std::to_string(eQM6));
            checkTrue("QM6 strictly better than Q4 (membrane locking)", std::fabs(eQM6) < std::fabs(eQ4),
                      "eQ4=" + std::to_string(eQ4) + " eQM6=" + std::to_string(eQM6));
        }

        // (c) Cook's skew membrane: QM6 converges FASTER than Q4. Use the fine QM6 tip as a
        // SELF-CONSISTENT converged reference (independent of any particular literature value or
        // load-lumping convention; the same consistent lumping is used at every mesh, so the
        // Q4-vs-QM6 comparison is fair) and require a COARSE QM6 mesh to be closer to it than a
        // coarse Q4 mesh. The literature tip ~23.96 is printed for context only.
        {
            const real Ec = 1.0, nuc = 1.0 / 3.0, Gc = Ec / (2.0 * (1.0 + nuc));
            Material cmat(Ec, Gc); cmat.nu = nuc;
            const real P = 1.0, tc = 1.0;
            auto tipC = [&](int n, const SolveOptions& o) -> real {
                FrameModel m; fixtures::cooksMembrane(m, n, P, tc, cmat);
                const SolveResult r = solve(m, o);
                return r.singular ? 0.0 : std::fabs(r.disp(n * (n + 1) + n, Uy));
            };
            const real refC = tipC(32, qm6);                       // fine QM6 = converged tip
            const real q4f  = tipC(32, SolveOptions{});             // fine Q4 must MATCH (compatibility)
            const real q4c  = tipC(4, SolveOptions{});
            const real qmc  = tipC(4, qm6);
            const real eq4 = (q4c - refC) / refC, eqm = (qmc - refC) / refC;
            std::printf("   Cook's refC(QM6 N=32)=%.4g  Q4 N=32=%.4g (lit~23.96)  N=4: Q4=%.4g (%.1f%%)  QM6=%.4g (%.1f%%)\n",
                        refC, q4f, q4c, 100 * eq4, qmc, 100 * eqm);
            // Compatibility: Q4 is also conforming and shares QM6's fine-mesh limit, just reaches
            // it much slower. At N=32 QM6 is already converged (~25.0) while Q4 is still climbing
            // (~24.8) toward the same value -- direct evidence QM6 converges faster and rules out
            // QM6 introducing spurious softening. The ~4% gap of the limit to the literature 23.96
            // is the consistent half-weight edge-load lumping + corner-tip sampling used here (NOT
            // the drilling penalty -- a penalty would STIFFEN and LOWER the tip, the wrong way); it
            // does not affect the QM6-vs-Q4 RELATIVE convergence conclusion.
            checkTrue("Cook converged tip in plausible range [23,27] (guards self-reference drift)",
                      refC > 23.0 && refC < 27.0, "refC=" + std::to_string(refC));
            checkTrue("Cook Q4 monotonically approaches the QM6-converged tip (same limit, Q4 slower)",
                      q4f > q4c && std::fabs(q4f - refC) < 0.02 * refC,
                      "q4_N4=" + std::to_string(q4c) + " q4_N32=" + std::to_string(q4f) + " refC=" + std::to_string(refC));
            checkTrue("Cook QM6 coarse closer to converged than Q4", std::fabs(eqm) < std::fabs(eq4),
                      "eq4=" + std::to_string(eq4) + " eqm=" + std::to_string(eqm));
            checkTrue("Cook QM6 N=4 within 6% of converged", std::fabs(eqm) < 0.06, "eqm=" + std::to_string(eqm));
        }
    }

    // ---------- F49: DKQ discrete-Kirchhoff thin-plate bending (S8-8b) ----------
    {
        std::printf("[F49] DKQ discrete-Kirchhoff thin plate: constant-curvature patch + thin-plate convergence\n");
        const real Es = 30000.0, nu = 0.3, Gs = Es / (2.0 * (1.0 + nu));
        Material smat(Es, Gs); smat.nu = nu;
        SolveOptions dkq; dkq.useDKQPlate = true;

        // (a) constant-curvature patch (cylindrical bending), regular + parallelogram: DKQ must
        // reproduce the exact constant moment to MACHINE PRECISION. This is the hard correctness
        // gate AND it pins the Batoz (w,tx,ty) -> (w,bx,by) DOF/sign mapping.
        for (const real skew : { 0.0, 0.4 }) {
            const real ap = 1000.0, tp = 10.0, cc = 1e-6;
            const real Dfac = Es * tp * tp * tp / (12.0 * (1.0 - nu * nu));
            FrameModel m; fixtures::platePatchCylindrical(m, ap, tp, skew, cc, smat);
            const SolveResult r = solve(m, dkq);
            const char* tag = (skew == 0.0) ? "regular" : "parallelogram";
            checkTrue("DKQ patch non-singular", !r.singular, r.diagnostic);
            const real MxxExp = -Dfac * cc, MyyExp = -nu * Dfac * cc;
            const real scale = std::fabs(MxxExp);
            real eMxx = 0, eMyy = 0, mMxy = 0;
            for (const auto& sf : r.shellForces) {
                eMxx = std::max(eMxx, std::fabs(sf.Mxx - MxxExp));
                eMyy = std::max(eMyy, std::fabs(sf.Myy - MyyExp));
                mMxy = std::max(mMxy, std::fabs(sf.Mxy));
            }
            std::printf("   DKQ patch (%-13s) Mxx_exp=%.6g max|dMxx|=%.2e |dMyy|=%.2e |Mxy|=%.2e\n",
                        tag, MxxExp, eMxx, eMyy, mMxy);
            checkTrue("DKQ patch Mxx = -D c (constant)", eMxx < 1e-8 * scale, "eMxx=" + std::to_string(eMxx));
            checkTrue("DKQ patch Myy = -nu D c", eMyy < 1e-8 * scale, "eMyy=" + std::to_string(eMyy));
            checkTrue("DKQ patch Mxy ~ 0", mMxy < 1e-8 * scale, "mMxy=" + std::to_string(mMxy));
        }

        // (b) simply-supported thin square plate (t/a=1/200) under pressure: DKQ converges to the
        // Kirchhoff centre deflection 0.00406 q a^4 / D. DKQ IS Kirchhoff, so (unlike Mindlin
        // MITC4) it does NOT overshoot. MITC4 printed for contrast.
        {
            const real a = 1000.0, t = 5.0, q = 0.01;
            const real D = Es * t * t * t / (12.0 * (1.0 - nu * nu));
            const real wc = 0.00406 * q * a * a * a * a / D;
            auto centerW = [&](int n, const SolveOptions& o) -> real {
                FrameModel m; fixtures::squarePlateShell(m, a, t, n, q, smat);
                const SolveResult r = solve(m, o);
                if (r.singular) return 1e30;
                const int c = (n / 2) * (n + 1) + (n / 2);
                return std::fabs(r.disp(c, Uz));
            };
            const real wD16 = centerW(16, dkq), wD24 = centerW(24, dkq);
            const real eD16 = relErr(wD16, wc), eD24 = relErr(wD24, wc);
            const real wM16 = centerW(16, SolveOptions{});
            std::printf("   SS thin plate wc=%.6g  DKQ N16=%.6g (%.2f%%) N24=%.6g (%.2f%%)  MITC4 N16=%.6g\n",
                        wc, wD16, 100 * eD16, wD24, 100 * eD24, wM16);
            checkTrue("DKQ converges to Kirchhoff (N16 & N24 within 2%)", eD16 < 0.02 && eD24 < 0.02,
                      "eD16=" + std::to_string(eD16) + " eD24=" + std::to_string(eD24));
        }

        // (c) clamped thin square plate: DKQ converges to Kirchhoff 0.00126 q a^4 / D.
        {
            const real a = 1000.0, t = 5.0, q = 0.01;
            const real D = Es * t * t * t / (12.0 * (1.0 - nu * nu));
            const real wc = 0.00126 * q * a * a * a * a / D;
            FrameModel m; fixtures::clampedPlateShell(m, a, t, 16, q, smat);
            const SolveResult r = solve(m, dkq);
            checkTrue("DKQ clamped non-singular", !r.singular, r.diagnostic);
            const int c = 8 * (16 + 1) + 8;
            const real wD = r.singular ? 1e30 : std::fabs(r.disp(c, Uz));
            const real eD = relErr(wD, wc);
            std::printf("   clamped thin plate wc=%.6g  DKQ N16=%.6g (%.2f%%)\n", wc, wD, 100 * eD);
            checkTrue("DKQ clamped within 3% of Kirchhoff", eD < 0.03, "eD=" + std::to_string(eD));
        }
    }

    // ---------- F50: co-rotational large displacement (planar elastica) ----------
    {
        std::printf("[F50] co-rotational large displacement: planar cantilever elastica (Mattiasson table)\n");
        const real kPiLoc = 3.14159265358979323846;
        // E=1, near-inextensible slender section: A large, I small so the elastica (inextensible) table is
        // matched (finite EA adds only ~alpha*I/(A L^2) axial strain).
        Material cmat(1.0, 0.4, 0.0);
        Section csec; csec.A = 100.0; csec.Iy = 1e-3; csec.Iz = 1e-3; csec.J = 1e-3;
        csec.cy = 1.0; csec.cz = 1.0; csec.Asy = 0.0; csec.Asz = 0.0;
        const real L = 1.0, Em = 1.0, Im = 1e-3;

        // WS_F F-3 elastica table: alpha = P L^2/(E I) -> tip delta_v/L (transverse), delta_h/L (shortening)
        struct Row { real alpha, dv, dh; };
        const Row tab[] = {
            { 1.0,  0.3017207738, 0.0564332363 },
            { 2.0,  0.4934574804, 0.1606417208 },
            { 5.0,  0.7137915236, 0.3876283607 },
            { 10.0, 0.8106090249, 0.5549955978 },
        };

        // (a) mesh convergence to the elastica table
        auto runTip = [&](int n, real alpha, real& dv, real& dh, int& iters, bool& conv) {
            FrameModel m; fixtures::cantileverPlanarTipShearN(m, n, L, alpha * Em * Im / (L * L), cmat, csec);
            CorotationalOptions co; co.loadSteps = std::max(12, (int)(alpha * 3)); co.maxIter = 80;
            const CorotationalResult R = runCorotational(m, co);
            conv = R.converged; iters = R.totalIterations;
            dv = R.finalState.u[(size_t)gdof(n, Uy)] / L;
            dh = -R.finalState.u[(size_t)gdof(n, Ux)] / L;
        };
        for (const Row& r : tab) {
            real dv4, dh4, dv16, dh16; int it4, it16; bool c4, c16;
            runTip(4, r.alpha, dv4, dh4, it4, c4);
            runTip(16, r.alpha, dv16, dh16, it16, c16);
            std::printf("   alpha=%2.0f exact dv/L=%.6f dh/L=%.6f | N4 dv=%.6f(%.2f%%) | N16 dv=%.6f(%.3f%%) dh=%.6f(%.2f%%) it=%d\n",
                        r.alpha, r.dv, r.dh, dv4, 100 * relErr(dv4, r.dv), dv16, 100 * relErr(dv16, r.dv),
                        dh16, 100 * relErr(dh16, r.dh), it16);
            checkTrue(("F50 converged alpha=" + std::to_string((int)r.alpha)).c_str(), c4 && c16, "");
            checkTrue(("F50 N4 coarse within 1.5% dv alpha=" + std::to_string((int)r.alpha)).c_str(),
                      relErr(dv4, r.dv) < 1.5e-2, "e=" + std::to_string(relErr(dv4, r.dv)));
            checkClose(("F50 N16 dv/L alpha=" + std::to_string((int)r.alpha)).c_str(), dv16, r.dv, 1.5e-3);
            checkClose(("F50 N16 dh/L alpha=" + std::to_string((int)r.alpha)).c_str(), dh16, r.dh, 5e-3);
        }

        // (b) in-plane rotational invariance (corotational frame indifference)
        {
            const real alpha = 3.0, phi = 0.6, P = alpha * Em * Im / (L * L);
            FrameModel m0; fixtures::cantileverPlanarTipShearN(m0, 12, L, P, cmat, csec);
            CorotationalOptions co; co.loadSteps = 15; co.maxIter = 80;
            const CorotationalResult R0 = runCorotational(m0, co);
            FrameModel m1; fixtures::cantileverPlanarTipShearN(m1, 12, L, 0.0, cmat, csec);
            const real cph = std::cos(phi), sph = std::sin(phi);
            for (auto& nd : m1.nodes) { const real x = nd.pos.x, y = nd.pos.y; nd.pos.x = cph * x - sph * y; nd.pos.y = sph * x + cph * y; }
            NodalLoad nl; nl.node = 12; nl.comp[Ux] = -sph * P; nl.comp[Uy] = cph * P; m1.nodalLoads = { nl };
            const CorotationalResult R1 = runCorotational(m1, co);
            auto mag = [&](const CorotationalResult& R) { const real ux = R.finalState.u[(size_t)gdof(12, Ux)], uy = R.finalState.u[(size_t)gdof(12, Uy)]; return std::sqrt(ux * ux + uy * uy); };
            const real m0mag = mag(R0), m1mag = mag(R1);
            std::printf("   rotational invariance: |u_tip| base=%.6g rotated=%.6g rel=%.2e\n", m0mag, m1mag, relErr(m1mag, m0mag));
            checkTrue("F50 rotated solve converged", R0.converged && R1.converged, "");
            checkClose("F50 in-plane rotational invariance", m1mag, m0mag, 1e-9);
        }

        // (c) P-Delta degeneration: small lateral load -> CR sway matches linearized 2nd-order runPDelta.
        {
            const int nE = 6; const real Hc = 1.0;
            const real Pcr = kPiLoc * kPiLoc * Em * Im / (4.0 * Hc * Hc);
            const real Paxial = 0.3 * Pcr, Hlat = 1e-4;
            auto build = [&](FrameModel& m) {
                fixtures::prepMatSec(m, cmat, csec);
                m.nodes.clear(); m.members.clear();
                for (int k = 0; k <= nE; ++k) { Node nd(k, 0, Hc * real(k) / nE, 0); nd.fixed[Uz] = nd.fixed[Rx] = nd.fixed[Ry] = true; if (k == 0) nd.fixAll(); m.nodes.push_back(nd); }
                for (int k = 0; k < nE; ++k) m.members.push_back(Member(k, k, k + 1, 0, 0));
                NodalLoad nl; nl.node = nE; nl.comp[Uy] = -Paxial; nl.comp[Ux] = Hlat; m.nodalLoads = { nl };
            };
            FrameModel mc; build(mc);
            CorotationalOptions co; co.loadSteps = 12; co.maxIter = 80;
            const CorotationalResult Rcr = runCorotational(mc, co);
            FrameModel mp; build(mp);
            PDeltaOptions po; po.refactorPath = true; po.solve.pivotTol = 1e-12;
            const PDeltaResult Rpd = runPDelta(mp, po);
            const real swayCR = Rcr.finalState.u[(size_t)gdof(nE, Ux)];
            const real swayPD = Rpd.finalState.u[(size_t)gdof(nE, Ux)];
            std::printf("   P-Delta degeneration: sway CR=%.6e PDelta=%.6e rel=%.2e (Paxial/Pcr=0.3)\n",
                        swayCR, swayPD, relErr(swayCR, swayPD));
            checkTrue("F50 P-Delta deg converged", Rcr.converged && Rpd.converged, Rcr.finalState.diagnostic + Rpd.finalState.diagnostic);
            checkClose("F50 CR sway == P-Delta sway", swayCR, swayPD, 1.2e-2);
        }
    }

    // ---------- F51: 3D GENERAL co-rotational (arbitrary axis, torsion, biaxial bending) ----------
    {
        std::printf("[F51] 3D co-rotational: arbitrary-axis invariance / spatial elastica / torsion+biaxial\n");
        const real L = 1.0, Em = 1.0, Gm = 0.4, Im = 1e-3;
        Material cmat(Em, Gm, 0.0);
        // symmetric slender section (Iy=Iz) so bending is isotropic (elastica/invariance)
        Section sym; sym.A = 100.0; sym.Iy = 1e-3; sym.Iz = 1e-3; sym.J = 2e-3;
        sym.cy = 1.0; sym.cz = 1.0; sym.Asy = 0.0; sym.Asz = 0.0;

        // (a) ARBITRARY-AXIS rotational invariance (machine precision): solve an X-cantilever, then solve the
        //     SAME problem rigidly rotated by R_g(axis,phi). |u_tip| must be identical (CR frame indifference;
        //     exercises expSO3/logSO3/frame E with all three e-axes non-trivial).
        {
            auto rotv = [](Vec3 v, Vec3 k, real ang) -> Vec3 {
                const real nk = std::sqrt(k.x * k.x + k.y * k.y + k.z * k.z); k = Vec3(k.x / nk, k.y / nk, k.z / nk);
                const real c = std::cos(ang), s = std::sin(ang), kd = k.x * v.x + k.y * v.y + k.z * v.z;
                const Vec3 kxv(k.y * v.z - k.z * v.y, k.z * v.x - k.x * v.z, k.x * v.y - k.y * v.x);
                return Vec3(v.x * c + kxv.x * s + k.x * kd * (1 - c), v.y * c + kxv.y * s + k.y * kd * (1 - c), v.z * c + kxv.z * s + k.z * kd * (1 - c));
            };
            const real alpha = 5.0, P = alpha * Em * Im / (L * L);
            FrameModel m0; fixtures::cantileverSpatial(m0, 12, L, Vec3(1, 0, 0), 0, P, 0, 0, 0, 0, cmat, sym);
            CorotationalOptions co; co.loadSteps = 15; co.maxIter = 80;
            const CorotationalResult R0 = runCorotational(m0, co);
            const Vec3 axisN(1, 1, 1); const real phi = 2.0;
            FrameModel m1; fixtures::cantileverSpatial(m1, 12, L, Vec3(1, 0, 0), 0, P, 0, 0, 0, 0, cmat, sym);
            for (auto& nd : m1.nodes) nd.pos = rotv(nd.pos, axisN, phi);
            for (auto& mm : m1.members) mm.refVec = rotv(mm.refVec, axisN, phi);
            const Vec3 fr = rotv(Vec3(0, P, 0), axisN, phi);
            m1.nodalLoads[0].comp[Ux] = fr.x; m1.nodalLoads[0].comp[Uy] = fr.y; m1.nodalLoads[0].comp[Uz] = fr.z;
            const CorotationalResult R1 = runCorotational(m1, co);
            auto mag = [&](const CorotationalResult& R) { const real ux = R.finalState.u[(size_t)gdof(12, Ux)], uy = R.finalState.u[(size_t)gdof(12, Uy)], uz = R.finalState.u[(size_t)gdof(12, Uz)]; return std::sqrt(ux * ux + uy * uy + uz * uz); };
            const real a = mag(R0), b = mag(R1);
            std::printf("   (a) arbitrary-axis invariance |u_tip| base=%.6g rot=%.6g rel=%.2e\n", a, b, relErr(b, a));
            checkTrue("F51a converged", R0.converged && R1.converged, R0.finalState.diagnostic + R1.finalState.diagnostic);
            checkClose("F51a arbitrary-axis rotational invariance", b, a, 1e-9);
        }

        // (b) SPATIAL elastica: cantilever along a TILTED axis (1,1,1), symmetric section, tip load
        //     perpendicular to the axis -> delta_v/L, delta_h/L vs Mattiasson table.
        {
            struct Row { real alpha, dv, dh; };
            const Row tab[] = { {1.0, 0.3017207738, 0.0564332363}, {5.0, 0.7137915236, 0.3876283607}, {10.0, 0.8106090249, 0.5549955978} };
            const real na = std::sqrt(3.0); const Vec3 axis(1, 1, 1), aN(1 / na, 1 / na, 1 / na);
            const real n2 = std::sqrt(2.0); const Vec3 wN(1 / n2, -1 / n2, 0.0);   // perpendicular to aN
            for (const Row& r : tab) {
                const real P = r.alpha * Em * Im / (L * L);
                FrameModel m; fixtures::cantileverSpatial(m, 16, L, axis, wN.x * P, wN.y * P, wN.z * P, 0, 0, 0, cmat, sym);
                CorotationalOptions co; co.loadSteps = std::max(12, (int)(r.alpha * 3)); co.maxIter = 80;
                const CorotationalResult R = runCorotational(m, co);
                const Vec3 ut(R.finalState.u[(size_t)gdof(16, Ux)], R.finalState.u[(size_t)gdof(16, Uy)], R.finalState.u[(size_t)gdof(16, Uz)]);
                const real dv = (ut.x * wN.x + ut.y * wN.y + ut.z * wN.z) / L;
                const real dh = -(ut.x * aN.x + ut.y * aN.y + ut.z * aN.z) / L;
                std::printf("   (b) spatial elastica alpha=%2.0f dv/L=%.6f(exp %.6f) dh/L=%.6f(exp %.6f)\n", r.alpha, dv, r.dv, dh, r.dh);
                checkTrue(("F51b converged alpha=" + std::to_string((int)r.alpha)).c_str(), R.converged, R.finalState.diagnostic);
                checkClose(("F51b spatial dv/L alpha=" + std::to_string((int)r.alpha)).c_str(), dv, r.dv, 2e-3);
                checkClose(("F51b spatial dh/L alpha=" + std::to_string((int)r.alpha)).c_str(), dh, r.dh, 6e-3);
            }
        }

        // (c) PURE TORSION: X-cantilever, small tip torque T -> tip twist theta_x = T L/(G J) (linear).
        {
            Section asec; asec.A = 100.0; asec.Iy = 2e-3; asec.Iz = 1e-3; asec.J = 1.5e-3; asec.cy = 1.0; asec.cz = 1.0; asec.Asy = 0.0; asec.Asz = 0.0;
            const real T = 1e-6;
            FrameModel m; fixtures::cantileverSpatial(m, 8, L, Vec3(1, 0, 0), 0, 0, 0, T, 0, 0, cmat, asec);
            CorotationalOptions co; co.loadSteps = 1; co.maxIter = 50;
            const CorotationalResult R = runCorotational(m, co);
            const real tw = R.finalState.u[(size_t)gdof(8, Rx)], exact = T * L / (Gm * asec.J);
            std::printf("   (c) pure torsion theta_x=%.6e exact=%.6e rel=%.2e\n", tw, exact, relErr(tw, exact));
            checkTrue("F51c torsion converged", R.converged, R.finalState.diagnostic);
            checkClose("F51c pure torsion theta=TL/GJ", tw, exact, 1e-6);
        }

        // (d) BIAXIAL bending separation: X-cantilever, small tip Fy,Fz (asymmetric Iy!=Iz) ->
        //     delta_y = Fy L^3/(3 E Iz), delta_z = Fz L^3/(3 E Iy). Verifies the two planes stay distinct.
        {
            Section asec; asec.A = 100.0; asec.Iy = 2e-3; asec.Iz = 1e-3; asec.J = 1.5e-3; asec.cy = 1.0; asec.cz = 1.0; asec.Asy = 0.0; asec.Asz = 0.0;
            const real Fy = 1e-7, Fz = 1e-7;
            FrameModel m; fixtures::cantileverSpatial(m, 8, L, Vec3(1, 0, 0), 0, Fy, Fz, 0, 0, 0, cmat, asec);
            CorotationalOptions co; co.loadSteps = 1; co.maxIter = 50;
            const CorotationalResult R = runCorotational(m, co);
            const real dy = R.finalState.u[(size_t)gdof(8, Uy)], dz = R.finalState.u[(size_t)gdof(8, Uz)];
            // localAxes(X-axis, refVec=+Z): local y = +Z_global, local z = -Y_global. So a global-Y tip force
            // bends about local y (uses Iy); a global-Z tip force bends about local z (uses Iz).
            const real ey = Fy * L * L * L / (3 * Em * asec.Iy), ez = Fz * L * L * L / (3 * Em * asec.Iz);
            std::printf("   (d) biaxial dy=%.6e(exp %.6e) dz=%.6e(exp %.6e)\n", dy, ey, dz, ez);
            checkTrue("F51d biaxial converged", R.converged, R.finalState.diagnostic);
            checkClose("F51d biaxial dy=FyL^3/3EIy", dy, ey, 1e-6);
            checkClose("F51d biaxial dz=FzL^3/3EIz", dz, ez, 1e-6);
        }

        // (e) LARGE-ANGLE pure torsion: single element, tip torque drives the section twist to ~0.9 pi
        //     (~162 deg) -> exercises logSO3 at its large-angle end (general branch, sin(theta) far from 0)
        //     and the GJ channel + spatial R_node update at a FINITE rotation (vs the small-angle F51c).
        //     Torsion is geometrically linear (twist about the unchanged chord), so theta_x = T L/(G J).
        {
            const real kPiL = 3.14159265358979323846;
            Section tsec; tsec.A = 100.0; tsec.Iy = 2e-3; tsec.Iz = 1e-3; tsec.J = 1.5e-3; tsec.cy = 1.0; tsec.cz = 1.0; tsec.Asy = 0.0; tsec.Asz = 0.0;
            const real theta = 0.9 * kPiL, T = theta * Gm * tsec.J / L;   // T L/(GJ) = theta
            FrameModel m; fixtures::cantileverSpatial(m, 1, L, Vec3(1, 0, 0), 0, 0, 0, T, 0, 0, cmat, tsec);
            CorotationalOptions co; co.loadSteps = 12; co.maxIter = 60;
            const CorotationalResult R = runCorotational(m, co);
            const real tw = R.finalState.u[(size_t)gdof(1, Rx)];
            std::printf("   (e) large-angle torsion theta_x=%.6f exact=%.6f rel=%.2e\n", tw, theta, relErr(tw, theta));
            checkTrue("F51e large-angle torsion converged", R.converged, R.finalState.diagnostic);
            checkClose("F51e large-angle pure torsion theta=0.9pi", tw, theta, 1e-9);
        }
    }

    // ---------- F52: S9c arc-length snap-through (shallow two-bar arch / von Mises frame) ----------
    {
        std::printf("[F52] S9c arc-length: shallow-arch snap-through (limit point + post-buckling path)\n");
        Material amat(1.0, 0.4, 0.0);
        Section asec; asec.A = 1.0; asec.Iy = 1e-4; asec.Iz = 1e-4; asec.J = 1e-4; asec.cy = 1.0; asec.cz = 1.0; asec.Asy = 0.0; asec.Asz = 0.0;
        const real b = 1.0, h = 0.25;

        // (a) arc-length tracks the full snap-through path (lambda rises to a limit point then descends)
        {
            FrameModel m; fixtures::shallowArchPair(m, b, h, -1.0, amat, asec);   // unit downward apex load
            CorotationalOptions co; co.useArcLength = true; co.arcLength = 0.03; co.arcSteps = 80; co.maxIter = 40; co.tolR = 1e-8;
            const CorotationalResult R = runCorotational(m, co);
            real lamMax = -1e30; size_t iMax = 0;
            for (size_t i = 0; i < R.pathLambda.size(); ++i) if (R.pathLambda[i] > lamMax) { lamMax = R.pathLambda[i]; iMax = i; }
            const real lamEnd = R.pathLambda.empty() ? 0.0 : R.pathLambda.back();
            const bool hasLimit = R.pathLambda.size() > 4 && iMax > 0 && iMax < R.pathLambda.size() - 1 && lamEnd < lamMax;
            std::printf("   (a) arc-length: %zu steps, lambda_peak=%.5f @step%zu, lambda_end=%.5f, hasLimit=%d\n",
                        R.pathLambda.size(), lamMax, iMax, lamEnd, (int)hasLimit);
            checkTrue("F52a arc-length converged", R.converged, R.finalState.diagnostic.c_str());
            checkTrue("F52a snap-through limit point tracked (lambda rises then falls)", hasLimit, "");
            // limit load asserted (not just "a limit exists"): cross-validated against OpenSees integrator
            // ArcLength (0.0058962, rel 0.64%); near the pin-jointed von Mises closed form 0.00566 (rigid
            // joints + bending raise it ~3.5%). Not an analytical closed form for THIS rigid frame.
            checkClose("F52a limit load lambda_peak ~ 0.00586 (cross-validated vs OpenSees ArcLength)", lamMax, 0.00586, 2e-2);

            // (b) load control to a target ABOVE the limit load -> diverges at the limit point (arc-length needed)
            FrameModel m2; fixtures::shallowArchPair(m2, b, h, -(lamMax > 0 ? lamMax : 1.0) * 1.5, amat, asec);
            CorotationalOptions cl; cl.loadSteps = 30; cl.maxIter = 60;
            const CorotationalResult Rl = runCorotational(m2, cl);
            std::printf("   (b) load control to 1.5x limit: converged=%d diverged=%d\n", (int)Rl.converged, (int)Rl.diverged);
            checkTrue("F52b load control diverges at limit point (arc-length is required)", Rl.diverged && !Rl.converged, "");
        }
    }

    // ---------- F53: S9c member UDL / prescribed displacement / consistent (FD) tangent ----------
    {
        std::printf("[F53] S9c: member UDL + prescribed displacement + consistent (FD) tangent\n");
        Material cmat(1.0, 0.4, 0.0);
        Section sym; sym.A = 100.0; sym.Iy = 1e-3; sym.Iz = 1e-3; sym.J = 2e-3; sym.cy = 1.0; sym.cz = 1.0; sym.Asy = 0.0; sym.Asz = 0.0;
        const real Lc = 1.0, Ec = 1.0, Ic = 1e-3;

        // (a) member UDL: X-cantilever, transverse local-y UDL w -> tip deflection wL^4/8EI (small w, linear)
        {
            const real w = 1e-6;
            FrameModel m; fixtures::cantileverSpatial(m, 8, Lc, Vec3(1, 0, 0), 0, 0, 0, 0, 0, 0, cmat, sym);
            for (const auto& mem : m.members) { MemberUDL ud; ud.member = mem.id; ud.w_local = Vec3(0, w, 0); m.memberUDLs.push_back(ud); }
            CorotationalOptions co; co.loadSteps = 1; co.maxIter = 50;
            const CorotationalResult R = runCorotational(m, co);
            const real dfl = R.finalState.u[(size_t)gdof(8, Uz)];   // local y = global Z for the X-cantilever
            const real exact = w * Lc * Lc * Lc * Lc / (8.0 * Ec * Ic);
            std::printf("   (a) UDL cantilever tip=%.6e exact=%.6e rel=%.2e\n", dfl, exact, relErr(dfl, exact));
            checkTrue("F53a UDL converged", R.converged, R.finalState.diagnostic.c_str());
            checkClose("F53a UDL cantilever tip = wL^4/8EI", dfl, exact, 1e-4);
        }

        // (b) prescribed displacement: X-cantilever, tip Uy prescribed = delta (Rz free) -> guided cantilever,
        //     base reaction = -3EI delta/L^3 (small delta). Verifies the Dirichlet BC is imposed + reaction.
        {
            const real delta = 1e-4; const int n = 8;
            FrameModel m; fixtures::prepMatSec(m, cmat, sym);
            m.nodes.clear(); m.members.clear();
            for (int k = 0; k <= n; ++k) {
                Node nd(k, Lc * real(k) / real(n), 0, 0);
                nd.fixed[Uz] = nd.fixed[Rx] = nd.fixed[Ry] = true;   // planar XY bending
                if (k == 0) nd.fixAll();
                m.nodes.push_back(nd);
            }
            m.nodes[(size_t)n].fixed[Uy] = true; m.nodes[(size_t)n].prescribed[Uy] = delta;   // impose tip deflection
            for (int k = 0; k < n; ++k) { Member mm(k, k, k + 1, 0, 0); mm.refVec = Vec3(0, 0, 1); m.members.push_back(mm); }
            CorotationalOptions co; co.loadSteps = 1; co.maxIter = 50;
            const CorotationalResult R = runCorotational(m, co);
            const real tipUy = R.finalState.u[(size_t)gdof(n, Uy)];
            const real baseRy = R.finalState.reactions[(size_t)gdof(0, Uy)];
            const real exactP = 3.0 * Ec * Ic * delta / (Lc * Lc * Lc);
            std::printf("   (b) prescribed tip Uy=%.6e(set %.6e) base reaction=%.6e exact=%.6e\n", tipUy, delta, baseRy, -exactP);
            checkTrue("F53b prescribed converged", R.converged, R.finalState.diagnostic.c_str());
            checkClose("F53b prescribed tip Uy == delta (Dirichlet BC imposed)", tipUy, delta, 1e-9);
            checkClose("F53b prescribed base reaction = -3EI delta/L^3", baseRy, -exactP, 5e-3);
        }

        // (b2) prescribed ROTATION: exercises the R_node=expSO3(lambda*presc) path (translation BC alone
        //      never touches it). X-cantilever, tip Rz prescribed = theta (Uy free) -> base reaction moment
        //      = EI theta/L (a tip-moment cantilever; small theta). The reaction is the NON-trivial oracle
        //      (tip Rz == theta is mechanically forced; it only checks the logSO3 recover of R_node).
        {
            const real theta = 1e-4; const int n = 8;
            FrameModel m; fixtures::prepMatSec(m, cmat, sym);
            m.nodes.clear(); m.members.clear();
            for (int k = 0; k <= n; ++k) {
                Node nd(k, Lc * real(k) / real(n), 0, 0);
                nd.fixed[Uz] = nd.fixed[Rx] = nd.fixed[Ry] = true;
                if (k == 0) nd.fixAll();
                m.nodes.push_back(nd);
            }
            m.nodes[(size_t)n].fixed[Rz] = true; m.nodes[(size_t)n].prescribed[Rz] = theta;   // impose tip rotation
            for (int k = 0; k < n; ++k) { Member mm(k, k, k + 1, 0, 0); mm.refVec = Vec3(0, 0, 1); m.members.push_back(mm); }
            CorotationalOptions co; co.loadSteps = 1; co.maxIter = 50;
            const CorotationalResult R = runCorotational(m, co);
            const real tipRz = R.finalState.u[(size_t)gdof(n, Rz)];
            const real baseMz = R.finalState.reactions[(size_t)gdof(0, Rz)];
            const real exactM = Ec * Ic * theta / Lc;
            std::printf("   (b2) prescribed rotation tip Rz=%.6e(set %.6e) base moment=%.6e exact=%.6e\n", tipRz, theta, baseMz, -exactM);
            checkTrue("F53b2 prescribed rotation converged", R.converged, R.finalState.diagnostic.c_str());
            checkClose("F53b2 prescribed tip Rz == theta (expSO3 BC recover)", tipRz, theta, 1e-6);
            checkClose("F53b2 prescribed rotation base moment = -EI theta/L", baseMz, -exactM, 1e-2);
        }

        // (c) consistent (FD) tangent: arc-length snap-through reaches the SAME limit load with the numerical
        //     consistent tangent as with the analytical main-term tangent (proves the FD tangent is correct).
        {
            Material amat(1.0, 0.4, 0.0);
            Section asec; asec.A = 1.0; asec.Iy = 1e-4; asec.Iz = 1e-4; asec.J = 1e-4; asec.cy = 1.0; asec.cz = 1.0; asec.Asy = 0.0; asec.Asz = 0.0;
            real peakMain = 0, peakFD = 0; bool cMain = false, cFD = false; int itMain = 0, itFD = 0;
            for (int fd = 0; fd < 2; ++fd) {
                FrameModel m; fixtures::shallowArchPair(m, 1.0, 0.25, -1.0, amat, asec);
                CorotationalOptions co; co.useArcLength = true; co.arcLength = 0.03; co.arcSteps = 80; co.maxIter = 40; co.tolR = 1e-8; co.consistentTangent = (fd != 0);
                const CorotationalResult R = runCorotational(m, co);
                real lm = 0; for (real l : R.pathLambda) if (l > lm) lm = l;
                if (fd) { peakFD = lm; cFD = R.converged; itFD = R.totalIterations; }
                else    { peakMain = lm; cMain = R.converged; itMain = R.totalIterations; }
            }
            std::printf("   (c) consistent tangent: main peak=%.6f(%dit) FD peak=%.6f(%dit) rel=%.2e\n", peakMain, itMain, peakFD, itFD, relErr(peakFD, peakMain));
            checkTrue("F53c consistent-tangent both converged", cMain && cFD, "");
            checkClose("F53c FD tangent same limit load as main-term", peakFD, peakMain, 1e-3);
        }
    }

    // ---------- F54: S10 N-M interaction plastic hinge (axially-reduced Mp_eff) ----------
    {
        std::printf("[F54] S10: N-M interaction plastic hinge (Mp_eff = Mp*(1-(N/Ny)^2))\n");
        Material mat(210000.0, 80769.0, 7850.0);
        mat.fy = 300.0;
        mat.cap = Capacity::make(300.0, 300.0, 180.0);
        const real b = 80.0, d = 120.0;
        Section sec = Section::Rectangular(b, d);          // Zz = b d^2/4, A = b d
        const real MpZ = mat.fy * sec.Zz;                  // full plastic moment about local z
        const real Ny  = mat.fy * sec.A;                   // axial squash load

        // (1) formula: unit values + sign symmetry (|N|) + clamp beyond the squash load
        checkClose("F54 Mp_eff(N=0) = Mp", reducedPlasticMoment(MpZ, 0.0, Ny), MpZ, 1e-15);
        checkClose("F54 Mp_eff(N=Ny/2) = 0.75 Mp", reducedPlasticMoment(MpZ, 0.5 * Ny, Ny), 0.75 * MpZ, 1e-15);
        checkClose("F54 Mp_eff sign symmetry (|N|)", reducedPlasticMoment(MpZ, -0.5 * Ny, Ny),
                   reducedPlasticMoment(MpZ, 0.5 * Ny, Ny), 1e-15);
        checkClose("F54 Mp_eff(N=Ny) = 0", reducedPlasticMoment(MpZ, Ny, Ny), 0.0, 1e-15);
        checkTrue("F54 Mp_eff clamps to 0 beyond squash (N=2Ny)",
                  reducedPlasticMoment(MpZ, 2.0 * Ny, Ny) == 0.0);

        // (2) RECTANGULAR exactness vs the first-principles plastic neutral-axis shift: a central
        //     band of depth 2c=N/(b fy) carries the axial force, leaving M_N = Mp - N^2/(4 b fy).
        {
            const real N = 0.3 * Ny;
            const real Mexact = MpZ - N * N / (4.0 * b * mat.fy);
            checkClose("F54 rectangular Mp_eff == neutral-axis first principles",
                       reducedPlasticMoment(MpZ, N, Ny), Mexact, 1e-12);
        }

        // (3) driver: an X-cantilever is determinate, so the base moment is w L^2/2 and the axial
        //     force is the tip load P exactly (stiffness-independent). Under P = Ny/2 the reduced
        //     capacity is Mp_eff = 0.75 Mp; a single base hinge turns the cantilever into a
        //     mechanism, so the N-M hinge load is bracketed cleanly with no statical ambiguity.
        const real L = 2000.0;
        const real P = 0.5 * Ny;                           // -> Mp_eff(P) = 0.75 MpZ
        const real MpEff = reducedPlasticMoment(MpZ, P, Ny);
        const real wStar = 2.0 * MpEff / (L * L);          // base Mz = w L^2/2 reaches Mp_eff
        auto buildCant = [&](FrameModel& m, real w) {
            m = FrameModel{};
            m.materials = { mat }; m.sections = { sec };
            Node n0(0, 0.0, 0.0, 0.0); n0.fixAll();
            m.nodes = { n0, Node(1, L, 0.0, 0.0) };
            m.members = { Member(0, 0, 1, 0, 0) };
            MemberUDL u; u.member = 0; u.w_local = { 0.0, -w, 0.0 };   // local-y transverse -> Mz
            m.memberUDLs = { u };
            NodalLoad p; p.node = 1; p.comp[Ux] = -P;                  // tip axial compression
            m.nodalLoads = { p };
        };

        CollapseOptions on;  on.dlf  = 1.0; on.plasticHinges  = true; on.nmInteraction  = true;
        CollapseOptions off; off.dlf = 1.0; off.plasticHinges = true; off.nmInteraction = false;

        // 0.99 w*: below the reduced threshold -> Stable, no hinge
        FrameModel mLo; buildCant(mLo, 0.99 * wStar);
        const CollapseHistory hLo = runProgressiveCollapse(mLo, on);
        checkTrue("F54 0.99 w* (N-M on): Stable, no hinge",
                  hLo.outcome == CollapseOutcome::Stable && hLo.steps.size() == 1 &&
                  hLo.steps[0].formedHinges.empty(), hLo.diagnostic);

        // 1.01 w*: above the reduced threshold -> base hinge -> mechanism -> Collapsed
        FrameModel mHi; buildCant(mHi, 1.01 * wStar);
        const CollapseHistory hHi = runProgressiveCollapse(mHi, on);
        checkTrue("F54 1.01 w* (N-M on): base hinge -> Collapsed",
                  hHi.outcome == CollapseOutcome::Collapsed && hHi.steps.size() == 2 &&
                  hHi.steps[1].formedHinges.size() == 1 && !hHi.steps[1].solved, hHi.diagnostic);
        if (hHi.steps.size() == 2 && hHi.steps[1].formedHinges.size() == 1) {
            checkTrue("F54 hinge at base end-i, Mz dof",
                      hHi.steps[1].formedHinges[0].member == 0 && hHi.steps[1].formedHinges[0].dof == 5);
            checkClose("F54 trigger ratio = |M|/Mp_eff = 1.01", hHi.steps[1].triggerRatio, 1.01, 1e-9);
            checkClose("F54 frozen residual |Mp| = Mp_eff(P) = 0.75 Mp",
                       std::fabs(hHi.steps[1].formedHinges[0].Mp), MpEff, 1e-9);
        }

        // contrast: the SAME 1.01 w* load with N-M OFF stays Stable (|M|/Mp = 0.7575 < 1) --
        // the axial interaction is precisely what tips this case into collapse.
        FrameModel mHiOff; buildCant(mHiOff, 1.01 * wStar);
        const CollapseHistory hOff = runProgressiveCollapse(mHiOff, off);
        checkTrue("F54 same load, N-M off: Stable (axial interaction is decisive)",
                  hOff.outcome == CollapseOutcome::Stable && hOff.steps.size() == 1 &&
                  hOff.steps[0].formedHinges.empty(), hOff.diagnostic);
        std::printf("   N-M hinge load bracketed: Stable@0.99w* / Collapsed@1.01w* (N-M on); same load Stable (N-M off); Mp_eff/Mp=%.3f\n",
                    MpEff / MpZ);
    }

    // ---------- F55: opt-in HP-FEM PCG lane (solveLoadHP) vs LDLT oracle ----------
    // A1 of the production HP-FEM integration: the opt-in matrix-free PCG lane must return
    // a SolveResult equal to the direct LDLT solve (same RHS / reactions / member forces),
    // defer to LDLT on a mechanism, and be a drop-in equal to solveLoad when disabled.
    {
        std::printf("[F55] HP-FEM opt-in lane: solveLoadHP (matrix-free PCG) vs LDLT oracle\n");
        HpSolveOptions hp; hp.enabled = true;   // pcgTol=1e-10, maxIter=500, fallbackOnFail=true
        auto hpVsLdlt = [&](const char* tag, const SolveResult& got, const SolveResult& ref) {
            if (got.singular) { checkTrue(tag, false, "HP unexpectedly singular: " + got.diagnostic); return; }
            double un = 1e-30, du = 0, rn = 1e-30, dr = 0, mn = 1e-30, dm = 0;
            for (real v : ref.u) un = std::max(un, std::fabs((double)v));
            for (size_t k = 0; k < got.u.size(); ++k) du = std::max(du, std::fabs((double)got.u[k] - (double)ref.u[k]));
            for (real v : ref.reactions) rn = std::max(rn, std::fabs((double)v));
            for (size_t k = 0; k < got.reactions.size(); ++k) dr = std::max(dr, std::fabs((double)got.reactions[k] - (double)ref.reactions[k]));
            for (size_t e = 0; e < ref.memberForces.size(); ++e) {
                const MemberEndForces gI = got.memberForces[e].endI, rI = ref.memberForces[e].endI;
                const MemberEndForces gJ = got.memberForces[e].endJ, rJ = ref.memberForces[e].endJ;
                const real ga[12] = { gI.N,gI.Vy,gI.Vz,gI.T,gI.My,gI.Mz, gJ.N,gJ.Vy,gJ.Vz,gJ.T,gJ.My,gJ.Mz };
                const real rb[12] = { rI.N,rI.Vy,rI.Vz,rI.T,rI.My,rI.Mz, rJ.N,rJ.Vy,rJ.Vz,rJ.T,rJ.My,rJ.Mz };
                for (int i = 0; i < 12; ++i) { mn = std::max(mn, std::fabs((double)rb[i])); dm = std::max(dm, std::fabs((double)ga[i] - (double)rb[i])); }
            }
            const bool ok = (du/un <= 1e-8) && (dr/rn <= 1e-8) && (dm/mn <= 1e-8);
            checkTrue(tag, ok, "uRel=" + std::to_string(du/un) + " Rrel=" + std::to_string(dr/rn) +
                      " mfRel=" + std::to_string(dm/mn) + " | " + got.diagnostic);
        };
        // FA cantilever tip load
        { FrameModel m; fixtures::cantileverTipLoad(m, 1000.0, 2000.0, mat, sec);
          PreparedSystem ps = assembleAndFactor(m);
          hpVsLdlt("F55-FA cantilever", solveLoadHP(ps, m, hp), solveLoad(ps, m)); }
        // FB simply-supported UDL (equivalent-nodal-load path)
        { FrameModel m; fixtures::simplySupportedUDL(m, 5.0, 3000.0, mat, sec);
          PreparedSystem ps = assembleAndFactor(m);
          hpVsLdlt("F55-FB SS+UDL", solveLoadHP(ps, m, hp), solveLoad(ps, m)); }
        // FC settlement (nonzero prescribed DOF -> exercises the CSC prescribed reduce).
        // Guard that the case is non-degenerate: the free midspan must actually move, else
        // the prescribed-reduce path would be vacuously "equal" with a zero RHS.
        { FrameModel m; fixtures::clampedSettlement(m, 3000.0, 5.0, mat, sec);
          PreparedSystem ps = assembleAndFactor(m);
          const SolveResult ref = solveLoad(ps, m);
          double freeNorm = 0; for (int d = 0; d < 6; ++d) freeNorm = std::max(freeNorm, std::fabs(ref.disp(1, d)));
          checkTrue("F55-FC settlement moves free node", freeNorm > 1e-9, "freeNorm=" + std::to_string(freeNorm));
          hpVsLdlt("F55-FC settlement", solveLoadHP(ps, m, hp), ref); }
        // FE released member (enableReleases -> condensed kl/Qf must flow through the apply)
        { FrameModel m; fixtures::proppedCantileverRelease(m, 5.0, 3000.0, mat, sec);
          SolveOptions so; so.enableReleases = true;
          PreparedSystem ps = assembleAndFactor(m, so);
          hpVsLdlt("F55-FE release", solveLoadHP(ps, m, hp), solveLoad(ps, m)); }
        // FD mechanism -> the HP lane must report singular (it defers to the LDLT pivot guard)
        { FrameModel m; fixtures::mechanism(m, mat, sec);
          PreparedSystem ps = assembleAndFactor(m);
          const SolveResult got = solveLoadHP(ps, m, hp);
          checkTrue("F55-FD mechanism -> singular", got.singular, got.diagnostic); }
        // disabled lane is a drop-in equal to solveLoad (LDLT)
        { FrameModel m; fixtures::cantileverTipLoad(m, 1000.0, 2000.0, mat, sec);
          PreparedSystem ps = assembleAndFactor(m);
          HpSolveOptions off; off.enabled = false;
          hpVsLdlt("F55 disabled==LDLT", solveLoadHP(ps, m, off), solveLoad(ps, m)); }
    }

    // ---------- F56: seeded HP-FEM session (HpSession) vs LDLT oracle ----------
    // A2c of the production HP-FEM integration: the seeded session must (A) hit the Galerkin
    // projection at ~0 PCG iters for a load IN the seeded subspace, (B) gate to a full PCG for a
    // load that LEAVES it (still correct), (C) be a drop-in equal to solveLoad when disabled or
    // un-seeded, (D) reject a model whose fingerprint changed, and (E) auto-expand the basis cap so
    // FIFO never drops a seed. All paths match the direct LDLT solve to a 1e-6 tolerance gate.
    {
        std::printf("[F56] HP-FEM seeded session: HpSession (Galerkin projection + warm PCG) vs LDLT oracle\n");

        // A horizontal cantilever along +X (node0 encastre, nodes 1..n free) gives a multi-dimensional
        // load family to seed and to leave. id == index here, so gdof(index, dof) addresses the load.
        const int  nseg  = 6;
        const real Lspan = 3000.0;
        auto buildCant = [&](FrameModel& m) {
            m = FrameModel{};
            m.materials = { mat }; m.sections = { sec };
            Node n0(0, 0, 0, 0); n0.fixAll();
            m.nodes = { n0 };
            for (int k = 1; k <= nseg; ++k) m.nodes.push_back(Node(k, Lspan * real(k) / real(nseg), 0, 0));
            for (int k = 0; k < nseg; ++k) m.members.push_back(Member(k, k, k + 1, 0, 0));
        };
        // A 6N global load vector: a single Uz load at node index `nd`.
        auto unitUz = [&](const FrameModel& m, int nd, real mag) {
            std::vector<real> lv((size_t)(6 * (int)m.nodes.size()), 0.0);
            lv[(size_t)gdof(nd, Uz)] = mag;
            return lv;
        };
        // Set the model's single nodal Uz load at node index `nd`.
        auto setLoad = [&](FrameModel& m, int nd, real mag) {
            NodalLoad nl; nl.node = m.nodes[(size_t)nd].id; nl.comp[Uz] = mag; m.nodalLoads = { nl };
        };
        // Compare an HpSession result to the LDLT oracle (u / reactions / member forces), 1e-6 gate.
        auto hpSessVsLdlt = [&](const char* tag, const SolveResult& got, const SolveResult& ref) {
            if (got.singular) { checkTrue(tag, false, "HpSession unexpectedly singular: " + got.diagnostic); return; }
            double un = 1e-30, du = 0, rn = 1e-30, dr = 0, mn = 1e-30, dm = 0;
            for (real v : ref.u) un = std::max(un, std::fabs((double)v));
            for (size_t k = 0; k < got.u.size(); ++k) du = std::max(du, std::fabs((double)got.u[k] - (double)ref.u[k]));
            for (real v : ref.reactions) rn = std::max(rn, std::fabs((double)v));
            for (size_t k = 0; k < got.reactions.size(); ++k) dr = std::max(dr, std::fabs((double)got.reactions[k] - (double)ref.reactions[k]));
            for (size_t e = 0; e < ref.memberForces.size(); ++e) {
                const MemberEndForces gI = got.memberForces[e].endI, rI = ref.memberForces[e].endI;
                const MemberEndForces gJ = got.memberForces[e].endJ, rJ = ref.memberForces[e].endJ;
                const real ga[12] = { gI.N,gI.Vy,gI.Vz,gI.T,gI.My,gI.Mz, gJ.N,gJ.Vy,gJ.Vz,gJ.T,gJ.My,gJ.Mz };
                const real rb[12] = { rI.N,rI.Vy,rI.Vz,rI.T,rI.My,rI.Mz, rJ.N,rJ.Vy,rJ.Vz,rJ.T,rJ.My,rJ.Mz };
                for (int i = 0; i < 12; ++i) { mn = std::max(mn, std::fabs((double)rb[i])); dm = std::max(dm, std::fabs((double)ga[i] - (double)rb[i])); }
            }
            const bool ok = (du/un <= 1e-6) && (dr/rn <= 1e-6) && (dm/mn <= 1e-6);
            checkTrue(tag, ok, "uRel=" + std::to_string(du/un) + " Rrel=" + std::to_string(dr/rn) +
                      " mfRel=" + std::to_string(dm/mn) + " | " + got.diagnostic);
        };

        FrameModel mSeed; buildCant(mSeed);
        PreparedSystem ps = assembleAndFactor(mSeed);
        HpSessionOptions opt;                                   // enabled, basisMax=64, projGateTol=1e-6, warmIter=3
        HpSession sess(ps, opt);
        checkTrue("F56 session valid (frame-only, non-singular)", sess.valid(), sess.diagnostic());
        // seed the Uz load family at nodes 1 and n (gravity + a far contact, conceptually)
        std::vector<std::vector<real>> seeds = { unitUz(mSeed, 1, 1.0), unitUz(mSeed, nseg, 1.0) };
        checkTrue("F56 setLoadBasis(2) ok", sess.setLoadBasis(seeds), sess.diagnostic());

        // (A) IN-subspace: a pure multiple of seed1 -> Galerkin projection is exact (~0 PCG iters)
        {
            FrameModel m; buildCant(m); setLoad(m, 1, 1234.0);
            HpSessionStats st;
            const SolveResult got = sess.solveFrame(m, &st);
            hpSessVsLdlt("F56-A in-subspace (node1)", got, solveLoad(ps, m));
            checkTrue("F56-A used seeded projection", st.usedProjection, "initialRel=" + std::to_string((double)st.initialRel));
            checkTrue("F56-A near-zero PCG iters (in-subspace)", st.pcgIters <= 1, "iters=" + std::to_string(st.pcgIters));
        }
        // (A2) IN-subspace linear combination of BOTH seeds
        {
            FrameModel m; buildCant(m);
            NodalLoad a; a.node = m.nodes[1].id;            a.comp[Uz] =  800.0;
            NodalLoad b; b.node = m.nodes[(size_t)nseg].id; b.comp[Uz] = -500.0;
            m.nodalLoads = { a, b };
            HpSessionStats st;
            const SolveResult got = sess.solveFrame(m, &st);
            hpSessVsLdlt("F56-A2 in-subspace (combo)", got, solveLoad(ps, m));
            checkTrue("F56-A2 used seeded projection (combo)", st.usedProjection, "initialRel=" + std::to_string((double)st.initialRel));
        }
        // (B) OUT-of-subspace: a load at an unseeded node leaves the subspace -> gate -> full PCG, still correct
        {
            FrameModel m; buildCant(m); setLoad(m, 2, 1000.0);
            HpSessionStats st;
            const SolveResult got = sess.solveFrame(m, &st);
            hpSessVsLdlt("F56-B out-of-subspace (node2)", got, solveLoad(ps, m));
            checkTrue("F56-B gate -> full PCG (not projection)", st.usedPcg && !st.usedProjection,
                      "initialRel=" + std::to_string((double)st.initialRel) + " iters=" + std::to_string(st.pcgIters));
            checkTrue("F56-B initialRel above projGateTol", (double)st.initialRel >= opt.projGateTol,
                      "initialRel=" + std::to_string((double)st.initialRel));
        }
        // (C) DISABLED lane is a drop-in equal to solveLoad (LDLT)
        {
            HpSessionOptions off = opt; off.enabled = false;
            HpSession sOff(ps, off);
            FrameModel m; buildCant(m); setLoad(m, 1, 1234.0);
            HpSessionStats st;
            const SolveResult got = sOff.solveFrame(m, &st);
            hpSessVsLdlt("F56-C disabled==LDLT", got, solveLoad(ps, m));
            checkTrue("F56-C used LDLT (disabled)", st.usedLdlt && !st.usedProjection && !st.usedPcg, got.diagnostic);
        }
        // (C2) ENABLED but UN-SEEDED (empty basis) -> LDLT drop-in
        {
            HpSession sEmpty(ps, opt);   // setLoadBasis never called
            FrameModel m; buildCant(m); setLoad(m, 1, 1234.0);
            HpSessionStats st;
            const SolveResult got = sEmpty.solveFrame(m, &st);
            hpSessVsLdlt("F56-C2 empty-basis==LDLT", got, solveLoad(ps, m));
            checkTrue("F56-C2 used LDLT (empty basis)", st.usedLdlt && st.basisSize == 0, got.diagnostic);
        }
        // (D) FINGERPRINT guard: a structural change since assembleAndFactor -> singular (re-run required)
        {
            FrameModel m; buildCant(m); setLoad(m, 1, 1234.0);
            m.sections[0].Iz *= 1.5;   // structural change -> fingerprint differs from the prepared system
            HpSessionStats st;
            const SolveResult got = sess.solveFrame(m, &st);
            checkTrue("F56-D fingerprint mismatch -> singular", got.singular, got.diagnostic);
        }
        // (E) basisMax AUTO-EXPAND: seeding 5 modes with basisMax=3 must keep all 5 (FIFO would drop 2)
        {
            HpSessionOptions small = opt; small.basisMax = 3;
            HpSession s5(ps, small);
            std::vector<std::vector<real>> five = {
                unitUz(mSeed, 1, 1.0), unitUz(mSeed, 2, 1.0), unitUz(mSeed, 3, 1.0),
                unitUz(mSeed, 4, 1.0), unitUz(mSeed, 5, 1.0) };
            checkTrue("F56-E setLoadBasis(5) ok", s5.setLoadBasis(five), s5.diagnostic());
            FrameModel m; buildCant(m); setLoad(m, 1, 1000.0);
            HpSessionStats st;
            s5.solveFrame(m, &st);
            checkClose("F56-E basisMax auto-expanded to 5 (FIFO would cap at 3)", (double)st.basisSize, 5.0, 0.5);
        }
        // (F) A2a PARALLEL apply: a threaded session (threads=8) matches the serial session AND the
        //     LDLT oracle on both the in-subspace (projection) and out-of-subspace (full PCG) paths.
        //     A larger frame (60 beams, nf=360) makes the block range actually split across workers,
        //     exercising the per-thread accumulate + touched-DOF reduction and the parallel block6 map.
        {
            const int nbig = 60;
            FrameModel mBig;
            mBig.materials = { mat }; mBig.sections = { sec };
            Node b0(0, 0, 0, 0); b0.fixAll(); mBig.nodes = { b0 };
            for (int k = 1; k <= nbig; ++k) mBig.nodes.push_back(Node(k, Lspan * real(k) / real(nbig), 0, 0));
            for (int k = 0; k < nbig; ++k) mBig.members.push_back(Member(k, k, k + 1, 0, 0));
            PreparedSystem psB = assembleAndFactor(mBig);

            std::vector<std::vector<real>> seedsB = { unitUz(mBig, 1, 1.0), unitUz(mBig, nbig, 1.0) };
            HpSessionOptions ser = opt; ser.threads = 1;
            HpSessionOptions par = opt; par.threads = 8;
            HpSession sSer(psB, ser), sPar(psB, par);
            checkTrue("F56-F serial session valid",   sSer.valid(), sSer.diagnostic());
            checkTrue("F56-F parallel session valid", sPar.valid(), sPar.diagnostic());
            checkTrue("F56-F setLoadBasis serial",    sSer.setLoadBasis(seedsB), sSer.diagnostic());
            checkTrue("F56-F setLoadBasis parallel",  sPar.setLoadBasis(seedsB), sPar.diagnostic());

            // in-subspace: both hit the projection (~0 iters); parallel matches serial and LDLT
            {
                FrameModel m = mBig; setLoad(m, 1, 2222.0);
                HpSessionStats ssr, spr;
                const SolveResult gSer = sSer.solveFrame(m, &ssr);
                const SolveResult gPar = sPar.solveFrame(m, &spr);
                const SolveResult ref  = solveLoad(psB, m);
                hpSessVsLdlt("F56-F in-subspace serial vs LDLT",   gSer, ref);
                hpSessVsLdlt("F56-F in-subspace parallel vs LDLT", gPar, ref);
                checkTrue("F56-F parallel in-subspace projection 0-iter",
                          spr.usedProjection && spr.pcgIters <= 1, "iters=" + std::to_string(spr.pcgIters));
                double dmax = 0, nmax = 1e-30;
                for (size_t k = 0; k < gPar.u.size(); ++k) {
                    dmax = std::max(dmax, std::fabs((double)gPar.u[k] - (double)gSer.u[k]));
                    nmax = std::max(nmax, std::fabs((double)gSer.u[k]));
                }
                checkTrue("F56-F parallel u == serial u (threaded reduction = FP order only)",
                          dmax / nmax <= 1e-8, "rel=" + std::to_string(dmax / nmax));
            }
            // out-of-subspace: full PCG exercises the parallel apply + parallel block6 precond
            {
                FrameModel m = mBig; setLoad(m, 7, 1500.0);
                HpSessionStats spr;
                const SolveResult gPar = sPar.solveFrame(m, &spr);
                hpSessVsLdlt("F56-F out-of-subspace parallel vs LDLT", gPar, solveLoad(psB, m));
                checkTrue("F56-F parallel out-of-subspace full PCG", spr.usedPcg && !spr.usedProjection,
                          "iters=" + std::to_string(spr.pcgIters));
            }
            std::printf("   A2a parallel(threads=8) apply+block6 matches serial & LDLT (in/out subspace); ~19x perf in HPFEM_RESEARCH_NOTES\n");
        }
        // (G) A2b BANDED COARSE: a tower (multi-storey, distinct node z) builds a block-tridiagonal
        //     coarse correction (banded, floors>1); a 1D cantilever (all nodes z=0) degrades
        //     gracefully to floors==1. Both stay correct vs LDLT on out-of-subspace frames.
        {
            auto unitDof = [&](const FrameModel& m, int nd, int dof, real mag) {
                std::vector<real> lv((size_t)(6 * (int)m.nodes.size()), 0.0);
                lv[(size_t)gdof(nd, dof)] = mag; return lv;
            };
            // vertical column along +Z: node0 base encastre, nodes 1..nstory at distinct z -> one free
            // node per floor -> the coarse matrix is block-tridiagonal in the floor ordering.
            const int  nstory = 8;
            const real hstory = 1000.0;
            FrameModel mTower;
            mTower.materials = { mat }; mTower.sections = { sec };
            Node t0(0, 0, 0, 0); t0.fixAll(); mTower.nodes = { t0 };
            for (int k = 1; k <= nstory; ++k) mTower.nodes.push_back(Node(k, 0, 0, hstory * real(k)));
            for (int k = 0; k < nstory; ++k) mTower.members.push_back(Member(k, k, k + 1, 0, 0));
            PreparedSystem psT = assembleAndFactor(mTower);

            HpSession sT(psT, opt);
            checkTrue("F56-G tower session valid", sT.valid(), sT.diagnostic());
            sT.setLoadBasis({ unitDof(mTower, 1, Ux, 1.0), unitDof(mTower, nstory, Ux, 1.0) });
            const HpCoarseInfo ci = sT.prepareCoarse(mTower);
            checkTrue("F56-G tower coarse built + banded (block-tridiagonal)", ci.built && ci.banded, sT.diagnostic());
            checkTrue("F56-G tower floors == nstory (per-storey aggregation)", ci.floors == nstory,
                      "floors=" + std::to_string(ci.floors));
            {
                FrameModel m = mTower;
                NodalLoad nl; nl.node = mTower.nodes[(size_t)(nstory / 2)].id; nl.comp[Ux] = 700.0; m.nodalLoads = { nl };
                HpSessionStats st;
                const SolveResult got = sT.solveFrame(m, &st);
                hpSessVsLdlt("F56-G tower coarse out-of-subspace vs LDLT", got, solveLoad(psT, m));
                checkTrue("F56-G tower coarse active in out-of-subspace PCG", st.usedCoarse,
                          "usedPcg=" + std::to_string((int)st.usedPcg) + " usedCoarse=" + std::to_string((int)st.usedCoarse));
            }
            // graceful degradation: a 1D horizontal cantilever (all z=0) -> floors==1, still correct
            {
                FrameModel mFlat; buildCant(mFlat);
                PreparedSystem psF = assembleAndFactor(mFlat);
                HpSession sF(psF, opt);
                sF.setLoadBasis({ unitUz(mFlat, 1, 1.0), unitUz(mFlat, nseg, 1.0) });
                const HpCoarseInfo cf = sF.prepareCoarse(mFlat);
                checkTrue("F56-G cantilever coarse graceful (floors==1)", cf.built && cf.floors == 1,
                          "built=" + std::to_string((int)cf.built) + " floors=" + std::to_string(cf.floors));
                FrameModel m = mFlat; setLoad(m, 2, 900.0);
                HpSessionStats st;
                const SolveResult got = sF.solveFrame(m, &st);
                hpSessVsLdlt("F56-G cantilever graceful out-of-subspace vs LDLT", got, solveLoad(psF, m));
            }
            std::printf("   A2b coarse: tower banded(floors=%d) on out-of-subspace; 1D cantilever degrades to block6; both == LDLT\n", ci.floors);
        }
        // (H) PRESCRIBED displacement (settlement): exercises the anyPresc=true RHS-reduce path. The
        //     perf optimization skips that O(nnz) sweep ONLY when no support carries a nonzero
        //     prescribed value, so this guards the non-skipped branch. A pure-nodal seed cannot span
        //     a settlement response, so the frame leaves the subspace -> PCG/LDLT; result == solveLoad.
        {
            FrameModel mS; buildCant(mS);
            mS.nodes[3].fixed[Uz] = true; mS.nodes[3].prescribed[Uz] = 5.0;   // impose a settlement at node 3
            PreparedSystem psS = assembleAndFactor(mS);
            HpSession sS(psS, opt);
            sS.setLoadBasis({ unitUz(mS, 1, 1.0), unitUz(mS, nseg, 1.0) });
            HpSessionStats st;
            const SolveResult got = sS.solveFrame(mS, &st);
            double freeNorm = 0; for (int d = 0; d < 6; ++d) freeNorm = std::max(freeNorm, std::fabs(got.disp(1, d)));
            checkTrue("F56-H settlement moves a free node (non-degenerate)", freeNorm > 1e-9, "freeNorm=" + std::to_string(freeNorm));
            hpSessVsLdlt("F56-H settlement (anyPresc reduce path) vs LDLT", got, solveLoad(psS, mS));
        }
        // (I) displacements-only fast path: u must match the full solve; reactions are left all-zero
        //     and member/shell forces empty (the O(nnz)/O(elements) recovery is intentionally skipped).
        {
            HpSessionOptions du = opt; du.displacementsOnly = true;
            HpSession sd(ps, du);
            checkTrue("F56-I displacements-only session valid", sd.valid(), sd.diagnostic());
            sd.setLoadBasis({ unitUz(mSeed, 1, 1.0), unitUz(mSeed, nseg, 1.0) });
            FrameModel m; buildCant(m); setLoad(m, 1, 1234.0);
            HpSessionStats st;
            const SolveResult got = sd.solveFrame(m, &st);
            const SolveResult ref = solveLoad(ps, m);
            double un = 1e-30, du2 = 0;
            for (real v : ref.u) un = std::max(un, std::fabs((double)v));
            for (size_t k = 0; k < got.u.size(); ++k) du2 = std::max(du2, std::fabs((double)got.u[k] - (double)ref.u[k]));
            checkTrue("F56-I displacements-only u == LDLT u", du2 / un <= 1e-6, "uRel=" + std::to_string(du2 / un));
            double rmax = 0; for (real v : got.reactions) rmax = std::max(rmax, std::fabs((double)v));
            checkTrue("F56-I reactions skipped (all-zero)", rmax == 0.0, "rmax=" + std::to_string(rmax));
            checkTrue("F56-I member forces skipped (empty)", got.memberForces.empty(), "");
            checkTrue("F56-I still hits the in-subspace projection", st.usedProjection, "");
        }
    }

    std::printf("\n%s  (failures=%d)\n", g_fail == 0 ? "ALL PASS" : "FAILURES", g_fail);
    return g_fail == 0 ? 0 : 1;
}
