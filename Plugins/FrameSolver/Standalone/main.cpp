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
#include "FrameCore/ResponseSpectrum.h"
#include "FrameCore/ModalDynamics.h"
#include "FrameCore/Connectivity.h"
#include "FrameCore/Collapse.h"
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

    std::printf("\n%s  (failures=%d)\n", g_fail == 0 ? "ALL PASS" : "FAILURES", g_fail);
    return g_fail == 0 ? 0 : 1;
}
