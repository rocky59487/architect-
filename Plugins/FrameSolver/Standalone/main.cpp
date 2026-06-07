// Standalone gate for FrameCore — builds fixtures, solves, compares to closed-form
// analytic solutions, prints PASS/FAIL, exits 0 iff all pass.
#include "FrameCore/FrameSolver.h"
#include "FrameCore/ElasticAllowable.h"
#include "FrameCore/Grillage.h"
#include "FrameCore/SelfWeight.h"
#include "FrameCore/Combination.h"
#include "FrameCore/InfluenceLine.h"
#include "FrameTestFixtures.h"

#include <vector>
#include <utility>

#include <cstdio>
#include <cmath>
#include <string>
#include <algorithm>

using namespace frame;

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
        mr.materials.reserve(1); mr.sections.reserve(1);
        mr.materials.push_back(mat); mr.sections.push_back(sec);
        const Material* pm = &mr.materials.back(); const Section* ps = &mr.sections.back();
        Node n0(0, 0, 0, 0); n0.fixAll();
        Node n1(1, rot(L,0,0,0), rot(L,0,0,1), rot(L,0,0,2));
        mr.nodes = { n0, n1 };
        Member mm(0, 0, 1, pm, ps);
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
        const real w4 = centerW(4), w8 = centerW(8), w16 = centerW(16);
        const real e4 = relErr(w4, wc), e8 = relErr(w8, wc), e16 = relErr(w16, wc);
        std::printf("   simply-supported square  plate w_c=%.6g mm\n", wc);
        std::printf("   N=4 w=%.5g (err %.2f%%)  N=8 w=%.5g (err %.2f%%)  N=16 w=%.5g (err %.2f%%)\n",
                    w4, 100 * e4, w8, 100 * e8, w16, 100 * e16);
        checkTrue("plate converges to Kirchhoff (N=16 within 2%)", e16 < 0.02,
                  "e16=" + std::to_string(e16));
        checkTrue("plate mesh-converging (e16 < e4)", e16 < e4,
                  "e4=" + std::to_string(e4) + " e16=" + std::to_string(e16));

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

    std::printf("\n%s  (failures=%d)\n", g_fail == 0 ? "ALL PASS" : "FAILURES", g_fail);
    return g_fail == 0 ? 0 : 1;
}
