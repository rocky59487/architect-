// Standalone gate for FrameCore — builds fixtures, solves, compares to closed-form
// analytic solutions, prints PASS/FAIL, exits 0 iff all pass.
#include "FrameCore/FrameSolver.h"
#include "FrameCore/ElasticAllowable.h"
#include "FrameCore/Grillage.h"
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

    std::printf("\n%s  (failures=%d)\n", g_fail == 0 ? "ALL PASS" : "FAILURES", g_fail);
    return g_fail == 0 ? 0 : 1;
}
