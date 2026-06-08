#include "FrameCore/BucklingAnalysis.h"
#include "FrameCore/Combination.h"
#include "FrameCore/FrameSolver.h"
#include "FrameCore/InfluenceLine.h"
#include "FrameCore/ModalAnalysis.h"
#include "FrameCore/ModalDynamics.h"
#include "FrameCore/ResponseSpectrum.h"
#include "FrameCore/SelfWeight.h"
#include "FrameCore/ElasticAllowable.h"
#include "FrameTestFixtures.h"
#include "MITC4ShellElement.h"   // Private seam: element + Eigen types for the element-level audit

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <Eigen/Eigenvalues>

using namespace frame;

#ifndef FRAMECORE_BUILD_SHA
#define FRAMECORE_BUILD_SHA "unknown"   // overridden by the build script via /D (git short SHA)
#endif

namespace {

constexpr real kPi = 3.14159265358979323846;

struct AuditRow {
    std::string area;
    std::string strengthened;
    std::string oracle;
    std::string metric;
    real value = 0.0;
    real tol = 0.0;
    bool ok = false;
};

std::vector<AuditRow> g_rows;

real relErr(real got, real expected) {
    return std::fabs(got - expected) / std::max<real>(1.0, std::fabs(expected));
}

std::string sci(real v) {
    std::ostringstream os;
    os << std::scientific << std::setprecision(3) << v;
    return os.str();
}

void addRow(const std::string& area, const std::string& strengthened,
            const std::string& oracle, const std::string& metric,
            real value, real tol, bool ok) {
    g_rows.push_back(AuditRow{ area, strengthened, oracle, metric, value, tol, ok });
}

std::array<real, 6> endVals(const MemberEndForces& f) {
    return { f.N, f.Vy, f.Vz, f.T, f.My, f.Mz };
}

std::array<real, 8> shellVals(const ShellElementForces& f) {
    return { f.Mxx, f.Myy, f.Mxy, f.Qx, f.Qy, f.Nxx, f.Nyy, f.Nxy };
}

real maxDiff(const SolveResult& a, const SolveResult& b) {
    real d = 0.0;
    const size_t nu = std::min(a.u.size(), b.u.size());
    const size_t nr = std::min(a.reactions.size(), b.reactions.size());
    for (size_t i = 0; i < nu; ++i) d = std::max(d, std::fabs(a.u[i] - b.u[i]));
    for (size_t i = 0; i < nr; ++i) d = std::max(d, std::fabs(a.reactions[i] - b.reactions[i]));
    const size_t nm = std::min(a.memberForces.size(), b.memberForces.size());
    for (size_t i = 0; i < nm; ++i) {
        const auto ai = endVals(a.memberForces[i].endI);
        const auto bi = endVals(b.memberForces[i].endI);
        const auto aj = endVals(a.memberForces[i].endJ);
        const auto bj = endVals(b.memberForces[i].endJ);
        for (size_t k = 0; k < ai.size(); ++k) {
            d = std::max(d, std::fabs(ai[k] - bi[k]));
            d = std::max(d, std::fabs(aj[k] - bj[k]));
        }
    }
    const size_t ns = std::min(a.shellForces.size(), b.shellForces.size());
    for (size_t i = 0; i < ns; ++i) {
        const auto av = shellVals(a.shellForces[i]);
        const auto bv = shellVals(b.shellForces[i]);
        for (size_t k = 0; k < av.size(); ++k) d = std::max(d, std::fabs(av[k] - bv[k]));
    }
    return d;
}

real maxAbs(const std::vector<real>& v) {
    real m = 0.0;
    for (real x : v) m = std::max(m, std::fabs(x));
    return m;
}

real peakDof(const ModalTimeHistory& th, int gd) {
    real p = 0.0;
    for (const auto& u : th.u) {
        if (gd >= 0 && static_cast<size_t>(gd) < u.size()) p = std::max(p, std::fabs(u[static_cast<size_t>(gd)]));
    }
    return p;
}

real sumOf(const std::vector<real>& v) {
    real s = 0.0;
    for (real x : v) s += x;
    return s;
}

void buildPureBendingSS(FrameModel& m, int n, real L, real rho, const Section& sec) {
    Material mat(210000.0, 80769.230769, rho);
    fixtures::simplySupportedBeamN(m, n, L, mat, sec);
    for (Node& nd : m.nodes) {
        nd.fixed[Ux] = true;
        nd.fixed[Uy] = true;
        nd.fixed[Rx] = true;
        nd.fixed[Rz] = true;
    }
}

void testShellCombinationEnvelope() {
    Material mat(30000.0, 11538.461538, 2500.0);
    mat.nu = 0.30;
    FrameModel m1, m2;
    fixtures::squarePlateShell(m1, 2000.0, 120.0, 4, 0.0020, mat);
    fixtures::squarePlateShell(m2, 2000.0, 120.0, 4, 0.0045, mat);
    const SolveResult r1 = solve(m1);
    const SolveResult r2 = solve(m2);
    const bool shellOk = !r1.singular && !r2.singular && !r1.shellForces.empty() && r1.shellForces.size() == r2.shellForces.size();
    addRow("Combination/envelope", "shell resultants included in post-processing",
           "two MITC4 pressure cases solve and expose shellForces",
           "shell case availability", shellOk ? 0.0 : 1.0, 0.0, shellOk);
    if (!shellOk) return;

    const SolveResult c = combine({ r1, r2 }, { 1.25, -0.40 });
    real combErr = 0.0;
    for (size_t s = 0; s < c.shellForces.size(); ++s) {
        const auto cv = shellVals(c.shellForces[s]);
        const auto a = shellVals(r1.shellForces[s]);
        const auto b = shellVals(r2.shellForces[s]);
        for (size_t k = 0; k < cv.size(); ++k) combErr = std::max(combErr, std::fabs(cv[k] - (1.25 * a[k] - 0.40 * b[k])));
    }
    addRow("Combination/envelope", "load combination now checks shell forces too",
           "component arithmetic: 1.25*caseA - 0.40*caseB",
           "max shell combine abs diff", combErr, 1e-9, combErr < 1e-9);

    const ResultEnvelope e = envelope({ r1, r2 });
    real envErr = (e.shellMax.size() == r1.shellForces.size() && e.shellMin.size() == r1.shellForces.size()) ? 0.0 : 1.0;
    if (envErr == 0.0) {
        for (size_t s = 0; s < r1.shellForces.size(); ++s) {
            const auto hi = shellVals(e.shellMax[s]);
            const auto lo = shellVals(e.shellMin[s]);
            const auto a = shellVals(r1.shellForces[s]);
            const auto b = shellVals(r2.shellForces[s]);
            for (size_t k = 0; k < hi.size(); ++k) {
                envErr = std::max(envErr, std::fabs(hi[k] - std::max(a[k], b[k])));
                envErr = std::max(envErr, std::fabs(lo[k] - std::min(a[k], b[k])));
            }
        }
    }
    addRow("Combination/envelope", "envelope covers shell resultants, not beam-only",
           "component max/min across two shell cases",
           "max shell envelope abs diff", envErr, 1e-9, envErr < 1e-9);
}

void testSlopedSelfWeight() {
    Material mat(200000.0, 76923.076923, 7850.0);
    Section sec = Section::Rectangular(220.0, 360.0);
    FrameModel m;
    m.materials.push_back(mat);   // index 0
    m.sections.push_back(sec);    // index 0
    Node n0(0, 0.0, 0.0, 0.0); n0.fixAll();
    Node n1(1, 3000.0, 0.0, 4000.0);
    m.nodes = { n0, n1 };
    m.members = { Member(0, 0, 1, 0, 0) };

    const real g = 9810.0;
    addSelfWeight(m, g);
    const real L = norm(n1.pos - n0.pos);
    const real wg = mat.rho * sec.A * g * 1.0e-12;
    const real expectedW = wg * L;

    real magErr = 1.0;
    if (m.memberUDLs.size() == 1) {
        const Vec3 w = m.memberUDLs[0].w_local;
        magErr = relErr(norm(w), wg);
    }
    addRow("Self-weight", "sloped member gravity rotation",
           "|R*w_global| must preserve distributed-load magnitude",
           "relative magnitude error", magErr, 1e-12, magErr < 1e-12);

    const SolveResult r = solve(m);
    real rz = 0.0;
    for (size_t i = Uz; i < r.reactions.size(); i += DOF_PER_NODE) rz += r.reactions[i];
    const real rzErr = relErr(rz, expectedW);
    addRow("Self-weight", "global equilibrium on an inclined cantilever",
           "sum Rz equals rho*A*g*L after local-axis load conversion",
           "relative vertical reaction error", rzErr, 1e-10, !r.singular && rzErr < 1e-10);
}

void testPreparedSystemReuse() {
    Material mat(200000.0, 76923.076923, 7850.0);
    Section sec = Section::Rectangular(250.0, 400.0);
    FrameModel m;
    fixtures::clampedSettlement(m, 6000.0, 1.0, mat, sec);
    PreparedSystem ps = assembleAndFactor(m);

    real worst = 0.0;
    for (real delta : { 0.25, 1.0, 2.0, -0.75 }) {
        m.nodes[2].prescribed[Uz] = -delta;
        const SolveResult reuse = solveLoad(ps, m);
        const SolveResult fresh = solve(m);
        worst = std::max(worst, maxDiff(reuse, fresh));
    }
    addRow("Prepared solve", "same factorization, multiple settlement values",
           "solveLoad(prepared) must match fresh solve for changed prescribed values",
           "worst abs result diff", worst, 1e-10, worst < 1e-10);
}

void testInfluenceLine() {
    Material mat(200000.0, 76923.076923, 7850.0);
    Section sec = Section::Rectangular(180.0, 300.0);
    const int n = 16;
    const real L = 8000.0;
    FrameModel m;
    fixtures::simplySupportedBeamN(m, n, L, mat, sec);
    PreparedSystem ps = assembleAndFactor(m);

    NodalLoad saved; saved.node = 3; saved.comp[Uz] = -123.0;
    m.nodalLoads = { saved };
    std::vector<NodeId> ids;
    for (int i = 0; i <= n; ++i) ids.push_back(i);
    const std::vector<real> il = reactionInfluenceLine(ps, m, ids, 0, Uz);

    real maxErr = 0.0;
    for (int i = 0; i <= n && static_cast<size_t>(i) < il.size(); ++i) {
        const real x = L * real(i) / real(n);
        const real exact = (L - x) / L;
        maxErr = std::max(maxErr, std::fabs(il[static_cast<size_t>(i)] - exact));
    }
    addRow("Influence line", "16-position moving unit load, not one midpoint",
           "Muller-Breslau/static equilibrium: R_A=(L-x)/L",
           "max ordinate abs error", maxErr, 1e-11, maxErr < 1e-11);

    const bool restored = m.nodalLoads.size() == 1 && m.nodalLoads[0].node == saved.node &&
                          std::fabs(m.nodalLoads[0].comp[Uz] - saved.comp[Uz]) < 1e-15;
    addRow("Influence line", "API mutation guard",
           "caller nodalLoads restored after marching load",
           "restoration flag", restored ? 0.0 : 1.0, 0.0, restored);
}

void testModalScaling() {
    Section sec = Section::Rectangular(120.0, 240.0);
    FrameModel m;
    buildPureBendingSS(m, 32, 6000.0, 7850.0, sec);
    PreparedSystem ps = assembleAndFactor(m);
    ModalResult mr = solveModal(ps, ModalOptions{ 6 });
    const bool modalOk = !mr.singular && mr.modes.size() >= 3;
    addRow("Modal", "pure vertical-bending modal family",
           "model must return at least first 3 bending modes",
           "modal availability", modalOk ? 0.0 : 1.0, 0.0, modalOk);
    if (!modalOk) return;

    const real r2 = mr.modes[1].omega / mr.modes[0].omega;
    const real r3 = mr.modes[2].omega / mr.modes[0].omega;
    const real ratioErr = std::max(relErr(r2, 4.0), relErr(r3, 9.0));
    addRow("Modal", "higher-mode ratios",
           "simply-supported beam omega_n/omega_1 = n^2",
           "max relative ratio error", ratioErr, 5e-3, ratioErr < 5e-3);

    FrameModel mdense;
    buildPureBendingSS(mdense, 32, 6000.0, 7850.0 * 4.0, sec);
    const ModalResult mrd = solveModal(assembleAndFactor(mdense), ModalOptions{ 1 });
    const real densityErr = (!mrd.singular && !mrd.modes.empty()) ? relErr(mr.modes[0].omega / mrd.modes[0].omega, 2.0) : 1.0;
    addRow("Modal", "density scaling",
           "rho x4 should make omega_1 half",
           "relative scaling error", densityErr, 1e-10, densityErr < 1e-10);

    FrameModel mlong;
    buildPureBendingSS(mlong, 32, 12000.0, 7850.0, sec);
    const ModalResult mrl = solveModal(assembleAndFactor(mlong), ModalOptions{ 1 });
    const real lengthErr = (!mrl.singular && !mrl.modes.empty()) ? relErr(mr.modes[0].omega / mrl.modes[0].omega, 4.0) : 1.0;
    addRow("Modal", "length scaling",
           "L x2 should make Euler-Bernoulli omega_1 quarter",
           "relative scaling error", lengthErr, 1e-10, lengthErr < 1e-10);
}

real bucklingPcr(real L, real pref, bool compression, bool& ok) {
    Material mat(210000.0, 80769.230769, 7850.0);
    Section sec = Section::Rectangular(160.0, 160.0);
    const int n = 20;
    FrameModel m;
    fixtures::simplySupportedBeamN(m, n, L, mat, sec);
    NodalLoad nl;
    nl.node = n;
    nl.comp[Ux] = compression ? -pref : pref;
    m.nodalLoads = { nl };
    PreparedSystem ps = assembleAndFactor(m);
    const BucklingResult br = solveBuckling(ps, m);
    ok = !br.singular;
    return ok ? br.criticalFactor * pref : 0.0;
}

void testBucklingScaling() {
    Material mat(210000.0, 80769.230769, 7850.0);
    Section sec = Section::Rectangular(160.0, 160.0);
    bool ok1 = false, ok2 = false, ok3 = false;
    const real L = 5000.0;
    const real p1 = bucklingPcr(L, 1000.0, true, ok1);
    const real p2 = bucklingPcr(L, 2000.0, true, ok2);
    const real pLong = bucklingPcr(2.0 * L, 1000.0, true, ok3);
    const real euler = kPi * kPi * mat.E * sec.Iz / (L * L);

    const real eulerErr = ok1 ? relErr(p1, euler) : 1.0;
    addRow("Buckling", "Euler oracle at a new length",
           "Pcr = pi^2 E I / L^2, pinned-pinned",
           "relative Pcr error", eulerErr, 1.5e-3, ok1 && eulerErr < 1.5e-3);

    const real loadScaleErr = (ok1 && ok2) ? relErr(p1 / p2, 1.0) : 1.0;
    addRow("Buckling", "reference load scaling",
           "lambda*Pref must be invariant when Pref changes",
           "relative Pcr scale error", loadScaleErr, 1e-10, loadScaleErr < 1e-10);

    const real lengthErr = (ok1 && ok3) ? relErr(p1 / pLong, 4.0) : 1.0;
    addRow("Buckling", "length scaling",
           "Pcr(L) / Pcr(2L) = 4",
           "relative length-scale error", lengthErr, 2e-3, lengthErr < 2e-3);

    bool tensionOk = false;
    (void)bucklingPcr(L, 1000.0, false, tensionOk);
    addRow("Buckling", "no-compression guard",
           "tension-only reference load should not return a positive buckling factor",
           "unexpected non-singular flag", tensionOk ? 1.0 : 0.0, 0.0, !tensionOk);
}

void testResponseSpectrum() {
    Material mat(210000.0, 80769.230769, 7850.0);
    Section sec = Section::Rectangular(90.0, 150.0);
    FrameModel m;
    fixtures::cantileverBeamN(m, 16, 4000.0, mat, sec);
    PreparedSystem ps = assembleAndFactor(m);
    const ModalResult mr4 = solveModal(ps, ModalOptions{ 4 });
    const ModalResult mr20 = solveModal(ps, ModalOptions{ 20 });
    Spectrum sp1; sp1.T = { 0.0, 10.0 }; sp1.Sa = { 9810.0, 9810.0 };
    Spectrum sp2; sp2.T = { 0.0, 10.0 }; sp2.Sa = { 19620.0, 19620.0 };

    const ResponseSpectrumResult r4 = solveResponseSpectrum(ps, mr4, sp1, Uz, SpectrumCombo::SRSS, 0.05);
    const ResponseSpectrumResult r20 = solveResponseSpectrum(ps, mr20, sp1, Uz, SpectrumCombo::SRSS, 0.05);
    const ResponseSpectrumResult r20x2 = solveResponseSpectrum(ps, mr20, sp2, Uz, SpectrumCombo::SRSS, 0.05);

    const real s4 = sumOf(r4.effMass);
    const real s20 = sumOf(r20.effMass);
    const bool massOk = !r4.singular && !r20.singular && s20 + 1e-8 >= s4 && s20 <= r20.totalMass * (1.0 + 1e-8);
    const real massMetric = (r20.totalMass > 0.0) ? s20 / r20.totalMass : 0.0;
    addRow("Response spectrum", "mode truncation and mass bound",
           "participating mass is monotone with more modes and <= total mass",
           "20-mode participating mass ratio", massMetric, 1.0, massOk);

    const real shearScaleErr = (!r20.singular && !r20x2.singular && r20.baseShear > 0.0)
        ? relErr(r20x2.baseShear / r20.baseShear, 2.0) : 1.0;
    const real dispScaleErr = (!r20.singular && !r20x2.singular && maxAbs(r20.u) > 0.0)
        ? relErr(maxAbs(r20x2.u) / maxAbs(r20.u), 2.0) : 1.0;
    const real scaleErr = std::max(shearScaleErr, dispScaleErr);
    addRow("Response spectrum", "Sa amplitude scaling",
           "doubling spectrum Sa doubles SRSS base shear and displacement",
           "max relative scaling error", scaleErr, 1e-12, scaleErr < 1e-12);

    Spectrum interp; interp.T = { 0.0, 1.0, 3.0 }; interp.Sa = { 10.0, 20.0, 60.0 };
    const real interpErr = std::fabs(interp.at(2.0) - 40.0);
    addRow("Response spectrum", "spectrum interpolation",
           "piecewise-linear Sa(T) at T=2.0 between 1.0 and 3.0",
           "absolute interpolation error", interpErr, 1e-14, interpErr < 1e-14);
}

void testModalDynamics() {
    Material mat(210000.0, 80769.230769, 7850.0);
    Section sec = Section::Rectangular(120.0, 240.0);
    const int n = 12;
    FrameModel m;
    fixtures::cantileverBeamN(m, n, 3500.0, mat, sec);
    NodalLoad nl; nl.node = n; nl.comp[Uz] = -1000.0; m.nodalLoads = { nl };
    PreparedSystem ps = assembleAndFactor(m);
    const ModalResult mr = solveModal(ps, ModalOptions{ 20 });
    const bool modalOk = !mr.singular && !mr.modes.empty();
    if (!modalOk) {
        addRow("Modal dynamics", "setup", "modal basis must be available", "setup flag", 1.0, 0.0, false);
        return;
    }

    const real T1 = 2.0 * kPi / mr.modes[0].omega;
    ModalDynamicsOptions opt; opt.zeta = 0.02; opt.dt = T1 / 250.0; opt.nSteps = 500;
    const ModalTimeHistory th1 = solveModalStepResponse(ps, m, mr, opt);
    FrameModel m2 = m;
    m2.nodalLoads[0].comp[Uz] = -2000.0;
    const ModalTimeHistory th2 = solveModalStepResponse(ps, m2, mr, opt);
    const int tipUz = gdof(n, Uz);
    const real p1 = peakDof(th1, tipUz);
    const real p2 = peakDof(th2, tipUz);
    const real loadScaleErr = (p1 > 0.0 && p2 > 0.0) ? relErr(p2 / p1, 2.0) : 1.0;
    addRow("Modal dynamics", "load scaling",
           "linear modal Newmark response doubles under doubled step load",
           "relative peak scaling error", loadScaleErr, 1e-12, loadScaleErr < 1e-12);

    ModalDynamicsOptions und = opt; und.zeta = 0.0;
    ModalDynamicsOptions damp = opt; damp.zeta = 0.25;
    const ModalTimeHistory thu = solveModalStepResponse(ps, m, mr, und);
    const ModalTimeHistory thd = solveModalStepResponse(ps, m, mr, damp);
    const real peakUnd = peakDof(thu, tipUz);
    const real peakDamp = peakDof(thd, tipUz);
    const bool dampingOk = peakDamp < peakUnd;
    addRow("Modal dynamics", "damping effect",
           "higher damping must reduce step-response peak over the same window",
           "peak damped / peak undamped", peakUnd > 0.0 ? peakDamp / peakUnd : 1.0, 1.0, dampingOk);

    ModalDynamicsOptions bad; bad.dt = 0.0; bad.nSteps = 10;
    const ModalTimeHistory badTh = solveModalStepResponse(ps, m, mr, bad);
    addRow("Modal dynamics", "invalid input guard",
           "dt <= 0 must return singular/bad modal-step input",
           "unexpected non-singular flag", badTh.singular ? 0.0 : 1.0, 0.0, badTh.singular);
}

void testAxialTorsionalModes() {
    // Closes the localMass12 axial/torsional coverage gap: F22 only checks BENDING modes.
    // A clamped-free rod restrained to a SINGLE dof family has the continuous-bar natural
    // frequency omega_1 = (pi/2L) * c, with c = sqrt(E/rho) (axial) or c = sqrt(G/rho)
    // (torsion of a CIRCULAR bar, where J = Ip = Iy+Iz so the torsional wave speed
    // sqrt(G*J/(rho*Ip)) reduces to sqrt(G/rho)). rho uses the kg/m^3 -> tonne/mm^3 bridge
    // (x1e-12), the same one the mass matrix applies. This validates the lumped/consistent
    // axial (ma) and torsional (mt = rho*Ip*L/6) terms of localMass12, never exercised by F22.
    const real E = 210000.0, nu = 0.30, rho = 7850.0;
    const real G = E / (2.0 * (1.0 + nu));
    const real L = 4000.0;
    const int n = 32;
    const Section sec = Section::Circular(80.0);
    const Material mat(E, G, rho);

    {   // axial: leave only Ux free
        FrameModel m;
        fixtures::cantileverBeamN(m, n, L, mat, sec);
        for (Node& nd : m.nodes) { nd.fixed[Uy] = nd.fixed[Uz] = nd.fixed[Rx] = nd.fixed[Ry] = nd.fixed[Rz] = true; }
        const ModalResult mr = solveModal(assembleAndFactor(m), ModalOptions{ 1 });
        const real wExact = (kPi / (2.0 * L)) * std::sqrt(E / (rho * 1e-12));
        const real err = (!mr.singular && !mr.modes.empty()) ? relErr(mr.modes[0].omega, wExact) : 1.0;
        addRow("Modal (axial)", "axial consistent mass (localMass12 ma term)",
               "clamped-free bar omega_1 = (pi/2L) sqrt(E/rho)",
               "relative omega_1 error", err, 5e-3, err < 5e-3);
    }
    {   // torsion: leave only Rx free (circular section -> J = Ip)
        FrameModel m;
        fixtures::cantileverBeamN(m, n, L, mat, sec);
        for (Node& nd : m.nodes) { nd.fixed[Ux] = nd.fixed[Uy] = nd.fixed[Uz] = nd.fixed[Ry] = nd.fixed[Rz] = true; }
        const ModalResult mr = solveModal(assembleAndFactor(m), ModalOptions{ 1 });
        const real wExact = (kPi / (2.0 * L)) * std::sqrt(G / (rho * 1e-12));
        const real err = (!mr.singular && !mr.modes.empty()) ? relErr(mr.modes[0].omega, wExact) : 1.0;
        addRow("Modal (torsion)", "torsional consistent mass (localMass12 mt = rho*Ip*L/6)",
               "clamped-free circular bar omega_1 = (pi/2L) sqrt(G/rho)",
               "relative omega_1 error", err, 5e-3, err < 5e-3);
    }
}

void testCQC() {
    // Closes the CQC-path coverage gap (every F24 fixture uses SRSS). The correlation
    // coefficient is re-implemented INDEPENDENTLY here and the engine's CQC base shear is
    // checked against sqrt(sum_ij rho_ij V_i V_j) with V_i = effMass_i * Sa(T_i). This runs
    // the engine's CQC double-loop (never touched by the SRSS gate); a wrong rho or a wrong
    // combination loop would mismatch. Also checks the bound SRSS <= CQC <= sum|V_i|.
    const Material mat(210000.0, 80769.230769, 7850.0);
    const Section sec = Section::Rectangular(90.0, 150.0);
    FrameModel m;
    fixtures::cantileverBeamN(m, 16, 4000.0, mat, sec);
    PreparedSystem ps = assembleAndFactor(m);
    const ModalResult mr = solveModal(ps, ModalOptions{ 8 });
    Spectrum sp; sp.T = { 0.0, 10.0 }; sp.Sa = { 9810.0, 9810.0 };
    const real zeta = 0.05;

    const ResponseSpectrumResult cqc  = solveResponseSpectrum(ps, mr, sp, Uz, SpectrumCombo::CQC, zeta);
    const ResponseSpectrumResult srss = solveResponseSpectrum(ps, mr, sp, Uz, SpectrumCombo::SRSS, zeta);
    if (cqc.singular || srss.singular || mr.modes.empty()) {
        addRow("Response spectrum (CQC)", "setup", "CQC result + modes available", "setup flag", 1.0, 0.0, false);
        return;
    }

    auto rhoIndep = [zeta](real wi, real wj) -> real {
        if (wi <= 0 || wj <= 0) return (wi == wj) ? 1.0 : 0.0;
        const real b = wi / wj;
        const real num = 8.0 * zeta * zeta * (1.0 + b) * std::pow(b, 1.5);
        const real den = (1.0 - b * b) * (1.0 - b * b) + 4.0 * zeta * zeta * b * (1.0 + b) * (1.0 + b);
        return den > 0 ? num / den : (wi == wj ? 1.0 : 0.0);
    };
    const int nm = static_cast<int>(mr.modes.size());
    std::vector<real> V(static_cast<size_t>(nm), 0.0);
    real srssSq = 0.0, absSum = 0.0;
    for (int i = 0; i < nm; ++i) {
        const real Ti = 2.0 * kPi / mr.modes[i].omega;
        V[static_cast<size_t>(i)] = cqc.effMass[static_cast<size_t>(i)] * sp.at(Ti);
        srssSq += V[i] * V[i];
        absSum += std::fabs(V[i]);
    }
    real cqcSq = 0.0;
    for (int i = 0; i < nm; ++i)
        for (int j = 0; j < nm; ++j)
            cqcSq += rhoIndep(mr.modes[i].omega, mr.modes[j].omega) * V[i] * V[j];
    const real expected = std::sqrt(std::max(cqcSq, 0.0));
    const real beErr = (expected > 0) ? relErr(cqc.baseShear, expected) : 1.0;
    addRow("Response spectrum (CQC)", "CQC double-loop combination path",
           "independent recompute sqrt(sum_ij rho_ij V_i V_j), V_i = effMass_i*Sa(T_i)",
           "relative base-shear error", beErr, 1e-8, beErr < 1e-8);

    const real srssBe = std::sqrt(std::max(srssSq, 0.0));
    const bool boundOk = cqc.baseShear >= srssBe * (1.0 - 1e-9) && cqc.baseShear <= absSum * (1.0 + 1e-9);
    addRow("Response spectrum (CQC)", "physical bound SRSS <= CQC <= sum|V|",
           "positively-correlated modes: cross terms raise CQC above SRSS",
           "CQC / SRSS base-shear ratio", srssBe > 0 ? cqc.baseShear / srssBe : 0.0, 1.0, boundOk);
}

void testRectBiaxialDC() {
    // Closes the rectangular biaxial-bending D/C gap (F9 only checks the circular resultant).
    // For a rectangle ElasticAllowable uses the worst-corner sum |My|/Wy + |Mz|/Wz, which is
    // the EXACT peak fibre stress at a corner AND a conservative upper bound on the SRSS-style
    // resultant sqrt((My/Wy)^2 + (Mz/Wz)^2). With N=0, sComp == sM (the bending stress).
    const Section sec = Section::Rectangular(220.0, 360.0);
    const Capacity cap = Capacity::make(300.0, 300.0, 180.0);
    const ElasticAllowable screen;

    {   // uniaxial: corner sum degenerates to |My|/Wy exactly
        MemberEndForces f; f.My = 5.0e7;
        const DemandResult d = screen.checkSection(f, sec, cap);
        const real expect = std::fabs(f.My) / sec.Wy();
        const real err = relErr(d.sComp, expect);
        addRow("D/C screen (rect)", "uniaxial bending exactness",
               "sComp (N=0) == |My|/Wy",
               "relative stress error", err, 1e-9, err < 1e-9);
    }
    {   // biaxial: exact worst-corner sum + conservatism vs the resultant
        MemberEndForces f; f.My = 5.0e7; f.Mz = 3.0e7;
        const DemandResult d = screen.checkSection(f, sec, cap);
        const real cornerSum = std::fabs(f.My) / sec.Wy() + std::fabs(f.Mz) / sec.Wz();
        const real err = relErr(d.sComp, cornerSum);
        addRow("D/C screen (rect)", "biaxial worst-corner sum",
               "sComp (N=0) == |My|/Wy + |Mz|/Wz",
               "relative stress error", err, 1e-9, err < 1e-9);

        const real sy = std::fabs(f.My) / sec.Wy(), sz = std::fabs(f.Mz) / sec.Wz();
        const real resultant = std::sqrt(sy * sy + sz * sz);
        const bool conservative = d.sComp >= resultant * (1.0 - 1e-12);
        addRow("D/C screen (rect)", "corner sum is a conservative upper bound",
               "|a|+|b| >= sqrt(a^2+b^2): rectangle screen >= resultant",
               "corner-sum / resultant ratio", resultant > 0 ? d.sComp / resultant : 0.0, 1.0, conservative);
    }
}

void testMITC4SoftMode() {
    // Pins down the disclosed MITC4 element-level trait (see MITC4ShellElement.h header):
    // the local 24x24 stiffness has EXACTLY 6 rigid-body zero modes PLUS one inherent
    // low-energy (near-zero, non-rigid) plate-bending mode, well separated from the first
    // true deformation mode. Built on a REGULAR square facet, so it is not distortion-induced.
    Material mat(30000.0, 0.0, 2500.0);
    mat.nu = 0.30;
    mat.G  = mat.E / (2.0 * (1.0 + mat.nu));
    const real a = 1000.0, t = 100.0;
    FrameModel m;
    m.materials.push_back(mat);   // index 0
    m.nodes.push_back(Node(0, 0, 0, 0));
    m.nodes.push_back(Node(1, a, 0, 0));
    m.nodes.push_back(Node(2, a, a, 0));
    m.nodes.push_back(Node(3, 0, a, 0));
    m.shells.push_back(ShellQuad(0, 0, 1, 2, 3, 0, t));

    MITC4ShellElement el(0);
    std::string why;
    if (!el.prepare(m, SolveOptions{}, why)) {
        addRow("MITC4 element", "local stiffness spectrum", "prepare() must succeed", "prepare flag", 1.0, 0.0, false);
        return;
    }
    const auto& K = el.localKForAudit();
    Eigen::MatrixXd Kd(24, 24);
    for (int i = 0; i < 24; ++i)
        for (int j = 0; j < 24; ++j) Kd(i, j) = K(i, j);
    const Eigen::VectorXd ev = Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(Kd).eigenvalues();  // ascending
    const real emax = ev(23);

    int nRigid = 0;
    for (int i = 0; i < 24; ++i) if (ev(i) / emax < 1e-8) ++nRigid;
    addRow("MITC4 element", "exactly 6 rigid-body zero modes",
           "local 24x24 K spectrum: count of eig/eigmax < 1e-8",
           "rigid-mode count", static_cast<real>(nRigid), 0.0, nRigid == 6);

    const real softRel = ev(6) / emax;                        // first non-rigid mode, relative
    const real sep = (ev(7) > 0) ? ev(6) / ev(7) : 1.0;       // soft mode vs first true deformation mode
    // Non-rigid (not one of the 6 zero modes) yet markedly softer than the first true mode.
    // The exact ratio depends on t/a (this thick 0.1 facet shows ~100x; thin plates show far
    // more), so the gate asserts the qualitative separation, not a geometry-specific magnitude.
    const bool softOk = softRel > 1e-12 && sep < 0.1;
    addRow("MITC4 element", "inherent low-energy plate mode (disclosed, monitored)",
           "mode-7 non-rigid (eig7/eigmax > 1e-12) yet >=10x softer than mode-8",
           "mode-7 / mode-8 ratio", sep, 1e-1, softOk);
}

void testSolveLoadFingerprint() {
    // B1 reuse-validity guard: solveLoad must ACCEPT changed nodal loads / prescribed values
    // (the interactive path) but REJECT a structurally changed model (geometry / UDL / support)
    // instead of silently reusing the stale factorization + baked distributed loads.
    Material mat(200000.0, 76923.076923, 7850.0);
    Section sec = Section::Rectangular(200.0, 300.0);
    FrameModel m;
    fixtures::cantileverBeamN(m, 4, 4000.0, mat, sec);
    NodalLoad nl; nl.node = 4; nl.comp[Uz] = -500.0; m.nodalLoads = { nl };
    PreparedSystem ps = assembleAndFactor(m);

    {   // (a) nodal-load change is allowed (excluded from the fingerprint)
        FrameModel m2 = m; m2.nodalLoads[0].comp[Uz] = -1000.0;
        const SolveResult r = solveLoad(ps, m2);
        addRow("solveLoad guard", "nodal-load change accepted",
               "fingerprint excludes nodal loads (interactive path)",
               "singular flag (want 0)", r.singular ? 1.0 : 0.0, 0.0, !r.singular);
    }
    {   // (b) prescribed-value (settlement) change is allowed
        FrameModel m2 = m; m2.nodes[0].prescribed[Uz] = -0.5;
        const SolveResult r = solveLoad(ps, m2);
        addRow("solveLoad guard", "prescribed-value change accepted",
               "fingerprint excludes prescribed values (settlement path)",
               "singular flag (want 0)", r.singular ? 1.0 : 0.0, 0.0, !r.singular);
    }
    {   // (c) geometry change must be rejected (stale factorization)
        FrameModel m2 = m; m2.nodes[2].pos.x += 100.0;
        const SolveResult r = solveLoad(ps, m2);
        addRow("solveLoad guard", "geometry change rejected",
               "stale factorization must not be reused silently",
               "singular flag (want 1)", r.singular ? 1.0 : 0.0, 0.0, r.singular);
    }
    {   // (d) distributed-load change must be rejected (baked UDL would be stale)
        FrameModel m2 = m; MemberUDL u; u.member = 0; u.w_local = { 0.0, -1.0, 0.0 }; m2.memberUDLs.push_back(u);
        const SolveResult r = solveLoad(ps, m2);
        addRow("solveLoad guard", "distributed-load change rejected",
               "baked UDL would be stale under factorization reuse",
               "singular flag (want 1)", r.singular ? 1.0 : 0.0, 0.0, r.singular);
    }
    {   // (e) material-VALUE change must be rejected: the factorization baked the OLD E into K
        // (this is the section/material gap PROJECT.txt P1 flagged — counts/geometry unchanged,
        // only the referenced Young's modulus differs, yet reuse would be a silent stale solve).
        FrameModel m2 = m; m2.materials[0].E *= 2.0;
        const SolveResult r = solveLoad(ps, m2);
        addRow("solveLoad guard", "material E change rejected",
               "fingerprint must hash material values (stale K otherwise)",
               "singular flag (want 1)", r.singular ? 1.0 : 0.0, 0.0, r.singular);
    }
    {   // (f) section-VALUE change must be rejected: changing Iz alters bending stiffness, so the
        // baked factorization is stale even though node/member/shell counts are identical.
        FrameModel m2 = m; m2.sections[0].Iz *= 2.0;
        const SolveResult r = solveLoad(ps, m2);
        addRow("solveLoad guard", "section Iz change rejected",
               "fingerprint must hash section values (stale K otherwise)",
               "singular flag (want 1)", r.singular ? 1.0 : 0.0, 0.0, r.singular);
    }
}

void testShellCornerMoments() {
    // D3 per-corner bending recovery: corner moments are EXACT on a constant-moment field (all
    // equal the centre) and REVEAL the spread that the single centre value averages out on a
    // varying field (so a design peak = max over corners, not the centre).
    const real Es = 30000.0, nu = 0.3;
    Material smat(Es, Es / (2.0 * (1.0 + nu))); smat.nu = nu;

    {   // (a) constant-curvature patch -> every corner moment equals the centre
        FrameModel m;
        fixtures::platePatchCylindrical(m, 1000.0, 10.0, 0.0, 1e-6, smat);
        const SolveResult r = solve(m);
        real scale = 1e-30, maxDev = 0.0;
        bool any = false;
        for (const auto& sf : r.shellForces) {
            any = true;
            scale = std::max(scale, std::fabs(sf.Mxx));
            for (int k = 0; k < 4; ++k) {
                maxDev = std::max(maxDev, std::fabs(sf.MxxC[k] - sf.Mxx));
                maxDev = std::max(maxDev, std::fabs(sf.MyyC[k] - sf.Myy));
                maxDev = std::max(maxDev, std::fabs(sf.MxyC[k] - sf.Mxy));
            }
        }
        const real rel = maxDev / scale;
        addRow("Shell corner moments", "constant-moment field: corners == centre",
               "patch test: |MxxC[k]-Mxx| ~ 0 on a constant-curvature field",
               "max corner-vs-centre rel dev", rel, 1e-8, any && rel < 1e-8);
    }
    {   // (b) varying field -> corners spread around the centre (info the centre averages out)
        FrameModel m;
        fixtures::squarePlateShell(m, 2000.0, 50.0, 6, 0.01, smat);
        const SolveResult r = solve(m);
        real maxSpread = 0.0;
        bool any = !r.singular && !r.shellForces.empty();
        for (const auto& sf : r.shellForces) {
            real hi = sf.MxxC[0], lo = sf.MxxC[0];
            for (int k = 1; k < 4; ++k) { hi = std::max(hi, sf.MxxC[k]); lo = std::min(lo, sf.MxxC[k]); }
            maxSpread = std::max(maxSpread, hi - lo);
        }
        addRow("Shell corner moments", "varying field: corners spread around centre",
               "non-constant moment field has corner spread > 0 (peak != centre average)",
               "max corner Mxx spread", maxSpread, 0.0, any && maxSpread > 0.0);
    }
}

void testSparseModal() {
    // E2 sparse subspace-iteration eigensolver: the opt-in sparse path must reproduce the dense
    // generalized-eigensolver frequencies (same M-normalized modes), on a model big enough that
    // it is a meaningful cross-check of the iteration.
    Section sec = Section::Rectangular(120.0, 240.0);
    FrameModel m;
    buildPureBendingSS(m, 40, 6000.0, 7850.0, sec);   // 41-node SS beam (vertical bending DOFs)
    PreparedSystem ps = assembleAndFactor(m);

    ModalOptions od; od.numModes = 6; od.useSparseSolver = false;
    ModalOptions os; os.numModes = 6; os.useSparseSolver = true;
    const ModalResult md  = solveModal(ps, od);
    const ModalResult msp = solveModal(ps, os);
    const bool ok = !md.singular && !msp.singular && md.modes.size() >= 6 && msp.modes.size() >= 6;
    if (!ok) {
        addRow("Sparse eigensolver", "setup", "dense + sparse modal both available", "flag", 1.0, 0.0, false);
        return;
    }
    real maxRel = 0.0;
    for (int i = 0; i < 6; ++i)
        maxRel = std::max(maxRel, relErr(msp.modes[i].omega, md.modes[i].omega));
    addRow("Sparse eigensolver", "subspace iteration matches the dense modal solve",
           "first 6 omega: sparse subspace (reusing LDLT) vs dense GeneralizedSelfAdjointEigenSolver",
           "max relative omega error", maxRel, 1e-6, maxRel < 1e-6);
}

void testElementRemoval() {
    // C1 element removal + C5 oracle. A propped cantilever (cantilever M0 + vertical prop M1) is
    // statically INDETERMINATE; deactivating the prop (Member::active=false) must drop it to a
    // DETERMINATE cantilever whose exact tip deflection is the closed form -PL^3/3EI. Also checks
    // the removal-by-flag invariant (active=false == physically omitting the member) and that
    // removing every member at a node yields a mechanism (LDLT singular).
    Material mat(210000.0, 80769.0, 7850.0);
    Section  sec = Section::Rectangular(100.0, 100.0);
    const real L = 2000.0, H = 1000.0, P = 1000.0;
    const real dExp = -P * L * L * L / (3.0 * mat.E * sec.Iy);

    auto buildPropped = [&](FrameModel& m) {
        m = FrameModel{};
        m.materials = { mat }; m.sections = { sec };
        Node n0(0, 0.0, 0.0, 0.0); n0.fixAll();
        Node n1(1,   L, 0.0, 0.0);
        Node n2(2,   L, 0.0,  -H); n2.fixAll();
        m.nodes = { n0, n1, n2 };
        m.members = { Member(0, 0, 1, 0, 0), Member(1, 2, 1, 0, 0) };
        NodalLoad p; p.node = 1; p.comp[Uz] = -P; m.nodalLoads = { p };
    };

    FrameModel mFull; buildPropped(mFull);
    const SolveResult rFull = solve(mFull);   // indeterminate baseline (prop present)

    FrameModel mCut = mFull; mCut.members[1].active = false;   // remove the prop
    const SolveResult rCut = solve(mCut);
    const real dErr = (!rCut.singular) ? relErr(rCut.disp(1, Uz), dExp) : 1.0;
    addRow("Element removal", "deactivated member drops indeterminate -> determinate",
           "remove prop from propped cantilever -> tip deflection = -PL^3/3EI (closed form)",
           "relative tip-deflection error", dErr, 1e-9, !rCut.singular && dErr < 1e-9);
    addRow("Element removal", "inactive member carries no force",
           "removed member's recovered end force must be exactly zero",
           "|N| of removed member", std::fabs(rCut.memberForces[1].endI.N), 1e-9,
           std::fabs(rCut.memberForces[1].endI.N) < 1e-9);

    FrameModel mOmit; buildPropped(mOmit); mOmit.members = { mOmit.members[0] };   // drop prop entirely
    const SolveResult rOmit = solve(mOmit);
    real duMax = 0.0;
    for (size_t k = 0; k < rCut.u.size() && k < rOmit.u.size(); ++k)
        duMax = std::max(duMax, std::fabs(rCut.u[k] - rOmit.u[k]));
    addRow("Element removal", "active=false is identical to omitting the member",
           "displacement field of (member inactive) == (member absent)",
           "max |du| between the two models", duMax, 1e-9, duMax < 1e-9);

    FrameModel mMech = mFull; mMech.members[0].active = false; mMech.members[1].active = false;
    const SolveResult rMech = solve(mMech);
    addRow("Element removal", "removing every member at a node -> mechanism",
           "isolated free node has no stiffness -> LDLT flags singular",
           "singular flag (want 1)", rMech.singular ? 1.0 : 0.0, 0.0, rMech.singular);
}

}  // namespace

int main() {
    std::cout << "# build " << FRAMECORE_BUILD_SHA << " | compiled " << __DATE__ << " " << __TIME__ << "\n";
    testShellCombinationEnvelope();
    testSlopedSelfWeight();
    testPreparedSystemReuse();
    testInfluenceLine();
    testModalScaling();
    testBucklingScaling();
    testResponseSpectrum();
    testModalDynamics();
    testAxialTorsionalModes();
    testCQC();
    testRectBiaxialDC();
    testMITC4SoftMode();
    testSolveLoadFingerprint();
    testShellCornerMoments();
    testSparseModal();
    testElementRemoval();

    int failures = 0;
    std::cout << "Linear-analysis deep audit (post F17-F25 strengthening)\n\n";
    std::cout << "| Major area | Strengthened coverage | Independent oracle / invariant | Metric | Result | Tol | Status |\n";
    std::cout << "|---|---|---|---|---:|---:|---|\n";
    for (const AuditRow& r : g_rows) {
        if (!r.ok) ++failures;
        std::cout << "| " << r.area
                  << " | " << r.strengthened
                  << " | " << r.oracle
                  << " | " << r.metric
                  << " | " << sci(r.value)
                  << " | " << sci(r.tol)
                  << " | " << (r.ok ? "PASS" : "FAIL") << " |\n";
    }
    std::cout << "\n" << (failures == 0 ? "PASS" : "FAIL")
              << " failures=" << failures
              << " checks=" << g_rows.size() << "\n";
    return failures == 0 ? 0 : 1;
}
