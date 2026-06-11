#include "FrameCore/BucklingAnalysis.h"
#include "FrameCore/PDeltaAnalysis.h"
#include "FrameCore/TensionOnly.h"
#include "FrameCore/SizeOpt.h"
#include "FrameCore/Topology.h"
#include "FrameCore/Reanalysis.h"
#include "FrameCore/Combination.h"
#include "FrameCore/FrameSolver.h"
#include "FrameCore/InfluenceLine.h"
#include "FrameCore/ModalAnalysis.h"
#include "FrameCore/ModalDynamics.h"
#include "FrameCore/ResponseSpectrum.h"
#include "FrameCore/SelfWeight.h"
#include "FrameCore/ElasticAllowable.h"
#include "FrameCore/Connectivity.h"
#include "FrameCore/Collapse.h"
#include "FrameTestFixtures.h"
#include "MITC4ShellElement.h"   // Private seam: element + Eigen types for the element-level audit
#include "PreparedSystemImpl.h"  // S2: assembled K + free map + prepared elements (mass)
#include "FragmentMomentum.h"    // S2: fragment momentum extraction (shared with the driver)
#include "FrameCore/DynamicCollapse.h"

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

    // S1: the sparse subspace-iteration path (forced via denseThreshold<=0) must reproduce the
    // dense GeneralizedSelfAdjointEigenSolver critical factor to tolerance (iterative, NOT bit-
    // identical). Rectangular section -> non-degenerate planes -> a clean smallest eigenvalue.
    {
        Section secR = Section::Rectangular(120.0, 200.0);
        FrameModel m;
        fixtures::simplySupportedBeamN(m, 20, L, mat, secR);
        NodalLoad nl; nl.node = 20; nl.comp[Ux] = -1000.0; m.nodalLoads = { nl };
        PreparedSystem ps = assembleAndFactor(m);
        const BucklingResult d = solveBuckling(ps, m);
        BucklingOptions opt; opt.denseThreshold = 0;             // force the sparse path
        const BucklingResult s = solveBuckling(ps, m, opt);
        const bool okSD = !d.singular && !s.singular;
        const real sdErr = okSD ? relErr(s.criticalFactor, d.criticalFactor) : 1.0;
        addRow("Buckling", "sparse path agrees with dense",
               "subspaceSmallest critical factor == dense eigensolver (forced sparse)",
               "relative sparse-vs-dense error", sdErr, 1e-6, okSD && sdErr < 1e-6);
    }
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

void testShellS8() {
    // S8: opt-in QM6 incompatible-mode membrane (8a) + DKQ discrete-Kirchhoff thin plate (8b).
    const real Es = 30000.0, nu = 0.3;
    Material smat(Es, Es / (2.0 * (1.0 + nu))); smat.nu = nu;
    SolveOptions qm6;  qm6.useIncompatibleMembrane = true;
    SolveOptions dkq;  dkq.useDKQPlate = true;

    // (a) QM6 still reproduces the distorted constant-strain membrane patch (centre-Jacobian).
    {
        const real a = 1000.0, t = 10.0, gx = 1e-4, f = Es / (1.0 - nu * nu);
        FrameModel m; fixtures::membranePatch(m, a, t, 0.4, gx, smat);
        const SolveResult r = solve(m, qm6);
        const real NxxExp = t * f * gx;
        real e = 0; for (const auto& sf : r.shellForces) e = std::max(e, std::fabs(sf.Nxx - NxxExp));
        const real rel = e / std::fabs(NxxExp);
        addRow("Shell S8 QM6", "QM6 passes the distorted membrane patch",
               "constant-strain patch: Nxx == t f gx to machine precision (centre-Jacobian)",
               "max |Nxx - exact| rel", rel, 1e-8, !r.singular && rel < 1e-8);
    }
    // (b) QM6 releases the in-plane membrane lock that cripples Q4 on a slender cantilever.
    {
        const real L = 100.0, H = 10.0, t = 1.0, P = 1.0;
        const real I = t * H * H * H / 12.0, dEB = P * L * L * L / (3.0 * Es * I);
        auto tip = [&](const SolveOptions& o) -> real {
            FrameModel m; fixtures::slenderMembraneCantilever(m, L, H, t, 4, 1, P, smat);
            const SolveResult r = solve(m, o);
            return r.singular ? 0.0 : std::fabs(r.disp(1 * 5 + 4, Uy));
        };
        const real dQ4 = tip(SolveOptions{}), dQM6 = tip(qm6);
        addRow("Shell S8 QM6", "Q4 membrane LOCKS on slender in-plane bending",
               "Q4 tip < 60% of Euler-Bernoulli (parasitic shear lock)",
               "dQ4 / dEB", dQ4 / dEB, 0.60, dQ4 > 0 && dQ4 < 0.60 * dEB);
        addRow("Shell S8 QM6", "QM6 RELEASES the membrane lock",
               "QM6 tip within 15% of Euler-Bernoulli AND strictly better than Q4",
               "|dQM6/dEB - 1|", std::fabs(dQM6 / dEB - 1.0), 0.15,
               std::fabs(dQM6 - dEB) < 0.15 * dEB && std::fabs(dQM6 - dEB) < std::fabs(dQ4 - dEB));
    }
    // (c) DKQ reproduces the distorted constant-curvature patch + is Kirchhoff (zero shear).
    {
        const real a = 1000.0, t = 10.0, cc = 1e-6;
        const real Dfac = Es * t * t * t / (12.0 * (1.0 - nu * nu));
        FrameModel m; fixtures::platePatchCylindrical(m, a, t, 0.4, cc, smat);
        const SolveResult r = solve(m, dkq);
        const real MxxExp = -Dfac * cc;
        real eM = 0, mQ = 0;
        for (const auto& sf : r.shellForces) {
            eM = std::max(eM, std::fabs(sf.Mxx - MxxExp));
            mQ = std::max(mQ, std::max(std::fabs(sf.Qx), std::fabs(sf.Qy)));
        }
        addRow("Shell S8 DKQ", "DKQ passes the distorted curvature patch",
               "constant-curvature patch: Mxx == -D c to machine precision",
               "max |Mxx - exact| rel", eM / std::fabs(MxxExp), 1e-8, !r.singular && eM < 1e-8 * std::fabs(MxxExp));
        addRow("Shell S8 DKQ", "DKQ is Kirchhoff (no transverse shear)",
               "recovered Qx, Qy are exactly zero", "max |Q|", mQ, 1e-12, mQ < 1e-12);
    }
    // (d) DKQ converges to the Kirchhoff thin-plate centre deflection.
    {
        const real a = 1000.0, t = 5.0, q = 0.01;
        const real D = Es * t * t * t / (12.0 * (1.0 - nu * nu));
        const real wc = 0.00406 * q * a * a * a * a / D;
        FrameModel m; fixtures::squarePlateShell(m, a, t, 16, q, smat);
        const SolveResult r = solve(m, dkq);
        const real w = r.singular ? 1e30 : std::fabs(r.disp(8 * 17 + 8, Uz));
        addRow("Shell S8 DKQ", "DKQ converges to the Kirchhoff thin-plate deflection",
               "SS thin plate centre w == 0.00406 q a^4 / D within 2%",
               "rel err vs Kirchhoff", relErr(w, wc), 2e-2, relErr(w, wc) < 2e-2);
    }
    // (e) opt-in orthogonality: on a FLAT plate the in-plane and out-of-plane DOF decouple, so
    // QM6 (membrane-only) leaves pure-bending Uz bit-for-bit, and DKQ (plate-only) leaves the
    // in-plane membrane Nxx bit-for-bit -- each opt-in touches ONLY its own block.
    {
        const real a = 1000.0, t = 10.0, q = 0.01;
        FrameModel m; fixtures::clampedPlateShell(m, a, t, 12, q, smat);
        const real wBase = std::fabs(solve(m).disp(6 * 13 + 6, Uz));
        const real wQM6  = std::fabs(solve(m, qm6).disp(6 * 13 + 6, Uz));
        addRow("Shell S8 default", "QM6 (membrane-only) leaves pure bending unchanged",
               "flat plate decouples in/out-of-plane -> QM6 does not touch Uz",
               "|wQM6 - wBase| rel", std::fabs(wQM6 - wBase) / std::max<real>(1e-30, wBase), 1e-12,
               wBase > 0 && std::fabs(wQM6 - wBase) < 1e-12 * wBase);
    }
    {
        const real a = 1000.0, t = 10.0, gx = 1e-4, f = Es / (1.0 - nu * nu);
        const real NxxExp = t * f * gx;
        FrameModel m; fixtures::membranePatch(m, a, t, 0.0, gx, smat);
        real nBase = 0, nDKQ = 0;
        for (const auto& sf : solve(m).shellForces)      nBase = std::max(nBase, std::fabs(sf.Nxx));
        for (const auto& sf : solve(m, dkq).shellForces) nDKQ  = std::max(nDKQ,  std::fabs(sf.Nxx));
        addRow("Shell S8 default", "DKQ (plate-only) leaves pure membrane unchanged",
               "membrane patch (out-of-plane restrained) -> DKQ does not touch Nxx",
               "|nDKQ - nBase| rel", std::fabs(nDKQ - nBase) / NxxExp, 1e-12,
               std::fabs(nDKQ - nBase) < 1e-12 * NxxExp);
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

void testShellRemoval() {
    // 3a shell element removal (ShellQuad::active), the facet mirror of C1. The strongest
    // invariant is the MODAL one: deactivating a facet must equal physically omitting it on
    // BOTH sides of the generalized eigenproblem (stiffness AND mass assembly skip it), so the
    // frequencies of the two models must agree to round-off. Plus the reuse fingerprint must
    // treat a flipped shell.active as a structural change (stale-factor rejection).
    const real Es = 30000.0, nu = 0.3;
    Material smat(Es, Es / (2.0 * (1.0 + nu)), 2500.0); smat.nu = nu;
    const real e = 500.0, t = 10.0;

    auto buildPlate = [&](FrameModel& m) {
        m = FrameModel{};
        m.materials = { smat };
        for (int j = 0; j <= 2; ++j)
            for (int i = 0; i <= 2; ++i) {
                Node n(j * 3 + i, i * e, j * e, 0.0);
                if (i == 0) n.fixAll();          // clamped edge x=0
                m.nodes.push_back(n);
            }
        m.shells = { ShellQuad(0, 0, 1, 4, 3, 0, t), ShellQuad(1, 1, 2, 5, 4, 0, t),
                     ShellQuad(2, 3, 4, 7, 6, 0, t), ShellQuad(3, 4, 5, 8, 7, 0, t) };
    };

    {   // (a) modal mirror: (q0 inactive) == (q0 absent) for the first 3 frequencies
        FrameModel mCut; buildPlate(mCut); mCut.shells[0].active = false;
        FrameModel mOmit; buildPlate(mOmit); mOmit.shells.erase(mOmit.shells.begin());
        const PreparedSystem pc = assembleAndFactor(mCut);
        const PreparedSystem po = assembleAndFactor(mOmit);
        const ModalResult rc = solveModal(pc, ModalOptions{ 3 });
        const ModalResult ro = solveModal(po, ModalOptions{ 3 });
        const bool ok = !rc.singular && !ro.singular && rc.modes.size() >= 3 && ro.modes.size() >= 3;
        real maxRel = ok ? 0.0 : 1.0;
        if (ok)
            for (int i = 0; i < 3; ++i)
                maxRel = std::max(maxRel, relErr(rc.modes[i].omega, ro.modes[i].omega));
        addRow("Shell removal", "inactive facet leaves K AND M (modal mirror of C1)",
               "first 3 omega of (facet inactive) == (facet absent): both assemblies skip it",
               "max relative omega diff", maxRel, 1e-14, ok && maxRel < 1e-14);
    }
    {   // (b) flipped shell.active must reject a stale factorization reuse
        FrameModel m; buildPlate(m);
        const PreparedSystem ps = assembleAndFactor(m);
        m.shells[0].active = false;
        const SolveResult r = solveLoad(ps, m);
        addRow("solveLoad guard", "shell active change rejected",
               "fingerprint must hash shell.active (silent stale solve otherwise)",
               "singular flag (want 1)", r.singular ? 1.0 : 0.0, 0.0, r.singular);
    }
}

void testConnectivity() {
    // 3b fragment mass properties. (a) COMPOSITE inertia oracle: the F29 gate proves single
    // pieces; here an L of two rods exercises the parallel-axis aggregation across pieces
    // (the part a single-piece test cannot see). Hand values for two L=1000 rods of mass m
    // meeting at a right angle: com=(750,250), Ixx=Iyy=5mL^2/24, Izz=5mL^2/12=Ixx+Iyy (planar
    // consistency), Ixy=-mL^2/8 (tensor matrix entry). (b) DETERMINISM: shuffling the caller's
    // members/shells vectors must leave every cluster bit-identical (id-sorted accumulation),
    // because the collapse driver's removal sequence must not depend on storage order.
    Material mat(210000.0, 80769.0, 7850.0);
    Section  sec = Section::Rectangular(100.0, 100.0);
    const real L = 1000.0;
    const real m = 7850.0 * 1e-12 * sec.A * L;   // 0.0785 tonne per rod

    {   // (a) L-shape composite inertia
        FrameModel mod; mod.materials = { mat }; mod.sections = { sec };
        mod.nodes = { Node(0, 0, 0, 0), Node(1, L, 0, 0), Node(2, L, L, 0) };
        mod.members = { Member(0, 0, 1, 0, 0), Member(1, 1, 2, 0, 0) };
        const ConnectivityResult c = analyzeConnectivity(mod);
        bool ok = c.valid && c.detached.size() == 1;
        real maxRel = 1.0;
        if (ok) {
            const FragmentCluster& f = c.detached[0];
            const real IxxE = 5.0 * m * L * L / 24.0, IzzE = 5.0 * m * L * L / 12.0, IxyE = -m * L * L / 8.0;
            maxRel = std::max({ relErr(f.inertia[0], IxxE), relErr(f.inertia[1], IxxE),
                                relErr(f.inertia[2], IzzE), relErr(f.inertia[3], IxyE),
                                relErr(f.com.x, 750.0), relErr(f.com.y, 250.0),
                                relErr(f.inertia[2], f.inertia[0] + f.inertia[1]) });
            ok = maxRel < 1e-12;
        }
        addRow("Connectivity", "composite inertia via parallel-axis aggregation",
               "L of two rods: com=(750,250), Ixx=Iyy=5mL^2/24, Izz=Ixx+Iyy, Ixy=-mL^2/8 (hand)",
               "max relative error vs hand values", maxRel, 1e-12, ok);
    }
    {   // (b) shuffle determinism (mixed members + shells, two detached clusters + one grounded)
        Material smat(30000.0, 11538.46, 2500.0); smat.nu = 0.3;
        auto build = [&](FrameModel& mod, bool shuffled) {
            mod = FrameModel{};
            mod.materials = { mat, smat }; mod.sections = { sec };
            Node g0(0, 0, 0, 0); g0.fixAll();
            mod.nodes = { g0, Node(1, 0, 0, 1000),                        // grounded tower
                          Node(2, 3000, 0, 0), Node(3, 4000, 0, 0), Node(4, 4000, 1000, 0),   // L fragment
                          Node(5, 8000, 0, 0), Node(6, 9000, 0, 0), Node(7, 9000, 1000, 0), Node(8, 8000, 1000, 0) };
            mod.members = { Member(0, 0, 1, 0, 0), Member(1, 2, 3, 0, 0), Member(2, 3, 4, 0, 0) };
            mod.shells  = { ShellQuad(0, 5, 6, 7, 8, 1, 10.0) };          // free plate fragment
            if (shuffled) {
                std::swap(mod.members[0], mod.members[2]);                // storage order changes,
                std::reverse(mod.nodes.begin() + 2, mod.nodes.end());     // ids untouched
            }
        };
        FrameModel mA; build(mA, false);
        FrameModel mB; build(mB, true);
        const ConnectivityResult cA = analyzeConnectivity(mA);
        const ConnectivityResult cB = analyzeConnectivity(mB);
        bool ok = cA.valid && cB.valid && cA.detached.size() == 2 && cB.detached.size() == 2 &&
                  cA.groundedComponents == 1 && cB.groundedComponents == 1;
        real maxAbs = ok ? 0.0 : 1.0;
        if (ok) {
            for (size_t k = 0; k < 2; ++k) {
                const FragmentCluster& a = cA.detached[k];
                const FragmentCluster& b = cB.detached[k];
                ok = ok && a.nodes == b.nodes && a.members == b.members && a.shells == b.shells;
                maxAbs = std::max(maxAbs, std::fabs(a.mass - b.mass));
                maxAbs = std::max({ maxAbs, std::fabs(a.com.x - b.com.x), std::fabs(a.com.y - b.com.y),
                                    std::fabs(a.com.z - b.com.z) });
                for (int q = 0; q < 6; ++q) maxAbs = std::max(maxAbs, std::fabs(a.inertia[q] - b.inertia[q]));
            }
            ok = ok && maxAbs == 0.0;
        }
        addRow("Connectivity", "output independent of caller storage order",
               "shuffled members/nodes vectors -> bit-identical clusters (id-sorted accumulation)",
               "max |field difference| (want exact 0)", maxAbs, 0.0, ok);
    }
}

void testCollapseDriver() {
    // 3c progressive-collapse driver invariants beyond the F30 closed-form sequence oracles:
    // (a) dlf scales FORCE loads only -- a settlement-driven model must be BIT-identical under
    //     any dlf (prescribed values are kinematic, not forces);
    // (b) the removal sequence must not depend on the caller's vector storage order;
    // (c) the driver's internal D/C screen must agree with the PUBLIC worstUtilization to the
    //     bit on the baseline step (two implementations of the same math must not drift).
    Material mat(210000.0, 80769.0, 7850.0);
    mat.cap = Capacity::make(300.0, 300.0, 180.0);
    Section sec = Section::Rectangular(100.0, 100.0);

    {   // (a) settlement-only model: dlf=1 vs dlf=2 bit-identical
        FrameModel m; m.materials = { mat }; m.sections = { sec };
        Node n0(0, 0, 0, 0); n0.fixAll();
        Node n1(1, 1000, 0, 0);
        Node n2(2, 2000, 0, 0); n2.fixAll(); n2.prescribed[Uz] = -5.0;   // settles
        m.nodes = { n0, n1, n2 };
        m.members = { Member(0, 0, 1, 0, 0), Member(1, 1, 2, 0, 0) };
        CollapseOptions c1; c1.dlf = 1.0;
        CollapseOptions c2; c2.dlf = 2.0;
        const CollapseHistory h1 = runProgressiveCollapse(m, c1);
        const CollapseHistory h2 = runProgressiveCollapse(m, c2);
        bool ok = !h1.steps.empty() && h1.steps.size() == h2.steps.size() && h1.steps[0].solved;
        real maxAbs = ok ? 0.0 : 1.0;
        if (ok) {
            maxAbs = std::fabs(h1.steps[0].maxDC - h2.steps[0].maxDC);
            for (size_t k = 0; k < h1.steps[0].u.size(); ++k)
                maxAbs = std::max(maxAbs, std::fabs(h1.steps[0].u[k] - h2.steps[0].u[k]));
            ok = maxAbs == 0.0;
        }
        addRow("Collapse driver", "dlf never scales prescribed displacements",
               "settlement-only model: dlf=1 vs dlf=2 -> bit-identical u and D/C",
               "max |difference| (want exact 0)", maxAbs, 0.0, ok);
    }
    {   // (b) storage-order determinism of the full removal sequence
        Material matProp(210000.0, 80769.0, 7850.0);
        matProp.cap = Capacity::make(5.0, 300.0, 180.0);
        Section secProp = Section::Circular(20.0);
        auto build = [&](FrameModel& m, bool shuffled) {
            m = FrameModel{};
            m.materials = { mat, matProp }; m.sections = { sec, secProp };
            Node n0(0, 0, 0, 0); n0.fixAll();
            Node n3(3, 3000, 0, -1000); n3.fixAll();
            m.nodes = { n0, Node(1, 1500, 0, 0), Node(2, 3000, 0, 0), n3 };
            Member prop(2, 3, 2, 1, 1);
            prop.release[4] = prop.release[5] = prop.release[10] = prop.release[11] = true;
            m.members = { Member(0, 0, 1, 0, 0), Member(1, 1, 2, 0, 0), prop };
            if (shuffled) std::swap(m.members[0], m.members[2]);
            NodalLoad p; p.node = 1; p.comp[Uz] = -40000.0; m.nodalLoads = { p };
        };
        CollapseOptions co; co.dlf = 1.0; co.solve.enableReleases = true;
        FrameModel mA; build(mA, false);
        FrameModel mB; build(mB, true);
        const CollapseHistory hA = runProgressiveCollapse(mA, co);
        const CollapseHistory hB = runProgressiveCollapse(mB, co);
        bool ok = hA.outcome == hB.outcome && hA.steps.size() == hB.steps.size();
        real maxAbs = ok ? 0.0 : 1.0;
        if (ok) {
            for (size_t s = 0; s < hA.steps.size(); ++s) {
                ok = ok && hA.steps[s].removedMembers == hB.steps[s].removedMembers &&
                     hA.steps[s].mode == hB.steps[s].mode;
                maxAbs = std::max(maxAbs, std::fabs(hA.steps[s].maxDC - hB.steps[s].maxDC));
                maxAbs = std::max(maxAbs, std::fabs(hA.steps[s].triggerRatio - hB.steps[s].triggerRatio));
            }
            ok = ok && maxAbs == 0.0;
        }
        addRow("Collapse driver", "removal sequence independent of storage order",
               "shuffled members vector -> identical removal ids, modes, and D/C bits",
               "max |D/C difference| (want exact 0)", maxAbs, 0.0, ok);
    }
    {   // (c) driver screen == public worstUtilization on the baseline step
        FrameModel m; m.materials = { mat }; m.sections = { sec };
        Node n0(0, 0, 0, 0); n0.fixAll();
        m.nodes = { n0, Node(1, 2000, 0, 0) };
        m.members = { Member(0, 0, 1, 0, 0) };
        NodalLoad p; p.node = 1; p.comp[Uz] = -1000.0; m.nodalLoads = { p };
        CollapseOptions co; co.dlf = 1.0;
        const CollapseHistory h = runProgressiveCollapse(m, co);
        const DemandSummary ds = worstUtilization(m, solve(m));
        const bool shape = !h.steps.empty() && h.steps[0].solved && ds.valid;
        const real diff = shape ? std::fabs(h.steps[0].maxDC - ds.maxDC) : 1.0;
        addRow("Collapse driver", "internal D/C screen mirrors worstUtilization",
               "baseline-step maxDC == public worstUtilization maxDC (no drift between the two)",
               "|maxDC difference| (want exact 0)", diff, 0.0, shape && diff == 0.0);
    }
}

void testShellFailureScreen() {
    // 3d shell von Mises screen. (a) FE-vs-handbook: simply-supported square plate under UDL,
    // Navier/Timoshenko centre moment Mx = My = 0.0479 q a^2 (nu = 0.3). With Mx == My the
    // surface von Mises equals the bending stress itself (sx = sy -> vM = sx), so the screen's
    // worst D/C must approach 6*0.0479*q*a^2 / t^2 / vm on a fine mesh (FE-grade tolerance:
    // discretization + corner sampling + Mindlin-vs-Kirchhoff gap). (b) DETERMINISM of the
    // merged member+shell screen: shuffling both vectors must not change the removal history.
    const real Es = 30000.0, nu = 0.3;
    Material smat(Es, Es / (2.0 * (1.0 + nu)), 2500.0); smat.nu = nu;
    smat.cap = Capacity::make(10.0, 10.0, 6.0);   // vm = 10

    {   // (a) Navier plate -- compare the CENTRE-ELEMENT screen to the handbook centre moment.
        // (The plate-WORST D/C deliberately lands elsewhere: the simply-supported corner twist
        // Mxy = 0.0325 q a^2 gives vM = sqrt(3)*6*Mxy/t^2 ~ 3.38 > centre 2.874, so the global
        // worst is corner-governed and mesh-sensitive; the centre element is the converged,
        // handbook-checkable sample.)
        const real a = 1000.0, t = 10.0, q = 0.01;
        const real riskExp = 6.0 * 0.0479 * q * a * a / (t * t) / smat.cap.vm;   // 2.874
        FrameModel m; fixtures::squarePlateShell(m, a, t, 16, q, smat);
        const SolveResult r = solve(m);
        real rel = 1.0;
        if (!r.singular) {
            size_t best = 0; real bestD = 1e30;                // element centroid nearest (a/2, a/2)
            for (size_t s = 0; s < m.shells.size(); ++s) {
                real cx = 0, cy = 0;
                for (int k = 0; k < 4; ++k) {
                    const Vec3& p = m.nodes[(size_t)m.nodeIndex(m.shells[s].n[k])].pos;
                    cx += 0.25 * p.x; cy += 0.25 * p.y;
                }
                const real dd = (cx - a / 2) * (cx - a / 2) + (cy - a / 2) * (cy - a / 2);
                if (dd < bestD) { bestD = dd; best = s; }
            }
            const ShellDemandResult d = checkShellSurface(r.shellForces[best], t, smat.cap);
            rel = relErr(d.risk, riskExp);
        }
        addRow("Shell vM screen", "centre-element D/C matches the plate handbook value",
               "SS square plate UDL: Mx=My=0.0479qa^2 at centre -> vM D/C = 6*0.0479*q*a^2/(t^2*vm)",
               "relative error vs Navier (centre element)", rel, 2e-2, rel < 2e-2);
    }
    {   // (b) merged screen determinism under storage shuffle (members AND shells)
        Material bmat(210000.0, 80769.0, 7850.0);
        bmat.cap = Capacity::make(300.0, 300.0, 180.0);
        Section bsec = Section::Rectangular(100.0, 100.0);
        auto build = [&](FrameModel& m, bool shuffled) {
            m = FrameModel{};
            m.materials = { smat, bmat }; m.sections = { bsec };
            const real e = 500.0, t = 10.0;
            Node n0(0, 0, 0, 0); n0.fixAll();
            Node n3(3, 0, e, 0); n3.fixAll();
            Node n6(6, 0, 2000, 0); n6.fixAll();               // independent grounded beam root
            m.nodes = { n0, Node(1, e, 0, 0), Node(2, 2 * e, 0, 0),
                        n3, Node(4, e, e, 0), Node(5, 2 * e, e, 0),
                        n6, Node(7, 1000, 2000, 0) };
            m.shells = { ShellQuad(0, 0, 1, 4, 3, 0, t), ShellQuad(1, 1, 2, 5, 4, 0, t) };
            m.members = { Member(0, 6, 7, 1, 0) };             // lightly loaded, never governs
            NodalLoad p2; p2.node = 2; p2.comp[Uz] = -500.0;
            NodalLoad p5; p5.node = 5; p5.comp[Uz] = -500.0;
            NodalLoad pb; pb.node = 7; pb.comp[Uz] = -100.0;
            m.nodalLoads = { p2, p5, pb };
            if (shuffled) {
                std::swap(m.shells[0], m.shells[1]);
                std::reverse(m.nodes.begin() + 1, m.nodes.end());
            }
        };
        CollapseOptions co; co.dlf = 1.0;
        FrameModel mA; build(mA, false);
        FrameModel mB; build(mB, true);
        const CollapseHistory hA = runProgressiveCollapse(mA, co);
        const CollapseHistory hB = runProgressiveCollapse(mB, co);
        // NODE reversal renumbers the DOFs, which legitimately changes the LDLT round-off
        // (~1e-11 here) -- so the contract is: the DISCRETE removal history (ids, modes,
        // outcome) must be identical, and the D/C values must agree to solver precision.
        // (Bit-exactness under MEMBER-only shuffles is separately gated in testCollapseDriver.)
        bool ok = hA.outcome == hB.outcome && hA.steps.size() == hB.steps.size();
        real maxAbs = ok ? 0.0 : 1.0;
        if (ok) {
            for (size_t s = 0; s < hA.steps.size(); ++s) {
                ok = ok && hA.steps[s].removedMembers == hB.steps[s].removedMembers &&
                     hA.steps[s].removedShells == hB.steps[s].removedShells &&
                     hA.steps[s].mode == hB.steps[s].mode;
                maxAbs = std::max(maxAbs, std::fabs(hA.steps[s].maxDC - hB.steps[s].maxDC));
            }
            ok = ok && maxAbs < 1e-9;
        }
        addRow("Shell vM screen", "merged member+shell screen: decisions survive renumbering",
               "shuffled shells + reversed nodes -> identical removal history; D/C to solver precision",
               "max |D/C difference|", maxAbs, 1e-9, ok);
    }
}

void testPlasticHingeMechanics() {
    // 4a hinge condensation edge cases beyond the F32 closed-form oracles:
    // (a) FOUR hinges on one member (both ends, both bending axes) is a legal condensation
    //     (each bending plane's released block is (EI/L)[[4,2],[2,4]], non-singular) while
    //     the both-end TORSION release stays rejected by the conditioning gate (regression);
    // (b) a hinge on an INACTIVE member is inert: the element is not assembled, so the model
    //     must solve bit-identically to one without the hinge record.
    Material mat(210000.0, 80769.0, 7850.0);
    Section sec = Section::Rectangular(100.0, 100.0);

    auto buildTwoSpan = [&](FrameModel& m) {
        m = FrameModel{};
        m.materials = { mat }; m.sections = { sec };
        Node n0(0, 0, 0, 0); n0.fixAll();
        Node n2(2, 4000, 0, 0); n2.fixAll();
        m.nodes = { n0, Node(1, 2000, 0, 0), n2 };
        m.members = { Member(0, 0, 1, 0, 0), Member(1, 1, 2, 0, 0) };
        NodalLoad p; p.node = 1; p.comp[Uz] = -1000.0; m.nodalLoads = { p };
    };

    {   // (a) four-hinge member solvable; both-end torsion release still gated
        FrameModel m; buildTwoSpan(m);
        for (int dof : { 4, 5, 10, 11 }) m.hinges.push_back(PlasticHinge{ 0, dof, 0.0 });
        const SolveResult r = solve(m);   // member1 still clamps node 1 -> stable
        bool ok = !r.singular;

        FrameModel mT; buildTwoSpan(mT);
        mT.members[0].release[3] = mT.members[0].release[9] = true;   // both torsional ends
        SolveOptions so; so.enableReleases = true;
        const SolveResult rT = solve(mT, so);
        ok = ok && rT.singular;           // the conditioning gate must still reject this
        addRow("Plastic hinge", "4-hinge member condenses; torsion gate intact",
               "both-end bi-axial bending releases are non-singular blocks; both-end torsion is not",
               "solvable & gated (want 1)", ok ? 1.0 : 0.0, 0.0, ok);
    }
    {   // (b) hinge on an inactive member is inert (bit-identical solve)
        FrameModel mA; buildTwoSpan(mA);
        mA.members.push_back(Member(7, 0, 2, 0, 0));   // an extra member that gets removed
        mA.members.back().active = false;
        FrameModel mB = mA;
        mB.hinges.push_back(PlasticHinge{ 7, 5, 1.0e7 });   // hinge on the removed member
        const SolveResult rA = solve(mA);
        const SolveResult rB = solve(mB);
        real duMax = (rA.singular || rB.singular) ? 1.0 : 0.0;
        for (size_t k = 0; k < rA.u.size() && k < rB.u.size(); ++k)
            duMax = std::max(duMax, std::fabs(rA.u[k] - rB.u[k]));
        addRow("Plastic hinge", "hinge on an inactive member is inert",
               "removed member is never assembled -> its hinge record changes nothing",
               "max |du| (want exact 0)", duMax, 0.0, duMax == 0.0);
    }
}

void testHingeDriverCompat() {
    // 4b driver-mode invariants:
    // (a) plasticHinges=false must behave EXACTLY as before 4b existed (bit-identical history
    //     on a hinge-capable model), and at a load where no event triggers in either mode the
    //     two modes must agree bit-for-bit;
    // (b) ductility is BENDING-only: a hinge-capable member failing in pure TENSION must still
    //     be removed brittly, so the hinge-mode history equals the brittle-mode history on the
    //     F30 hanging-chain fixture even with fy set.
    Material mat(210000.0, 80769.0, 7850.0);
    mat.fy = 300.0;
    mat.cap = Capacity::make(300.0, 300.0, 180.0);
    Section sec = Section::Rectangular(100.0, 100.0);

    {   // (a) sub-threshold fixed-fixed beam: both modes Stable with identical bits
        FrameModel m; m.materials = { mat }; m.sections = { sec };
        Node n0(0, 0, 0, 0); n0.fixAll();
        Node n2(2, 4000, 0, 0); n2.fixAll();
        m.nodes = { n0, Node(1, 2000, 0, 0), n2 };
        m.members = { Member(0, 0, 1, 0, 0), Member(1, 1, 2, 0, 0) };
        MemberUDL u0; u0.member = 0; u0.w_local = { 0.0, -30.0, 0.0 };   // D/C = 0.8 < 1
        MemberUDL u1; u1.member = 1; u1.w_local = { 0.0, -30.0, 0.0 };
        m.memberUDLs = { u0, u1 };
        CollapseOptions cb; cb.dlf = 1.0;
        CollapseOptions ch = cb; ch.plasticHinges = true;
        const CollapseHistory hb = runProgressiveCollapse(m, cb);
        const CollapseHistory hh = runProgressiveCollapse(m, ch);
        const bool shape = hb.outcome == CollapseOutcome::Stable && hh.outcome == CollapseOutcome::Stable &&
                           hb.steps.size() == 1 && hh.steps.size() == 1;
        const real d = shape ? std::fabs(hb.steps[0].maxDC - hh.steps[0].maxDC) : 1.0;
        addRow("Hinge driver", "hinge mode is a strict superset below every threshold",
               "sub-threshold beam: brittle vs hinge mode -> identical Stable step (bit-exact D/C)",
               "|maxDC difference| (want exact 0)", d, 0.0, shape && d == 0.0);
    }
    {   // (b) pure tension stays brittle in hinge mode (ductility is bending-only)
        Section secThin = Section::Rectangular(20.0, 20.0);
        FrameModel m; m.materials = { mat }; m.sections = { sec, secThin };
        Node n0(0, 0, 0, 0); n0.fixAll();
        m.nodes = { n0, Node(1, 0, 0, -1000), Node(2, 0, 0, -2000), Node(3, 0, 0, -3000) };
        m.members = { Member(0, 0, 1, 0, 0), Member(1, 1, 2, 0, 1), Member(2, 2, 3, 0, 0) };
        NodalLoad p; p.node = 3; p.comp[Uz] = -150000.0; m.nodalLoads = { p };
        CollapseOptions cb; cb.dlf = 1.0;
        CollapseOptions ch = cb; ch.plasticHinges = true;
        const CollapseHistory hb = runProgressiveCollapse(m, cb);
        const CollapseHistory hh = runProgressiveCollapse(m, ch);
        bool ok = hb.outcome == hh.outcome && hb.steps.size() == hh.steps.size();
        if (ok)
            for (size_t s = 0; s < hb.steps.size(); ++s)
                ok = ok && hb.steps[s].removedMembers == hh.steps[s].removedMembers &&
                     hb.steps[s].mode == hh.steps[s].mode && hh.steps[s].formedHinges.empty();
        addRow("Hinge driver", "pure tension stays brittle under hinge mode",
               "hanging chain with fy set: hinge-mode history == brittle history, zero hinges",
               "histories identical (want 1)", ok ? 1.0 : 0.0, 0.0, ok);
    }
}

void testSafetyAndMargin() {
    // C3 safety factor + C4 criticality (pivot) margin. A cantilever (root moment PL) has the
    // closed-form worst utilization D/C = (PL/W)/cap.bend and safety factor 1/(D/C). pivotMargin
    // = min/max LDLT pivot ratio: exactly 1 for a single DOF, scale-invariant, 0 for a mechanism.
    Section  sec = Section::Rectangular(100.0, 100.0);
    Material mat(210000.0, 80769.0, 7850.0);
    mat.cap = Capacity::make(300.0, 300.0, 180.0);
    const real L = 2000.0, P = 1000.0;
    const real dcExact = (P * L / sec.Wz()) / mat.cap.bend;

    FrameModel m;
    m.materials = { mat }; m.sections = { sec };
    Node n0(0, 0, 0, 0); n0.fixAll();
    Node n1(1, L, 0, 0);
    m.nodes = { n0, n1 };
    m.members = { Member(0, 0, 1, 0, 0) };
    NodalLoad p; p.node = 1; p.comp[Uz] = -P; m.nodalLoads = { p };
    const SolveResult r = solve(m);
    const DemandSummary ds = worstUtilization(m, r);
    addRow("Safety factor", "structural worst Demand/Capacity (C3)",
           "cantilever closed form D/C = (PL/W)/cap.bend = 0.04",
           "relative D/C error", relErr(ds.maxDC, dcExact), 1e-9, ds.valid && relErr(ds.maxDC, dcExact) < 1e-9);
    addRow("Safety factor", "safety factor = 1/maxDC (C3)",
           "elastic load multiplier to first allowable-stress failure",
           "relative SF error", relErr(ds.safetyFactor, 1.0 / dcExact), 1e-9, relErr(ds.safetyFactor, 1.0 / dcExact) < 1e-9);

    FrameModel ma; ma.materials = { mat }; ma.sections = { sec };
    Node a0(0, 0, 0, 0); a0.fixAll();
    Node a1(1, 1000, 0, 0);
    a1.fixed[Uy] = a1.fixed[Uz] = a1.fixed[Rx] = a1.fixed[Ry] = a1.fixed[Rz] = true;  // free: Ux only
    ma.nodes = { a0, a1 };
    ma.members = { Member(0, 0, 1, 0, 0) };
    NodalLoad pa; pa.node = 1; pa.comp[Ux] = 1000.0; ma.nodalLoads = { pa };
    const real m1 = solve(ma).pivotMargin;
    addRow("Criticality margin", "single-DOF pivotMargin = 1 (C4)",
           "one free DOF -> one LDLT pivot -> min/max ratio is exactly 1",
           "|pivotMargin - 1|", std::fabs(m1 - 1.0), 1e-12, std::fabs(m1 - 1.0) < 1e-12);

    FrameModel mE = m; mE.materials[0].E *= 1000.0; mE.materials[0].G *= 1000.0;
    const real mScaled = solve(mE).pivotMargin;
    addRow("Criticality margin", "pivotMargin scale-invariant (C4)",
           "min/max pivot RATIO unchanged when all stiffness x1000",
           "relative margin change", relErr(mScaled, r.pivotMargin), 1e-12, relErr(mScaled, r.pivotMargin) < 1e-12);

    FrameModel mS = m; mS.members[0].active = false;
    const SolveResult rS = solve(mS);
    addRow("Criticality margin", "mechanism -> pivotMargin 0 (C4)",
           "removing the only member -> singular -> margin floored at 0",
           "singular & margin==0 (want 1)", (rS.singular && rS.pivotMargin == 0.0) ? 1.0 : 0.0, 0.0,
           rS.singular && rS.pivotMargin == 0.0);
}

void testReanalysis() {
    Material mat(210000.0, 80769.230769, 7850.0);
    Section  sec = Section::Rectangular(120.0, 120.0);

    // Propped cantilever: a horizontal beam (n0 fixed - n1) plus a prop (n2 fixed - n1). Statically
    // indeterminate; removing the prop leaves a stable cantilever -> one exact Tier-1 removal.
    auto build = [&]() {
        FrameModel m;
        m.materials.push_back(mat);
        m.sections.push_back(sec);
        Node n0(0, 0.0,    0, 0.0);      n0.fixAll();
        Node n1(1, 3000.0, 0, 0.0);
        Node n2(2, 3000.0, 0, -2000.0);  n2.fixAll();
        m.nodes = { n0, n1, n2 };
        Member beam(0, 0, 1, 0, 0); beam.refVec = Vec3(0, 0, 1);
        Member prop(1, 2, 1, 0, 0); prop.refVec = Vec3(1, 0, 0);
        m.members = { beam, prop };
        NodalLoad nl; nl.node = 1; nl.comp[Uz] = -5000.0; m.nodalLoads = { nl };
        return m;
    };
    auto relU = [](const SolveResult& a, const SolveResult& b) -> real {
        real num = 0, den = 1e-30;
        const size_t n = std::min(a.u.size(), b.u.size());
        for (size_t i = 0; i < n; ++i) { num = std::max(num, std::fabs(a.u[i] - b.u[i])); den = std::max(den, std::fabs(b.u[i])); }
        return num / den;
    };

    FrameModel base = build();
    ReSolveSession s(base);
    const SolveResult fr0 = solve(base);

    // (1) Tier-1 Woodbury removal == fresh assembleAndFactor+solveLoad
    s.setMemberActive(1, false);
    ReanalysisStats st1;
    const SolveResult re1 = s.solve(&st1);
    FrameModel w1 = build(); w1.members[1].active = false;
    const real e1 = relU(re1, solve(w1));
    addRow("Reanalysis", "Tier-1 Woodbury == fresh",
           "remove a member; low-rank update on the baseline factor vs a fresh assembleAndFactor+solveLoad",
           "relative u error", e1, 1e-9, st1.tier == 1 && !re1.singular && e1 < 1e-9);

    // (2) remove + restore returns to baseline (the ladder's signed columns cancel)
    s.setMemberActive(1, true);
    const real eR = relU(s.solve(), fr0);
    addRow("Reanalysis", "remove+restore returns to baseline",
           "removing then restoring a member cancels in the ladder to factorization round-off",
           "relative u drift", eR, 1e-11, eR < 1e-11);

    // (3) capacitance mechanism detection (no connectivity needed)
    FrameModel chain;
    chain.materials.push_back(mat);
    chain.sections.push_back(sec);
    Node c0(0, 0, 0, 0);       c0.fixAll();
    Node c1(1, 0, 0, 1500.0);
    Node c2(2, 0, 0, 3000.0);
    chain.nodes = { c0, c1, c2 };
    chain.members = { Member(0, 0, 1, 0, 0), Member(1, 1, 2, 0, 0) };
    NodalLoad nl; nl.node = 2; nl.comp[Ux] = 1000.0; chain.nodalLoads = { nl };
    ReSolveSession cs(chain);
    cs.setMemberActive(0, false);
    ReanalysisStats sm;
    const SolveResult rm = cs.solve(&sm);
    addRow("Reanalysis", "capacitance mechanism detection",
           "removing the base of a 2-member chain -> capacitance singular -> mechanism flag + singular",
           "flagged & singular (want 1)", (sm.mechanism && rm.singular) ? 1.0 : 0.0, 0.0,
           sm.mechanism && rm.singular);

    // (4) Tier-2 stale-LDLT PCG: a small maxRank forces one removal (rank 6) past Tier-1 into Tier-2.
    {
        ReanalysisOptions o2; o2.maxRank = 5;
        FrameModel b2 = build();
        ReSolveSession s2(b2, o2);
        s2.setMemberActive(1, false);
        ReanalysisStats t2;
        const SolveResult r2 = s2.solve(&t2);
        FrameModel w2 = build(); w2.members[1].active = false;
        const real e2 = relU(r2, solve(w2));
        addRow("Reanalysis", "Tier-2 stale-LDLT PCG == fresh (tolerance)",
               "rank>maxRank -> CG preconditioned by the stale baseline factor; tolerance-level (not bit-identical)",
               "relative u error", e2, 1e-8, t2.tier == 2 && !r2.singular && e2 < 1e-8);
    }
}

// ---------------------------------------------------------------------------
// S2 dynamic collapse: cross-event inheritance == full-system Newmark, momentum
// handoff closure (imposed rigid motion), and full-basis exactness.
// ---------------------------------------------------------------------------
SpMat dcReduceFF(const SpMat& K, const std::vector<int>& fmap, int nf) {
    std::vector<Triplet> t;
    for (int c = 0; c < K.outerSize(); ++c)
        for (SpMat::InnerIterator it(K, c); it; ++it) {
            const int r = it.row();
            if (fmap[(size_t)r] >= 0 && fmap[(size_t)c] >= 0) t.emplace_back(fmap[(size_t)r], fmap[(size_t)c], it.value());
        }
    SpMat R(nf, nf); R.setFromTriplets(t.begin(), t.end()); return R;
}
SpMat dcMassFF(const PreparedSystem::Impl& S) {
    std::vector<Triplet> mt; for (const auto& el : S.elems) el->assembleMass(mt);
    SpMat M(S.N, S.N); M.setFromTriplets(mt.begin(), mt.end());
    return dcReduceFF(M, S.fmap, S.nf);
}
VecX dcReduceVec(const VecX& F, const std::vector<int>& fmap, int nf) {
    VecX f = VecX::Zero(nf);
    for (int g = 0; g < (int)fmap.size(); ++g) if (fmap[(size_t)g] >= 0) f(fmap[(size_t)g]) = F(g);
    return f;
}

void testDynamicCollapse() {
    // ---- Check 1: cross-event modal inheritance == full-system Newmark (full basis) ----
    {
        const real E = 200000.0, Gs = 80000.0, rho = 7850.0, W = 3000.0, Hh = 3000.0, P = 50000.0;
        Section sec = Section::Rectangular(200.0, 200.0);
        Material mat(E, Gs, rho); mat.cap = Capacity::make(1e9, 1e9, 1e9);
        FrameModel m; m.materials = { mat }; m.sections = { sec };
        Node n0(0, 0, 0, 0); n0.fixAll();
        Node n1(1, W, 0, 0); n1.fixAll();
        Node n2(2, 0, 0, Hh); Node n3(3, W, 0, Hh);
        m.nodes   = { n0, n1, n2, n3 };
        m.members = { Member(0, 0, 2, 0, 0), Member(1, 1, 3, 0, 0), Member(2, 2, 3, 0, 0), Member(3, 0, 3, 0, 0) };
        PreparedSystem ps1 = assembleAndFactor(m); const auto& S1 = *ps1.impl;
        const SpMat K1 = dcReduceFF(S1.K, S1.fmap, S1.nf), M1 = dcMassFF(S1);
        VecX Fg = VecX::Zero(S1.N); Fg(gdof(2, Ux)) = P;
        const VecX F1 = dcReduceVec(Fg, S1.fmap, S1.nf);
        FrameModel m2 = m; m2.members[3].active = false;             // remove brace -> stays connected
        PreparedSystem ps2 = assembleAndFactor(m2); const auto& S2 = *ps2.impl;
        const SpMat K2 = dcReduceFF(S2.K, S2.fmap, S2.nf), M2 = dcMassFF(S2);
        const VecX F2 = dcReduceVec(Fg, S2.fmap, S2.nf);

        const MatX K1d = MatX(K1), M1d = MatX(M1), K2d = MatX(K2), M2d = MatX(M2);
        Eigen::GeneralizedSelfAdjointEigenSolver<MatX> g1(K1d, M1d), g2(K2d, M2d);
        const VecX w1 = g1.eigenvalues(), w2 = g2.eigenvalues();
        const MatX P1 = g1.eigenvectors(), P2 = g2.eigenvectors();
        const real T1 = 2 * kPi / std::sqrt(w1(0)); const real dt = T1 / 200.0;
        const int nSteps = 400, eventStep = 160;
        const real a0 = 4.0 / (dt * dt), a2 = 4.0 / dt;

        LDLTSolver h1; h1.compute(SpMat(K1 + a0 * M1));
        LDLTSolver h2; h2.compute(SpMat(K2 + a0 * M2));
        LDLTSolver Ml1; Ml1.compute(M1); LDLTSolver Ml2; Ml2.compute(M2);
        auto physStep = [&](const LDLTSolver& hh, const SpMat& M, const SpMat& K, const VecX& F, VecX& u, VecX& v, VecX& a) {
            (void)K; const VecX Fhat = F + M * (a0 * u + a2 * v + a);
            const VecX un = hh.solve(Fhat);
            const VecX an = a0 * (un - u) - a2 * v - a;
            v = v + (dt / 2.0) * (a + an); u = un; a = an;
        };
        VecX u = VecX::Zero(S1.nf), v = VecX::Zero(S1.nf), a = Ml1.solve(F1);
        std::vector<VecX> traceA; traceA.reserve((size_t)nSteps);
        for (int s = 0; s < nSteps; ++s) {
            if (s == eventStep) a = Ml2.solve(F2 - K2 * u);
            if (s < eventStep) physStep(h1, M1, K1, F1, u, v, a); else physStep(h2, M2, K2, F2, u, v, a);
            traceA.push_back(u);
        }
        auto modalStep = [&](VecX& q, VecX& qd, VecX& qdd, const VecX& f, const VecX& w2c) {
            for (int i = 0; i < (int)q.size(); ++i) {
                const real kh = w2c(i) + a0;
                const real qn = (f(i) + a0 * q(i) + a2 * qd(i) + qdd(i)) / kh;
                const real qddn = a0 * (qn - q(i)) - a2 * qd(i) - qdd(i);
                qd(i) += (dt / 2.0) * (qdd(i) + qddn); q(i) = qn; qdd(i) = qddn;
            }
        };
        const VecX f1 = P1.transpose() * F1, f2 = P2.transpose() * F2;
        VecX q = VecX::Zero(S1.nf), qd = VecX::Zero(S1.nf), qdd = f1;
        real maxErr = 0, maxRef = 0;
        for (int s = 0; s < nSteps; ++s) {
            if (s == eventStep) {
                const VecX ue = P1 * q, ve = P1 * qd;
                q = P2.transpose() * (M2 * ue); qd = P2.transpose() * (M2 * ve); qdd = f2 - w2.cwiseProduct(q);
            }
            modalStep(q, qd, qdd, (s < eventStep ? f1 : f2), (s < eventStep ? w1 : w2));
            const VecX ub = (s < eventStep ? P1 : P2) * q;
            maxErr = std::max(maxErr, (ub - traceA[(size_t)s]).cwiseAbs().maxCoeff());
            maxRef = std::max(maxRef, traceA[(size_t)s].cwiseAbs().maxCoeff());
        }
        const real rel = maxErr / std::max<real>(1e-30, maxRef);
        addRow("Dynamic collapse", "cross-event modal inheritance == full-system Newmark",
               "modal projection q'=Phi'^T M' u (full basis) vs avg-accel Newmark with state carry",
               "relMax u", rel, 1e-8, rel < 1e-8);
    }

    // ---- Check 2 & 3: fragment momentum closure under imposed rigid motion ----
    {
        FrameModel c;
        Material cm(200000.0, 76923.07692307692, 7850.0); cm.cap = Capacity::make(1e9, 1e9, 1e9);
        Section cs = Section::Rectangular(150.0, 150.0);
        c.materials = { cm }; c.sections = { cs };
        for (int k = 0; k <= 3; ++k) { Node n(k, 0, 0, 1000.0 * k); if (k == 0) n.fixAll(); c.nodes.push_back(n); }
        for (int k = 0; k < 3; ++k) { Member mm(k, k, k + 1, 0, 0); mm.refVec = Vec3(1, 0, 0); c.members.push_back(mm); }
        FrameModel work = c; work.members[1].active = false;
        const ConnectivityResult cr = analyzeConnectivity(work);
        const FragmentCluster& fc = cr.detached[0];   // {n2, n3, m2}

        const real v0 = 7.3;
        VecX vT = VecX::Zero(c.dofCount());
        for (NodeId nid : fc.nodes) vT(gdof(work.nodeIndex(nid), Ux)) = v0;
        Vec3 p, L; fragmentMomentum(fc, work, vT, p, L);
        const real px_exp = fc.mass * v0;
        const bool ok2 = relErr(p.x, px_exp) < 1e-10 && std::fabs(p.y) < 1e-8 * px_exp && std::fabs(p.z) < 1e-8 * px_exp;
        addRow("Dynamic collapse", "fragment linear momentum closure",
               "imposed rigid translation: p = m*v0 (FE consistent mass), transverse p ~ 0",
               "p_x vs m*v0", relErr(p.x, px_exp), 1e-10, ok2);

        const real w0 = 0.31;
        VecX vR = VecX::Zero(c.dofCount());
        for (NodeId nid : fc.nodes) {
            const int gi = work.nodeIndex(nid);
            const Vec3 rel = work.nodes[(size_t)gi].pos - fc.com;
            const Vec3 tr = cross(Vec3(0, 1, 0), rel);
            vR(gdof(gi, Ux)) = w0 * tr.x; vR(gdof(gi, Uy)) = w0 * tr.y; vR(gdof(gi, Uz)) = w0 * tr.z;
            vR(gdof(gi, Ry)) = w0;
        }
        Vec3 p3, L3; fragmentMomentum(fc, work, vR, p3, L3);
        const real Ly_exp = fc.inertia[1] * w0;
        addRow("Dynamic collapse", "fragment transverse angular momentum",
               "imposed rigid rotation about com (axis y): L_y = I_yy*w0 (slender-rod cluster form)",
               "L_y vs I_yy*w0", relErr(L3.y, Ly_exp), 1e-2, relErr(L3.y, Ly_exp) < 1e-2);
    }

    // ---- Check 4: full-basis inheritance is exact (zero per-event truncation residual) ----
    {
        const real E = 200000.0, Gs = 80000.0, rho = 7850.0, W = 3000.0, Hh = 3000.0, P = 50000.0;
        Section sec = Section::Rectangular(200.0, 200.0);
        Material strong(E, Gs, rho); strong.cap = Capacity::make(1e9, 1e9, 1e9);
        Material braceMat(E, Gs, rho); braceMat.cap = Capacity::make(0.5, 0.5, 1e9);
        FrameModel m; m.materials = { strong, braceMat }; m.sections = { sec };
        Node n0(0, 0, 0, 0); n0.fixAll();
        Node n1(1, W, 0, 0); n1.fixAll();
        Node n2(2, 0, 0, Hh); Node n3(3, W, 0, Hh);
        m.nodes   = { n0, n1, n2, n3 };
        m.members = { Member(0, 0, 2, 0, 0), Member(1, 1, 3, 0, 0), Member(2, 2, 3, 0, 0), Member(3, 0, 3, 1, 0) };
        NodalLoad nl; nl.node = 2; nl.comp[Ux] = P; m.nodalLoads = { nl };
        DynCollapseOptions opt; opt.dt = 1e-5; opt.maxTime = 40 * opt.dt; opt.basisSize = 200;
        opt.screenEvery = 1; opt.frameStride = 5; opt.removeThreshold = 1.0; opt.maxEvents = 3;
        const DynCollapseHistory h = runDynamicCollapse(m, opt);
        real maxTrunc = 0; for (const auto& ev : h.events) maxTrunc = std::max(maxTrunc, ev.truncationResidual);
        addRow("Dynamic collapse", "full-basis inheritance is exact",
               "runDynamicCollapse basisSize >> nf: per-event projection residual ~ 0",
               "max truncationResidual", maxTrunc, 1e-9, !h.events.empty() && maxTrunc < 1e-9);
    }
}

void testPDelta() {
    // Independent P-Delta checks (S3). Cantilever column under axial P + lateral H; Euler P_cr.
    const real Lc = 6000.0, sideC = 200.0, Hlat = 1000.0, Ec = 200000.0;
    const real Ic = sideC * sideC * sideC * sideC / 12.0;
    const real Pcr = kPi * kPi * Ec * Ic / (4.0 * Lc * Lc);
    Section  pcs = Section::Rectangular(sideC, sideC);
    Material pcm(Ec, Ec / (2.0 * (1.0 + 0.3)), 7850.0); pcm.cap = Capacity::make(1e9, 1e9, 1e9);
    const int nE = 8;

    // ---- Check 1: frozen pseudo-load reuse path == K_T refactor reference (one fixed point) ----
    {
        const real P = 0.5 * Pcr;
        FrameModel m; fixtures::pdeltaColumn(m, nE, Lc, P, Hlat, pcm, pcs);
        PDeltaOptions of;  of.refactorPath = false; of.maxIter = 5000; of.tolU = 1e-13;
        PDeltaOptions orf; orf.refactorPath = true;
        const PDeltaResult rf = runPDelta(m, of);
        const PDeltaResult rr = runPDelta(m, orf);
        const real rel = relErr(rf.finalState.disp(nE, Ux), rr.finalState.disp(nE, Ux));
        addRow("P-Delta", "frozen-reuse pseudo-load iteration vs K_T refactor",
               "two independent second-order paths share one fixed point (cantilever column, P/Pcr=0.5)",
               "frozen vs reference tip sway", rel, 1e-10, rf.converged && rr.converged && rel < 1e-10);
    }

    // ---- Check 2: P = 0 degeneracy -> the second-order state IS the linear state, bit-for-bit ----
    {
        FrameModel m; fixtures::pdeltaColumn(m, nE, Lc, 0.0, Hlat, pcm, pcs);
        const PDeltaResult rf = runPDelta(m, PDeltaOptions{});
        const SolveResult lin = solve(m);
        real maxAbs = 0;
        for (size_t g = 0; g < rf.finalState.u.size() && g < lin.u.size(); ++g)
            maxAbs = std::max(maxAbs, std::fabs(rf.finalState.u[g] - lin.u[g]));
        addRow("P-Delta", "P = 0 degeneracy",
               "no compression -> empty Kg -> runPDelta returns the linear solve bit-for-bit (0 iters)",
               "max |u_pdelta - u_linear|", maxAbs, 0.0, rf.converged && rf.iterations == 0 && maxAbs == 0.0);
    }

    // ---- Check 3: beyond the Euler load both paths flag divergence (no silent wrong answer) ----
    {
        const real P = 1.05 * Pcr;
        FrameModel m; fixtures::pdeltaColumn(m, nE, Lc, P, Hlat, pcm, pcs);
        PDeltaOptions of;  of.refactorPath = false; of.maxIter = 5000; of.tolU = 1e-13;
        PDeltaOptions orf; orf.refactorPath = true;
        const PDeltaResult rf = runPDelta(m, of);
        const PDeltaResult rr = runPDelta(m, orf);
        const bool both = rf.diverged && rr.diverged;
        addRow("P-Delta", "divergence past P_cr",
               "frozen sliding-window detector + reference K_T-not-PD both report instability",
               "both paths diverged (1=yes)", both ? 1.0 : 0.0, 0.5, both);
    }
}

void testTensionOnly() {
    Material tmat(200000.0, 76923.07692307692, 7850.0); tmat.cap = Capacity::make(1e9, 1e9, 1e9);
    Section stocky = Section::Rectangular(250.0, 250.0);
    Section tbrace = Section::Rectangular(60.0, 60.0);

    // ---- Check 1: converged tension-only state == model with the compressed brace OMITTED ----
    {
        FrameModel m; fixtures::xBracedPortal(m, 5.0e4, 0.0, tmat, stocky, tbrace);
        const TensionOnlyResult R = runTensionOnly(m);
        const MemberId slackId = R.slack.empty() ? -1 : R.slack[0];
        FrameModel ref; fixtures::xBracedPortal(ref, 5.0e4, 0.0, tmat, stocky, tbrace);
        for (size_t e = 0; e < ref.members.size(); ++e)
            if (ref.members[e].id == slackId) { ref.members.erase(ref.members.begin() + (long)e); break; }
        const SolveResult rr = solve(ref);
        real maxDiff = 0, maxU = 0;
        for (size_t g = 0; g < R.finalState.u.size() && g < rr.u.size(); ++g) {
            maxDiff = std::max(maxDiff, std::fabs(R.finalState.u[g] - rr.u[g]));
            maxU    = std::max(maxU, std::fabs(rr.u[g]));
        }
        const real rel = maxDiff / std::max<real>(1e-30, maxU);
        addRow("Tension-only", "ReSolve eliminator converged state == omitted-member model",
               "compressed tension-only brace dropped -> displacements equal omitting it (ReSolve round-off)",
               "rel max |u_TO - u_omitted|", rel, 1e-10, R.converged && R.slack.size() == 1 && rel < 1e-10);
    }

    // ---- Check 2: tensionOnly is in the solveLoad reuse fingerprint (a flip rejects a stale factor) ----
    {
        FrameModel m; fixtures::xBracedPortal(m, 5.0e4, 0.0, tmat, stocky, tbrace);
        PreparedSystem ps = assembleAndFactor(m);
        m.members[3].tensionOnly = !m.members[3].tensionOnly;
        const SolveResult stale = solveLoad(ps, m);
        addRow("Tension-only", "tensionOnly flag in solveLoad fingerprint",
               "flipping Member::tensionOnly rejects a reused factorization (no silent stale solve)",
               "stale solve flagged singular (1=yes)", stale.singular ? 1.0 : 0.0, 0.5, stale.singular);
    }
}

void testSizeOpt() {
    // ---- Check 1: statically DETERMINATE -> stress-ratio FSD = exact minimum-weight optimum ----
    // A single axial column: the member force is area-independent, so A converges to |N|/sigma_a
    // in one effective step (Haftka & Gurdal: determinate FSD is provably the lightest design).
    {
        const real E = 210000.0, sigA = 250.0, P = 5.0e5, h = 3000.0;
        Material mat(E, E / 2.6, 0.0); mat.cap = Capacity::make(sigA, sigA, sigA);
        Section sec = Section::Rectangular(80.0, 80.0);
        FrameModel m; fixtures::axialColumn(m, P, h, mat, sec);
        SizeOptOptions o; o.maxIter = 50; o.dcTol = 1e-12; o.Amin = 1.0;
        const SizeOptResult R = runSizeOptimization(m, o);
        const real Aexact = P / sigA;                          // sigma = P/A = sigma_a  ->  A = P/sigma_a
        const real rel = std::fabs(R.finalAreas[0] - Aexact) / Aexact;
        addRow("Size optimization (FSD)", "determinate stress-ratio = exact optimum",
               "single axial column: A -> |N|/sigma_allow (determinate FSD provably minimum-weight)",
               "rel |A_final - P/sigma_a|", rel, 1e-12, R.converged && rel < 1e-12);
    }

    // ---- Check 2: MULTI-CASE envelope -- the worst case sizes the member; D/C <= 1 under BOTH ----
    {
        const real E = 210000.0, sigA = 250.0, h = 3000.0, P1 = 3.0e5, P2 = 6.0e5;
        Material mat(E, E / 2.6, 0.0); mat.cap = Capacity::make(sigA, sigA, sigA);
        Section sec = Section::Rectangular(80.0, 80.0);
        FrameModel m; fixtures::axialColumn(m, P1, h, mat, sec);   // geometry; loads come from cases
        SizeOptOptions o; o.maxIter = 50; o.dcTol = 1e-12; o.Amin = 1.0;
        SizeOptLoadCase cA, cB;
        { NodalLoad nl; nl.node = 1; nl.comp[Uz] = -P1; cA.nodalLoads = { nl }; }
        { NodalLoad nl; nl.node = 1; nl.comp[Uz] = -P2; cB.nodalLoads = { nl }; }
        o.cases = { cA, cB };
        const SizeOptResult R = runSizeOptimization(m, o);
        const real Aexact = P2 / sigA;                            // driven by the larger case
        const bool sizedByMax = std::fabs(R.finalAreas[0] - Aexact) / Aexact < 1e-10;
        const bool safeBoth   = R.finalDC[0] <= 1.0 + 1e-9;       // envelope D/C at the cap, never above
        addRow("Size optimization (FSD)", "multi-case envelope sizing (worst case governs)",
               "two axial cases: area = max(P)/sigma_a; envelope D/C <= 1 under BOTH cases",
               "final envelope D/C", R.finalDC[0], 1.0 + 1e-9, R.converged && sizedByMax && safeBoth);
    }

    // ---- Check 3: DISCRETE section table -- round-up (conservative) + finite termination ----
    {
        const real E = 210000.0, sigA = 250.0, h = 3000.0, P = 5.0e5;
        Material mat(E, E / 2.6, 0.0); mat.cap = Capacity::make(sigA, sigA, sigA);
        Section sec = Section::Rectangular(80.0, 80.0);
        FrameModel m; fixtures::axialColumn(m, P, h, mat, sec);
        SizeOptOptions o; o.maxIter = 50; o.Amin = 1.0;
        o.sectionTable = { 500.0, 1000.0, 1500.0, 2200.0, 3000.0, 5000.0 };  // ascending mm^2
        const real Acont = P / sigA;                              // 2000 mm^2 continuous optimum
        const SizeOptResult R = runSizeOptimization(m, o);
        const bool terminated   = R.converged || R.cycled;        // never hangs
        const bool isTableVal   = std::fabs(R.finalAreas[0] - 2200.0) < 1e-9;  // round-up of 2000
        const bool conservative = R.finalAreas[0] >= Acont;       // round-UP is never unsafe
        addRow("Size optimization (FSD)", "discrete section table: round-up + oscillation guard",
               "coarse area table -> snap UP to the smallest sufficient section; finite termination",
               "final area (round-up 2000 -> 2200)", R.finalAreas[0], 2200.0,
               terminated && isTableVal && conservative);
    }
}

// S7 BESO topology optimization + N2 collapse robustness.
void testBESO() {
    const real E = 210000.0, G = 80769.0;
    Material mat(E, G, 7850.0); mat.cap = Capacity::make(300.0, 300.0, 180.0);
    const Section sec = Section::Rectangular(80.0, 80.0);

    // ---- Check 1: sensitivity = element strain energy -> energy balance sum(alpha) == 1/2 F.u ----
    // A portal frame (axial + bending + shear engaged) independent of F45's L-frame: the BESO
    // sensitivity is the exact element strain energy, so it must sum to half the external work.
    {
        FrameModel m; m.materials.push_back(mat); m.sections.push_back(sec);
        m.nodes.push_back(Node(0, 0, 0, 0));      m.nodes[0].fixAll();
        m.nodes.push_back(Node(1, 0, 0, 2000));
        m.nodes.push_back(Node(2, 2000, 0, 2000));
        m.nodes.push_back(Node(3, 2000, 0, 0));   m.nodes[3].fixAll();
        { Member a(0, 0, 1, 0, 0); a.refVec = { 1, 0, 0 }; m.members.push_back(a); }   // column
        { Member b(1, 1, 2, 0, 0); b.refVec = { 0, 0, 1 }; m.members.push_back(b); }   // beam
        { Member c(2, 2, 3, 0, 0); c.refVec = { 1, 0, 0 }; m.members.push_back(c); }   // column
        NodalLoad l1; l1.node = 1; l1.comp[Ux] = 10000.0; m.nodalLoads.push_back(l1);
        NodalLoad l2; l2.node = 2; l2.comp[Uz] = -6000.0; m.nodalLoads.push_back(l2);
        const SolveResult r = solve(m);
        real sumA = 0;
        for (size_t e = 0; e < m.members.size(); ++e) sumA += memberStrainEnergy(m, r, (int)e);
        real C = 0;
        for (const NodalLoad& L : m.nodalLoads) {
            const int ni = m.nodeIndex(L.node);
            for (int d = 0; d < 6; ++d) C += L.comp[d] * r.u[(size_t)gdof(ni, d)];
        }
        const real rel = std::fabs(sumA - 0.5 * C) / std::max<real>(1.0, std::fabs(0.5 * C));
        addRow("Topology (BESO)", "sensitivity = exact element strain energy (all components)",
               "portal frame: sum_e alpha_e == 1/2 F.u (energy balance, independent of F45 L-frame)",
               "rel |sum a - 1/2 F.u|", rel, 1e-9, !r.singular && rel < 1e-9);
    }

    // ---- Check 2: N2 robustness — robust BESO topology survives every single-member removal + locks
    // (unconstrained sister topology does NOT). Colinear main chain g+m with a small-section backup b.
    {
        FrameModel m; m.materials.push_back(mat);
        m.sections.push_back(sec);
        m.sections.push_back(Section::Rectangular(50.0, 50.0));
        m.nodes.push_back(Node(0, 0, 0, 0));      m.nodes[0].fixAll();
        m.nodes.push_back(Node(1, 0, 0, 2000));   m.nodes[1].fixAll();
        m.nodes.push_back(Node(2, 1500, 0, 0));
        m.nodes.push_back(Node(3, 3000, 0, 1000));
        m.members.push_back(Member(0, 0, 2, 0, 0));   // g (main)
        m.members.push_back(Member(1, 2, 3, 0, 0));   // m (main)
        m.members.push_back(Member(2, 1, 3, 0, 1));   // b (small backup)
        NodalLoad l; l.node = 3; l.comp[Uz] = -2000.0; m.nodalLoads.push_back(l);
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
        BESOOptions ou; ou.targetVolFrac = 0.6; ou.evolRate = 0.34; ou.maxIter = 20;
        const BESOResult Ru = runBESO(m, ou);
        BESOOptions orb = ou; orb.redundancyCheckEvery = 1; orb.redundancySamples = 0;
        orb.redundancy.dlf = 1.0; orb.redundancy.removeThreshold = 1.0;
        const BESOResult Rr = runBESO(m, orb);
        const bool uncFragile = fragile(Ru.finalActive);
        const bool robSafe    = !fragile(Rr.finalActive) && !Rr.protectedMembers.empty();
        addRow("Topology (BESO N2)", "collapse-robustness constraint rolls back + locks bars",
               "unconstrained topology collapses on a single removal; robust one survives all + locks",
               "robust-safe AND unc-fragile (1=yes)", (uncFragile && robSafe) ? 1.0 : 0.0, 0.0,
               uncFragile && robSafe);
    }

    // ---- Check 3: mechanism guard — BESO never stops on a singular topology (FRESH-confirmed) ----
    {
        const int NX = 5, NZ = 3; const real SP = 500.0;
        auto nid = [&](int i, int k) { return k * (NX + 1) + i; };
        FrameModel m; m.materials.push_back(mat); m.sections.push_back(sec);
        for (int k = 0; k <= NZ; ++k)
            for (int i = 0; i <= NX; ++i) {
                Node n(nid(i, k), i * SP, 0.0, k * SP);
                if (i == 0) n.fixAll();
                else { n.fixed[Uy] = true; n.fixed[Rx] = true; n.fixed[Rz] = true; }
                m.nodes.push_back(n);
            }
        auto add = [&](int a, int b) {
            Member mem((int)m.members.size(), a, b, 0, 0);
            const Vec3 d = m.nodes[(size_t)b].pos - m.nodes[(size_t)a].pos;
            mem.refVec = (std::fabs(d.z) > std::fabs(d.x)) ? Vec3{ 1, 0, 0 } : Vec3{ 0, 0, 1 };
            m.members.push_back(mem);
        };
        for (int k = 0; k <= NZ; ++k) for (int i = 0; i < NX; ++i) add(nid(i, k), nid(i + 1, k));
        for (int k = 0; k < NZ; ++k) for (int i = 0; i <= NX; ++i) add(nid(i, k), nid(i, k + 1));
        for (int k = 0; k < NZ; ++k) for (int i = 0; i < NX; ++i) { add(nid(i, k), nid(i + 1, k + 1)); add(nid(i + 1, k), nid(i, k + 1)); }
        NodalLoad l; l.node = nid(NX, NZ / 2); l.comp[Uz] = -1.0e4; m.nodalLoads.push_back(l);
        BESOOptions o; o.targetVolFrac = 0.45; o.evolRate = 0.06; o.maxIter = 80;
        const BESOResult R = runBESO(m, o);
        FrameModel mb = m;
        for (size_t e = 0; e < mb.members.size(); ++e) mb.members[e].active = R.bestActive[e] != 0;
        const SolveResult rb = solve(mb);
        addRow("Topology (BESO)", "mechanism guard (capacitance fast + fresh confirm)",
               "best topology factors under a FRESH solve -> BESO never stops on a singular structure",
               "best topology singular flag (0=ok)", rb.singular ? 1.0 : 0.0, 0.0,
               !rb.singular && !R.finalActive.empty());
    }
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
    testShellS8();
    testSparseModal();
    testElementRemoval();
    testShellRemoval();
    testConnectivity();
    testCollapseDriver();
    testShellFailureScreen();
    testPlasticHingeMechanics();
    testHingeDriverCompat();
    testSafetyAndMargin();
    testReanalysis();
    testDynamicCollapse();
    testPDelta();
    testTensionOnly();
    testSizeOpt();
    testBESO();

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
