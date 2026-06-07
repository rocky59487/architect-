// Stage-1 automation tests (mirror standalone F17/F18): load-case superposition (the
// combination primitive) and self-weight derived from Material.rho with the kg/m^3 ->
// N-mm-tonne-s unit bridge. Same solver path + same fixtures as the standalone gate.
#include "CoreMinimal.h"
#include "Misc/AutomationTest.h"
#include "FrameCore/FrameSolver.h"
#include "FrameCore/SelfWeight.h"
#include "FrameCore/Combination.h"
#include "FrameCore/InfluenceLine.h"
#include "FrameCore/ModalAnalysis.h"
#include "FrameCore/Material.h"
#include "FrameCore/Section.h"
#include "FrameTestFixtures.h"

#include <cmath>

#if WITH_DEV_AUTOMATION_TESTS

// ---- F17 mirror: superposition identity ----
IMPLEMENT_SIMPLE_AUTOMATION_TEST(FFrameCoreLoadSuperpositionTest,
	"FrameCore.Load.Superposition",
	EAutomationTestFlags_ApplicationContextMask | EAutomationTestFlags::SmokeFilter)
bool FFrameCoreLoadSuperpositionTest::RunTest(const FString&)
{
	using namespace frame;
	Section sec = Section::Rectangular(100.0, 100.0);
	Material mat(210000.0, 80769.0, 7850.0);

	auto run = [&](double Pz, double My) -> SolveResult {
		FrameModel m; fixtures::cantileverBare(m, 2000.0, mat, sec);
		NodalLoad nl; nl.node = 1; nl.comp[Uz] = Pz; nl.comp[Ry] = My;
		m.nodalLoads = { nl };
		return solve(m);
	};
	const SolveResult rAB = run(1000.0, 2.0e6);
	const SolveResult rC  = combine({ run(1000.0, 0.0), run(0.0, 2.0e6) }, { 1.0, 1.0 });

	double us = 1e-30, du = 0;
	for (size_t k = 0; k < rAB.u.size(); ++k) { us = FMath::Max(us, FMath::Abs(rAB.u[k])); du = FMath::Max(du, FMath::Abs(rC.u[k] - rAB.u[k])); }
	TestTrue(TEXT("u(A)+u(B) == u(A+B) to machine precision"), du < 1e-10 * us);
	return true;
}

// ---- F18 mirror: self-weight + unit bridge ----
IMPLEMENT_SIMPLE_AUTOMATION_TEST(FFrameCoreLoadSelfWeightTest,
	"FrameCore.Load.SelfWeight",
	EAutomationTestFlags_ApplicationContextMask | EAutomationTestFlags::SmokeFilter)
bool FFrameCoreLoadSelfWeightTest::RunTest(const FString&)
{
	using namespace frame;
	Section sec = Section::Rectangular(100.0, 100.0);
	Material mat(210000.0, 80769.0, 7850.0);
	const double g = 9810.0, E = 210000.0;
	const double wg = mat.rho * sec.A * g * 1.0e-12;   // N/mm

	// cantilever: fixed-end moment wL^2/2, reaction wL.
	{
		const double L = 2000.0;
		FrameModel m; fixtures::cantileverBare(m, L, mat, sec);
		addSelfWeight(m, g);
		const SolveResult r = solve(m);
		TestFalse(TEXT("cantilever SW non-singular"), r.singular);
		const auto& mf = r.memberForces[0].endI;
		const double Mroot = FMath::Sqrt(mf.My * mf.My + mf.Mz * mf.Mz);
		TestTrue(TEXT("fixed-end moment = wL^2/2"), FMath::Abs(Mroot - wg * L * L / 2.0) < 1e-6 * wg * L * L);
		TestTrue(TEXT("reaction = wL"), FMath::Abs(FMath::Abs(r.reaction(0, Uz)) - wg * L) < 1e-6 * wg * L);
	}
	// shell self-weight (body load) == equivalent rho*t*g pressure.
	{
		const double a = 1000.0, t = 10.0, nu = 0.3, Es = 30000.0;
		Material smat(Es, Es / (2.0 * (1.0 + nu)), 2400.0); smat.nu = nu;
		const double p = smat.rho * t * g * 1.0e-12;
		const int n = 12, c = (n / 2) * (n + 1) + (n / 2);
		FrameModel mp; fixtures::squarePlateShell(mp, a, t, n, p, smat);
		const SolveResult rp = solve(mp);
		FrameModel mw; fixtures::squarePlateShell(mw, a, t, n, 0.0, smat);
		addSelfWeight(mw, g);
		const SolveResult rw = solve(mw);
		const double wP = FMath::Abs(rp.disp(c, Uz)), wW = FMath::Abs(rw.disp(c, Uz));
		TestTrue(TEXT("shell self-weight == rho*t*g pressure"), FMath::Abs(wW - wP) < 1e-9 * FMath::Max(wP, 1e-30));
	}
	(void)E;
	return true;
}

// ---- F19 mirror (a): factorize-once-solve-many reuse ----
IMPLEMENT_SIMPLE_AUTOMATION_TEST(FFrameCoreSolverFactorizeOnceTest,
	"FrameCore.Solver.FactorizeOnce",
	EAutomationTestFlags_ApplicationContextMask | EAutomationTestFlags::SmokeFilter)
bool FFrameCoreSolverFactorizeOnceTest::RunTest(const FString&)
{
	using namespace frame;
	Section sec = Section::Rectangular(100.0, 100.0);
	Material mat(210000.0, 80769.0, 7850.0);

	FrameModel m; fixtures::cantileverBare(m, 2000.0, mat, sec);
	PreparedSystem ps = assembleAndFactor(m);
	auto maxDiff = [](const SolveResult& a, const SolveResult& b) {
		double d = 0; for (size_t k = 0; k < a.u.size(); ++k) d = FMath::Max(d, FMath::Abs(a.u[k] - b.u[k])); return d;
	};
	m.nodalLoads = { [] { NodalLoad n; n.node = 1; n.comp[Uz] = 1000.0; return n; }() };
	TestTrue(TEXT("solveLoad(A) == fresh solve(A)"), maxDiff(solveLoad(ps, m), solve(m)) < 1e-12);
	m.nodalLoads = { [] { NodalLoad n; n.node = 1; n.comp[Uy] = 700.0; n.comp[Ry] = 1.0e6; return n; }() };
	TestTrue(TEXT("solveLoad(B) reuses same factorization == fresh"), maxDiff(solveLoad(ps, m), solve(m)) < 1e-12);
	return true;
}

// ---- F19 mirror (b): prescribed support settlement ----
IMPLEMENT_SIMPLE_AUTOMATION_TEST(FFrameCoreSolverSettlementTest,
	"FrameCore.Solver.Settlement",
	EAutomationTestFlags_ApplicationContextMask | EAutomationTestFlags::SmokeFilter)
bool FFrameCoreSolverSettlementTest::RunTest(const FString&)
{
	using namespace frame;
	Section sec = Section::Rectangular(100.0, 100.0);
	Material mat(210000.0, 80769.0, 7850.0);
	const double L = 2000.0, delta = 1.0, E = 210000.0;
	FrameModel m; fixtures::clampedSettlement(m, L, delta, mat, sec);
	const SolveResult r = solve(m);
	TestFalse(TEXT("settlement non-singular"), r.singular);
	const auto& mf = r.memberForces[0].endI;
	const double Mend = FMath::Sqrt(mf.My * mf.My + mf.Mz * mf.Mz);
	TestTrue(TEXT("end moment 6EI d/L^2"),   FMath::Abs(Mend - 6.0 * E * sec.Iz * delta / (L * L)) < 1e-6 * 6.0 * E * sec.Iz * delta / (L * L));
	TestTrue(TEXT("reaction 12EI d/L^3"), FMath::Abs(FMath::Abs(r.reaction(0, Uz)) - 12.0 * E * sec.Iz * delta / (L * L * L)) < 1e-6 * 12.0 * E * sec.Iz * delta / (L * L * L));
	return true;
}

// ---- F20 mirror: live-load pattern loading + envelope ----
IMPLEMENT_SIMPLE_AUTOMATION_TEST(FFrameCoreLoadPatternEnvelopeTest,
	"FrameCore.Load.PatternEnvelope",
	EAutomationTestFlags_ApplicationContextMask | EAutomationTestFlags::SmokeFilter)
bool FFrameCoreLoadPatternEnvelopeTest::RunTest(const FString&)
{
	using namespace frame;
	Section sec = Section::Rectangular(100.0, 100.0);
	Material mat(210000.0, 80769.0, 7850.0);
	const double L = 3000.0, w = 5.0;
	auto pattern = [&](bool left, bool right) -> SolveResult {
		FrameModel m; fixtures::twoSpanContinuous(m, L, mat, sec);
		auto addUDL = [&](int mem) { MemberUDL u; u.member = mem; u.w_local = { 0, -w, 0 }; m.memberUDLs.push_back(u); };
		if (left)  { addUDL(0); addUDL(1); }
		if (right) { addUDL(2); addUDL(3); }
		return solve(m);
	};
	const SolveResult rL = pattern(true, false), rR = pattern(false, true), rAll = pattern(true, true);
	auto Mres = [](const MemberEndForces& f) { return FMath::Sqrt(f.My * f.My + f.Mz * f.Mz); };
	TestTrue(TEXT("support moment (full) = wL^2/8"),
		FMath::Abs(Mres(rAll.memberForces[1].endJ) - w * L * L / 8.0) < 1e-4 * w * L * L / 8.0);
	const ResultEnvelope env = envelope({ rL, rR, rAll });
	const double MBworst = FMath::Sqrt(env.endJMin[1].My * env.endJMin[1].My + env.endJMin[1].Mz * env.endJMin[1].Mz);
	TestTrue(TEXT("envelope worst support moment = wL^2/8"), FMath::Abs(MBworst - w * L * L / 8.0) < 1e-4 * w * L * L / 8.0);
	TestTrue(TEXT("single-span pattern bends its span more than full load"),
		Mres(rL.memberForces[0].endJ) > Mres(rAll.memberForces[0].endJ) * 1.05);
	return true;
}

// ---- F21 mirror: influence line + Muller-Breslau cross-check ----
IMPLEMENT_SIMPLE_AUTOMATION_TEST(FFrameCoreLoadInfluenceLineTest,
	"FrameCore.Load.InfluenceLine",
	EAutomationTestFlags_ApplicationContextMask | EAutomationTestFlags::SmokeFilter)
bool FFrameCoreLoadInfluenceLineTest::RunTest(const FString&)
{
	using namespace frame;
	Section sec = Section::Rectangular(100.0, 100.0);
	Material mat(210000.0, 80769.0, 7850.0);
	const int n = 8; const double L = 4000.0;
	FrameModel m; fixtures::simplySupportedBeamN(m, n, L, mat, sec);
	PreparedSystem ps = assembleAndFactor(m);
	std::vector<NodeId> loadNodes; for (int i = 0; i <= n; ++i) loadNodes.push_back(i);

	const std::vector<double> ilR = reactionInfluenceLine(ps, m, loadNodes, 0, Uz);
	m.nodalLoads.clear(); m.nodes[0].prescribed[Uz] = 1.0;
	const SolveResult rMB = solveLoad(ps, m);
	m.nodes[0].prescribed[Uz] = 0.0;
	double eAna = 0, eMB = 0;
	for (int i = 0; i <= n; ++i) {
		const double x = L * i / n, exact = (L - x) / L;
		eAna = FMath::Max(eAna, FMath::Abs(ilR[i] - exact));
		eMB  = FMath::Max(eMB, FMath::Abs(ilR[i] - rMB.disp(i, Uz)));
	}
	TestTrue(TEXT("reaction IL == (L-x)/L"), eAna < 1e-9);
	TestTrue(TEXT("reaction IL == Muller-Breslau deflected shape"), eMB < 1e-9);
	return true;
}

// ---- F22 mirror: modal analysis (natural frequencies) ----
IMPLEMENT_SIMPLE_AUTOMATION_TEST(FFrameCoreModalBeamTest,
	"FrameCore.Modal.Beam",
	EAutomationTestFlags_ApplicationContextMask | EAutomationTestFlags::SmokeFilter)
bool FFrameCoreModalBeamTest::RunTest(const FString&)
{
	using namespace frame;
	Section sec = Section::Rectangular(100.0, 100.0);
	Material mat(210000.0, 80769.0, 7850.0);
	const double E = 210000.0, kPi = 3.14159265358979323846, rhoC = mat.rho * 1.0e-12;

	{   // simply-supported fundamental omega1 = (pi/L)^2 sqrt(EI/rhoA)
		const int n = 12; const double L = 4000.0;
		FrameModel m; fixtures::simplySupportedBeamN(m, n, L, mat, sec);
		PreparedSystem ps = assembleAndFactor(m);
		ModalOptions mo; mo.numModes = 3;
		const ModalResult mr = solveModal(ps, mo);
		TestFalse(TEXT("SS modal non-singular"), mr.singular);
		TestTrue(TEXT("has modes"), mr.modes.size() >= 1);
		const double w1ex = kPi * kPi / (L * L) * FMath::Sqrt(E * sec.Iz / (rhoC * sec.A));
		TestTrue(TEXT("SS omega1 within 1%"), FMath::Abs(mr.modes[0].omega - w1ex) < 0.01 * w1ex);
	}
	{   // cantilever fundamental omega1 = 1.875^2 sqrt(EI/rhoAL^4)
		const int n = 12; const double L = 3000.0;
		FrameModel m; fixtures::cantileverBeamN(m, n, L, mat, sec);
		PreparedSystem ps = assembleAndFactor(m);
		const ModalResult mr = solveModal(ps, ModalOptions{});
		const double b1 = 1.8751040687;
		const double w1ex = b1 * b1 * FMath::Sqrt(E * sec.Iz / (rhoC * sec.A * L * L * L * L));
		TestTrue(TEXT("cantilever omega1 within 1%"), FMath::Abs(mr.modes[0].omega - w1ex) < 0.01 * w1ex);
	}
	return true;
}

#endif // WITH_DEV_AUTOMATION_TESTS
