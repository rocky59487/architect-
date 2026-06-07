// Stage-1 automation tests (mirror standalone F17/F18): load-case superposition (the
// combination primitive) and self-weight derived from Material.rho with the kg/m^3 ->
// N-mm-tonne-s unit bridge. Same solver path + same fixtures as the standalone gate.
#include "CoreMinimal.h"
#include "Misc/AutomationTest.h"
#include "FrameCore/FrameSolver.h"
#include "FrameCore/SelfWeight.h"
#include "FrameCore/Combination.h"
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

#endif // WITH_DEV_AUTOMATION_TESTS
