// M6 upgrades (mirrors of the standalone F10-F12 + an IElement-refactor regression):
//  * Element.IElementRegression : cantilever + simply-supported still hit the analytic
//    solutions after the solver was routed through the IElement seam (zero regression).
//  * Timoshenko.DeepCantilever   : shear-flexible tip deflection = PL^3/3EI + PL/GAs;
//    deep beam shows shear, slender converges back to Euler-Bernoulli.
//  * Section.Circular            : I = pi r^4/4, J = 2I, and the round-section biaxial
//    bending uses the resultant sqrt(My^2+Mz^2)/W (< the rectangular corner sum).
//  * Curved.ArchConvergence      : a faceted quarter-circle cantilever converges to the
//    closed-form (pi/4) P R^3 / (E I) as the segment count grows.
#include "CoreMinimal.h"
#include "Misc/AutomationTest.h"
#include "FrameCore/FrameSolver.h"
#include "FrameCore/Section.h"
#include "FrameCore/Material.h"
#include "FrameCore/ElasticAllowable.h"
#include "FrameTestFixtures.h"

#include <cmath>

#if WITH_DEV_AUTOMATION_TESTS

namespace
{
	frame::Material SteelMat()
	{
		frame::Material m(210000.0, 80769.0, 7850.0);
		m.cap = frame::Capacity::make(300.0, 300.0, 180.0);
		return m;
	}
}

IMPLEMENT_SIMPLE_AUTOMATION_TEST(FFrameCoreIElementRegressionTest,
	"FrameCore.Element.IElementRegression",
	EAutomationTestFlags_ApplicationContextMask | EAutomationTestFlags::SmokeFilter)
bool FFrameCoreIElementRegressionTest::RunTest(const FString&)
{
	using namespace frame;
	const Section sec = Section::Rectangular(100.0, 100.0);
	const Material mat = SteelMat();

	// cantilever F1 through the IElement path
	{
		const real P = 1000.0, L = 2000.0;
		FrameModel m; fixtures::cantileverTipLoad(m, P, L, mat, sec);
		const SolveResult r = solve(m);
		TestTrue(TEXT("cantilever not singular"), !r.singular);
		TestEqual(TEXT("element recovery filled member forces"), (int32)r.memberForces.size(), 1);
		const real dExp = P * L * L * L / (3.0 * 210000.0 * sec.Iz);
		TestTrue(TEXT("tip deflection PL^3/3EI"),
			FMath::Abs(std::fabs(r.disp(1, Uz)) - dExp) < 1e-6 * dExp);
		const real Mroot = std::sqrt(r.memberForces[0].endI.My * r.memberForces[0].endI.My
			+ r.memberForces[0].endI.Mz * r.memberForces[0].endI.Mz);
		TestTrue(TEXT("root moment PL"), FMath::Abs(Mroot - P * L) < 1e-6 * (P * L));
	}
	// simply-supported F2 through the IElement path
	{
		const real w = 5.0, L = 3000.0;
		FrameModel m; fixtures::simplySupportedUDL(m, w, L, mat, sec);
		const SolveResult r = solve(m);
		TestTrue(TEXT("ss not singular"), !r.singular);
		const real dExp = 5.0 * w * L * L * L * L / (384.0 * 210000.0 * sec.Iz);
		TestTrue(TEXT("midspan deflection 5wL^4/384EI"),
			FMath::Abs(std::fabs(r.disp(1, Uz)) - dExp) < 1e-6 * dExp);
		TestTrue(TEXT("midspan moment wL^2/8"),
			FMath::Abs(std::fabs(r.memberForces[0].endJ.Mz) - w * L * L / 8.0) < 1e-5 * (w * L * L / 8.0));
	}
	return true;
}

IMPLEMENT_SIMPLE_AUTOMATION_TEST(FFrameCoreTimoshenkoTest,
	"FrameCore.Timoshenko.DeepCantilever",
	EAutomationTestFlags_ApplicationContextMask | EAutomationTestFlags::SmokeFilter)
bool FFrameCoreTimoshenkoTest::RunTest(const FString&)
{
	using namespace frame;
	const Section sec = Section::Rectangular(100.0, 100.0);
	const Material mat = SteelMat();
	const real E = 210000.0, G = 80769.0;

	const real P = 1000.0, L = 2000.0;
	FrameModel m; fixtures::cantileverTipLoad(m, P, L, mat, sec);
	SolveOptions opt; opt.useTimoshenko = true;
	const SolveResult r = solve(m, opt);
	const real dEB    = P * L * L * L / (3.0 * E * sec.Iz);
	const real dShear = P * L / (G * sec.Asy);
	TestTrue(TEXT("not singular"), !r.singular);
	TestTrue(TEXT("tip = PL^3/3EI + PL/GAs"),
		FMath::Abs(std::fabs(r.disp(1, Uz)) - (dEB + dShear)) < 1e-6 * (dEB + dShear));

	auto ratio = [&](real Ls) -> real
	{
		FrameModel ms; fixtures::cantileverTipLoad(ms, P, Ls, mat, sec);
		SolveOptions o; o.useTimoshenko = true;
		const real dt = std::fabs(solve(ms, o).disp(1, Uz));
		const real de = P * Ls * Ls * Ls / (3.0 * E * sec.Iz);
		return dt / de;
	};
	const real rDeep = ratio(300.0);   // L/d = 3
	const real rSlen = ratio(6000.0);  // L/d = 60
	TestTrue(TEXT("deep-beam shear > 5%"), rDeep > 1.05);
	TestTrue(TEXT("slender converges to EB"), FMath::Abs(rSlen - 1.0) < 1e-3);
	TestTrue(TEXT("shear effect shrinks with slenderness"), rSlen < rDeep);
	return true;
}

IMPLEMENT_SIMPLE_AUTOMATION_TEST(FFrameCoreCircularSectionTest,
	"FrameCore.Section.Circular",
	EAutomationTestFlags_ApplicationContextMask | EAutomationTestFlags::SmokeFilter)
bool FFrameCoreCircularSectionTest::RunTest(const FString&)
{
	using namespace frame;
	const real rad = 50.0;
	const Section cir = Section::Circular(rad);
	const real Iexp = kPi * rad * rad * rad * rad / 4.0;
	TestTrue(TEXT("I = pi r^4/4"), FMath::Abs(cir.Iz - Iexp) < 1e-9 * Iexp);
	TestTrue(TEXT("Iy == Iz"), FMath::Abs(cir.Iy - cir.Iz) < 1e-12 * Iexp);
	TestTrue(TEXT("J = 2I"), FMath::Abs(cir.J - 2.0 * Iexp) < 1e-9 * Iexp);

	const Material mat = SteelMat();
	const real E = 210000.0, P = 1000.0, L = 2000.0;
	FrameModel m; fixtures::cantileverTipLoad(m, P, L, mat, cir);
	const SolveResult r = solve(m);
	const real dExp = P * L * L * L / (3.0 * E * Iexp);
	TestTrue(TEXT("circular cantilever tip"), FMath::Abs(std::fabs(r.disp(1, Uz)) - dExp) < 1e-6 * dExp);

	ElasticAllowable strength;
	const Capacity cap = Capacity::make(250.0, 250.0, 150.0);
	MemberEndForces f; f.My = 1.0e6; f.Mz = 1.0e6;
	const DemandResult d = strength.checkSection(f, cir, cap);
	const real W = cir.Wz();
	TestTrue(TEXT("biaxial resultant sqrt2*M/W"),
		FMath::Abs(d.sComp - std::sqrt(2.0) * 1.0e6 / W) < 1e-6 * (std::sqrt(2.0) * 1.0e6 / W));
	TestTrue(TEXT("round resultant < rectangular corner sum"), d.sComp < 2.0e6 / W - 1e-9);
	return true;
}

IMPLEMENT_SIMPLE_AUTOMATION_TEST(FFrameCoreArchConvergenceTest,
	"FrameCore.Curved.ArchConvergence",
	EAutomationTestFlags_ApplicationContextMask | EAutomationTestFlags::SmokeFilter)
bool FFrameCoreArchConvergenceTest::RunTest(const FString&)
{
	using namespace frame;
	const Section sec = Section::Rectangular(100.0, 100.0);   // square: Iy == Iz
	const Material mat = SteelMat();
	const real E = 210000.0, R = 2000.0, P = 1000.0;
	const real dExp = (kPi / 4.0) * P * R * R * R / (E * sec.Iz);

	const int Ns[4] = { 8, 16, 32, 64 };
	real prevErr = 1e30, lastVal = 0; bool monotone = true, nonsing = true;
	for (int i = 0; i < 4; ++i)
	{
		FrameModel m; fixtures::circularArchCantilever(m, R, Ns[i], P, mat, sec);
		const SolveResult r = solve(m);
		if (r.singular) nonsing = false;
		const real d   = std::fabs(r.disp(Ns[i], Uy));
		const real err = std::fabs(d - dExp) / dExp;
		if (err > prevErr) monotone = false;
		prevErr = err; lastVal = d;
	}
	TestTrue(TEXT("arch solves (non-singular)"), nonsing);
	TestTrue(TEXT("arch error shrinks as N grows"), monotone);
	TestTrue(TEXT("arch tip (N=64) -> (pi/4)PR^3/EI"), FMath::Abs(lastVal - dExp) < 1.5e-2 * dExp);
	return true;
}

#endif // WITH_DEV_AUTOMATION_TESTS
