// MITC4 flat-shell automation tests (mirror standalone F13-F16). Same solver path,
// same fixtures as the standalone gate, so the UE module exercises the shell element
// identically. Covers: plate bending (patch + convergence + no shear locking),
// membrane + drilling (patch + coplanar non-singular + 3D rotation invariance), and the
// two curved-shell benchmarks (Scordelis-Lo roof, pinched cylinder).
#include "CoreMinimal.h"
#include "Misc/AutomationTest.h"
#include "FrameCore/FrameSolver.h"
#include "FrameCore/Material.h"
#include "FrameTestFixtures.h"

#include <cmath>

#if WITH_DEV_AUTOMATION_TESTS

namespace {
constexpr double kShellPi = 3.14159265358979323846;
inline frame::Material shellMat()
{
	const double E = 30000.0, nu = 0.3;
	frame::Material m(E, E / (2.0 * (1.0 + nu)));
	m.nu = nu;
	return m;
}
}

// ---- F13 mirror: plate bending (patch reproduction + Kirchhoff convergence + no lock) ----
IMPLEMENT_SIMPLE_AUTOMATION_TEST(FFrameCoreShellPlateBendingTest,
	"FrameCore.Shell.PlateBending",
	EAutomationTestFlags_ApplicationContextMask | EAutomationTestFlags::SmokeFilter)
bool FFrameCoreShellPlateBendingTest::RunTest(const FString&)
{
	using namespace frame;
	const Material smat = shellMat();
	const double E = 30000.0, nu = 0.3;

	// constant-curvature patch (parallelogram mesh) -> exact constant moment.
	{
		const double a = 1000.0, t = 10.0, c = 1e-6;
		const double Dfac = E * t * t * t / (12.0 * (1.0 - nu * nu));
		FrameModel m; fixtures::platePatchCylindrical(m, a, t, 0.4, c, smat);
		const SolveResult r = solve(m);
		TestFalse(TEXT("patch non-singular"), r.singular);
		const double MxxExp = -Dfac * c, scale = FMath::Abs(MxxExp);
		double eMxx = 0;
		for (const auto& sf : r.shellForces) eMxx = FMath::Max(eMxx, FMath::Abs(sf.Mxx - MxxExp));
		TestTrue(TEXT("patch reproduces constant moment"), eMxx < 1e-8 * scale);
	}

	// simply-supported square plate -> Kirchhoff w_c = 0.00406 q a^4 / D, within 2% at N=16.
	{
		const double a = 1000.0, t = 10.0, q = 0.01;
		const double D = E * t * t * t / (12.0 * (1.0 - nu * nu));
		const double wc = 0.00406 * q * a * a * a * a / D;
		FrameModel m; fixtures::squarePlateShell(m, a, t, 16, q, smat);
		const SolveResult r = solve(m);
		TestFalse(TEXT("square plate non-singular"), r.singular);
		const double w = FMath::Abs(r.disp(8 * 17 + 8, Uz));
		TestTrue(TEXT("square plate within 2% of Kirchhoff"), FMath::Abs(w - wc) < 0.02 * wc);
	}

	// thin plate (t/a = 0.001) must NOT shear-lock.
	{
		const double a = 1000.0, t = 1.0, q = 0.01;
		const double D = E * t * t * t / (12.0 * (1.0 - nu * nu));
		const double wc = 0.00406 * q * a * a * a * a / D;
		FrameModel m; fixtures::squarePlateShell(m, a, t, 16, q, smat);
		const SolveResult r = solve(m);
		const double w = FMath::Abs(r.disp(8 * 17 + 8, Uz));
		TestTrue(TEXT("thin plate not shear-locked"), FMath::Abs(w - wc) < 0.03 * wc);
	}
	return true;
}

// ---- F14 mirror: membrane + drilling + 3D facet (patch + coplanar non-singular + rotation) ----
IMPLEMENT_SIMPLE_AUTOMATION_TEST(FFrameCoreShellMembraneTest,
	"FrameCore.Shell.MembraneDrilling",
	EAutomationTestFlags_ApplicationContextMask | EAutomationTestFlags::SmokeFilter)
bool FFrameCoreShellMembraneTest::RunTest(const FString&)
{
	using namespace frame;
	const Material smat = shellMat();
	const double E = 30000.0, nu = 0.3;

	// membrane constant-strain patch (parallelogram) -> constant N.
	{
		const double a = 1000.0, t = 10.0, gx = 1e-4, f = E / (1.0 - nu * nu);
		FrameModel m; fixtures::membranePatch(m, a, t, 0.4, gx, smat);
		const SolveResult r = solve(m);
		TestFalse(TEXT("membrane patch non-singular"), r.singular);
		const double NxxExp = t * f * gx, scale = FMath::Abs(NxxExp);
		double eNxx = 0;
		for (const auto& sf : r.shellForces) eNxx = FMath::Max(eNxx, FMath::Abs(sf.Nxx - NxxExp));
		TestTrue(TEXT("membrane patch reproduces constant N"), eNxx < 1e-8 * scale);
	}

	// drilling gate: flat clamped plate with interior in-plane + drilling free -> non-singular.
	{
		const double a = 1000.0, t = 10.0, q = 0.01;
		const double D = E * t * t * t / (12.0 * (1.0 - nu * nu));
		const double wcC = 0.00126 * q * a * a * a * a / D;
		FrameModel m; fixtures::clampedPlateShell(m, a, t, 16, q, smat);
		const SolveResult r = solve(m);
		TestFalse(TEXT("flat clamped shell non-singular (drilling gate)"), r.singular);
		const double w = FMath::Abs(r.disp(8 * 17 + 8, Uz));
		TestTrue(TEXT("clamped plate within 3% of theory"), FMath::Abs(w - wcC) < 0.03 * wcC);
	}

	// rotation invariance: |u_centre| preserved under an arbitrary 3D rigid rotation.
	{
		const double a = 1000.0, t = 10.0, q = 0.01;
		const int n = 12, c = (n / 2) * (n + 1) + (n / 2);
		FrameModel m0; fixtures::clampedPlateShell(m0, a, t, n, q, smat);
		const SolveResult r0 = solve(m0);
		const double d0 = std::sqrt(r0.disp(c, Ux) * r0.disp(c, Ux) +
			r0.disp(c, Uy) * r0.disp(c, Uy) + r0.disp(c, Uz) * r0.disp(c, Uz));
		const double an = std::sqrt(14.0), kx = 1 / an, ky = 2 / an, kz = 3 / an, th = 0.9;
		const double ct = std::cos(th), st = std::sin(th), vt = 1 - ct;
		const double R[3][3] = {
			{ ct + kx * kx * vt,      kx * ky * vt - kz * st, kx * kz * vt + ky * st },
			{ ky * kx * vt + kz * st, ct + ky * ky * vt,      ky * kz * vt - kx * st },
			{ kz * kx * vt - ky * st, kz * ky * vt + kx * st, ct + kz * kz * vt }
		};
		FrameModel m1; fixtures::clampedPlateShell(m1, a, t, n, q, smat);
		fixtures::rotateModelRigid(m1, R);
		const SolveResult r1 = solve(m1);
		TestFalse(TEXT("rotated shell non-singular"), r1.singular);
		const double d1 = std::sqrt(r1.disp(c, Ux) * r1.disp(c, Ux) +
			r1.disp(c, Uy) * r1.disp(c, Uy) + r1.disp(c, Uz) * r1.disp(c, Uz));
		TestTrue(TEXT("|u_centre| invariant under 3D rotation"), FMath::Abs(d1 - d0) < 1e-9 * d0);
	}
	return true;
}

// ---- F15 mirror: Scordelis-Lo roof (curved membrane-dominated shell) ----
IMPLEMENT_SIMPLE_AUTOMATION_TEST(FFrameCoreShellScordelisLoTest,
	"FrameCore.Shell.ScordelisLo",
	EAutomationTestFlags_ApplicationContextMask | EAutomationTestFlags::SmokeFilter)
bool FFrameCoreShellScordelisLoTest::RunTest(const FString&)
{
	using namespace frame;
	const double R = 25.0, L = 50.0, phi0 = 40.0 * kShellPi / 180.0, t = 0.25, g = 90.0;
	const double Er = 4.32e8, Gr = Er / 2.0, ref = 0.3024;
	Material rmat(Er, Gr); rmat.nu = 0.0;
	auto edgeW = [&](int n) -> double {
		FrameModel m; fixtures::scordelisLoRoof(m, R, L, phi0, t, g, n, n, rmat);
		const SolveResult r = solve(m);
		if (r.singular) return 1e30;
		return FMath::Abs(r.disp(n * (n + 1) + n, Uz));
	};
	const double w8 = edgeW(8), w24 = edgeW(24);
	TestTrue(TEXT("roof non-singular"), w8 < 1e29 && w24 < 1e29);
	TestTrue(TEXT("roof within 3% of reference at N=24"), FMath::Abs(w24 - ref) < 0.03 * ref);
	TestTrue(TEXT("roof mesh-converging"), FMath::Abs(w24 - ref) < FMath::Abs(w8 - ref));
	return true;
}

// ---- F16 mirror: pinched cylinder (hard inextensional benchmark) ----
IMPLEMENT_SIMPLE_AUTOMATION_TEST(FFrameCoreShellPinchedCylinderTest,
	"FrameCore.Shell.PinchedCylinder",
	EAutomationTestFlags_ApplicationContextMask | EAutomationTestFlags::SmokeFilter)
bool FFrameCoreShellPinchedCylinderTest::RunTest(const FString&)
{
	using namespace frame;
	const double R = 300.0, L = 600.0, t = 3.0, P = 1.0;
	const double Ec = 3.0e6, nuc = 0.3, Gc = Ec / (2.0 * (1.0 + nuc)), ref = 1.8248e-5;
	Material cmat(Ec, Gc); cmat.nu = nuc;
	auto loadW = [&](int n) -> double {
		FrameModel m; fixtures::pinchedCylinder(m, R, L, t, P, n, n, cmat);
		const SolveResult r = solve(m);
		if (r.singular) return 0.0;
		return FMath::Abs(r.disp(n * (n + 1) + 0, Ux));
	};
	const double w8 = loadW(8), w16 = loadW(16), w32 = loadW(32);
	TestTrue(TEXT("cylinder converging upward"), w32 > w16 && w16 > w8);
	TestTrue(TEXT("cylinder N=32 >= 90% of reference"), w32 >= 0.90 * ref);
	TestTrue(TEXT("cylinder N=32 <= 105% of reference"), w32 <= 1.05 * ref);
	return true;
}

#endif // WITH_DEV_AUTOMATION_TESTS
