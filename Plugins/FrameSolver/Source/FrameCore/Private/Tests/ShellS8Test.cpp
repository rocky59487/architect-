// S8 shell automation tests (mirror standalone F48-F49). Same solver path / fixtures as the
// standalone gate, so the UE module exercises the opt-in QM6 incompatible-mode membrane (8a)
// and the DKQ discrete-Kirchhoff thin-plate bending (8b) identically. Both default OFF, so the
// existing MITC4 shell gate is untouched (verified by the default-path checks in F13-F16).
// NOTE: do NOT name any local constant IN / OUT here -- those are Windows SAL macros pulled in
// through CoreMinimal.h and would not compile.
#include "CoreMinimal.h"
#include "Misc/AutomationTest.h"
#include "FrameCore/FrameSolver.h"
#include "FrameCore/SolveOptions.h"
#include "FrameCore/Material.h"
#include "FrameTestFixtures.h"

#include <cmath>

#if WITH_DEV_AUTOMATION_TESTS

namespace {
inline frame::Material s8ShellMat()
{
	const double E = 30000.0, nu = 0.3;
	frame::Material m(E, E / (2.0 * (1.0 + nu)));
	m.nu = nu;
	return m;
}
}

// ---- F48 mirror: QM6 incompatible-mode membrane (opt-in) ----
IMPLEMENT_SIMPLE_AUTOMATION_TEST(FFrameCoreShellIncompatibleMembraneTest,
	"FrameCore.Shell.IncompatibleMembrane",
	EAutomationTestFlags_ApplicationContextMask | EAutomationTestFlags::SmokeFilter)
bool FFrameCoreShellIncompatibleMembraneTest::RunTest(const FString&)
{
	using namespace frame;
	const Material smat = s8ShellMat();
	const double E = 30000.0, nu = 0.3;
	SolveOptions qm6; qm6.useIncompatibleMembrane = true;

	// (a) QM6 still passes the distorted constant-strain membrane patch (centre-Jacobian
	// correction keeps integral(B_inc)=0), so it reproduces the exact constant N.
	{
		const double a = 1000.0, t = 10.0, gx = 1e-4, f = E / (1.0 - nu * nu);
		FrameModel m; fixtures::membranePatch(m, a, t, 0.4, gx, smat);
		const SolveResult r = solve(m, qm6);
		TestFalse(TEXT("QM6 patch non-singular"), r.singular);
		const double NxxExp = t * f * gx, scale = FMath::Abs(NxxExp);
		double eNxx = 0;
		for (const auto& sf : r.shellForces) eNxx = FMath::Max(eNxx, FMath::Abs(sf.Nxx - NxxExp));
		TestTrue(TEXT("QM6 patch reproduces constant N"), eNxx < 1e-8 * scale);
	}

	// (b) slender in-plane cantilever (4x1): Q4 LOCKS (<60% of Euler-Bernoulli), QM6 RELEASES
	// (<15% from EB) and is strictly better than Q4.
	{
		const double L = 100.0, H = 10.0, t = 1.0, P = 1.0;
		const int nx = 4, ny = 1;
		const double I = t * H * H * H / 12.0, dEB = P * L * L * L / (3.0 * E * I);
		auto tip = [&](const SolveOptions& o) -> double {
			FrameModel m; fixtures::slenderMembraneCantilever(m, L, H, t, nx, ny, P, smat);
			const SolveResult r = solve(m, o);
			return r.singular ? 0.0 : FMath::Abs(r.disp(ny * (nx + 1) + nx, Uy));
		};
		const double dQ4 = tip(SolveOptions{}), dQM6 = tip(qm6);
		TestTrue(TEXT("Q4 membrane locks (<60% of Euler-Bernoulli)"), dQ4 < 0.60 * dEB);
		TestTrue(TEXT("QM6 releases lock (<15% from Euler-Bernoulli)"), FMath::Abs(dQM6 - dEB) < 0.15 * dEB);
		TestTrue(TEXT("QM6 strictly better than Q4"), FMath::Abs(dQM6 - dEB) < FMath::Abs(dQ4 - dEB));
	}
	return true;
}

// ---- F49 mirror: DKQ discrete-Kirchhoff thin-plate bending (opt-in) ----
IMPLEMENT_SIMPLE_AUTOMATION_TEST(FFrameCoreShellDKQPlateTest,
	"FrameCore.Shell.DKQPlate",
	EAutomationTestFlags_ApplicationContextMask | EAutomationTestFlags::SmokeFilter)
bool FFrameCoreShellDKQPlateTest::RunTest(const FString&)
{
	using namespace frame;
	const Material smat = s8ShellMat();
	const double E = 30000.0, nu = 0.3;
	SolveOptions dkq; dkq.useDKQPlate = true;

	// (a) DKQ constant-curvature patch (parallelogram) -> exact constant moment; Kirchhoff
	// means recovered transverse shears are exactly zero.
	{
		const double a = 1000.0, t = 10.0, c = 1e-6;
		const double Dfac = E * t * t * t / (12.0 * (1.0 - nu * nu));
		FrameModel m; fixtures::platePatchCylindrical(m, a, t, 0.4, c, smat);
		const SolveResult r = solve(m, dkq);
		TestFalse(TEXT("DKQ patch non-singular"), r.singular);
		const double MxxExp = -Dfac * c, scale = FMath::Abs(MxxExp);
		double eMxx = 0, mQ = 0;
		for (const auto& sf : r.shellForces) {
			eMxx = FMath::Max(eMxx, FMath::Abs(sf.Mxx - MxxExp));
			mQ = FMath::Max(mQ, FMath::Max(FMath::Abs(sf.Qx), FMath::Abs(sf.Qy)));
		}
		TestTrue(TEXT("DKQ patch reproduces constant moment"), eMxx < 1e-8 * scale);
		TestTrue(TEXT("DKQ Kirchhoff: recovered transverse shear is zero"), mQ < 1e-12);
	}

	// (b) simply-supported thin square plate (t/a=1/200): DKQ converges to the Kirchhoff centre
	// deflection 0.00406 q a^4 / D within 2% at N=16.
	{
		const double a = 1000.0, t = 5.0, q = 0.01;
		const double D = E * t * t * t / (12.0 * (1.0 - nu * nu));
		const double wc = 0.00406 * q * a * a * a * a / D;
		FrameModel m; fixtures::squarePlateShell(m, a, t, 16, q, smat);
		const SolveResult r = solve(m, dkq);
		TestFalse(TEXT("DKQ plate non-singular"), r.singular);
		const double w = FMath::Abs(r.disp(8 * 17 + 8, Uz));
		TestTrue(TEXT("DKQ within 2% of Kirchhoff"), FMath::Abs(w - wc) < 0.02 * wc);
	}
	return true;
}

#endif // WITH_DEV_AUTOMATION_TESTS
