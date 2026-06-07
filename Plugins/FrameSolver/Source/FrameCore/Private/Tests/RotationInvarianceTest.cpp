// M7 #7 verification framework (mirrors standalone F15): the solver must be rotationally
// EQUIVARIANT. Rotate the whole model (node positions, loads, and the member ref vector)
// by R; an encastre support is rotation-invariant, so the tip displacement must come out
// as exactly R times the un-rotated displacement. A strong anti-bug invariant for the
// coordinate-transform / assembly path.
#include "CoreMinimal.h"
#include "Misc/AutomationTest.h"
#include "FrameCore/FrameSolver.h"
#include "FrameCore/Section.h"
#include "FrameCore/Material.h"
#include "FrameTestFixtures.h"

#include <cmath>

#if WITH_DEV_AUTOMATION_TESTS

IMPLEMENT_SIMPLE_AUTOMATION_TEST(FFrameCoreRotationInvarianceTest,
	"FrameCore.Verify.RotationInvariance",
	EAutomationTestFlags_ApplicationContextMask | EAutomationTestFlags::SmokeFilter)
bool FFrameCoreRotationInvarianceTest::RunTest(const FString&)
{
	using namespace frame;
	const Section sec = Section::Rectangular(100.0, 100.0);
	Material mat(210000.0, 80769.0, 7850.0);
	mat.cap = Capacity::make(300.0, 300.0, 180.0);

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
	// material + section are referenced by index 0 (mat at materials[0], sec at sections[0])
	Node n0(0, 0, 0, 0); n0.fixAll();
	Node n1(1, rot(L,0,0,0), rot(L,0,0,1), rot(L,0,0,2));
	mr.nodes = { n0, n1 };
	Member mm(0, 0, 1, 0, 0);
	mm.refVec = Vec3(rot(0,0,1,0), rot(0,0,1,1), rot(0,0,1,2));
	mr.members = { mm };
	NodalLoad nl; nl.node = 1;
	nl.comp[Ux] = rot(0,0,P,0); nl.comp[Uy] = rot(0,0,P,1); nl.comp[Uz] = rot(0,0,P,2);
	mr.nodalLoads = { nl };
	const SolveResult rr = solve(mr);

	auto sq = [](real v) { return v * v; };
	const real ex = rot(u0[0], u0[1], u0[2], 0), ey = rot(u0[0], u0[1], u0[2], 1), ez = rot(u0[0], u0[1], u0[2], 2);
	const real du  = std::sqrt(sq(rr.disp(1,Ux) - ex) + sq(rr.disp(1,Uy) - ey) + sq(rr.disp(1,Uz) - ez));
	const real nrm = std::sqrt(sq(u0[0]) + sq(u0[1]) + sq(u0[2]));

	TestTrue(TEXT("rotated solve non-singular"), !rr.singular);
	TestTrue(TEXT("tip |u| preserved under rotation"),
		FMath::Abs(std::sqrt(sq(rr.disp(1,Ux)) + sq(rr.disp(1,Uy)) + sq(rr.disp(1,Uz))) - nrm) < 1e-9 * nrm);
	TestTrue(TEXT("u_rotated == R u_base (equivariance)"), du < 1e-6 * nrm);
	return true;
}

#endif // WITH_DEV_AUTOMATION_TESTS
