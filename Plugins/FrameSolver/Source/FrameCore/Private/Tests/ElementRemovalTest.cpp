// UE automation mirror of the standalone F26 — element removal / collapse foundation (C1 + C5).
// Re-proves the SAME solver path the standalone gate validates: a propped cantilever drops to a
// statically determinate cantilever when the prop member is deactivated (Member::active=false),
// recovering the exact closed-form tip deflection -PL^3/3EI; deactivation is identical to omitting
// the member; and removing every member at a node is flagged as a mechanism (singular).
#include "CoreMinimal.h"
#include "Misc/AutomationTest.h"
#include "FrameCore/FrameSolver.h"
#include "FrameTestFixtures.h"

#if WITH_DEV_AUTOMATION_TESTS

IMPLEMENT_SIMPLE_AUTOMATION_TEST(FFrameCoreElementRemovalTest,
    "FrameCore.Member.ElementRemoval",
    EAutomationTestFlags_ApplicationContextMask | EAutomationTestFlags::SmokeFilter)

bool FFrameCoreElementRemovalTest::RunTest(const FString& /*Parameters*/)
{
    using namespace frame;
    Section  sec = Section::Rectangular(100.0, 100.0);
    Material mat(210000.0, 80769.0, 7850.0);
    const real L = 2000.0, H = 1000.0, P = 1000.0;
    const double dExp = -P * L * L * L / (3.0 * mat.E * sec.Iz);   // determinate cantilever tip deflection

    auto buildPropped = [&](FrameModel& m)
    {
        m = FrameModel{};
        m.materials = { mat }; m.sections = { sec };
        Node n0(0, 0.0, 0.0, 0.0); n0.fixAll();   // encastre base
        Node n1(1,   L, 0.0, 0.0);                // free tip
        Node n2(2,   L, 0.0,  -H); n2.fixAll();   // prop foot below the tip
        m.nodes = { n0, n1, n2 };
        m.members = { Member(0, 0, 1, 0, 0),      // cantilever beam
                      Member(1, 2, 1, 0, 0) };    // vertical prop
        NodalLoad p; p.node = 1; p.comp[Uz] = -P; m.nodalLoads = { p };
    };

    FrameModel mFull; buildPropped(mFull);
    SolveResult rFull = solve(mFull);
    TestFalse(TEXT("propped cantilever not singular"), rFull.singular);

    // remove the prop -> determinate cantilever -> tip deflection = -PL^3/3EI
    FrameModel mCut = mFull; mCut.members[1].active = false;
    SolveResult rCut = solve(mCut);
    TestFalse(TEXT("not singular after removing prop"), rCut.singular);
    TestTrue(TEXT("tip deflection = -PL^3/3EI"),
             FMath::Abs(rCut.disp(1, Uz) - dExp) <= 1e-6 * FMath::Abs(dExp));
    TestTrue(TEXT("removed member force is zero"),
             FMath::Abs(rCut.memberForces[1].endI.N) <= 1e-9);

    // active=false MUST be identical to physically omitting the member
    FrameModel mOmit; buildPropped(mOmit); mOmit.members = { mOmit.members[0] };
    SolveResult rOmit = solve(mOmit);
    double duMax = 0.0;
    for (int k = 0; k < (int)rCut.u.size() && k < (int)rOmit.u.size(); ++k)
        duMax = FMath::Max(duMax, FMath::Abs(rCut.u[k] - rOmit.u[k]));
    TestTrue(TEXT("active=false identical to omitted member"), duMax <= 1e-9);

    // remove every member at the tip -> isolated node -> mechanism
    FrameModel mMech = mFull; mMech.members[0].active = false; mMech.members[1].active = false;
    SolveResult rMech = solve(mMech);
    TestTrue(TEXT("isolated node -> mechanism (singular)"), rMech.singular);
    return true;
}

#endif // WITH_DEV_AUTOMATION_TESTS
