// UE automation mirror of the standalone F33 — the event-to-event hinge driver (collapse
// stage 4b). Re-proves the SAME classic plastic-collapse bracket the standalone gate
// validates on the fixed-fixed beam (Mp = fy*Zz): at 0.98 w* (w* = 16 Mp/L^2) the run ends
// Stable after exactly two support hinges; at 1.02 w* the third (midspan) hinge forms and
// the run goes singular -> Collapsed. The analytic collapse load is bracketed to +/-2%.
#include "CoreMinimal.h"
#include "Misc/AutomationTest.h"
#include "FrameCore/Collapse.h"

#if WITH_DEV_AUTOMATION_TESTS

IMPLEMENT_SIMPLE_AUTOMATION_TEST(FFrameCoreProgressiveHingeTest,
    "FrameCore.Collapse.ProgressiveHinge",
    EAutomationTestFlags_ApplicationContextMask | EAutomationTestFlags::SmokeFilter)

bool FFrameCoreProgressiveHingeTest::RunTest(const FString& /*Parameters*/)
{
    using namespace frame;
    Material mat(210000.0, 80769.0, 7850.0);
    mat.fy = 300.0;
    mat.cap = Capacity::make(300.0, 300.0, 180.0);
    Section sec = Section::Rectangular(100.0, 100.0);    // Zz = 250000 -> Mp = 7.5e7
    const real Lh = 2000.0, Ltot = 2.0 * Lh;
    const real Mp = mat.fy * sec.Zz;
    const real wStar = 16.0 * Mp / (Ltot * Ltot);        // 75 N/mm

    auto buildBeam = [&](FrameModel& m, real w)
    {
        m = FrameModel{};
        m.materials = { mat }; m.sections = { sec };
        Node n0(0, 0.0, 0.0, 0.0);  n0.fixAll();
        Node n2(2, Ltot, 0.0, 0.0); n2.fixAll();
        m.nodes = { n0, Node(1, Lh, 0.0, 0.0), n2 };
        m.members = { Member(0, 0, 1, 0, 0), Member(1, 1, 2, 0, 0) };
        MemberUDL u0; u0.member = 0; u0.w_local = { 0.0, -w, 0.0 };
        MemberUDL u1; u1.member = 1; u1.w_local = { 0.0, -w, 0.0 };
        m.memberUDLs = { u0, u1 };
    };

    CollapseOptions co; co.dlf = 1.0; co.plasticHinges = true;

    FrameModel mA; buildBeam(mA, 0.98 * wStar);
    const CollapseHistory hA = runProgressiveCollapse(mA, co);
    TestTrue(TEXT("0.98 w*: Stable after exactly two hinges"),
             hA.outcome == CollapseOutcome::Stable && hA.steps.size() == 3 &&
             hA.steps[1].formedHinges.size() == 1 && hA.steps[2].formedHinges.size() == 1);
    if (hA.steps.size() == 3)
    {
        TestTrue(TEXT("hinge order: support i of member 0, then support j of member 1"),
                 hA.steps[1].formedHinges[0].member == 0 && hA.steps[1].formedHinges[0].dof == 5 &&
                 hA.steps[2].formedHinges[0].member == 1 && hA.steps[2].formedHinges[0].dof == 11);
        TestTrue(TEXT("hinge 1 ratio = 0.98*16/12"),
                 FMath::Abs(hA.steps[1].triggerRatio - 0.98 * 16.0 / 12.0) <= 1e-9);
        TestTrue(TEXT("Stable with reported D/C above 1 (ductile bending waits at Mp)"),
                 hA.steps[2].maxDC > 1.0);
    }

    FrameModel mB; buildBeam(mB, 1.02 * wStar);
    const CollapseHistory hB = runProgressiveCollapse(mB, co);
    TestTrue(TEXT("1.02 w*: third hinge -> mechanism -> Collapsed"),
             hB.outcome == CollapseOutcome::Collapsed && hB.steps.size() == 4 &&
             hB.steps[3].formedHinges.size() == 1 && !hB.steps[3].solved);
    if (hB.steps.size() == 4)
        TestTrue(TEXT("hinge 3 ratio = 1.02*2 - 1"),
                 FMath::Abs(hB.steps[3].triggerRatio - (1.02 * 2.0 - 1.0)) <= 1e-9);
    return true;
}

#endif // WITH_DEV_AUTOMATION_TESTS
