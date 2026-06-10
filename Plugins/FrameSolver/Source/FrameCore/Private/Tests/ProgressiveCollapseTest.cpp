// UE automation mirror of the standalone F30 — the C2 progressive-collapse driver
// (collapse stage 3c). Re-proves the SAME closed-form sequence the standalone gate
// validates on the hanging-chain fixture: the thin link's D/C = 150000/(400*300) = 1.25
// exactly, it is removed as the governing member (Tension), the unsupported lower half
// detaches as debris with closed-form mass properties and takes its load with it, and the
// surviving link reads exactly zero force -> Stable in exactly two steps.
#include "CoreMinimal.h"
#include "Misc/AutomationTest.h"
#include "FrameCore/Collapse.h"

#if WITH_DEV_AUTOMATION_TESTS

IMPLEMENT_SIMPLE_AUTOMATION_TEST(FFrameCoreProgressiveCollapseTest,
    "FrameCore.Collapse.ProgressiveDriver",
    EAutomationTestFlags_ApplicationContextMask | EAutomationTestFlags::SmokeFilter)

bool FFrameCoreProgressiveCollapseTest::RunTest(const FString& /*Parameters*/)
{
    using namespace frame;
    Material mat(210000.0, 80769.0, 7850.0);
    mat.cap = Capacity::make(300.0, 300.0, 180.0);
    Section sec     = Section::Rectangular(100.0, 100.0);   // A = 1e4
    Section secThin = Section::Rectangular(20.0, 20.0);     // A = 400
    const real Lc = 1000.0, P = 150000.0;

    FrameModel m;
    m.materials = { mat };
    m.sections  = { sec, secThin };
    Node n0(0, 0.0, 0.0, 0.0); n0.fixAll();                  // hung from the top
    m.nodes = { n0, Node(1, 0.0, 0.0, -Lc), Node(2, 0.0, 0.0, -2.0 * Lc), Node(3, 0.0, 0.0, -3.0 * Lc) };
    m.members = { Member(0, 0, 1, 0, 0),
                  Member(1, 1, 2, 0, 1),                     // thin link governs: D/C = 1.25
                  Member(2, 2, 3, 0, 0) };
    NodalLoad p; p.node = 3; p.comp[Uz] = -P; m.nodalLoads = { p };

    CollapseOptions co; co.dlf = 1.0;
    const CollapseHistory h = runProgressiveCollapse(m, co);

    TestTrue(TEXT("outcome Stable in exactly 2 steps"),
             h.outcome == CollapseOutcome::Stable && h.steps.size() == 2);
    if (h.steps.size() != 2) return true;   // remaining asserts would index out of range

    TestTrue(TEXT("step0 maxDC = 1.25 (closed form)"),
             FMath::Abs(h.steps[0].maxDC - 1.25) <= 1e-9);
    const CollapseStep& s1 = h.steps[1];
    TestTrue(TEXT("step1 removes the thin link in Tension"),
             s1.removedMembers == std::vector<MemberId>({ 1 }) && s1.mode == FailMode::Tension);
    TestTrue(TEXT("step1 detaches the lower half as debris"),
             s1.detached.size() == 1 &&
             s1.detached[0].nodes == std::vector<NodeId>({ 2, 3 }) &&
             s1.detached[0].members == std::vector<MemberId>({ 2 }));
    if (s1.detached.size() == 1)
    {
        const double mExp = 7850.0 * 1e-12 * sec.A * Lc;     // 0.0785 tonne
        TestTrue(TEXT("debris mass = rho*A*L"), FMath::Abs(s1.detached[0].mass - mExp) <= 1e-12);
        TestTrue(TEXT("debris com z = -2500"), FMath::Abs(s1.detached[0].com.z + 2500.0) <= 1e-9);
    }
    TestTrue(TEXT("load leaves with the debris: surviving link reads D/C = 0"),
             s1.maxDC == 0.0 && s1.solved);

    // dlf is a pure load scale on a linear model
    CollapseOptions co2; co2.dlf = 2.0;
    const CollapseHistory h2 = runProgressiveCollapse(m, co2);
    TestTrue(TEXT("dlf=2 doubles the baseline D/C exactly"),
             !h2.steps.empty() && h2.steps[0].maxDC == 2.5);
    return true;
}

#endif // WITH_DEV_AUTOMATION_TESTS
