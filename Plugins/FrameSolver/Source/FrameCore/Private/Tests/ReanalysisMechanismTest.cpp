// UE automation mirror of the standalone F35 mechanism case — removing the base member of a
// 2-member cantilever chain makes the Woodbury capacitance singular, so ReSolveSession flags a
// mechanism (singular), matching a fresh assembleAndFactor on the same removed set.
#include "CoreMinimal.h"
#include "Misc/AutomationTest.h"
#include "FrameCore/FrameSolver.h"
#include "FrameCore/Reanalysis.h"

#if WITH_DEV_AUTOMATION_TESTS

IMPLEMENT_SIMPLE_AUTOMATION_TEST(FFrameCoreReanalysisMechanismTest,
    "FrameCore.Reanalysis.MechanismDetection",
    EAutomationTestFlags_ApplicationContextMask | EAutomationTestFlags::SmokeFilter)

bool FFrameCoreReanalysisMechanismTest::RunTest(const FString& /*Parameters*/)
{
    using namespace frame;
    Material mat(210000.0, 80769.0, 7850.0);
    Section  sec = Section::Rectangular(120.0, 120.0);

    FrameModel chain;
    chain.materials = { mat }; chain.sections = { sec };
    Node c0(0, 0, 0, 0);       c0.fixAll();
    Node c1(1, 0, 0, 1500.0);
    Node c2(2, 0, 0, 3000.0);
    chain.nodes = { c0, c1, c2 };
    chain.members = { Member(0, 0, 1, 0, 0), Member(1, 1, 2, 0, 0) };
    NodalLoad nl; nl.node = 2; nl.comp[Ux] = 1000.0; chain.nodalLoads = { nl };

    ReSolveSession s(chain);
    TestTrue(TEXT("baseline valid"), s.valid());
    s.setMemberActive(0, false);                 // remove the base -> top floats
    ReanalysisStats st;
    const SolveResult re = s.solve(&st);
    TestTrue(TEXT("capacitance mechanism flagged"), st.mechanism);
    TestTrue(TEXT("ReSolve singular"), re.singular);

    FrameModel w = chain; w.members[0].active = false;
    TestTrue(TEXT("fresh singular too"), solve(w).singular);
    return true;
}

#endif // WITH_DEV_AUTOMATION_TESTS
