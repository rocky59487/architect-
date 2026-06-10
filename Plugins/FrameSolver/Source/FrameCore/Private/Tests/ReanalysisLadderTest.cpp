// UE automation mirror of the standalone F35 (Tier-1) — the ReSolve Woodbury ladder reproduces a
// fresh assembleAndFactor+solveLoad after a member removal. Propped cantilever -> remove the prop ->
// statically determinate cantilever; the low-rank update on the baseline factor must match a fresh
// solve to factorization round-off, and report Tier-1.
#include "CoreMinimal.h"
#include "Misc/AutomationTest.h"
#include "FrameCore/FrameSolver.h"
#include "FrameCore/Reanalysis.h"

#if WITH_DEV_AUTOMATION_TESTS

IMPLEMENT_SIMPLE_AUTOMATION_TEST(FFrameCoreReanalysisLadderTest,
    "FrameCore.Reanalysis.LadderAgreesFresh",
    EAutomationTestFlags_ApplicationContextMask | EAutomationTestFlags::SmokeFilter)

bool FFrameCoreReanalysisLadderTest::RunTest(const FString& /*Parameters*/)
{
    using namespace frame;
    Material mat(210000.0, 80769.0, 7850.0);
    Section  sec = Section::Rectangular(120.0, 120.0);
    const real L = 3000.0, H = 2000.0, P = 5000.0;

    auto build = [&](FrameModel& m)
    {
        m = FrameModel{};
        m.materials = { mat }; m.sections = { sec };
        Node n0(0, 0.0, 0.0, 0.0);  n0.fixAll();
        Node n1(1,   L, 0.0, 0.0);
        Node n2(2,   L, 0.0,  -H);  n2.fixAll();
        m.nodes = { n0, n1, n2 };
        m.members = { Member(0, 0, 1, 0, 0), Member(1, 2, 1, 0, 0) };
        NodalLoad p; p.node = 1; p.comp[Uz] = -P; m.nodalLoads = { p };
    };

    FrameModel base; build(base);
    ReSolveSession s(base);
    TestTrue(TEXT("session valid"), s.valid());

    s.setMemberActive(1, false);                 // remove the prop
    ReanalysisStats st;
    const SolveResult re = s.solve(&st);
    TestTrue(TEXT("Tier-1 selected"), st.tier == 1);
    TestFalse(TEXT("ReSolve not singular"), re.singular);

    FrameModel w; build(w); w.members[1].active = false;
    const SolveResult fr = solve(w);
    double num = 0.0, den = 1e-30;
    for (int k = 0; k < (int)re.u.size() && k < (int)fr.u.size(); ++k)
    {
        num = FMath::Max(num, FMath::Abs(re.u[k] - fr.u[k]));
        den = FMath::Max(den, FMath::Abs(fr.u[k]));
    }
    TestTrue(TEXT("ReSolve == fresh (rel <= 1e-9)"), num / den <= 1e-9);
    return true;
}

#endif // WITH_DEV_AUTOMATION_TESTS
