// UE automation mirror of the standalone F42 — tension-only (cable/brace eliminator). The active-
// set driver runs its inner re-solves through a ReSolveSession; under a lateral load on the X-braced
// portal one diagonal compresses, drops out as a tension-only member, and the converged state equals
// the model with that brace omitted entirely. Also gates the tensionOnly reuse fingerprint.
#include "CoreMinimal.h"
#include "Misc/AutomationTest.h"
#include "FrameCore/FrameSolver.h"
#include "FrameCore/TensionOnly.h"
#include "FrameTestFixtures.h"

#if WITH_DEV_AUTOMATION_TESTS

IMPLEMENT_SIMPLE_AUTOMATION_TEST(FFrameCoreTensionOnlyTest,
    "FrameCore.TensionOnly.Eliminator",
    EAutomationTestFlags_ApplicationContextMask | EAutomationTestFlags::SmokeFilter)

bool FFrameCoreTensionOnlyTest::RunTest(const FString& /*Parameters*/)
{
    using namespace frame;
    Material tmat(200000.0, 76923.07692307692, 7850.0); tmat.cap = Capacity::make(1e9, 1e9, 1e9);
    Section stocky = Section::Rectangular(250.0, 250.0);
    Section tbrace = Section::Rectangular(60.0, 60.0);

    FrameModel m; fixtures::xBracedPortal(m, 5.0e4, 0.0, tmat, stocky, tbrace);
    const TensionOnlyResult R = runTensionOnly(m);
    TestTrue(TEXT("converged"), R.converged);
    TestTrue(TEXT("converged in <=3 iters"), R.iterations <= 3);
    TestEqual(TEXT("exactly one brace slack"), (int)R.slack.size(), 1);

    const MemberId slackId = R.slack.empty() ? (MemberId)-1 : R.slack[0];
    const MemberId tenId   = (slackId == 3) ? 4 : 3;     // the surviving diagonal
    const double Nten = R.finalState.memberForces[(size_t)tenId].endI.N;
    TestTrue(TEXT("surviving brace in tension (N<0)"), Nten < 0.0);

    // oracle: the model with the slack (compressed) brace omitted from the member list entirely
    FrameModel ref; fixtures::xBracedPortal(ref, 5.0e4, 0.0, tmat, stocky, tbrace);
    for (size_t e = 0; e < ref.members.size(); ++e)
        if (ref.members[e].id == slackId) { ref.members.erase(ref.members.begin() + (long)e); break; }
    const SolveResult rr = solve(ref);
    double maxDiff = 0, maxU = 0;
    for (size_t g = 0; g < R.finalState.u.size() && g < rr.u.size(); ++g) {
        maxDiff = FMath::Max(maxDiff, FMath::Abs(R.finalState.u[g] - rr.u[g]));
        maxU    = FMath::Max(maxU, FMath::Abs(rr.u[g]));
    }
    TestTrue(TEXT("converged == omitted-brace (rel <= 1e-10)"), maxDiff <= 1e-10 * FMath::Max(maxU, 1e-30));

    // fingerprint: flipping tensionOnly rejects a stale factor
    FrameModel m2; fixtures::xBracedPortal(m2, 5.0e4, 0.0, tmat, stocky, tbrace);
    PreparedSystem ps2 = assembleAndFactor(m2);
    m2.members[3].tensionOnly = !m2.members[3].tensionOnly;
    const SolveResult stale = solveLoad(ps2, m2);
    TestTrue(TEXT("tensionOnly in fingerprint (stale rejected)"), stale.singular);
    return true;
}

#endif // WITH_DEV_AUTOMATION_TESTS
