// UE automation mirror of the standalone F40 — P-Delta second-order amplification. The frozen
// pseudo-load reuse path and the K_T refactor reference reach the same second-order tip sway, the
// reference matches the beam-column stability-function closed form, P=0 is the linear solve
// bit-for-bit, and past the Euler load both paths flag divergence (no silent wrong answer).
#include "CoreMinimal.h"
#include "Misc/AutomationTest.h"
#include "FrameCore/FrameSolver.h"
#include "FrameCore/PDeltaAnalysis.h"
#include "FrameTestFixtures.h"

#if WITH_DEV_AUTOMATION_TESTS

IMPLEMENT_SIMPLE_AUTOMATION_TEST(FFrameCorePDeltaTest,
    "FrameCore.PDelta.AmplificationOracle",
    EAutomationTestFlags_ApplicationContextMask | EAutomationTestFlags::SmokeFilter)

bool FFrameCorePDeltaTest::RunTest(const FString& /*Parameters*/)
{
    using namespace frame;
    const real   Lc = 6000.0, sideC = 200.0, Hlat = 1000.0, Ec = 200000.0;
    const real   Ic = sideC * sideC * sideC * sideC / 12.0;
    const double kPiC = 3.14159265358979323846;
    const real   Pcr = kPiC * kPiC * Ec * Ic / (4.0 * Lc * Lc);
    Section  pcs = Section::Rectangular(sideC, sideC);
    Material pcm(Ec, Ec / (2.0 * (1.0 + 0.3)), 7850.0); pcm.cap = Capacity::make(1e9, 1e9, 1e9);
    const int nE = 8;

    // f = 0.3: frozen == reference, and reference == beam-column exact
    {
        const real P = 0.3 * Pcr;
        FrameModel m; fixtures::pdeltaColumn(m, nE, Lc, P, Hlat, pcm, pcs);
        PDeltaOptions of;  of.refactorPath = false; of.maxIter = 5000; of.tolU = 1e-13;
        PDeltaOptions orf; orf.refactorPath = true;
        const PDeltaResult rf = runPDelta(m, of);
        const PDeltaResult rr = runPDelta(m, orf);
        TestTrue(TEXT("frozen converged"),    rf.converged);
        TestTrue(TEXT("reference converged"), rr.converged);
        const double tf = rf.finalState.disp(nE, Ux), tr = rr.finalState.disp(nE, Ux);
        TestTrue(TEXT("frozen == reference (rel <= 1e-10)"), FMath::Abs(tf - tr) <= 1e-10 * FMath::Abs(tr));
        const double kk = FMath::Sqrt(P / (Ec * Ic));
        const double tipExact = Hlat * (FMath::Tan(kk * Lc) - kk * Lc) / (P * kk);
        TestTrue(TEXT("reference == beam-column exact (rel <= 1e-3)"),
                 FMath::Abs(tr - tipExact) <= 1e-3 * FMath::Abs(tipExact));
    }

    // P = 0: the second-order state is the linear solve, bit-for-bit, in zero iterations
    {
        FrameModel m; fixtures::pdeltaColumn(m, nE, Lc, 0.0, Hlat, pcm, pcs);
        const PDeltaResult rf = runPDelta(m, PDeltaOptions{});
        const SolveResult lin = solve(m);
        TestTrue(TEXT("P=0 zero iterations"), rf.iterations == 0);
        TestTrue(TEXT("P=0 bit-identical to linear"), rf.finalState.disp(nE, Ux) == lin.disp(nE, Ux));
    }

    // f = 1.05 (beyond Euler): both paths must report divergence
    {
        const real P = 1.05 * Pcr;
        FrameModel m; fixtures::pdeltaColumn(m, nE, Lc, P, Hlat, pcm, pcs);
        PDeltaOptions of;  of.refactorPath = false; of.maxIter = 5000; of.tolU = 1e-13;
        PDeltaOptions orf; orf.refactorPath = true;
        const PDeltaResult rf = runPDelta(m, of);
        const PDeltaResult rr = runPDelta(m, orf);
        TestTrue(TEXT("reference diverged (K_T not PD)"), rr.diverged);
        TestTrue(TEXT("frozen diverged (sliding window)"), rf.diverged);
    }

    return true;
}

#endif // WITH_DEV_AUTOMATION_TESTS
