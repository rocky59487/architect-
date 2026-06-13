// UE automation mirror of the standalone F34 — the sparse subspace-iteration buckling path agrees
// with the dense eigensolver and the Euler analytic load. Forces the sparse path (denseThreshold<=0)
// so this exercises subspaceSmallest, not the dense fallback. A rectangular section keeps the two
// bending planes non-degenerate so the lowest eigenvalue is clean.
#include "CoreMinimal.h"
#include "Misc/AutomationTest.h"
#include "FrameCore/FrameSolver.h"
#include "FrameCore/BucklingAnalysis.h"
#include "FrameTestFixtures.h"

#if WITH_DEV_AUTOMATION_TESTS

IMPLEMENT_SIMPLE_AUTOMATION_TEST(FFrameCoreSparseBucklingTest,
    "FrameCore.Buckling.SparseAgreesDense",
    EAutomationTestFlags_ApplicationContextMask | EAutomationTestFlags::SmokeFilter)

bool FFrameCoreSparseBucklingTest::RunTest(const FString& /*Parameters*/)
{
    using namespace frame;
    Material mat(210000.0, 80769.0, 7850.0);
    Section  sec = Section::Rectangular(60.0, 100.0);     // Iy != Iz -> non-degenerate modes
    const double Ibuck = FMath::Min(sec.Iy, sec.Iz);
    const real   L = 3000.0, Pref = 1000.0;
    FrameModel m; fixtures::simplySupportedBeamN(m, 10, L, mat, sec);
    NodalLoad nl; nl.node = 10; nl.comp[Ux] = -Pref; m.nodalLoads = { nl };
    PreparedSystem ps = assembleAndFactor(m);

    const BucklingResult dense = solveBuckling(ps, m);
    BucklingOptions opt; opt.denseThreshold = 0;          // force the sparse path
    const BucklingResult sparse = solveBuckling(ps, m, opt);

    TestFalse(TEXT("dense non-singular"),  dense.singular);
    TestFalse(TEXT("sparse non-singular"), sparse.singular);
    TestTrue(TEXT("sparse == dense (rel <= 1e-6)"),
             FMath::Abs(sparse.criticalFactor - dense.criticalFactor) <= 1e-6 * FMath::Abs(dense.criticalFactor));
    const double Pcr   = sparse.criticalFactor * Pref;
    const double PcrEx = kPi * kPi * mat.E * Ibuck / (L * L);
    TestTrue(TEXT("sparse Pcr == Euler (rel <= 1e-4)"), FMath::Abs(Pcr - PcrEx) <= 1e-4 * PcrEx);
    return true;
}

#endif // WITH_DEV_AUTOMATION_TESTS
