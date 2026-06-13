// UE automation mirror of the standalone F54 -- the S10 N-M interaction plastic hinge. A
// hinge-capable member's plastic moment is reduced under axial force by Mp_eff(N) =
// Mp*(1-(N/Ny)^2) (CollapseOptions::nmInteraction). The fixture is the F54 determinate
// X-cantilever: base moment w L^2/2 and axial force = tip load P, both stiffness-independent,
// so the N-M hinge load is bracketed cleanly. Under P = Ny/2 the reduced capacity is 0.75 Mp;
// 1.01 w* tips it into a base-hinge mechanism (Collapsed) while the SAME load with N-M off
// stays Stable -- the axial interaction is precisely what causes collapse.
#include "CoreMinimal.h"
#include "Misc/AutomationTest.h"
#include "FrameCore/Collapse.h"
#include "FrameCore/NMInteraction.h"

#if WITH_DEV_AUTOMATION_TESTS

IMPLEMENT_SIMPLE_AUTOMATION_TEST(FFrameCoreNMInteractionTest,
    "FrameCore.Collapse.NMInteraction",
    EAutomationTestFlags_ApplicationContextMask | EAutomationTestFlags::SmokeFilter)

bool FFrameCoreNMInteractionTest::RunTest(const FString& /*Parameters*/)
{
    using namespace frame;
    Material mat(210000.0, 80769.0, 7850.0);
    mat.fy = 300.0;
    mat.cap = Capacity::make(300.0, 300.0, 180.0);
    const real b = 80.0, d = 120.0;
    Section sec = Section::Rectangular(b, d);
    const real MpZ = mat.fy * sec.Zz;
    const real Ny  = mat.fy * sec.A;

    // formula spot checks (machine precision)
    TestTrue(TEXT("Mp_eff(0)=Mp"), FMath::Abs(reducedPlasticMoment(MpZ, 0.0, Ny) - MpZ) <= 1e-9 * MpZ);
    TestTrue(TEXT("Mp_eff(Ny/2)=0.75Mp"),
             FMath::Abs(reducedPlasticMoment(MpZ, 0.5 * Ny, Ny) - 0.75 * MpZ) <= 1e-9 * MpZ);
    TestTrue(TEXT("Mp_eff clamps to 0 beyond squash"), reducedPlasticMoment(MpZ, 2.0 * Ny, Ny) == 0.0);

    const real L = 2000.0;
    const real P = 0.5 * Ny;                              // Mp_eff(P) = 0.75 MpZ
    const real MpEff = reducedPlasticMoment(MpZ, P, Ny);
    const real wStar = 2.0 * MpEff / (L * L);

    auto buildCant = [&](FrameModel& m, real w)
    {
        m = FrameModel{};
        m.materials = { mat }; m.sections = { sec };
        Node n0(0, 0.0, 0.0, 0.0); n0.fixAll();
        m.nodes = { n0, Node(1, L, 0.0, 0.0) };
        m.members = { Member(0, 0, 1, 0, 0) };
        MemberUDL u; u.member = 0; u.w_local = { 0.0, -w, 0.0 };
        m.memberUDLs = { u };
        NodalLoad p; p.node = 1; p.comp[Ux] = -P;
        m.nodalLoads = { p };
    };

    CollapseOptions on;  on.dlf  = 1.0; on.plasticHinges  = true; on.nmInteraction  = true;
    CollapseOptions off; off.dlf = 1.0; off.plasticHinges = true; off.nmInteraction = false;

    FrameModel mLo; buildCant(mLo, 0.99 * wStar);
    const CollapseHistory hLo = runProgressiveCollapse(mLo, on);
    TestTrue(TEXT("0.99 w* (N-M on): Stable, no hinge"),
             hLo.outcome == CollapseOutcome::Stable && hLo.steps.size() == 1 &&
             hLo.steps[0].formedHinges.empty());

    FrameModel mHi; buildCant(mHi, 1.01 * wStar);
    const CollapseHistory hHi = runProgressiveCollapse(mHi, on);
    TestTrue(TEXT("1.01 w* (N-M on): base hinge -> Collapsed"),
             hHi.outcome == CollapseOutcome::Collapsed && hHi.steps.size() == 2 &&
             hHi.steps[1].formedHinges.size() == 1 && !hHi.steps[1].solved);
    if (hHi.steps.size() == 2 && hHi.steps[1].formedHinges.size() == 1)
    {
        TestTrue(TEXT("hinge at base end-i, Mz dof"),
                 hHi.steps[1].formedHinges[0].member == 0 && hHi.steps[1].formedHinges[0].dof == 5);
        TestTrue(TEXT("trigger ratio = |M|/Mp_eff = 1.01"),
                 FMath::Abs(hHi.steps[1].triggerRatio - 1.01) <= 1e-9);
        TestTrue(TEXT("frozen residual |Mp| = Mp_eff(P)"),
                 FMath::Abs(FMath::Abs(hHi.steps[1].formedHinges[0].Mp) - MpEff) <= 1e-6 * MpEff);
    }

    FrameModel mHiOff; buildCant(mHiOff, 1.01 * wStar);
    const CollapseHistory hOff = runProgressiveCollapse(mHiOff, off);
    TestTrue(TEXT("same load, N-M off: Stable (axial interaction is decisive)"),
             hOff.outcome == CollapseOutcome::Stable && hOff.steps.size() == 1 &&
             hOff.steps[0].formedHinges.empty());
    return true;
}

#endif // WITH_DEV_AUTOMATION_TESTS
