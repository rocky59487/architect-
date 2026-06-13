// UE automation mirror of the standalone F39 — analytic SDOF inheritance. An axial chain's mid
// member trips, the upper part detaches, and the retained piece (node0 encastre - member0 - node1)
// is a fixed-free axial bar (one analytic SDOF). The inheritance must hand node1 its pre-event
// static displacement exactly (full basis).
#include "CoreMinimal.h"
#include "Misc/AutomationTest.h"
#include "FrameCore/FrameSolver.h"
#include "FrameCore/DynamicCollapse.h"

#if WITH_DEV_AUTOMATION_TESTS

IMPLEMENT_SIMPLE_AUTOMATION_TEST(FFrameCoreDynCollapseOutcomesTest,
    "FrameCore.DynamicCollapse.Outcomes",
    EAutomationTestFlags_ApplicationContextMask | EAutomationTestFlags::SmokeFilter)

bool FFrameCoreDynCollapseOutcomesTest::RunTest(const FString& /*Parameters*/)
{
    using namespace frame;
    const real E = 200000.0, Gs = 80000.0, rho = 7850.0, L = 1000.0, P = 10000.0;
    Section sec = Section::Rectangular(100.0, 100.0);
    Material strong(E, Gs, rho); strong.cap = Capacity::make(1e9, 1e9, 1e9);
    Material weak(E, Gs, rho);   weak.cap   = Capacity::make(0.5, 0.5, 1e9);
    FrameModel m; m.materials = { strong, weak }; m.sections = { sec };
    Node n0(0, 0, 0, 0); n0.fixAll();
    Node n1(1, 0, 0, L); Node n2(2, 0, 0, 2 * L); Node n3(3, 0, 0, 3 * L);
    m.nodes   = { n0, n1, n2, n3 };
    m.members = { Member(0, 0, 1, 0, 0), Member(1, 1, 2, 1, 0), Member(2, 2, 3, 0, 0) };
    NodalLoad nl; nl.node = 3; nl.comp[Uz] = P; m.nodalLoads = { nl };

    const real rho_b = rho * 1e-12;
    const real omega = FMath::Sqrt(3.0 * E / (rho_b * L * L));
    const real T = 2.0 * kPi / omega;
    const real u1_static = P * L / (E * sec.A);

    DynCollapseOptions o; o.dt = T / 200.0; o.maxTime = 80 * o.dt; o.basisSize = 200; o.useRitzVectors = false;
    o.screenEvery = 1; o.frameStride = 1; o.removeThreshold = 1.0; o.maxEvents = 5;
    const DynCollapseHistory h = runDynamicCollapse(m, o);

    TestTrue(TEXT("exactly one event"), h.events.size() == 1);
    if (h.events.size() == 1) {
        TestTrue(TEXT("member1 removed"), h.events[0].removedMembers.size() == 1 && h.events[0].removedMembers[0] == 1);
        TestTrue(TEXT("upper part detached"), !h.events[0].detached.empty());
        const real te = h.events[0].t; int fi = -1;
        for (size_t k = 0; k < h.frames.size(); ++k) if (FMath::Abs(h.frames[k].t - te) < 0.5 * o.dt) fi = (int)k;
        TestTrue(TEXT("post-event frame exists"), fi >= 0);
        if (fi >= 0) {
            const real u1 = h.frames[(size_t)fi].u[(size_t)gdof(1, Uz)];
            TestTrue(TEXT("inherited u1 == static (analytic)"), FMath::Abs(u1 - u1_static) / FMath::Abs(u1_static) < 1e-9);
        }
    }
    return true;
}

#endif // WITH_DEV_AUTOMATION_TESTS
