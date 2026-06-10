// UE automation mirror of the standalone F31 — shell surface von Mises screen + driver
// shell removal (collapse stage 3d). Re-proves the SAME oracles the standalone gate
// validates: the closed-form synthetic screen (membrane vM = sqrt(10)/vm, face symmetry,
// corner governing), and the driver condemning the root facet of a cantilever strip via
// ShellVonMises with the unsupported tip facet detaching as closed-form debris.
#include "CoreMinimal.h"
#include "Misc/AutomationTest.h"
#include "FrameCore/Collapse.h"
#include "FrameCore/ElasticAllowable.h"

#if WITH_DEV_AUTOMATION_TESTS

IMPLEMENT_SIMPLE_AUTOMATION_TEST(FFrameCoreShellFailureTest,
    "FrameCore.Collapse.ShellFailure",
    EAutomationTestFlags_ApplicationContextMask | EAutomationTestFlags::SmokeFilter)

bool FFrameCoreShellFailureTest::RunTest(const FString& /*Parameters*/)
{
    using namespace frame;
    const real Es = 30000.0, nu = 0.3, Gs = Es / (2.0 * (1.0 + nu));
    Material smat(Es, Gs, 2500.0); smat.nu = nu;
    smat.cap = Capacity::make(10.0, 10.0, 6.0);   // vm = 10

    {   // pure-function screen oracles (no solve)
        const real t = 10.0;
        ShellElementForces f{};                    // membrane only: sigma = (3,1,1)
        f.Nxx = 30.0; f.Nyy = 10.0; f.Nxy = 10.0;
        const ShellDemandResult d1 = checkShellSurface(f, t, smat.cap);
        TestTrue(TEXT("membrane vM = sqrt(10)/vm"),
                 FMath::Abs(d1.risk - FMath::Sqrt(10.0) / 10.0) <= 1e-12);
        TestTrue(TEXT("constant field: centre governs"), d1.corner == -1);

        ShellElementForces cnr{};                  // a single hot corner governs
        cnr.MxxC[2] = 50.0;                        // 6M/t^2 = 3 at corner 2 only
        const ShellDemandResult d3 = checkShellSurface(cnr, t, smat.cap);
        TestTrue(TEXT("corner peak vM = 0.3 at corner 2"),
                 FMath::Abs(d3.risk - 0.3) <= 1e-12 && d3.corner == 2);
    }

    {   // driver shell removal: root facet condemned, tip facet detaches, Collapsed
        const real e = 500.0, t = 10.0;
        FrameModel m;
        m.materials = { smat };
        Node n0(0, 0.0, 0.0, 0.0); n0.fixAll();
        Node n3(3, 0.0, e, 0.0);   n3.fixAll();
        m.nodes = { n0, Node(1, e, 0.0, 0.0), Node(2, 2.0 * e, 0.0, 0.0),
                    n3, Node(4, e, e, 0.0),   Node(5, 2.0 * e, e, 0.0) };
        m.shells = { ShellQuad(0, 0, 1, 4, 3, 0, t),       // root facet
                     ShellQuad(1, 1, 2, 5, 4, 0, t) };     // tip facet
        NodalLoad p2; p2.node = 2; p2.comp[Uz] = -500.0;
        NodalLoad p5; p5.node = 5; p5.comp[Uz] = -500.0;
        m.nodalLoads = { p2, p5 };

        CollapseOptions co; co.dlf = 1.0;
        const CollapseHistory h = runProgressiveCollapse(m, co);
        TestTrue(TEXT("outcome Collapsed in exactly 2 steps"),
                 h.outcome == CollapseOutcome::Collapsed && h.steps.size() == 2);
        if (h.steps.size() != 2) return true;

        const CollapseStep& s1 = h.steps[1];
        TestTrue(TEXT("step1 removes the root facet via ShellVonMises"),
                 s1.removedShells == std::vector<int>({ 0 }) && s1.removedMembers.empty() &&
                 s1.mode == FailMode::ShellVonMises);
        TestTrue(TEXT("tip facet detaches as debris"),
                 s1.detached.size() == 1 && s1.detached[0].shells == std::vector<int>({ 1 }));
        if (s1.detached.size() == 1)
        {
            const double mExp = 2500.0 * 1e-12 * t * e * e;   // 0.00625 tonne
            TestTrue(TEXT("debris mass = rho*t*a^2"), FMath::Abs(s1.detached[0].mass - mExp) <= 1e-12);
            TestTrue(TEXT("debris com x = 750"), FMath::Abs(s1.detached[0].com.x - 750.0) <= 1e-9);
        }
    }
    return true;
}

#endif // WITH_DEV_AUTOMATION_TESTS
