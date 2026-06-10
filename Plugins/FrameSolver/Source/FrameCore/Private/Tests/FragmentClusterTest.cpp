// UE automation mirror of the standalone F29 — connectivity + debris mass properties
// (collapse stage 3b, the Chaos-handoff data path). Re-proves the SAME closed forms the
// standalone gate validates: a free rod's mass/com/inertia (rho*A*L, mL^2/12), the
// 45-degree rotation pinning the product-of-inertia SIGN convention (inertia[] stores
// tensor MATRIX entries, Ixy = -integral(xy dm)), and the twin-tower graph classification
// (grounded vs detached vs loose nodes).
#include "CoreMinimal.h"
#include "Misc/AutomationTest.h"
#include "FrameCore/Connectivity.h"

#if WITH_DEV_AUTOMATION_TESTS

IMPLEMENT_SIMPLE_AUTOMATION_TEST(FFrameCoreFragmentClusterTest,
    "FrameCore.Collapse.FragmentCluster",
    EAutomationTestFlags_ApplicationContextMask | EAutomationTestFlags::SmokeFilter)

bool FFrameCoreFragmentClusterTest::RunTest(const FString& /*Parameters*/)
{
    using namespace frame;
    Section  sec = Section::Rectangular(100.0, 100.0);
    Material mat(210000.0, 80769.0, 7850.0);
    const real L = 2000.0;
    const double mRod = 7850.0 * 1e-12 * sec.A * L;   // 0.157 tonne

    auto relOk = [](double got, double exp, double tol)
    {
        return FMath::Abs(got - exp) <= tol * FMath::Max(1.0, FMath::Abs(exp));
    };

    {   // free rod along +x
        FrameModel m;
        m.materials = { mat }; m.sections = { sec };
        m.nodes = { Node(0, 0.0, 0.0, 0.0), Node(1, L, 0.0, 0.0) };   // no supports
        m.members = { Member(0, 0, 1, 0, 0) };
        ConnectivityResult c = analyzeConnectivity(m);
        const bool shape = c.valid && c.detached.size() == 1 && c.groundedComponents == 0;
        TestTrue(TEXT("free rod: 1 detached, 0 grounded"), shape);
        if (!shape) return true;   // the remaining asserts would index an empty vector
        const FragmentCluster& f = c.detached[0];
        TestTrue(TEXT("rod mass = rho*A*L"), relOk(f.mass, mRod, 1e-12));
        TestTrue(TEXT("rod com x = L/2"), relOk(f.com.x, L / 2.0, 1e-12));
        TestTrue(TEXT("rod Iyy = mL^2/12"), relOk(f.inertia[1], mRod * L * L / 12.0, 1e-9));
        TestTrue(TEXT("rod Ixx ~ 0"), FMath::Abs(f.inertia[0]) <= 1e-9);
    }

    {   // 45-degree rod: product-of-inertia sign convention (tensor matrix entry)
        const real h = L / 2.0 * FMath::Sqrt(2.0);
        FrameModel m;
        m.materials = { mat }; m.sections = { sec };
        m.nodes = { Node(0, 0.0, 0.0, 0.0), Node(1, h, h, 0.0) };
        m.members = { Member(0, 0, 1, 0, 0) };
        ConnectivityResult c = analyzeConnectivity(m);
        TestTrue(TEXT("rotated rod: 1 detached"), c.valid && c.detached.size() == 1);
        const FragmentCluster& f = c.detached[0];
        TestTrue(TEXT("rot rod Ixx = mL^2/24"), relOk(f.inertia[0], mRod * L * L / 24.0, 1e-9));
        TestTrue(TEXT("rot rod Ixy = -mL^2/24"), relOk(f.inertia[3], -mRod * L * L / 24.0, 1e-9));
    }

    {   // twin towers + bridge: grounded / detached / loose classification
        const real Ht = 1000.0;
        FrameModel m;
        m.materials = { mat }; m.sections = { sec };
        Node a0(0, 0.0, 0.0, 0.0); a0.fixAll();
        m.nodes = { a0, Node(1, 0.0, 0.0, Ht), Node(2, 5000.0, 0.0, 0.0), Node(3, 5000.0, 0.0, Ht) };
        m.members = { Member(0, 0, 1, 0, 0), Member(1, 2, 3, 0, 0), Member(2, 1, 3, 0, 0) };

        FrameModel mCut = m; mCut.members[2].active = false;   // cut bridge, tower B unfooted
        ConnectivityResult c = analyzeConnectivity(mCut);
        TestTrue(TEXT("unfooted tower detaches"), c.valid && c.groundedComponents == 1 && c.detached.size() == 1);
        const FragmentCluster& f = c.detached[0];
        TestTrue(TEXT("cluster = tower B"),
                 f.nodes == std::vector<NodeId>({ 2, 3 }) && f.members == std::vector<MemberId>({ 1 }));
        TestTrue(TEXT("cluster mass = rho*A*H"), relOk(f.mass, 7850.0 * 1e-12 * sec.A * Ht, 1e-12));

        FrameModel mBare = mCut; mBare.members[1].active = false;   // kill tower B too
        ConnectivityResult c2 = analyzeConnectivity(mBare);
        TestTrue(TEXT("bare free nodes -> looseNodes"),
                 c2.detached.empty() && c2.looseNodes == std::vector<NodeId>({ 2, 3 }));
    }
    return true;
}

#endif // WITH_DEV_AUTOMATION_TESTS
