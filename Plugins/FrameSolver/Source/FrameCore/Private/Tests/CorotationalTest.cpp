// UE automation mirror of the standalone F50 -- S9 planar co-rotational large displacement. A
// transverse end-loaded cantilever's tip deflection matches the independent elastica shooting table
// (Bisshopp-Drucker / Mattiasson), an in-plane rigid rotation of the whole model leaves the tip
// displacement magnitude invariant (co-rotational frame indifference -> zero spurious rigid-rotation
// force), and at small displacement the CR sway degenerates to the linearized P-Delta amplification.
#include "CoreMinimal.h"
#include "Misc/AutomationTest.h"
#include "FrameCore/FrameSolver.h"
#include "FrameCore/CorotationalAnalysis.h"
#include "FrameCore/PDeltaAnalysis.h"
#include "FrameTestFixtures.h"

#if WITH_DEV_AUTOMATION_TESTS

IMPLEMENT_SIMPLE_AUTOMATION_TEST(FFrameCoreCorotationalTest,
    "FrameCore.Corotational.ElasticaCantilever",
    EAutomationTestFlags_ApplicationContextMask | EAutomationTestFlags::SmokeFilter)

bool FFrameCoreCorotationalTest::RunTest(const FString& /*Parameters*/)
{
    using namespace frame;
    const double kPiC = 3.14159265358979323846;

    // E=1, near-inextensible slender section (A large, I small) so the inextensible elastica table holds.
    Material cmat(1.0, 0.4, 0.0); cmat.cap = Capacity::make(1e9, 1e9, 1e9);
    Section csec; csec.A = 100.0; csec.Iy = 1e-3; csec.Iz = 1e-3; csec.J = 1e-3;
    csec.cy = 1.0; csec.cz = 1.0; csec.Asy = 0.0; csec.Asz = 0.0;
    const real Lb = 1.0, Eb = 1.0, Ib = 1e-3;

    // (a) elastica table (WS_F F-3): alpha -> tip dv/L. N=16 mesh, within 1.5e-3.
    struct Row { real alpha, dv; };
    const Row tab[] = { { 1.0, 0.3017207738 }, { 5.0, 0.7137915236 }, { 10.0, 0.8106090249 } };
    for (const Row& r : tab)
    {
        FrameModel m; fixtures::cantileverPlanarTipShearN(m, 16, Lb, r.alpha * Eb * Ib / (Lb * Lb), cmat, csec);
        CorotationalOptions co; co.loadSteps = FMath::Max(12, (int)(r.alpha * 3)); co.maxIter = 80;
        const CorotationalResult R = runCorotational(m, co);
        const double dv = R.finalState.u[(size_t)gdof(16, Uy)] / Lb;
        TestTrue(TEXT("elastica converged"), R.converged);
        TestTrue(TEXT("elastica tip dv/L matches shooting table (rel<1.5e-3)"),
                 FMath::Abs(dv - r.dv) <= 1.5e-3 * FMath::Abs(r.dv));
    }

    // (b) in-plane rigid-rotation frame indifference: rotate model+load by phi -> |u_tip| unchanged.
    {
        const real alpha = 3.0, phi = 0.6, P = alpha * Eb * Ib / (Lb * Lb);
        FrameModel m0; fixtures::cantileverPlanarTipShearN(m0, 12, Lb, P, cmat, csec);
        CorotationalOptions co; co.loadSteps = 15; co.maxIter = 80;
        const CorotationalResult R0 = runCorotational(m0, co);
        FrameModel m1; fixtures::cantileverPlanarTipShearN(m1, 12, Lb, 0.0, cmat, csec);
        const real cph = FMath::Cos(phi), sph = FMath::Sin(phi);
        for (auto& nd : m1.nodes) { const real x = nd.pos.x, y = nd.pos.y; nd.pos.x = cph * x - sph * y; nd.pos.y = sph * x + cph * y; }
        NodalLoad nl; nl.node = 12; nl.comp[Ux] = -sph * P; nl.comp[Uy] = cph * P; m1.nodalLoads = { nl };
        const CorotationalResult R1 = runCorotational(m1, co);
        auto mag = [](const CorotationalResult& R) { const double ux = R.finalState.u[(size_t)gdof(12, Ux)], uy = R.finalState.u[(size_t)gdof(12, Uy)]; return FMath::Sqrt(ux * ux + uy * uy); };
        const double m0mag = mag(R0), m1mag = mag(R1);
        TestTrue(TEXT("rotated solve converged"), R0.converged && R1.converged);
        TestTrue(TEXT("in-plane rotational invariance (|u_tip| preserved, rel<1e-9)"),
                 FMath::Abs(m1mag - m0mag) <= 1e-9 * FMath::Abs(m0mag) + 1e-12);
    }

    // (c) P-Delta degeneration: small lateral load -> CR sway matches linearized runPDelta.
    {
        const int nE = 6; const real Hc = 1.0;
        const real Pcr = kPiC * kPiC * Eb * Ib / (4.0 * Hc * Hc);
        const real Paxial = 0.3 * Pcr, Hlat = 1e-4;
        auto build = [&](FrameModel& m) {
            fixtures::prepMatSec(m, cmat, csec); m.nodes.clear(); m.members.clear();
            for (int k = 0; k <= nE; ++k) { Node nd(k, 0, Hc * real(k) / nE, 0); nd.fixed[Uz] = nd.fixed[Rx] = nd.fixed[Ry] = true; if (k == 0) nd.fixAll(); m.nodes.push_back(nd); }
            for (int k = 0; k < nE; ++k) m.members.push_back(Member(k, k, k + 1, 0, 0));
            NodalLoad nl; nl.node = nE; nl.comp[Uy] = -Paxial; nl.comp[Ux] = Hlat; m.nodalLoads = { nl };
        };
        FrameModel mc; build(mc); CorotationalOptions co; co.loadSteps = 12; co.maxIter = 80;
        const CorotationalResult Rcr = runCorotational(mc, co);
        FrameModel mp; build(mp); PDeltaOptions po; po.refactorPath = true;
        const PDeltaResult Rpd = runPDelta(mp, po);
        const double swayCR = Rcr.finalState.u[(size_t)gdof(nE, Ux)];
        const double swayPD = Rpd.finalState.u[(size_t)gdof(nE, Ux)];
        TestTrue(TEXT("P-Delta degeneration converged"), Rcr.converged && Rpd.converged);
        TestTrue(TEXT("CR sway == P-Delta sway (small-disp, rel<1.5e-2)"),
                 FMath::Abs(swayCR - swayPD) <= 1.5e-2 * FMath::Abs(swayPD));
    }

    return true;
}

#endif // WITH_DEV_AUTOMATION_TESTS
