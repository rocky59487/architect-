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

// UE automation mirror of the standalone F51 -- S9b 3D GENERAL co-rotational. A cantilever along a TILTED
// axis (genuinely 3D: all three frame axes non-aligned) reproduces the elastica shooting table; an
// arbitrary-axis rigid rotation of the whole problem leaves the tip-displacement magnitude invariant (3D
// frame indifference, exercising expSO3/logSO3); and a tip torque gives the linear twist T L/(G J).
IMPLEMENT_SIMPLE_AUTOMATION_TEST(FFrameCoreCorotational3DTest,
    "FrameCore.Corotational.SpatialElastica",
    EAutomationTestFlags_ApplicationContextMask | EAutomationTestFlags::SmokeFilter)

bool FFrameCoreCorotational3DTest::RunTest(const FString& /*Parameters*/)
{
    using namespace frame;
    Material cmat(1.0, 0.4, 0.0); cmat.cap = Capacity::make(1e9, 1e9, 1e9);
    Section csec; csec.A = 100.0; csec.Iy = 1e-3; csec.Iz = 1e-3; csec.J = 2e-3;
    csec.cy = 1.0; csec.cz = 1.0; csec.Asy = 0.0; csec.Asz = 0.0;
    const real Lb = 1.0, Eb = 1.0, Ib = 1e-3;

    auto rotv = [](Vec3 v, Vec3 k, real ang) -> Vec3 {
        const real nk = FMath::Sqrt(k.x * k.x + k.y * k.y + k.z * k.z); k = Vec3(k.x / nk, k.y / nk, k.z / nk);
        const real c = FMath::Cos(ang), s = FMath::Sin(ang), kd = k.x * v.x + k.y * v.y + k.z * v.z;
        const Vec3 kxv(k.y * v.z - k.z * v.y, k.z * v.x - k.x * v.z, k.x * v.y - k.y * v.x);
        return Vec3(v.x * c + kxv.x * s + k.x * kd * (1 - c), v.y * c + kxv.y * s + k.y * kd * (1 - c), v.z * c + kxv.z * s + k.z * kd * (1 - c));
    };

    // (a) spatial elastica along a TILTED axis (1,1,1), symmetric section, perpendicular tip load.
    {
        const real na = FMath::Sqrt(3.0); const Vec3 axis(1, 1, 1);
        const real n2 = FMath::Sqrt(2.0); const Vec3 wN(1 / n2, -1 / n2, 0.0);
        struct Row { real alpha, dv; };
        const Row tab[] = { { 1.0, 0.3017207738 }, { 5.0, 0.7137915236 }, { 10.0, 0.8106090249 } };
        for (const Row& r : tab) {
            const real P = r.alpha * Eb * Ib / (Lb * Lb);
            FrameModel m; fixtures::cantileverSpatial(m, 16, Lb, axis, wN.x * P, wN.y * P, wN.z * P, 0, 0, 0, cmat, csec);
            CorotationalOptions co; co.loadSteps = FMath::Max(12, (int)(r.alpha * 3)); co.maxIter = 80;
            const CorotationalResult R = runCorotational(m, co);
            const Vec3 ut(R.finalState.u[(size_t)gdof(16, Ux)], R.finalState.u[(size_t)gdof(16, Uy)], R.finalState.u[(size_t)gdof(16, Uz)]);
            const real dv = (ut.x * wN.x + ut.y * wN.y + ut.z * wN.z) / Lb;
            TestTrue(TEXT("spatial elastica converged"), R.converged);
            TestTrue(TEXT("spatial elastica tip dv/L vs shooting table (rel<2e-3)"), FMath::Abs(dv - r.dv) <= 2e-3 * FMath::Abs(r.dv));
        }
    }

    // (b) arbitrary-axis rigid rotation invariance (machine precision).
    {
        const real P = 5.0 * Eb * Ib / (Lb * Lb);
        CorotationalOptions co; co.loadSteps = 15; co.maxIter = 80;
        FrameModel m0; fixtures::cantileverSpatial(m0, 12, Lb, Vec3(1, 0, 0), 0, P, 0, 0, 0, 0, cmat, csec);
        const CorotationalResult R0 = runCorotational(m0, co);
        const Vec3 axisN(1, 1, 1); const real phi = 2.0;
        FrameModel m1; fixtures::cantileverSpatial(m1, 12, Lb, Vec3(1, 0, 0), 0, P, 0, 0, 0, 0, cmat, csec);
        for (auto& nd : m1.nodes) nd.pos = rotv(nd.pos, axisN, phi);
        for (auto& mm : m1.members) mm.refVec = rotv(mm.refVec, axisN, phi);
        const Vec3 fr = rotv(Vec3(0, P, 0), axisN, phi);
        m1.nodalLoads[0].comp[Ux] = fr.x; m1.nodalLoads[0].comp[Uy] = fr.y; m1.nodalLoads[0].comp[Uz] = fr.z;
        const CorotationalResult R1 = runCorotational(m1, co);
        auto mag = [](const CorotationalResult& R) { const real ux = R.finalState.u[(size_t)gdof(12, Ux)], uy = R.finalState.u[(size_t)gdof(12, Uy)], uz = R.finalState.u[(size_t)gdof(12, Uz)]; return FMath::Sqrt(ux * ux + uy * uy + uz * uz); };
        TestTrue(TEXT("3D rotated solve converged"), R0.converged && R1.converged);
        TestTrue(TEXT("arbitrary-axis rotational invariance (rel<1e-9)"), FMath::Abs(mag(R1) - mag(R0)) <= 1e-9 * FMath::Abs(mag(R0)) + 1e-12);
    }

    // (c) pure torsion: tip torque -> twist T L/(G J).
    {
        Section tsec; tsec.A = 100.0; tsec.Iy = 2e-3; tsec.Iz = 1e-3; tsec.J = 1.5e-3; tsec.cy = 1.0; tsec.cz = 1.0; tsec.Asy = 0.0; tsec.Asz = 0.0;
        const real T = 1e-6;
        FrameModel m; fixtures::cantileverSpatial(m, 8, Lb, Vec3(1, 0, 0), 0, 0, 0, T, 0, 0, cmat, tsec);
        const CorotationalResult R = runCorotational(m, CorotationalOptions{});
        const real tw = R.finalState.u[(size_t)gdof(8, Rx)], ex = T * Lb / (0.4 * tsec.J);
        TestTrue(TEXT("torsion converged"), R.converged);
        TestTrue(TEXT("pure torsion theta=TL/GJ (rel<1e-6)"), FMath::Abs(tw - ex) <= 1e-6 * FMath::Abs(ex));
    }

    return true;
}

#endif // WITH_DEV_AUTOMATION_TESTS
