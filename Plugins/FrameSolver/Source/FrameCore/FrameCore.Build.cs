using UnrealBuildTool;

// FrameCore — engine-agnostic structural solver. The Private *.cpp files use Eigen
// only through Private/FrameEigen.h; Public headers stay Eigen-free (POD + std).
// NOTE: this module is authored drop-in-ready but is NOT compiled in milestone 1
// (no host .uproject yet). The standalone harness under Standalone/ is the gate.
public class FrameCore : ModuleRules
{
    public FrameCore(ReadOnlyTargetRules Target) : base(Target)
    {
        PCHUsage    = PCHUsageMode.UseExplicitOrSharedPCHs;
        CppStandard = CppStandardVersion.Cpp20;   // UE 5.7 / VS2026 no longer allow C++17
        bEnableExceptions = true;   // Eigen may throw / assert

        // Each analysis .cpp keeps its OWN anonymous-namespace numeric helpers (reduceFF, kPi, ...).
        // Those names repeat across translation units BY DESIGN — they are file-local. A unity/jumbo
        // build merges several .cpp into one TU, which turns those file-local helpers into COLLIDING
        // namespace-scope symbols (redefinition + C4459 shadow-as-error), and which files land in the
        // same blob depends on file count/order — so adding a new analysis file can silently detonate
        // a pre-existing latent collision. Disable unity for this pure-numeric module: the helpers
        // stay clean, the build is order-independent, and the compile-time cost here is negligible.
        bUseUnity = false;

        PublicDependencyModuleNames.AddRange(new string[] { "Core" });

        // UE-bundled Eigen 3.4.0 (header-only, MPL2). Private so consumers of the
        // FrameCore public API never see Eigen.
        AddEngineThirdPartyPrivateStaticDependencies(Target, "Eigen");

        // Flip FrameEigen.h to the UE-guarded include path (THIRD_PARTY_INCLUDES_*).
        PublicDefinitions.Add("FRAMECORE_UE=1");
        PublicDefinitions.Add("EIGEN_MPL2_ONLY");
    }
}
