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

        // Shared helpers are hoisted to header-only inlines (FrameTypes.h, PreparedSystemImpl.h,
        // FrameEigen.h, MemberQuery.h, CollapseSupport.h, ElementFactory.h) — no anonymous-
        // namespace name collisions remain across TUs, so unity build is safe.
        bUseUnity = true;

        PublicDependencyModuleNames.AddRange(new string[] { "Core" });

        // UE-bundled Eigen 3.4.0 (header-only, MPL2). Private so consumers of the
        // FrameCore public API never see Eigen.
        AddEngineThirdPartyPrivateStaticDependencies(Target, "Eigen");

        // Flip FrameEigen.h to the UE-guarded include path (THIRD_PARTY_INCLUDES_*).
        PublicDefinitions.Add("FRAMECORE_UE=1");
        PublicDefinitions.Add("EIGEN_MPL2_ONLY");
    }
}
