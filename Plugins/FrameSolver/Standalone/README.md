# FrameCore standalone tools

This directory builds the UE-free FrameCore executables used for validation.

## Tools

- `build.bat` -> `frametest.exe`: standalone analytic and benchmark gate. It covers
  beam-column closed forms, releases, Timoshenko, grillage, MITC4 shell plate tests,
  Scordelis-Lo roof, and pinched cylinder.
- `build_cli.bat` -> `frame_cli.exe`: stdin/stdout solver driver used by Python audits
  and OpenSees comparisons.
- `build_perf.bat` -> `frame_perf.exe`: in-process performance benchmark.

## Quick start

From the repository root:

```bat
Plugins\FrameSolver\Standalone\build.bat
```

Expected result:

```text
ALL PASS  (failures=0)
```

The scripts locate Visual Studio through `vswhere` and compile against the UE-bundled
Eigen headers. The public FrameCore API remains UE-free; only the build scripts know the
local UE install path.
