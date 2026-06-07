# FrameSolver plugin

FrameSolver contains `FrameCore`, the engine-agnostic structural mechanics core used by
the ArchSim UE host project and by the standalone validation tools.

## Current scope

- 3-D linear-elastic direct-stiffness solver.
- Euler-Bernoulli beam-column elements, with optional Timoshenko shear flexibility.
- Member end releases through static condensation.
- MITC4 Reissner-Mindlin flat-shell facets for plate/shell analysis.
- Legacy grillage plate idealization as a cheap beam-grid approximation.
- Elastic allowable-stress D/C screening.

Out of scope: geometric nonlinearity, dynamics/modal analysis, plastic collapse, and RC
fiber-section ultimate-strength design.

## Layout

- `Source/FrameCore/`: pure C++17 core. Public headers are POD/std-only; Eigen stays in
  private implementation files.
- `Source/FrameCore/Private/Tests/`: UE automation tests for `FrameCore.*`.
- `Standalone/`: console gates and CLI drivers used by the Python audits.

## Verification

Fast standalone gate:

```bat
Standalone\build.bat
```

Full gate from the repository root:

```powershell
powershell -ExecutionPolicy Bypass -File Scripts\run_gate.ps1 -RequireOpenSees
```

The full gate runs standalone analytic/benchmark fixtures, UE headless automation, and
OpenSees cross-validation. Use `-RequireOpenSees` in CI so a missing OpenSeesPy install
cannot silently skip the external comparison.
