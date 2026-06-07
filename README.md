# FrameCore — Structural Mechanics Engine

A self-contained **C++17 + Eigen** 3-D linear-elastic **beam-column finite-element engine**.
It is the structural core of an "architect simulator" graduation project: you describe a
frame (nodes, members, sections, supports, loads), and it returns nodal displacements,
member end forces, support reactions, a **mechanism / instability** verdict, and an elastic
**demand/capacity (D/C)** screen. It also computes **continuous surfaces directly** with a
**MITC4 Reissner–Mindlin flat-shell element** (membrane + plate bending + drilling), for true
plate/shell analysis; a legacy **grillage idealization** (a woven beam grid) is retained as a
cheap approximation alongside it.

The engine is deliberately small, **rigorously validated**, and **engine-agnostic**: the
public API uses only plain C++/POD types (no UE, no Eigen leakage), so the same source
compiles as a standalone console gate *and* as an Unreal Engine module.

> **Status:** clean, audited baseline. Verification gate is green
> (`standalone ALL PASS` · `16 UE automation tests` · `OpenSees strict cross-validation PASS`,
> including same-element agreement with OpenSees `ShellMITC4` to ~1e-10).

---

## What it computes

| Capability | Notes |
|---|---|
| 3-D linear-elastic direct-stiffness solve | 12×12 Euler–Bernoulli beam-column; sparse assembly; `SimplicialLDLT` |
| **Timoshenko** shear-flexible element | optional (`SolveOptions.useTimoshenko`); reduces to Euler–Bernoulli as slenderness grows |
| **End releases** / static indeterminacy | per-member `release[12]`; static condensation of the released sub-block (stiffness **and** fixed-end forces) |
| Sections | rectangular and **circular** (A, Iy, Iz, J, shear areas, extreme-fibre distances) |
| Loads | nodal forces/moments (global) + member UDL (local) |
| **Mechanism / instability detection** | from the LDLᵀ factorization (near-zero / negative pivots), **not** from connectivity — refuses to report forces on an unstable model |
| Member-end force / reaction recovery | local `{N,Vy,Vz,T,My,Mz}` at both ends; global reactions `R = K·u − F` |
| **Elastic D/C screen** (`ElasticAllowable`) | combined axial + biaxial bending + transverse shear (peak-factored) + torsion vs allowable capacities; reports the governing failure mode |
| **MITC4 flat-shell element** | 4-node Reissner–Mindlin facet: plane-stress membrane + plate bending with **MITC4 assumed transverse shear** (no shear locking) + Hughes–Brezzi **drilling**, assembled as 24 DOF and rotated into 3-D; recovers `{Mxx,Myy,Mxy,Qx,Qy,Nxx,Nyy,Nxy}` |
| **Grillage plate idealization** | a simply-supported isotropic plate → ν-inflated woven beam grid; a cheaper beam-grid approximation kept alongside the shell |

### Scope boundaries (read this — the engine is honest about what it is *not*)

- The D/C check is an **elastic / allowable-stress screen**, **not** RC ultimate strength.
  Transverse shear is screened on the **peak** stress (1.5·V/A rectangle, 4⁄3·V/A circle);
  rectangular **torsion** uses a conservative diagonal-corner heuristic (circular torsion is
  exact, `T·r/J`).
- The grillage is a **beam-grid approximation** of a plate: center deflection lands within
  ~2 % of Kirchhoff plate theory and is mesh-stable, but transverse bending moments are
  **over-estimated** (the grillage analogy trades the Poisson cross-moment for extra twist).
  It is an engineering idealization, not an exact plate solver.
- The **MITC4 shell** is a **4-node bilinear flat facet**: curved surfaces are approximated by
  flat panels, so there is a faceting error that vanishes under mesh refinement (the curved
  benchmarks report the convergence). It reproduces constant curvature/strain **exactly** on
  regular/parallelogram meshes; general (non-parallelogram) quadrilaterals show an O(h) patch
  residual that converges away. The drilling DOF is a **Hughes–Brezzi** treatment (it makes a
  coplanar shell non-singular and vanishes in constant-strain states, so it does not pollute
  the patch tests). It is **linear-elastic, small-deformation**.
- **No** dynamics, modal analysis, geometric nonlinearity (P-Δ / buckling), RC fiber-section
  precision, or plastic-collapse pushover. (Earlier experimental layers for these were removed
  in the cleanup to keep a small, fully-verified core.)

---

## Build & test

**Standalone gate (fastest — seconds):** compiles FrameCore + the analytic-oracle fixtures
and runs them.

```bat
Plugins\FrameSolver\Standalone\build.bat
```
Expected: each `[PASS] Fn …` then `ALL PASS (failures=0)`, exit 0. (Needs Visual Studio with
the C++ toolset; the script locates it via `vswhere`.)

**One-click full gate** (standalone + UE headless automation + OpenSees cross-validation):

```powershell
powershell -ExecutionPolicy Bypass -File Scripts\run_gate.ps1
# CI: -RequireOpenSees makes a missing openseespy a hard failure instead of a soft skip
```

**Unreal Engine** (the engine as a UE module): open `ArchSim.uproject`, or headless:

```bat
Engine\Build\BatchFiles\Build.bat ArchSimEditor Win64 Development -project=...\ArchSim.uproject
Engine\Binaries\Win64\UnrealEditor-Cmd.exe ...\ArchSim.uproject -ExecCmds="Automation RunTests FrameCore; Quit" -unattended -nullrhi -nopause
```

---

## How it is validated (the point of the project)

Every capability is anchored to an **independent oracle**, not just a self-consistent re-run:

- **Closed-form analytic solutions** — cantilever `PL³/3EI` & root moment `PL`; simply-supported
  UDL `5wL⁴/384EI`, `wL²/8`, `wL/2`; propped-cantilever-via-release `5wL/8`, `3wL/8`, `wL²/8`;
  Timoshenko `PL³/3EI + PL/GAₛ`; circular `I = πr⁴/4`; quarter-circle arch convergence;
  grillage center deflection vs Kirchhoff plate theory.
- **Shell oracles** — the MITC4 **patch test** (constant curvature *and* constant membrane
  strain reproduced to machine precision on regular & parallelogram meshes); simply-supported
  and clamped square plates vs Kirchhoff theory (~0.1–0.3 %); a thin-plate **no-shear-locking**
  check; and the two MacNeal–Harder curved-shell benchmarks — **Scordelis-Lo roof** (0.83 % at
  N=24) and **pinched cylinder** (98.8 % of the `1.8248e-5` reference at N=32, converging from
  below). A flat-plate-with-free-drilling case is the **coplanar non-singular** gate.
- **Rotation equivariance** — rotate the whole model by an arbitrary `R`; displacements must
  transform as `R·u` (catches transform / off-diagonal errors a norm-only check would miss).
  The shell variant rotates a clamped plate into an arbitrary 3-D orientation and checks the
  centre displacement magnitude is preserved (exercises the facet 3-D rotation).
- **An independent dense Gaussian solver** inside the gate (not Eigen) cross-checks results
  where it matters.
- **OpenSees** offline cross-validation (`Tools/opensees_compare.py`) — strict `1e-8` agreement
  by default for the beam models, `--relaxed` for cross-platform float drift. The MITC4 shell is
  cross-checked against OpenSees' **own `ShellMITC4`** (the same element): node displacements
  agree to **~1e-10**. OpenSees is a **validation tool only**; never shipped or linked.

See **[`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md)** for the data model, solve pipeline, sign
/ unit / DOF conventions, and the element abstraction.

---

## Minimal use (C++)

```cpp
#include "FrameCore/FrameSolver.h"
#include "FrameCore/ElasticAllowable.h"
using namespace frame;

Material mat(210000.0, 80769.0, 7850.0);            // E, G (MPa), rho
mat.cap = Capacity::make(300.0, 300.0, 180.0);      // allowable comp/tens/shear (MPa)
Section sec = Section::Rectangular(100.0, 100.0);   // b, d (mm)

FrameModel m;
m.materials = { mat };  m.sections = { sec };
const Material* pm = &m.materials.back();  const Section* ps = &m.sections.back();
Node n0(0, 0,0,0);  n0.fixAll();                    // encastre base
Node n1(1, 2000,0,0);                               // 2 m cantilever
m.nodes = { n0, n1 };
m.members = { Member(0, 0, 1, pm, ps) };
NodalLoad p;  p.node = 1;  p.comp[Uz] = -1000.0;    // 1 kN tip load
m.nodalLoads = { p };

SolveResult r = solve(m);                           // SolveOptions optional (releases / Timoshenko)
if (!r.singular) {
    double tip = r.disp(1, Uz);                     // = -PL^3/3EI
    DemandResult d = ElasticAllowable{}.checkSection(r.memberForces[0].endI, sec, mat.cap);
    // d.risk (D/C), d.mode (governing failure mode)
}
```

> **Pointer-lifetime rule:** `Member` holds `const Material*`/`const Section*` into the model's
> own vectors. `reserve()` those vectors to their final size *before* pushing, and capture the
> pointers only after all pushes — otherwise a reallocation dangles them.

---

## Repository layout

```
ArchSim.uproject                       UE host project (engine-as-module shell)
Plugins/FrameSolver/
  Source/FrameCore/                     the engine (pure C++17 + Eigen, UE-agnostic)
    Public/FrameCore/*.h                public API: FrameModel, Node, Material, Section,
                                        Member, Load, FrameSolver, SolveOptions/Result,
                                        ISectionStrength, ElasticAllowable, Grillage, ...
    Private/*.cpp                       implementation (+ Private/FrameEigen.h: the single
                                        Eigen include, dual-build guarded)
    Private/Tests/*.cpp                 UE automation tests (mirror the standalone oracles)
  Standalone/                           the console gate: build.bat (frametest), build_cli.bat
                                        (frame_cli, the OpenSees driver), main.cpp fixtures
Scripts/run_gate.ps1                    one-click standalone + UE + OpenSees gate
Tools/                                  engine-validation tools (all drive frame_cli.exe):
                                        opensees_compare, independent_precision_audit,
                                        complex_structure_benchmark, grillage_curve_audit
docs/ARCHITECTURE.md                    data model, solve pipeline, conventions
```

---

## Roadmap

The continuous-surface goal is met: a **MITC4 Reissner–Mindlin flat-shell element** sits behind
the existing `IElement` seam, so the engine computes plates and shells directly (cross-validated
against OpenSees `ShellMITC4`). Possible next steps, in rough order of value: a **prescribed
support-settlement** oracle (the math path exists but is not yet independently checked); a
warping-aware / higher-order shell to cut the flat-facet error on strongly curved surfaces; and,
further out, geometric nonlinearity (P-Δ / buckling), dynamics / modal analysis, or an RC
fiber-section precision layer — each gated by its own independent oracle before being claimed.

## License / use

Graduation-project code. FrameCore depends only on Eigen (MPL2, header-only). OpenSees is used
for offline validation only and is **not** redistributed or linked into the engine.
