# FrameCore — Architecture

FrameCore is a 3-D linear-elastic direct-stiffness FEM with beam-column **and** MITC4 flat-shell
elements. This document covers the conventions, data model, solve pipeline, element abstraction
(beam + shell), strength screen, the grillage idealization, and the validation framework. See
[`../README.md`](../README.md) for the high-level overview and scope boundaries.

---

## 1. Conventions (the authoritative contract)

- **Units:** N · mm · MPa (MPa = N/mm²). All checks use relative error, so the unit choice
  does not affect tolerances — but loads, sections and capacities must all be in this system.
- **Sign:** compression is **positive** for axial force `N`. Member end forces are in **local**
  coordinates; reactions are in **global** coordinates.
- **DOF order:** per node `[Ux, Uy, Uz, Rx, Ry, Rz]`. Global DOF index = `6·nodeIndex + localDof`
  (`gdof(node,d)`). An element's 12 DOFs are ordered `[node-i 6][node-j 6]`. The coordinate
  transform and the fixed-end-force vector use the **same** ordering — a mismatch would put a
  moment where a translation belongs. This ordering is declared once in `FrameTypes.h`.
- **Local axes:** local *x* runs i→j. The local *x*–*y* plane contains the member's `refVec`;
  local *z* completes the right-handed triad. A member parallel to its default `refVec`
  (a vertical column vs `+Z`) is degenerate, so `memberLocalAxes` falls back to a `+Y`
  reference. `Iz` governs bending in the local *x*–*y* plane, `Iy` the *x*–*z* plane.

---

## 2. Data model (`Public/FrameCore/`)

| Type | Holds |
|---|---|
| `Node` | `id`, position, `fixed[6]` (boundary mask), `prescribed[6]` (imposed support displacements); `fixAll()` |
| `Material` | `E`, `G`, **`nu`** (Poisson, used by the shell constitutive; beams ignore it), density, `Capacity{comp,tens,shear,bend=min(comp,tens),tors=shear}` (allowable stresses) |
| `Section` | `A, Iy, Iz, J, cy, cz, Asy, Asz, shape`; section moduli `Wy()=Iy/cy`, `Wz()=Iz/cz`; factories `Rectangular(b,d)` and `Circular(r)` (Iy=Iz=πr⁴/4, J=2I); the factories also set the Timoshenko shear areas `Asy/Asz` |
| `Member` | `id`, end node ids `i,j`, `const Material*`, `const Section*`, `refVec`, `release[12]` |
| `ShellQuad` | `id`, 4 corner node ids `n[4]` (CCW about +normal), `const Material*` (needs `nu`), thickness `t` — a MITC4 flat-shell facet |
| `NodalLoad` | node id + 6-component global force/moment |
| `MemberUDL` | member id + local distributed load `w_local` (force/length) |
| `ShellPressure` | shell id + transverse pressure `p` (along the facet normal) |
| `ShellElementForces` | per-shell stress resultants in the facet frame: `{Mxx,Myy,Mxy,Qx,Qy,Nxx,Nyy,Nxy}` (in `SolveResult.shellForces`) |
| `FrameModel` | vectors of the above (`members` **and** `shells` are parallel element sources) + `materials`/`sections` storage (keeps the pointers alive) + `validate()` + `dofCount()` |

**Pointer-lifetime invariant.** `Member` captures raw `const Material*`/`const Section*` into
`FrameModel::materials`/`sections`. Those vectors must be `reserve()`-d to their final size and
fully populated *before* any `Member` captures a pointer; otherwise a `push_back` reallocation
dangles every captured pointer. The fixtures and builders all follow "reserve → push → capture".

**`validate()`** rejects: no nodes, *no members **and** no shells*, a member referencing a missing
node, identical or coincident endpoints, null material/section, non-positive `E/G` or `A/Iy/Iz/J`;
for shells — a missing/duplicate corner node, null material, `nu ∉ [0,0.5)`, non-positive
thickness, a degenerate (zero-area) quad; and — after the audit — **loads (nodal, UDL, pressure)
referencing a missing node / member / shell** (these would otherwise be silently dropped by the
solver, yielding a quietly under-loaded "successful" solve).

---

## 3. Element abstraction (`IElement` seam)

The solver iterates over an `IElement` interface, not a hard-coded beam, so different element
types drop in behind the same seam (the solver assembles, applies BCs, factorizes, and recovers
the same way for every element — only `localDof()` differs):

```
struct IElement {
  int  localDof();                                   // 12 for a beam, 24 for a shell facet
  bool prepare(model, opts, why);                    // build local k + loads; condense releases
  void assemble(triplets);                           // scatter Tᵀ k_local T into the global K
  void addEquivalentNodalLoads(globalF);             // scatter element loads into F
  void recover(u, result);                           // element internal forces -> SolveResult
};
```

**`BeamColumnElement`** (12 DOF). Local stiffness `localStiffness12` is the textbook 12×12: axial
`EA/L`; bending `12EI/L³, 6EI/L², 4EI/L, 2EI/L` in both planes (the *y*-block carries sign flips
relative to the *z*-block); torsion `GJ/L`. The transform is `T = diag(R,R,R,R)` and
`k_global = Tᵀ k_local T`. The optional Timoshenko variant scales the bending block by `1/(1+Φ)`
with `Φ = 12EI/(G·Aₛ·L²)` per plane and reduces to Euler–Bernoulli as `Φ → 0`.

**`MITC4ShellElement`** (24 DOF = 4 nodes × 6). A flat Reissner–Mindlin facet built from three
blocks, each integrated 2×2 Gauss and mapped into the node DOFs `[Ux,Uy,Uz,Rx,Ry,Rz]`:

- **Plate bending** (fiber rotations `w, βx, βy`): `Kb = ∫ Bbᵀ Db Bb`, with the MITC4
  **assumed natural transverse shear** (Bathe–Dvorkin) — covariant shears sampled at the four
  edge-midpoint *tying points* and interpolated, which defeats shear locking in the thin limit.
  Fiber rotations map to the nodal rotations by `βx = Ry`, `βy = −Rx`, `w = Uz`.
- **Membrane** (plane stress `u, v`): `Km = ∫ t Bmᵀ Dm Bm`.
- **Drilling** (`Rz`): a **Hughes–Brezzi** penalty `∫ γ (θz − ½(v,x − u,y))²`, `γ = G·t`. It gives
  the in-plane rotation real stiffness — removing the coplanar `Rz` zero-energy mode that the
  LDLᵀ mechanism check would otherwise (correctly) flag singular — yet vanishes in constant-strain
  states, so it does not pollute the patch tests.

The facet's local frame is built from the corner geometry (`e1,e2` in-plane, `n` the averaged
normal); `T = diag(R,…)` (eight 3×3 blocks) rotates the 24-DOF local stiffness into 3-D, exactly
as for the beam. `recover` returns the **element-centre** stress resultants
`{Mxx,Myy,Mxy,Qx,Qy,Nxx,Nyy,Nxy}` (a single Gauss point at the centroid — the average for a
linearly-varying field, **not** nodal/peak values; no nodal extrapolation is performed).

**Known element-level trait (disclosed, monitored).** The local 24×24 stiffness carries the 6
true rigid-body zero modes **plus one inherent low-energy (near-zero, non-rigid) plate-bending
mode** — a standard MITC4 trait, *not* distortion-induced and *not* the drilling DOF. Adjacent
elements constrain it on assembly, so the drilling gate and every curved-shell benchmark stay
non-singular; it is documented rather than hidden, and pinned by an element-level spectrum oracle
in `linear_deep_audit` (asserts exactly 6 zero modes + the soft mode ≥10× softer than the first
true deformation mode).

**End releases.** With `SolveOptions.enableReleases`, a member's released local DOFs `c` are
statically condensed out of the retained set `r`: `k* = k_rr − k_rc k_cc⁻¹ k_cr` **and**
`Qf* = Qf_r − k_rc k_cc⁻¹ Qf_c` (condensing only `k` and not the fixed-end-force vector `Qf` is
the classic bug — a loaded member would report a phantom moment at the hinge). A singular
released sub-block (e.g. both torsional ends released → a free mechanism) is caught by a
rank-revealing factorization and reported as singular with a diagnostic, never as `NaN`.

---

## 4. Solve pipeline (`FrameSolver::solve`)

```
validate(model)                      -> singular + diagnostic on failure
build IElement per member + shell    BeamColumnElement per member, MITC4ShellElement per shell
assemble global F                    nodal loads + element equivalent/fixed-end loads (scattered by Tᵀ)
assemble global K                    triplets -> SimplicialLDLT-friendly sparse matrix
reduce to free DOFs                  drop fixed rows/cols; move prescribed (≠0) terms to the RHS
SimplicialLDLT.compute(K_ff)
mechanism detection                  factorization info != Success, OR a pivot |D_k| < eps·max|D|,
                                     OR a negative pivot (indefinite) -> SINGULAR + which DOF;
                                     driven by the factorization, NOT by graph connectivity
solve K_ff · u_f = F_f               (skipped / flagged if singular)
backfill u                           free + prescribed values
reactions R = K · u − F              (constrained rows carry the reactions)
recover member end forces            q_local = k_local·(T·u_e) + Qf*
```

Mechanism detection coming from the LDLᵀ factorization (rather than a connectivity graph) is a
deliberate rule: a member can be topologically connected yet form a kinematic mechanism, and
only the stiffness factorization sees that. NaN coordinates and near-zero-length members both
return `SINGULAR` rather than silent garbage.

---

## 4b. Analysis suite (factorize-once + the analysis modules)

`solve(model, opts)` is now a thin wrapper over a two-phase split that is the basis of every
analysis beyond one-shot static:

```
PreparedSystem assembleAndFactor(model, opts)   // EXPENSIVE, once: build K, reduce, LDLᵀ factor,
                                                //   mechanism check, bake distributed loads + geometry
SolveResult    solveLoad(prepared, model)       // CHEAP, many: rebuild RHS (nodal + prescribed) -> back-substitute
```

`PreparedSystem` is an opaque PIMPL (`Public/FrameCore/FrameSolver.h`); its Eigen-carrying body
(`Impl`: `K`, the free-DOF map, the LDLᵀ factorization, the prepared elements) lives in
`Private/PreparedSystemImpl.h`, shared by the analysis modules. Reusing the factorization makes
interactive (UE5) re-solves and multi-load-case work near-free; it is valid while geometry,
topology, support FLAGS and distributed loads are unchanged (nodal loads and prescribed VALUES
may vary — that is the interactive / settlement path). `solve()` stays **bit-for-bit identical**
to the previous monolithic solver (regression-checked, F19).

Two new `IElement` hooks feed the dynamic/stability analyses (default-empty, so existing elements
are unaffected): `assembleMass` (consistent mass `localMass12` / shell `shellMass24`) and
`assembleGeometric` (geometric stiffness `localGeometric12` from the prior solve's axial forces).
The analysis modules (each a free function + POD result, **no `solve()` flag bloat**):

| Module | Entry point | Math |
|---|---|---|
| Load combination / envelope | `combine`, `envelope` (`Combination.h`) | linear superposition; component-wise max/min |
| Self-weight | `addSelfWeight` (`SelfWeight.h`) | `w=ρgA` / `p=ρgt`; unit bridge ρ·1e-12 |
| Influence line | `reactionInfluenceLine` (`InfluenceLine.h`) | unit load marched, reusing the factorization |
| Modal | `solveModal` (`ModalAnalysis.h`) | generalized eigenproblem `Kφ=ω²Mφ` (dense) |
| Buckling | `solveBuckling` (`BucklingAnalysis.h`) | `(-Kg_ff)φ = γ K_ff φ`, λ_cr = 1/γ_max |
| Response spectrum | `solveResponseSpectrum` (`ResponseSpectrum.h`) | modal participation + SRSS/CQC |
| Transient | `solveModalStepResponse` (`ModalDynamics.h`) | Newmark-β per modal coordinate |

Units for mass/self-weight: the engine is consistent **N-mm-tonne-s**, so `Material.rho` (kg/m³)
is bridged by `×1e-12` (→ tonne/mm³). The self-weight and modal oracles validate this conversion
implicitly — it only matches `wL²/8` (F18) and `ωₙ=(nπ/L)²√(EI/ρA)` (F22) when the factor is right.

## 5. Strength screen (`ISectionStrength` → `ElasticAllowable`)

`checkSection(endForces, section, capacity) -> {risk, mode, sComp, sTens, tau, sTor}`:

- axial `sN = N/A` (compression-positive);
- biaxial bending `sM` — **rectangle:** conservative corner sum `|My|/Wy + |Mz|/Wz`;
  **circle:** exact resultant `√(My²+Mz²)/W`;
- combined `sComp = max(sM + sN, 0)`, `sTens = max(sM − sN, 0)`;
- transverse shear `tau = k·√(Vy²+Vz²)/A` with the **peak factor** `k = 1.5` (rectangle) or
  `4/3` (circle) — screening on the peak, not the average, so a shear-controlled member is not
  under-checked;
- torsion `sTor = |T|·c/J` — **circle:** exact `T·r/J` (`c = r`); **rectangle:** conservative
  diagonal corner `c = hypot(cy,cz)` (true St-Venant rectangular torsion needs warping, out of
  scope for an elastic screen);
- `risk` = the largest demand/capacity ratio across the five modes; `mode` = its argmax.

This is an **elastic / allowable-stress screen**, not RC ultimate strength.

---

## 6. Grillage idealization (`Grillage.h`) — the continuous-surface approximation

A simply-supported isotropic rectangular plate is idealized as a grid of longitudinal +
transverse beams (Hambly). Each strip of tributary width `w`, thickness `t` is **ν-inflated**
so the equivalent orthotropic plate matches the isotropic plate's `Dx = Dy = 2H = D`:

```
bending  I = w·t³ / [12 (1 − ν²)]      (E·I per width = D)
torsion  J = w·t³ / [6  (1 − ν)]       (G·J per width = D, with G = E/[2(1+ν)])
```

(The plain Hambly recipe `t³/12, t³/6` is ~26 % too flexible at ν=0.3; the inflations remove
that.) Only out-of-plane action is physical, so in-plane DOFs `(Ux,Uy,Rz)` are locked at every
node. Node index `node(i,j) = j·(nx+1) + i`. The uniform pressure is lumped to consistent nodal
loads, so total applied load equals `q·a·b` exactly (a load-conservation oracle). **Accuracy:**
center deflection within ~2 % of plate theory and mesh-stable; transverse moments
over-estimated. It is kept as a cheap approximation alongside the MITC4 shell (§3), which solves
plates/shells directly and converges to the exact solution.

---

## 7. Validation framework

- **`Private/FrameTestFixtures.h`** — pure-`frame` fixture builders (beams: cantilever, simply
  supported, mechanism, axial column, propped-cantilever-via-release, torsion-release mechanism,
  circular arch; shells: square/clamped plate, plate & membrane patch, Scordelis-Lo roof, pinched
  cylinder, rigid model rotation). Shared by both the standalone gate and the UE automation tests,
  so a green standalone run exercises the *same* solver path UE compiles.
- **`Standalone/main.cpp`** — F1…F25 fixtures vs closed-form oracles, benchmark references, and
  patch tests (see README §validation), printing `[PASS]/[FAIL]` and `ALL PASS (failures=n)`.
  F13–F16 cover the MITC4 shell; F17–F25 cover the linear-analysis suite (superposition + self-
  weight, factorize-once + settlement, pattern envelope, influence lines, modal, buckling,
  response spectrum, transient).
- **`Standalone/frame_cli.cpp`** — a stdin/stdout solver driver used by the Python validation
  tools (parses a model incl. `SMAT/SHELL/SPRESS` + node `prescribed` + `EIGEN`, solves, prints
  displacements + member/shell forces + `FREQ`).
- **`Private/Tests/*.cpp`** — UE automation mirrors (`FrameCore.*`), **26** tests (4 `Shell.*`,
  the load/solver/modal/buckling/response-spectrum/dynamics suite), same oracles.
- **`Tools/`** — `opensees_compare.py` (OpenSees cross-validation: beams strict 1e-8; prescribed
  settlement vs `sp()` to 0; the MITC4 shell vs OpenSees' own `ShellMITC4` to ~1e-10 on the
  flat/tilted plates gated here, ~1e-7–1e-8 on skewed/warped meshes in `shell_mitc4_deep_audit.py`;
  natural frequencies vs `eigen -cMass` to ~1e-11; `--relaxed` for cross-platform drift) — note
  these are *measured* agreements; the gate tolerances are looser on purpose (shell 1e-7, modal
  1e-4) to leave float/library headroom,
  `independent_precision_audit.py`, `complex_structure_benchmark.py`, `grillage_curve_audit.py`
  — all black-box the engine through `frame_cli.exe`.
- **`Scripts/run_gate.ps1`** — runs all three legs and prints a combined verdict + exit code
  (`-RequireOpenSees` to make a missing OpenSees a hard failure for CI).

---

## 8. Dual-build (standalone ⇄ Unreal)

The same `Private/*.cpp` compile in both targets. `Private/FrameEigen.h` is the **single** Eigen
include site: under `FRAMECORE_UE` it wraps Eigen in UE's third-party include guards; standalone
it includes Eigen plainly. The **public** API never exposes Eigen or UE types (POD boundary), and
cross-DLL symbols are tagged with `FRAMECORE_API` (expands to the UE export macro under UE, empty
standalone). Pure-core translation units never include `CoreMinimal.h`; only the UE module file
and the UE test files do.
