# FrameCore — Architecture

FrameCore is a 3-D linear-elastic direct-stiffness beam-column FEM. This document covers the
conventions, data model, solve pipeline, element abstraction, strength screen, the grillage
idealization, and the validation framework. See [`../README.md`](../README.md) for the
high-level overview and scope boundaries.

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
| `Material` | `E`, `G`, density, `Capacity{comp,tens,shear,bend=min(comp,tens),tors=shear}` (allowable stresses) |
| `Section` | `A, Iy, Iz, J, cy, cz, Asy, Asz, shape`; section moduli `Wy()=Iy/cy`, `Wz()=Iz/cz`; factories `Rectangular(b,d)` and `Circular(r)` (Iy=Iz=πr⁴/4, J=2I); the factories also set the Timoshenko shear areas `Asy/Asz` |
| `Member` | `id`, end node ids `i,j`, `const Material*`, `const Section*`, `refVec`, `release[12]` |
| `NodalLoad` | node id + 6-component global force/moment |
| `MemberUDL` | member id + local distributed load `w_local` (force/length) |
| `FrameModel` | vectors of the above + `materials`/`sections` storage (keeps the pointers alive) + `validate()` + `dofCount()` |

**Pointer-lifetime invariant.** `Member` captures raw `const Material*`/`const Section*` into
`FrameModel::materials`/`sections`. Those vectors must be `reserve()`-d to their final size and
fully populated *before* any `Member` captures a pointer; otherwise a `push_back` reallocation
dangles every captured pointer. The fixtures and builders all follow "reserve → push → capture".

**`validate()`** rejects: no nodes/members, a member referencing a missing node, identical or
coincident endpoints, null material/section, non-positive `E/G` or `A/Iy/Iz/J`, and — after the
audit — **loads referencing a missing node or member** (these would otherwise be silently
dropped by the solver, yielding a quietly under-loaded "successful" solve).

---

## 3. Element abstraction (`IElement` seam)

The solver iterates over an `IElement` interface, not a hard-coded beam, so new element types
(e.g. future shells) drop in behind the same seam:

```
struct IElement {
  int  dofCount();
  void globalDofs(model, out);
  void assembleStiffness(model, opts, triplets);   // scatter k_global into the global K
  void fixedEndForces(model, opts, globalF);        // scatter element loads (FEF) into F
  bool condense(opts, why);                          // e.g. release static condensation
  void recoverEndForces(model, u, out);              // local {N,Vy,Vz,T,My,Mz} at both ends
};
```

`BeamColumnElement` is the only implementation today. Its local stiffness `localStiffness12`
is the textbook 12×12: axial `EA/L`; bending `12EI/L³, 6EI/L², 4EI/L, 2EI/L` in both planes
(the *y*-block carries sign flips relative to the *z*-block); torsion `GJ/L`. The transform is
`T = diag(R,R,R,R)` and `k_global = Tᵀ k_local T`. The optional Timoshenko variant scales the
bending block by `1/(1+Φ)` with `Φ = 12EI/(G·Aₛ·L²)` per plane and reduces to Euler–Bernoulli
as `Φ → 0`.

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
build IElement per member
assemble global F                    nodal loads + element fixed-end forces (scattered by Tᵀ)
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
over-estimated. This is the engine's bridge toward true shell elements.

---

## 7. Validation framework

- **`Private/FrameTestFixtures.h`** — pure-`frame` fixture builders (cantilever, simply
  supported, mechanism, axial column, propped-cantilever-via-release, torsion-release
  mechanism, circular arch). Shared by both the standalone gate and the UE automation tests, so
  a green standalone run exercises the *same* solver path UE compiles.
- **`Standalone/main.cpp`** — F1…F12 fixtures vs closed-form oracles (see README §validation),
  printing `[PASS]/[FAIL]` and `ALL PASS (failures=n)`.
- **`Standalone/frame_cli.cpp`** — a stdin/stdout solver driver used by the Python validation
  tools (it parses a model, solves, prints displacements + member forces).
- **`Private/Tests/*.cpp`** — UE automation mirrors (`FrameCore.*`), 12 tests, same oracles.
- **`Tools/`** — `opensees_compare.py` (OpenSees cross-validation, strict 1e-8 / `--relaxed`),
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
