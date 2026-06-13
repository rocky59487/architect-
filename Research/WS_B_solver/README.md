# FrameCore HP-FEM Solver Research Lane

A research-only investigation into breaking the sparse-direct (LDLT) factorization /
memory wall for repeated structural FEM solves, **without adding dependencies, without
changing FrameCore's public API, and without lowering the default solver's accuracy
contract**. The exact assembled `SimplicialLDLT` always remains the oracle and fallback.

> **Status:** research lane. It does **not** change the engine's default `solve()` and is
> **not** part of the five-leg engine gate. Every claim below is backed by an independent
> oracle (analytic / assembled-matrix / LDLT-relative) and a tolerance gate. Numbers are
> from an AMD Ryzen 9 8940HX (16 cores / 32 threads, 64 MB L3 across two CCDs).

---

## 1. Problem

A direct sparse `LDLT` factorization is robust and amortizes beautifully across many
right-hand sides — but its factorization time and fill-in grow super-linearly, and for
large gameplay-sized systems the factor wall dominates (measured: 18.7k DOF → 1.7 s;
48.4k DOF → 20.7 s; the existing scale ladder reaches >500 s and >10 GB beyond ~190k DOF).
A game engine re-solves a fixed (or slowly changing) structure every frame under
**low-dimensional** loads (gravity + a few contacts). The question: can a matrix-free
iterative lane beat a *reused* prefactored LDLT in the regimes that matter?

## 2. Two-lane architecture

* **Lane A — exact oracle / fallback (unchanged):** the assembled `SimplicialLDLT`. Used
  for small systems, singularity / mechanism detection, prescribed-DOF RHS, reactions, and
  any case where an iterative residual or pivot guard is inconclusive.
* **Lane B — HP iterative (this work):** a matrix-free element-by-element preconditioned
  conjugate gradient with
  * a persistent thread pool that runs both the element apply (per-thread touched-DOF
    reduction) and a parallel `6×6` block-Jacobi preconditioner map over the same warm
    threads;
  * a compact geometric coarse correction whose floor-aggregated matrix is **exactly
    block-tridiagonal**, solved by a hand-rolled block-Thomas `LDLᵀ` (`--coarseSolve
    banded`) or the engine's own sparse `SimplicialLDLT` (`--coarseSolve sparse`);
  * multi-RHS recycling: an A-orthogonal basis of recent solutions gives a Galerkin
    initial guess; a **seeded load-basis** path solves a few load modes first and freezes
    the response basis for low-dimensional parametric loads;
  * a projection-residual gate (`maxCombinedInitialRel`) and the LDLT oracle as the safety
    net.

## 3. Headline result — seeded HP for game-engine loads

For low-dimensional loads, after a one-time seed setup (measured `seedMs ≈ 1.3 s`, five
load-mode PCG solves) the seeded basis gives a near-perfect initial guess
(`maxCombinedInitialRel ≈ 1.1e-10`), so **each frame converges in ≈ 1 PCG iteration**:

| per-frame cost | reused LDLT back-substitution | seeded HP solve | accuracy vs LDLT |
| --- | --- | --- | --- |
| 18.7k DOF | ≈ 10.9 ms | ≈ 0.73 ms (**≈ 15×**) | `2e-12` |

Over a session of many frames the realized batch speedup over a *reused* LDLT climbs as the
one-time seed setup amortizes (single-session sweep, `--noSerialBaseline`; ~±20 % run-to-run
noise at low frame counts — see §4 for the pinned-statistics methodology):

| frames (RHS) | realized batch speedup vs reused LDLT |
| --- | --- |
| 200 | 2.37× |
| 1000 | 5.65× |
| 2000 | 8.04× |
| 4000 | 10.21× |

This is the production-relevant lane: for the load pattern games actually have, the *realized*
batch speedup over a reused prefactored LDLT reaches **10.2× at 4000 frames and keeps climbing
toward the ~15× per-frame ratio** (the asymptote, not a guaranteed batch number), at
effectively-exact accuracy.

## 4. Validated micro-architecture wins (pinned, 7-rep interleaved)

`--parallelPrecond` + `--coarseSolve banded` cut the preconditioner apply ~24 % and improve
every PCG-based path. `factor_bypass_batch_speedup`, median [min..max]:

| config | OLD (dense, serial precond) | NEW (parallel precond + banded) |
| --- | --- | --- |
| seeded 32 RHS | 1.211 [1.114..1.334] | **1.374 [1.347..1.578]** |
| non-seeded 8 RHS | 0.740 [0.657..1.017] | 0.826 [0.729..1.005] |

The seeded-NEW minimum (1.347) exceeds the seeded-OLD maximum (1.334) — the distributions do
not overlap. All keep `maxCombinedVsLdlt ≈ 2e-12`. (`factor_bypass_batch_speedup < 1.0` means
the HP batch is *slower* than a reused prefactored LDLT for that frame count; for the
non-seeded path the win is relative — NEW vs OLD — not absolute. Beating reused LDLT in
absolute terms for arbitrary loads is not the goal; see §6.)

## 5. Scaling — the factor-bypass case grows with size

| nf | LDLT factor | HP setup | factor-bypass setup speedup | apply bandwidth |
| --- | --- | --- | --- | --- |
| 18 720 | 1.71 s | 98.6 ms | 17.3× | 87 GB/s |
| 34 560 | 8.78 s | 285 ms | 30.8× | 96 GB/s |
| 48 384 | 20.68 s | 523 ms | 39.6× | 87 GB/s |

The LDLT factor explodes super-linearly (2.58× DOF → 12.1× factor time, the fill-in wall);
HP setup grows ~linearly, so the bigger the problem, the more the HP lane wins by avoiding
the factor. Apply bandwidth is a flat ~87–96 GB/s: the element data (≤30 MB) is **L3-resident**
(64 MB L3), so the apply is L3-bandwidth / ILP bound, not DRAM bound, at these sizes.

## 6. Honest negative / neutral results (and the principle behind them)

These are documented so they are not retried blindly. Each is **correct** (`vsLdlt ≈ 2e-12`)
but not a net win for this matrix-free, already-well-preconditioned solver:

* **Spectral deflation** (eigCG harvest + Type-1 Deflated CG): iteration reduction *saturates*
  at ~1.41× even with k = 24..48 and a Lanczos window of 600 — the preconditioned spectrum is
  spread, with no small isolated low-eigenspace to deflate. The dense per-iteration projection
  is memory-bound and roughly cancels the saving (best `speedupMs ≈ 0.98`).
* **Finer coarse** (sparse `SimplicialLDLT` path, 2×2 → 8×8): cuts iterations strongly
  (756 → 321) but the coarse *solve* cost grows with size via fill-in (precond 76 → 946 ms);
  the banded path shows the same trend (886 → 565 iters at 3×3/4×4). `combinedMs` is minimized
  at **2×2-banded** (~207 ms). A finer coarse wants recursion (a multigrid V-cycle), not a
  direct factor.
* **Symmetric apply** (store the lower triangle only): correct but neutral at these sizes —
  the element data is L3-resident, so halving it does not speed a non-DRAM-bound apply, and the
  symmetric loop's cross-update hurts ILP.
* **Seed recycling between load modes:** neutral — the load modes are too distinct.

**Unifying principle.** The element apply is cheap (L3-resident, ~0.14 ms/iter at 16T).
Therefore any method that adds per-iteration work of order one matvec (deflation, finer direct
coarse, multiplicative multigrid, Chebyshev smoothing) must cut iterations by *more* than its
added cost to pay — and a well-preconditioned spread spectrum only yields ~1.4–2× iteration
reductions, which the added cost cancels. The robust wins come from (1) making each existing
iteration cheaper *without* adding matvecs, and (2) exploiting structure *across* solves
(seeded recycling, factor-bypass).

**Fundamental limit.** A reused LDLT back-substitution (~11 ms; observed 10.5–12.4 ms across
configs) is ~18× cheaper than a full non-seeded PCG solve (~207 ms); "beat reused LDLT for
arbitrary many-RHS" is not achievable by iteration at this size, and is the wrong target. The HP lane's genuine niches are (a) few RHS / first
solve [factor avoidance], (b) low-dimensional parametric loads [seeded, the game-engine case],
and (c) very large problems where the LDLT factor is prohibitive.

## 7. Correctness & validation

* `exp_full_apply_oracle` (4/4): element-by-element full `K·x`, prescribed-DOF RHS,
  equivalent-load reactions and tower prescribed reactions vs the assembled `solveLoad` oracle.
* `exp_mechanism_guard_oracle` (4/4): under-constrained / release-condensation / isolated-node
  singular cases the HP lane must defer to the LDLT pivot guard.
* In-solver gates per run: `applyRel ≤ 1e-10` (matrix-free vs assembled), `maxCombinedVsLdlt ≤
  1e-6` and `maxCombinedTrueRel ≤ 1e-8` (every HP solution vs the LDLT solution and true
  residual), banded one-time self-check (`Kc · solve(probe) ≈ probe`, dense fallback if it
  fails), seeded projection-residual gate. Probes are mixed-frequency + sign-alternating.
* Independent adversarial review (read-only agent): no critical bugs; the thread-pool protocol,
  block-Thomas factorization, Type-1 DCG projection, and Rayleigh-Ritz were verified correct,
  and the deflation negative result is honestly supported.

## 8. Build & reproduce (Windows, from the repo root `ArchSim/`)

```bat
:: build the research experiments (reuses prebuilt FrameCore objs)
Research\build_research.bat -skipcore exp_parallel_pcg

:: correctness oracles
Research\bin\exp_full_apply_oracle.exe
Research\bin\exp_mechanism_guard_oracle.exe
```

Headline seeded asymptotic (the ~15×/frame result), and the validated micro-arch win:

```bash
# seeded HP vs reused LDLT, climbing batch speedup (bash)
EXE=Research/bin/exp_parallel_pcg.exe
SEED="--preset xxl --threads 16 --precond block6_coarse --coarseBinsX 2 --coarseBinsY 2 \
      --basisMax 8 --seedLoadBasis --parallelPrecond --coarseSolve banded --noSerialBaseline"
for n in 200 1000 2000 4000; do $EXE $SEED --rhs $n; done

# stable OLD-vs-NEW comparison (interleaved, median of 7 reps)
bash Research/WS_B_solver/bench_pinned.sh 7 factorBypassBatchSpeedup
```

Key flags on `exp_parallel_pcg` (all default **off** = original behaviour, so prior artifacts
reproduce): `--parallelPrecond`, `--coarseSolve banded|sparse`, `--deflation k
[--deflationWindow m]`, `--symApply`, `--seedLoadBasis`, `--noSerialBaseline`,
`--towerNx/Ny/Stories`.

## 9. File guide

| file | role |
| --- | --- |
| `exp_parallel_pcg.cpp` | the HP solver: matrix-free PCG, parallel block6 precond, banded/sparse coarse, deflation, seeded recycling, symmetric apply |
| `exp_full_apply_oracle.cpp` | full `K·x` / prescribed / reaction oracle vs assembled `solveLoad` |
| `exp_mechanism_guard_oracle.cpp` | singular / mechanism / fallback oracle |
| `run_hp_tuning.py` | benchmark / autotuning harness (parses the experiment metrics) |
| `bench_pinned.sh` | repeated / interleaved benchmark for stable OLD-vs-NEW statistics |
| `HPFEM_RESEARCH_NOTES.md` | detailed running research log (all measurements, dated) |
| `research_common.h` (`../common/`) | shared tower benchmark model + timing + reduced-system helpers |

## 10. References (algorithmic ideas, not copied code)

* Stathopoulos & Orginos, *Computing and deflating eigenvalues while solving multiple RHS with
  CG* (eigCG), SIAM J. Sci. Comput. 2010 — <https://www.cs.wm.edu/~andreas/publications/eigCG.pdf>
* Saad, Yeung, Erhel, Guyomarc'h, *A deflated version of the CG algorithm*, 2000 —
  <https://inria.hal.science/inria-00523686/document>
* Parks, de Sturler et al., *Recycling Krylov subspaces for sequences of linear systems*, 2006.

All algorithms here are implemented from published descriptions in pure C++17 + Eigen; no
external solver library is linked.
