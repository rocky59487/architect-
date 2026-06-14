# FrameCore HP-FEM Solver Research Notes

Date: 2026-06-13

Scope: research-only experiments for breaking the sparse direct LDLT memory wall
without changing the public FrameCore API, adding dependencies, or lowering the
default solver accuracy contract.

## Feasibility Judgment

The direct sparse LDLT path is robust and should remain the default oracle, but
it does not scale as the primary path for large gameplay-sized FEM systems. The
research results support a two-lane architecture:

1. Keep assembled LDLT as the exact/fallback path for gates, singularity checks,
   prescribed DOF handling, reactions, and small systems.
2. Add an opt-in high-performance iterative lane for repeated large solves:
   matrix-free element apply, tuned preconditioners, persistent thread-local
   reduction, and multi-RHS recycling.

The strongest near-term wins are not from scalar CSR-to-BSR replacement alone.
They come from reducing data motion, reusing structure across solves, and
amortizing thread/local-buffer setup.

## Current Evidence

All numbers below are research tower `xxl` unless noted. They are single-run
measurements and should be treated as directional until repeated on pinned cores.

| Experiment | Result |
| --- | --- |
| Direct LDLT scaling | Existing scale ladder reaches fill-in wall: 186k free DOF around 522s, 389k free DOF around 3229s and >10GB. |
| Matrix-free apply correctness | `applyRel` around `2e-16` versus assembled `Kff*x`. |
| PCG baseline | Jacobi: 1249 iterations. |
| 6x6 block Jacobi | 1225 iterations; small benefit only. |
| Floor coarse + block6 | 1107 iterations; coarse space is the real contributor. |
| 2x2 per-floor coarse + block6 | 967 iterations, but dense coarse cost is close to the break-even point. |
| 3x3/4x4 coarse | Fewer iterations, worse wall time due to dense coarse cache footprint. |
| Real FrameCore BSR6 apply | Around 1.0-1.1x versus Eigen sparse on real `Kff`; synthetic regular BSR was about 2.4x, so real sparsity/layout dominates. |
| Thread-local persistent apply | Small case is slower; `xxl` sees about 1.7x at 4-8 threads after removing per-apply thread spawn. |
| Persistent-thread PCG | Full Jacobi-PCG on `xxl` keeps about 1.5-1.6x at 4 threads; 8 threads regresses, showing memory/reduction wall. |
| Parallel PCG + coarse | `block6_coarse` still reduces iterations, but its serial dense preconditioner cost eats part of apply parallelism. |
| Multi-RHS recycling | `basisMax=8` gives about 1.24-1.26x PCG time speedup for correlated 12-RHS sequence, with LDLT-relative error around `1e-12`. |
| Combined HP solver loop | `block6_coarse + persistent 4-thread apply + basisMax=8 recycling` on `xxl` gave about `1.8x` average solve-time speedup over serial zero-start with the same preconditioner, while staying around `2e-12` vs LDLT. |
| Sparse touched reduction | Persistent pool now clears/reduces only per-thread touched DOFs. On `xxl`, `8T/1x1/basis8` improved average combined solve from about `298.851 ms` to `255.488 ms`; apply time dropped to `140.639 ms`. |
| Factor-bypass first solve | Non-seeded `xxl` direct `assembleAndFactor` was around `1700 ms`; factor-free HP setup was around `97 ms`, so setup/factor bypass remains about `17x`. With `16T/2x2/basis8`, first-solve bypass was about `4.8x` without seed and about `1.0x` with seed cost included. |
| Factor-bypass short batch | Non-seeded 8-RHS autotune reached `0.932x` at `16T + 2x2 coarse + basis8`; it still cannot beat reused LDLT for short repeated loads. |
| Seeded load-basis batch | For parametric low-dimensional loads, solving 5 load-mode responses first and freezing the A-orthogonal basis changes the regime: 32 RHS reached `1.228x` factor-bypass batch speedup; 200 RHS reached `2.244x`, with `maxCombinedVsLdlt` around `2e-12`. This is valid only when the runtime loads remain in the seeded load subspace. |
| Seed projection gate | The seeded 32-RHS artifact reports `maxCombinedInitialRel=1.031e-10`; one PCG correction brings true residual to about `4.8e-11`. Production must keep this residual gate and fall back when it grows. |
| PCG timing breakdown | Best non-seeded 8-RHS autotune (`16T/2x2/basis8`) averaged `223.521 ms`: element apply `103.379 ms`, preconditioner `90.988 ms`, other Krylov work `29.155 ms`. Apply and preconditioner are both material. |
| Line additive Schwarz | Reduces iteration strongly (`xxl` floor-lines about 1249 -> 733) but current dense local solves make wall time slower. Keep as preconditioner research baseline, not production candidate yet. |
| Full-global apply oracle | Element-by-element full `K*x` now reproduces prescribed RHS, equivalent-load reactions, prescribed-value reuse, and tower prescribed reactions; 4/4 smoke cases pass. |
| Mechanism/fallback oracle | F3 under-constrained, F7 release-condensation, isolated-node removal, and a cantilever control all pass; HP lane must preserve this guard before PCG. |

Representative artifacts:

- `Research/out/hp_tuning/codex_hpfem_integrated_xxl/results.jsonl`
- `Research/out/hp_tuning/codex_hpfem_integrated_xxl/summary.csv`
- `Research/out/hp_tuning/codex_hpfem_integrated_xxl/summary.md`
- `Research/out/hp_tuning/codex_hpfem_factor_bypass_xxl/results.jsonl`
- `Research/out/hp_tuning/codex_hpfem_factor_bypass_xxl/summary.md`
- `Research/out/hp_tuning/codex_hpfem_touched_threads_xxl/summary.md`
- `Research/out/hp_tuning/codex_hpfem_batch_autotune_xxl/summary.md`
- `Research/out/hp_tuning/codex_hpfem_seed_basis_xxl/summary.md`
- `Research/out/hp_tuning/codex_hpfem_seed_basis_long_xxl/summary.md`

## Recommended Architecture

### Lane A: Exact Oracle/Fallback

- Existing `assembleAndFactor()` + `solveLoad()` remains authoritative.
- Use it for:
  - small systems;
  - singular/mechanism confirmation;
  - prescribed displacement RHS columns;
  - reactions and exact audit paths;
  - any case where iterative residual or pivot guards are inconclusive.

### Lane B: HP Iterative Solver

Core components:

- Element-by-element `LinearOperator` storing fixed-size member blocks.
- Persistent thread-local apply: partition element blocks once, reuse local `y`
  buffers, reduce in deterministic thread order.
- Preconditioner stack:
  - diagonal or 6x6 block Jacobi as baseline;
  - compact coarse correction constrained by dense coarse size;
  - line/column additive Schwarz for frame topology, only after reducing local
    solve cost;
  - next: geometric/algebraic multilevel once the coarse space is proven.
- Multi-RHS recycling:
  - maintain A-orthogonal basis of recent converged solutions;
  - project new RHS through Galerkin initial guess;
  - always finish with PCG residual and LDLT/assembled oracle validation in research.

The current best production-shaped candidate is now measured as a combined
research loop, not just separate microbenchmarks:

`matrix-free apply + persistent 4-thread pool + compact block/coarse preconditioner + RHS recycling for repeated solves`.

The first-solve factor-bypass target is now met for the `xxl` research tower:
avoid global LDLT, build a factor-free HP operator/preconditioner/thread pool,
and use LDLT only as the research oracle. The repeated-RHS target is not met
for arbitrary short RHS sequences because LDLT amortizes its expensive factor
very well; the HP lane must cut PCG iterations and preconditioner cost further
before it can dominate those cases.

For low-dimensional parametric gameplay loads, a seeded response-basis path is
now the stronger architecture: solve a small set of load modes through PCG,
freeze the A-orthogonal response basis, then project each frame's RHS and finish
with a residual gate. This preserves precision on the tested load family and
beats reused LDLT once enough frames amortize the seed solves.

The next optimization target is the preconditioner apply path and iteration
count. If the preconditioner remains mostly serial or too weak, it limits the
end-to-end benefit of parallel element apply.

## Core Risks

- Bit-identical gates: do not replace default `solve()` until all exact tests are
  protected by fallback.
- Prescribed DOFs and reactions require full global `K*u`, not only `Kff*x`.
  A research oracle now covers this, but production adoption still needs it wired
  into engine-facing gates.
- Iterative residual alone is not a mechanism detector; LDLT pivot/fallback must
  remain the guard. A research mechanism oracle now covers representative global
  pivot and element-condensation singular cases.
- Coarse preconditioners reduce iterations but can become cache-hostile quickly.
- First-solve factor bypass is solved on `xxl`; arbitrary short repeated RHS
  batches remain hard. Seeded response bases beat reused LDLT only when the
  load space is known and low-dimensional.
- Seeded load-basis acceleration must be guarded by projection residual. If a
  frame load leaves the seeded subspace, fall back to ordinary PCG or LDLT.
- Thread-level reduction changes floating-point summation order; research path
  currently validates at tolerance level, not bit identity.
- Additive Schwarz must handle local semidefinite subdomains by pivot guard and
  skip/fallback.
- Line-Schwarz local dense solve is currently compute/cache expensive; the next
  version needs cheaper subspace solves, batching, or a lower-rank line model.

## Benchmark Metrics

Required per run:

- correctness: `applyRel`, `pcgVsLdlt`, true residual, full-apply prescribed/reaction oracle;
- solver: PCG iterations, PCG ms, setup ms, apply ms, preconditioner apply ms;
- memory: element block bytes, coarse matrix bytes, local thread buffers, peak RSS if available;
- hardware: thread count, pool vs spawn apply, bandwidth proxy, cache-sensitive coarse size;
- tuning: preconditioner kind, coarse bins, basis size, RHS count, thread count.

## Next Research Route

1. Add a projection-residual gate for seeded load-basis solves, then route
   out-of-subspace loads to ordinary PCG/fallback.
2. Make non-seeded `factor_bypass_batch_speedup > 1.0` for short arbitrary RHS
   by reducing both iteration count and preconditioner cost.
3. Use the PCG timing breakdown to reduce the dominant element-apply term and
   the secondary preconditioner term together; optimizing only one is unlikely
   to push repeated batches past LDLT reuse.
4. Replace dense coarse inverse apply with a cheaper hierarchy or batched small
   coarse solve, then retest `1x1/2x2/3x3` coarse spaces under cache counters.
5. Optimize line-Schwarz local solve cost before considering it for production.
6. Run repeated/pinned benchmarks, then agentic autotuning over:
   `precond`, `coarseBins`, `basisMax`, `threads`, `rhsBatch`, and apply mode.

## 2026-06-13 Session 2: parallel precond + banded coarse + deflation plan

Hardware confirmed: AMD Ryzen 9 8940HX, 16 physical cores / 32 threads, L2 16MB
(1MB/core), L3 64MB across two CCDs (~32MB each). Cross-CCD L3 penalty explains
why 16 threads is optimal and 32 regresses. The 2.88MB dense coarse matvec blows
the 1MB/core L2; the banded coarse factor (~115KB-220KB, streamed) fits L2.

### Single-run noise is large — interleave or repeat

The previously-reported non-seeded 8-RHS `0.932x` was a lucky quiet run. Repeated
runs on the same exe gave `0.745-0.797x` when the machine was busy and `1.013x`
when quiet. **The factor_bypass_batch_speedup at 8 RHS is dominated by the one-time
LDLT factor (~1860ms) and is too noisy/factor-dominated to resolve a per-solve
optimization.** Always compare baseline vs candidate back-to-back in the SAME
session, and prefer `combinedMsAvg` (per-solve) over the batch ratio for tuning.

### Parallel block6 precond + banded block-tridiagonal coarse (DONE, validated)

`exp_parallel_pcg.cpp` gained two opt-in flags (default = old behavior):
- `--parallelPrecond`: the persistent ThreadApplyPool now also runs the block6
  Jacobi map z=inv6(.)r as a disjoint-output parallel job (no reduction). The
  serial PCG baseline still uses serial precond for an honest comparison.
- `--coarseSolve banded`: replaces the dense `coarseInv * rc` matvec with a
  precomputed block-Thomas LDL^T factor of the floor-aggregated coarse matrix.
  The coarse matrix IS exactly block-tridiagonal in floor order
  (`bandedOffBandRel=0.000e+00`; frame elements only couple same/adjacent floors).
  A one-time self-check (`bandedResidual=1.5e-11`) and a tridiagonality guard
  gate the path; either failing falls back to the dense inverse.

Measured on `xxl` (nf=18720, coarseDofs=576, 16T, 2x2, basis8, 8 RHS), quiet
same-session back-to-back:
- baseline serial-precond + dense coarse: combinedMsAvg ~239ms (apply ~111,
  precond ~97), iters 739.75.
- parallel-precond + banded: combinedMsAvg ~212ms (apply ~106, **precond ~73**),
  iters 758.38. `maxCombinedVsLdlt=2.1e-12` (correctness preserved).

So precond apply dropped ~24% and combined per-solve dropped ~11%. The banded
factor is a marginally weaker preconditioner (solves Kc to 1.5e-11 vs dense
~1e-14), costing ~2.5% more iterations — a benign, honest trade (final answer
still 2e-12 vs LDLT). The two optimizations COMPOUND with deflation: each cheaper
iteration multiplies the iteration reduction below.

### Next big lever: deflation / recycled PCG for ARBITRARY RHS (in progress)

Seeded load-basis only helps when loads lie in a known subspace. The RHS-INDEPENDENT
win is spectral deflation: the slowest-converging error modes are the smallest
eigenvectors of M^{-1}K and do not depend on b. Plan (from a literature scout):
1. Harvest k~30 smallest Ritz vectors W via **eigCG** (Stathopoulos & Orginos 2010)
   during the first solve(s): the CG alpha/beta/sigma sequence IS the Lanczos
   tridiagonal T; Rayleigh-Ritz on T gives approximate eigenvectors, thick-restart
   bounds the window to m~2k.
2. For every subsequent RHS run **Type-1 deflated PCG** (Saad/Yeung/Erhel/Guyomarc'h
   2000): precompute `E=W^T K W` and `KW=K W`; corrected init `x0 = W E^{-1} W^T b`;
   project the residual `r -= W E^{-1} (W^T r)` every few iterations to fight
   orthogonality drift.
   Predicted: iters 740 -> ~230 (kappa_eff = lambda_n/lambda_{k+1}), net ~2.5x after
   the O(kn)/iter overhead, RHS-independent. Pure Eigen, no new dependency.
3. Gates: Ritz residual `||K w - lambda w||/lambda < 1e-6` before trusting W;
   fall back to plain PCG if deflated iters exceed 1.3x plain; rebuild W if K changes.
References: eigCG https://www.cs.wm.edu/~andreas/publications/eigCG.pdf ;
Deflated-CG https://inria.hal.science/inria-00523686/document .

### Session 2 RESULTS (measured)

New opt-in flags on `exp_parallel_pcg.cpp` (all default off = old behavior, so existing
artifacts reproduce): `--parallelPrecond`, `--coarseSolve banded`, `--deflation k
[--deflationWindow m]`, `--symApply`, `--towerNx/Ny/Stories`. A reusable repeated /
interleaved benchmark lives in `bench_pinned.sh` (single runs are too noisy to compare).

CONFIRMED WINS (parallel block6 precond + banded 2x2 coarse), pinned 7-rep interleaved,
`factor_bypass_batch_speedup` median [min..max]:
- non-seeded 8 RHS: `0.740 [0.657..1.017]` -> `0.826 [0.729..1.005]` (+12%, NEW>OLD 7/7).
- seeded 32 RHS:    `1.211 [1.114..1.334]` -> `1.374 [1.347..1.578]` (+13%, NEW>OLD 7/7;
  the NEW minimum 1.347 exceeds the OLD maximum 1.334 — distributions do not overlap).
- seeded 200 RHS (single, same session): `2.233 -> 2.443`.
All keep `maxCombinedVsLdlt ~2e-12`. The precond apply dropped ~24% (97->73ms); these
gains flow into every PCG-based path, including the game-relevant seeded lane.

HARDWARE / SCALING (the strongest argument for the HP lane), single runs, 16T:
| nf | elem data | LDLT factor | HP setup | factorBypassSetup | apply BW |
| 18720 | 11.8MB | 1709ms | 98.6ms | 17.3x | 87 GB/s |
| 34560 | 21.6MB | 8782ms | 285ms  | 30.8x | 96 GB/s |
| 48384 | 30.1MB | 20683ms| 523ms  | 39.6x | 87 GB/s |
- The LDLT factor explodes superlinearly (2.58x DOF -> 12.1x factor time: the fill-in
  wall), while HP setup grows ~linearly, so the factor-bypass SETUP speedup GROWS with
  size (17x -> 40x). This is the core scaling case for HP at game-engine sizes.
- Element-apply bandwidth is a flat ~87-96 GB/s across all three sizes: the element data
  (<=30MB) is L3-RESIDENT (64MB L3), so the apply is L3-bandwidth/ILP bound, NOT DRAM
  bound, at these sizes. Apply parallelism saturates near 16 threads (cross-CCD wall).

HONEST NEGATIVE / NEUTRAL RESULTS (investigated, not net wins here — kept as gated
opt-in experiments, documented so they are not re-tried blindly):
- Spectral deflation (`--deflation`): correct (maxDeflVsLdlt 2e-12) but NOT a net win.
  Iteration reduction SATURATES at ~1.41x even with k=24..48 and window=600 -> the
  block6+coarse-preconditioned spectrum is spread, with no small isolated low-eigenspace
  to deflate. Worse, the dense per-iteration projection `z - W E^{-1}(KW)^T z` streams
  2*(n x k) and is memory-bound (~0.07ms/iter at k=16, more at larger k), roughly
  cancelling the iteration savings -> best `speedupMs ~0.98` (break-even). The literature's
  ~2.5x assumed isolated small eigenvalues AND an apply that is cheap relative to the dense
  projection; neither holds for this matrix-free, already-well-preconditioned solver.
- Finer banded coarse (3x3/4x4): cuts iterations (886->565) but the banded coarse solve
  cost grows ~fb^2 (floor block 24->96), so combinedMs has a MINIMUM at 2x2 (~207ms).
- Seed-basis recycling between load modes: neutral — the 5 load modes are too distinct for
  the A-orthogonal basis of earlier seeds to shortcut later seeds (seedIters ~unchanged).
- Symmetric apply (`--symApply`): correct (vsLdlt 5.9e-12) but neutral at these sizes —
  halving element storage (78 vs 144) does not speed the apply because the data is
  L3-resident (not DRAM-bound) and the symmetric loop's `ye[j]+=` cross-update hurts ILP,
  cancelling the memory saving. It DOES force exact symmetry (iters 758->742). Expected to
  help only once element data exceeds the 64MB L3 (>~77k elements / >~100k DOF).

FUNDAMENTAL LIMIT (why "non-seeded arbitrary many-RHS beats reused LDLT" is the wrong
goal at this scale): a reused sparse-LDLT back-substitution is ~11ms (10.5-12.4ms observed)
vs ~207ms for a full non-seeded PCG solve — a ~18x gap that no preconditioner/iteration
improvement closes. Direct
factorization amortizes a single factor beautifully over many cheap solves. The HP lane's
genuine niches are therefore: (a) FEW RHS / first solve [factor avoidance, 4-40x setup],
(b) LOW-DIMENSIONAL parametric load families [seeded, 1.4-2.4x batch], and (c) VERY LARGE
problems where the LDLT factor is prohibitive in time/memory [the scaling table above].
For a real game engine the loads ARE low-dimensional (gravity + a few contacts), so the
seeded lane — now ~14% faster — is the production-relevant path, not arbitrary RHS.

GAME-ENGINE HEADLINE — seeded HP asymptotic vs reused LDLT (the production-relevant
result). A real game loop is thousands of frames with LOW-DIMENSIONAL loads (gravity +
a few contacts). After a one-time seed setup (~1.3s, 5 load-mode solves), the seeded
basis gives a near-perfect initial guess (`maxCombinedInitialRel ~1.1e-10`) so each
frame converges in ~1 PCG iteration: `combinedMsAvg ~0.73ms` vs a reused sparse-LDLT
back-substitution `ldltSolveMsAvg ~10.9ms` — a ~15x PER-FRAME asymptote at `vsLdlt 2e-12`
(effectively exact). The realized batch speedup over reused LDLT climbs with frame count
as the seed setup amortizes:
| RHS (frames) | factor_bypass_batch_speedup |
| 200  | 2.37x |
| 1000 | 5.65x |
| 2000 | 8.04x |
| 4000 | 10.21x  (-> ~15x asymptote) |
This is the strongest case for the HP lane as a game-engine structural solver: for the
load pattern games actually have, it is ~10-15x faster than re-using a prefactored LDLT,
not the ~1x of arbitrary RHS. (`--noSerialBaseline` skips the per-RHS serial reference so
large-batch sweeps are tractable.)

SPARSE COARSE (`--coarseSolve sparse`, Eigen SimplicialLDLT on the sparse coarse matrix):
correct (vsLdlt 2e-12) and lets fine coarse spaces be built, but NOT a win. Finer coarse
does cut PCG iterations strongly (2x2:756 -> 4x4:568 -> 8x8:321) but the coarse SOLVE cost
grows with coarse size via fill-in (precond 76ms -> 182ms -> 946ms), and SimplicialLDLT is
slightly WORSE than the hand-rolled block-Thomas for this exactly-block-tridiagonal coarse
(4x4: sparse 182ms vs banded 127ms). combinedMs stays minimized at 2x2-banded (~207ms).
The lesson: a finer coarse wants to be solved by RECURSION (multigrid V-cycle), not a
direct factor — the direct coarse solve cost outruns the iteration savings.

UNIFYING SYNTHESIS — per-iteration economics of a matrix-free L3-resident solver. The
element apply is cheap (L3-resident, ~87 GB/s, ~0.14ms/iter at 16T). Therefore ANY method
that adds per-iteration work of order one matvec (spectral deflation's dense W projection,
multiplicative/finer direct coarse, Chebyshev smoothing) must cut iterations by MORE than
its added cost to pay — and this already-well-preconditioned spread spectrum only yields
~1.4-2x iteration reductions, which the added cost cancels. The robust wins therefore come
from (1) making each existing iteration cheaper without adding matvecs (parallel block6
precond + banded coarse — done, ~24% precond), and (2) exploiting structure across solves
(seeded recycling for low-dim loads -> ~15x/frame; factor-bypass at scale -> 17-40x setup).
This is why "beat reused LDLT for arbitrary many-RHS" is the wrong target at this size and
the seeded/scale lanes are the real product.

ADVERSARIAL CORRECTNESS REVIEW (independent agent, read-only): no CRITICAL bugs. The
ThreadApplyPool condition-variable protocol (spurious/lost-wakeup safe, disjoint precond
writes race-free), the block-Thomas LDL^T (Schur complement + forward/diag/back), the
Type-1 DCG projection and `x0=W E^{-1}W^T b`, the K-inner-product Rayleigh-Ritz, and the
LDLT-referenced `ok` gates were all verified correct. The deflation negative result is
honestly supported (maxDeflVsLdlt 2e-12 proves the deflated solution equals LDLT). Minor
notes addressed: the apply/banded self-check probes are now mixed-frequency + sign-
alternating instead of a single smooth sine.

## 2026-06-14 Session 3: real game-engine point-load validation

The headline seeded result above was measured on a SYNTHETIC load family: makeRhsFamily's
every frame is a smooth sinusoidal combination of the SAME 5 building modes
(base/windX/windY/liveZ/torsion) that the seed step solves, so every frame is in-subspace
BY CONSTRUCTION. That is circular -- it proves the basis reproduces the basis. This session
validates the seeded lane on a REAL game-engine load family instead -- constant gravity plus
a few MOVING sparse contact point loads -- and checks whether the projection-residual gate
detects loads that leave the seeded subspace.

### Method (makeGameLoadFamily, opt-in `--gameLoad`, default off)

* loadModes (seed) = { gravity } U { the -Z unit response of each in-candidate contact
  node }. sequence (per frame) = gravity + `active` moving contact point loads. Contact
  nodes/strengths come from a stateless Knuth hash (reproducible; no rand/time).
* Math: a load that only pushes on node set P has its solution in span{K^-1 e_p : p in P}.
  So a frame is in-subspace iff all its contacts sit on seeded nodes. `--gameOutFraction p`
  routes a fraction p of frames onto NON-seeded nodes (out-candidates, disjoint from the
  seeded set) to exercise the gate. Three regimes: R1 p=0, R2 p=0.5, R3 p=1.
* basisMax auto-raised to seedModes (=candidates+1) so RecycleBasis's FIFO eviction cannot
  silently drop early seeds (basisAccepted==seedModes==13, basisRejected=0 verified).

### Results (xxl nf=18720, 16T, parallel precond + banded coarse; small for correctness)

R1 (contacts all seeded): each frame converges in **0 PCG iter** (initialRel 9.66e-11 <
tol), solution vs LDLT 2.1e-12. A dense second-path cross-check (`--gameVerify`) confirms
the basis is K-orthonormal (||VtKV - I|| = 3.4e-15) and that the cheap Euclidean
initialGuess equals the Galerkin-optimal V(VtKV)^-1 Vt b (crossRel 2.4e-15).

R3 (contacts all non-seeded): initialRel 4.1e-2 -- EIGHT orders of magnitude above R1's
9.7e-11 -- and iterations return to baseline (~335 vs serial 349). The solution is still
correct (vs LDLT 4.1e-12): a bad initial guess costs only speed, never accuracy. This is
the gate's justification -- a 1e-6 threshold on initialRel cleanly separates in/out.
(initialRel is 4e-2, not 1.0, because gravity dominates ||b|| and is recovered exactly; the
residual is the un-projected contact part -- a real, honest physical dilution.)

R2 (half out): in-frames 9.6e-11 / 0 iter, out-frames 3.8e-2 / ~335 iter -- a clean split.

### HONEST per-frame cost: projection is NOT free (adversarial-review correction)

pcg's combinedMs does NOT time the Galerkin projection that basis.initialGuess does every
frame (it lives in projectMs). The meaningful per-frame cost is combinedMsAvg + projectMsAvg.
With that correction (new outputs `perFrameMsWithProj` / `perFrameSpeedupWithProj`):

| family | combinedMsAvg | projectMsAvg | per-frame (with proj) | per-frame speedup vs reused LDLT |
| --- | --- | --- | --- | --- |
| game (k=17) | 0.339 | 0.212 | 0.550 | **19.5x** |
| synthetic (k=5) | 0.473 | 0.071 | 0.544 | 19.7x |

So the real-load lane is ~19x per frame -- the SAME as synthetic, not faster. The game
lane's 0-iter advantage (vs synthetic's ~0.44 iter) is cancelled by its larger k=17
projection (more seeded contacts). The naive combinedMsAvg-only ratio (~31x for game)
OVERSTATES the win; the same projection cost applies to the synthetic headline above too.

### Seed-dimension cost and the batch crossover

Real loads need one seed per contact candidate, so seedMs is ~4500 ms (16 contacts) vs
~1300 ms (5 synthetic modes). The one-time seed amortizes more slowly: factor_bypass_batch_
speedup at 4000 frames is ~7.4x (game) vs ~10-13x (synthetic, single-run noisy). Per-frame
they are equal, so the game lane overtakes only after the extra seed cost amortizes --
crossover ~24000 frames (~3247 ms extra seed / ~0.13 ms per-frame saving). For a real
multi-thousand-frame session both sit at the ~19x per-frame asymptote; the batch number
just trails until the larger seed pays off.

### Validation & honesty

* Third-party numpy check (`validate_subspace.py`, no FrameCore/Eigen): independently
  re-derives R1 in-subspace residual ~1e-14, R3 out-of-subspace ~1.1, and VtKV=I (so the
  Euclidean projection equals the Galerkin one). Methodology confirmed in a separate impl.
* Regression: with `--gameLoad` off the synthetic path is unchanged (seedCount=5, vsLdlt
  2e-12). All new flags default off.
* Adversarial read-only review raised three points, all addressed: (1) projection cost now
  reported (above); (2) out/in-candidate overlap guarded with a loud warning when
  nFree < 2*cand; (3) `--gameAddHorizontal` deliberately leaves the -Z seed span (a
  controlled perturbation) -- documented, and the 0-iter claim is the pure -Z default only.

### Takeaway for production (task A)

Real low-dimensional game loads DO ride the seeded lane (~19x/frame, gate reliable), but
the win EQUALS the synthetic estimate rather than exceeding it, and the seeded subspace must
cover the actual contact set. The production fallback should key on the projection-residual
gate: in vs out initialRel differ by 8-9 orders, so a 1e-6 threshold separates them with
wide margin; out-of-subspace frames fall back to ordinary PCG (still correct, baseline cost)
or LDLT.
