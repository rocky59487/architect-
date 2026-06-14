#!/usr/bin/env python3
# validate_subspace.py -- research-only THIRD-PARTY adversarial check (numpy).
#
# Independently re-derives, in pure numpy (no FrameCore, no Eigen), the linear-algebra
# claim that the HP-FEM seeded-recycling lane relies on, so the C++ result in
# exp_parallel_pcg.cpp is cross-checked by a separate implementation:
#
#   The seeded "load basis" is the K-orthonormalized span of the unit point-load
#   responses {K^-1 e_p : p in P}. The cheap Euclidean projection x0 = V V^T b that
#   exp_parallel_pcg's RecycleBasis::initialGuess uses equals the Galerkin-optimal
#   x0 = V (V^T K V)^-1 V^T b ONLY when V is K-orthonormal (V^T K V = I). Therefore:
#     * R1 (in-subspace): any load b = sum a_p e_p with p in P projects to the EXACT
#       solution  -> initial relative residual ~ machine epsilon  -> ~1 PCG iter.
#     * R3 (out-of-subspace): a load e_q with q not in P leaves the span -> the same
#       projection leaves an O(1) residual -> the projection-residual gate must fire.
#
# This is the *methodology* check (the pure linear algebra). The C++ experiment then
# checks the same identities on the real assembled FrameCore K via --gameVerify.
#
# Run:  python validate_subspace.py        (exit 0 = all assertions passed)

import numpy as np


def k_orthonormal_basis(K, responses):
    """K-inner-product modified Gram-Schmidt, byte-faithful to RecycleBasis::add
    (exp_parallel_pcg.cpp:954). Returns V whose columns satisfy V^T K V = I, plus the
    accepted-count so we can detect near-linear-dependence rejection (the C++
    basisAccepted < seedCount case)."""
    V = []           # list of K-orthonormal columns v_i
    AV = []          # cached A v_i = K v_i (the C++ av[i])
    accepted = 0
    for x in responses:
        q = x.astype(float).copy()
        Aq = K @ q
        for vi, avi in zip(V, AV):
            c = vi @ Aq          # v_i^T (K q): the C++ uses v[i].dot(Aq)
            q = q - c * vi
            Aq = Aq - c * avi
        normA2 = q @ Aq          # q^T K q
        scale_ref = max(1e-300, np.linalg.norm(x))
        if not (normA2 > 0) or not np.isfinite(normA2) or np.linalg.norm(q) <= 1e-12 * scale_ref:
            continue             # rejected (matches the C++ rejection guard)
        inv = 1.0 / np.sqrt(normA2)
        V.append(q * inv)
        AV.append(Aq * inv)
        accepted += 1
    return (np.array(V).T if V else np.zeros((K.shape[0], 0))), accepted


def initial_guess(V, b):
    """Euclidean projection x0 = sum_i v_i (v_i . b) = V V^T b -- exactly what
    RecycleBasis::initialGuess (exp_parallel_pcg.cpp:945) computes."""
    return V @ (V.T @ b)


def rel_resid(K, x0, b):
    return np.linalg.norm(b - K @ x0) / max(1e-300, np.linalg.norm(b))


def main():
    rng = np.random.default_rng(0)
    n = 60

    # Symmetric positive-definite K with a non-trivial condition number, standing in
    # for the reduced free-free stiffness K_ff. The shift keeps it well clear of
    # singular while leaving a spread spectrum (the realistic, hard-to-deflate case).
    A = rng.standard_normal((n, n))
    K = A.T @ A + 0.5 * np.eye(n)
    assert np.allclose(K, K.T), "K must be symmetric"
    assert np.min(np.linalg.eigvalsh(K)) > 0, "K must be SPD"

    # Seedable contact set P (in-candidates): a handful of distinct DOFs carry a unit
    # point load. e_p is the canonical basis vector; its response is K^-1 e_p.
    in_dofs = [3, 7, 11, 19, 23, 31, 37, 41]
    out_dofs = [5, 17, 50]                       # NOT seeded -> must leave the span
    assert not (set(in_dofs) & set(out_dofs)), "in/out DOF sets must be disjoint"

    responses = [np.linalg.solve(K, np.eye(n)[:, p]) for p in in_dofs]
    V, accepted = k_orthonormal_basis(K, responses)

    # --- Claim 0: the basis really is K-orthonormal (the equivalence precondition) ---
    off = np.linalg.norm(V.T @ K @ V - np.eye(V.shape[1]))
    assert accepted == len(in_dofs), f"basis lost vectors: accepted={accepted} < {len(in_dofs)}"
    assert off < 1e-10, f"V^T K V != I : ||.|| = {off:.3e}"

    # --- Claim 1: Euclidean projection equals the Galerkin-optimal projection -------
    # x0_galerkin = V (V^T K V)^-1 V^T b ; with V^T K V = I this reduces to V V^T b.
    b_probe = rng.standard_normal(n)
    VtKV = V.T @ K @ V
    x0_galerkin = V @ np.linalg.solve(VtKV, V.T @ b_probe)
    x0_fast = initial_guess(V, b_probe)
    cross = np.linalg.norm(x0_galerkin - x0_fast) / max(1e-300, np.linalg.norm(x0_galerkin))
    assert cross < 1e-10, f"Euclidean != Galerkin projection: rel = {cross:.3e}"

    # --- Claim 2 (R1, in-subspace): any combination of seeded point loads is exact --
    worst_in = 0.0
    for _ in range(200):
        alpha = rng.standard_normal(len(in_dofs))
        b = np.zeros(n)
        for a, p in zip(alpha, in_dofs):
            b[p] += a
        worst_in = max(worst_in, rel_resid(K, initial_guess(V, b), b))
    assert worst_in < 1e-10, f"R1 in-subspace residual too large: {worst_in:.3e}"

    # --- Claim 3 (R3, out-of-subspace): a non-seeded point load leaves an O(1) gap --
    best_out = np.inf
    for q in out_dofs:
        b = np.eye(n)[:, q]
        best_out = min(best_out, rel_resid(K, initial_guess(V, b), b))
    assert best_out > 1e-2, f"R3 out-of-subspace residual unexpectedly small: {best_out:.3e}"

    # --- Claim 4 (mixed frame): in-part stays projected, out-part drives the residual.
    # A frame that adds even one out-of-subspace contact must show a large residual:
    # this is exactly what the projection-residual gate keys on.
    b_mixed = np.zeros(n)
    for a, p in zip(rng.standard_normal(len(in_dofs)), in_dofs):
        b_mixed[p] += a
    b_mixed[out_dofs[0]] += 1.0
    mixed = rel_resid(K, initial_guess(V, b_mixed), b_mixed)
    assert mixed > 1e-3, f"mixed-frame residual should be large: {mixed:.3e}"

    print(f"[validate_subspace] n={n} seeded={len(in_dofs)} basisAccepted={accepted} "
          f"offDiagVtKV={off:.2e} crossEuclidVsGalerkin={cross:.2e}")
    print(f"[validate_subspace] R1_inSubspaceMaxRel={worst_in:.2e} "
          f"R3_outSubspaceMinRel={best_out:.3f} mixedFrameRel={mixed:.3f}")
    print("[validate_subspace] all assertions passed")


if __name__ == "__main__":
    main()
