#!/usr/bin/env python
"""
Focused audit for the grillage plate idealization and segmented curved-member
approximation. Models are generated independently in Python and solved through
frame_cli.exe, complementing the C++ standalone F12/F17 tests.
"""
from __future__ import annotations

import math

import independent_precision_audit as audit


def strip_section(width: float, t: float, nu: float) -> dict:
    i = t**3 / (12.0 * (1.0 - nu * nu))
    c = t**3 / (6.0 * (1.0 - nu))
    return dict(A=width * t, Iy=width * i, Iz=width * i, J=width * c,
                cy=t / 2.0, cz=t / 2.0, Asy=0.0, Asz=0.0)


def grid_node(i: int, j: int, nx: int) -> int:
    return j * (nx + 1) + i


def grillage_model(n: int) -> tuple[dict, float]:
    E, nu = 30000.0, 0.3
    G = E / (2.0 * (1.0 + nu))
    a = b = 4000.0
    t = 250.0
    q = 0.025
    dx, dy = a / n, b / n
    sections = [
        strip_section(dy, t, nu),
        strip_section(dy / 2.0, t, nu),
        strip_section(dx, t, nu),
        strip_section(dx / 2.0, t, nu),
    ]
    nodes = []
    for j in range(n + 1):
        for i in range(n + 1):
            on_edge = i == 0 or i == n or j == 0 or j == n
            nodes.append(dict(id=grid_node(i, j, n), x=i * dx, y=j * dy, z=0.0,
                              fix=[1, 1, 1 if on_edge else 0, 0, 0, 1]))
    members = []
    mid = 0
    for j in range(n + 1):
        sec = 1 if j in (0, n) else 0
        for i in range(n):
            members.append(dict(id=mid, i=grid_node(i, j, n), j=grid_node(i + 1, j, n),
                                mat=0, sec=sec, refvec=(0.0, 0.0, 1.0)))
            mid += 1
    for i in range(n + 1):
        sec = 3 if i in (0, n) else 2
        for j in range(n):
            members.append(dict(id=mid, i=grid_node(i, j, n), j=grid_node(i, j + 1, n),
                                mat=0, sec=sec, refvec=(0.0, 0.0, 1.0)))
            mid += 1
    loads = []
    for j in range(n + 1):
        wy = dy / 2.0 if j in (0, n) else dy
        for i in range(n + 1):
            wx = dx / 2.0 if i in (0, n) else dx
            loads.append(dict(node=grid_node(i, j, n), comp=[0.0, 0.0, -q * wx * wy, 0.0, 0.0, 0.0]))
    D = E * t**3 / (12.0 * (1.0 - nu * nu))
    exact = 0.00406 * q * a**4 / D
    return dict(materials=[dict(E=E, G=G, rho=2400.0)], sections=sections, nodes=nodes,
                members=members, nloads=loads, udls=[]), exact


def arch_model(nseg: int) -> tuple[dict, float]:
    E, G = 210000.0, 80769.0
    R = 2000.0
    P = 1000.0
    sec = audit.rect_section(100.0, 100.0)
    nodes = []
    for k in range(nseg + 1):
        th = (math.pi / 2.0) * (k / nseg)
        nodes.append(dict(id=k, x=R * math.cos(th), y=R * math.sin(th), z=0.0,
                          fix=[1, 1, 1, 1, 1, 1] if k == 0 else [0, 0, 0, 0, 0, 0]))
    members = [dict(id=k, i=k, j=k + 1, mat=0, sec=0, refvec=(0.0, 0.0, 1.0)) for k in range(nseg)]
    exact = (math.pi / 4.0) * P * R**3 / (E * sec["Iz"])
    return dict(materials=[dict(E=E, G=G, rho=7850.0)], sections=[sec], nodes=nodes,
                members=members, nloads=[dict(node=nseg, comp=[0.0, P, 0.0, 0.0, 0.0, 0.0])],
                udls=[]), exact


def main() -> int:
    print("=" * 86)
    print("Grillage and segmented-curve audit")
    print("=" * 86)
    failures = 0

    print("\n[grillage square plate, simply supported, uniform pressure]")
    prev = None
    tail = []
    for n in [4, 6, 8, 10, 12, 16, 20]:
        model, exact = grillage_model(n)
        out = audit.run_cli(model)
        center = abs(out["disp"][grid_node(n // 2, n // 2, n)][audit.Uz])
        err = abs(center - exact) / exact
        change = 0.0 if prev is None else abs(center - prev) / max(abs(prev), 1e-30)
        print(f"  N={n:2d}  w_center={center:.8g}  err={100*err:6.3f}%  change={100*change:6.3f}%")
        if err > 0.05:
            failures += 1
        if n >= 12:
            tail.append(center)
        prev = center
    if len(tail) >= 2 and abs(tail[-1] - tail[-2]) / abs(tail[-2]) > 0.02:
        failures += 1

    print("\n[segmented quarter-circle arch cantilever]")
    prev_err = math.inf
    best_err = math.inf
    final_err = math.inf
    for n in [4, 8, 16, 32, 64, 128]:
        model, exact = arch_model(n)
        out = audit.run_cli(model)
        tip = abs(out["disp"][n][audit.Uy])
        err = abs(tip - exact) / exact
        ratio = prev_err / err if math.isfinite(prev_err) and err > 0 else float("inf")
        print(f"  N={n:3d}  tip={tip:.8g}  err={100*err:8.5f}%  errRatio={ratio:7.3f}")
        if n <= 64 and n >= 8 and err > prev_err:
            failures += 1
        best_err = min(best_err, err)
        final_err = err
        prev_err = err
    if final_err > 5.0e-4:
        failures += 1
    if final_err > best_err:
        print(f"  note: error bottoms before the finest mesh, then approaches a small plateau "
              f"({100*final_err:.5f}% vs bending-only oracle).")

    print("\n" + ("PASS" if failures == 0 else f"FAIL failures={failures}"))
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
