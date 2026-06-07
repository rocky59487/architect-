#!/usr/bin/env python3
"""
Independent MITC4 shell audit.

This script intentionally builds shell models in Python and drives the standalone
frame_cli.exe as a black box. It does not call the C++ test fixtures. OpenSeesPy
is used only for additional cross-validation when it is installed.
"""

from __future__ import annotations

import copy
import math
import os
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CLI = ROOT / "Plugins" / "FrameSolver" / "Standalone" / "frame_cli.exe"
BUILD_CLI = ROOT / "Plugins" / "FrameSolver" / "Standalone" / "build_cli.bat"

UX, UY, UZ, RX, RY, RZ = range(6)
E0 = 30000.0
NU0 = 0.3
G0 = E0 / (2.0 * (1.0 + NU0))


class Audit:
    def __init__(self) -> None:
        self.failures = 0
        self.warnings = 0

    def check(self, name: str, ok: bool, detail: str = "") -> None:
        print(f"{'PASS' if ok else 'FAIL'}  {name:<62} {detail}")
        if not ok:
            self.failures += 1

    def warn(self, name: str, detail: str = "") -> None:
        print(f"WARN  {name:<62} {detail}")
        self.warnings += 1

    def report(self, name: str, detail: str = "") -> None:
        print(f"INFO  {name:<62} {detail}")


def ensure_cli() -> None:
    if CLI.exists():
        return
    rc = subprocess.run(["cmd", "/c", str(BUILD_CLI)], cwd=ROOT).returncode
    if rc != 0 or not CLI.exists():
        raise RuntimeError("could not build frame_cli.exe")


def mat(E: float = E0, nu: float = NU0) -> dict:
    return {"E": E, "nu": nu, "G": E / (2.0 * (1.0 + nu))}


def vsub(a, b):
    return [a[i] - b[i] for i in range(3)]


def cross(a, b):
    return [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]


def norm(v) -> float:
    return math.sqrt(sum(x * x for x in v))


def quad_area(nodes_by_id: dict[int, dict], ids: list[int]) -> float:
    p = [[nodes_by_id[i]["x"], nodes_by_id[i]["y"], nodes_by_id[i]["z"]] for i in ids]
    return 0.5 * norm(cross(vsub(p[2], p[0]), vsub(p[3], p[1])))


def run_cli(model: dict) -> dict:
    ensure_cli()
    lines: list[str] = []
    for m in model["mats"]:
        lines.append(f"SMAT {m['E']} {m['nu']} {m['G']}")
    for n in model["nodes"]:
        f = n["fix"]
        lines.append(
            "NODE {id} {x:.17g} {y:.17g} {z:.17g} {f0} {f1} {f2} {f3} {f4} {f5}".format(
                id=n["id"], x=n["x"], y=n["y"], z=n["z"],
                f0=f[0], f1=f[1], f2=f[2], f3=f[3], f4=f[4], f5=f[5]
            )
        )
    for s in model["shells"]:
        n = s["n"]
        lines.append(f"SHELL {s['id']} {n[0]} {n[1]} {n[2]} {n[3]} {s['mat']} {s['t']}")
    for l in model.get("nloads", []):
        c = l["comp"]
        lines.append(f"NLOAD {l['node']} {c[0]} {c[1]} {c[2]} {c[3]} {c[4]} {c[5]}")
    for p in model.get("spress", []):
        lines.append(f"SPRESS {p['shell']} {p['p']}")
    lines.append("END")
    p = subprocess.run(
        [str(CLI)],
        input="\n".join(lines) + "\n",
        capture_output=True,
        text=True,
        cwd=ROOT,
    )
    if p.returncode != 0:
        raise RuntimeError(p.stderr or "frame_cli failed")

    out = {"singular": None, "disp": {}, "rf": {}, "sf": {}, "stdout": p.stdout}
    for ln in p.stdout.splitlines():
        t = ln.split()
        if not t:
            continue
        if t[0] == "SINGULAR":
            out["singular"] = int(t[1])
        elif t[0] == "DISP":
            out["disp"][int(t[1])] = [float(x) for x in t[2:8]]
        elif t[0] == "RF":
            out["rf"][int(t[1])] = [float(x) for x in t[2:8]]
        elif t[0] == "SF":
            out["sf"][int(t[1])] = [float(x) for x in t[2:10]]
    return out


def plate_D(E: float, nu: float, t: float) -> float:
    return E * t ** 3 / (12.0 * (1.0 - nu * nu))


def square_plate(
    n: int,
    support: str,
    a: float = 1000.0,
    t: float = 10.0,
    q: float = 0.01,
    reverse: bool = False,
    pressure: bool = True,
    lumped: bool = False,
    perturb: float = 0.0,
    id_base: int = 1000,
    shell_base: int = 5000,
) -> tuple[dict, int]:
    h = a / n
    nodes: list[dict] = []
    shells: list[dict] = []
    spress: list[dict] = []
    nload_by_node: dict[int, float] = {}

    def nid(i: int, j: int) -> int:
        return id_base + j * (n + 1) + i

    for j in range(n + 1):
        for i in range(n + 1):
            x = i * h
            y = j * h
            if perturb and 0 < i < n and 0 < j < n:
                x += perturb * h * math.sin(math.pi * i / n) * math.sin(2.0 * math.pi * j / n)
                y += perturb * h * math.sin(2.0 * math.pi * i / n) * math.sin(math.pi * j / n)
            edge = i == 0 or i == n or j == 0 or j == n
            fix = [0, 0, 0, 0, 0, 0]
            if support == "simple":
                fix[UX] = fix[UY] = fix[RZ] = 1
                if edge:
                    fix[UZ] = 1
            elif support == "clamped":
                if edge:
                    fix = [1, 1, 1, 1, 1, 1]
            else:
                raise ValueError(support)
            nodes.append({"id": nid(i, j), "x": x, "y": y, "z": 0.0, "fix": fix})

    nodes_by_id = {nd["id"]: nd for nd in nodes}
    sid = 0
    for j in range(n):
        for i in range(n):
            ccw = [nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)]
            conn = [ccw[0], ccw[3], ccw[2], ccw[1]] if reverse else ccw
            shell_id = shell_base + sid
            shells.append({"id": shell_id, "n": conn, "mat": 0, "t": t})
            if pressure:
                # Base orientation normal is +Z. Reversed orientation normal is -Z.
                # Choose p so both cases load global -Z.
                spress.append({"shell": shell_id, "p": q if reverse else -q})
            if lumped:
                area = quad_area(nodes_by_id, ccw)
                for node_id in ccw:
                    nload_by_node[node_id] = nload_by_node.get(node_id, 0.0) - q * area * 0.25
            sid += 1

    nloads = [
        {"node": node_id, "comp": [0.0, 0.0, fz, 0.0, 0.0, 0.0]}
        for node_id, fz in sorted(nload_by_node.items())
    ]
    model = {"mats": [mat()], "nodes": nodes, "shells": shells, "spress": spress, "nloads": nloads}
    return model, nid(n // 2, n // 2)


def rotate_xyz(R: list[list[float]], v: list[float]) -> list[float]:
    return [sum(R[i][j] * v[j] for j in range(3)) for i in range(3)]


def rotate_model(model: dict, R: list[list[float]]) -> dict:
    m = copy.deepcopy(model)
    for nd in m["nodes"]:
        x = rotate_xyz(R, [nd["x"], nd["y"], nd["z"]])
        nd["x"], nd["y"], nd["z"] = x
    return m


def rodrigues(axis: tuple[float, float, float], theta: float) -> list[list[float]]:
    ax = list(axis)
    an = norm(ax)
    kx, ky, kz = [x / an for x in ax]
    ct, st, vt = math.cos(theta), math.sin(theta), 1.0 - math.cos(theta)
    return [
        [ct + kx * kx * vt, kx * ky * vt - kz * st, kx * kz * vt + ky * st],
        [ky * kx * vt + kz * st, ct + ky * ky * vt, ky * kz * vt - kx * st],
        [kz * kx * vt - ky * st, kz * ky * vt + kx * st, ct + kz * kz * vt],
    ]


def max_disp_rel(a: dict[int, list[float]], b: dict[int, list[float]], dofs=range(6)) -> tuple[float, float, float]:
    max_abs = 0.0
    scale = 0.0
    for node_id, va in a.items():
        vb = b[node_id]
        max_abs = max(max_abs, max(abs(va[d] - vb[d]) for d in dofs))
        scale = max(scale, max(abs(va[d]) for d in dofs), max(abs(vb[d]) for d in dofs))
    return max_abs / max(scale, 1e-30), max_abs, scale


def cantilever_shell_model() -> dict:
    nx, ny = 7, 5
    Lx, Ly, t = 1200.0, 700.0, 12.0
    tilt = math.radians(23.0)
    ct, st = math.cos(tilt), math.sin(tilt)
    base = 2000

    def nid(i, j):
        return base + j * (nx + 1) + i

    nodes = []
    for j in range(ny + 1):
        for i in range(nx + 1):
            x0 = Lx * i / nx
            y0 = Ly * j / ny + 0.16 * (Lx / nx) * (i / nx) * math.sin(math.pi * j / ny)
            z0 = 0.0
            x = ct * x0 + st * z0
            z = -st * x0 + ct * z0
            fix = [1, 1, 1, 1, 1, 1] if i == 0 else [0, 0, 0, 0, 0, 0]
            nodes.append({"id": nid(i, j), "x": x, "y": y0, "z": z, "fix": fix})

    shells = []
    sid = 0
    for j in range(ny):
        for i in range(nx):
            shells.append({
                "id": 8000 + sid,
                "n": [nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)],
                "mat": 0,
                "t": t,
            })
            sid += 1
    nloads = []
    for j in range(ny + 1):
        nloads.append({"node": nid(nx, j), "comp": [-12.0, 3.0, -150.0 / (ny + 1), 0.0, 0.0, 0.0]})
    return {"mats": [mat()], "nodes": nodes, "shells": shells, "spress": [], "nloads": nloads}


def run_opensees_shell(model: dict) -> dict | None:
    try:
        import openseespy.opensees as ops
    except Exception:
        return None

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    for n in model["nodes"]:
        ops.node(n["id"], float(n["x"]), float(n["y"]), float(n["z"]))
        ops.fix(n["id"], *[int(x) for x in n["fix"]])
    for s in model["shells"]:
        m = model["mats"][s["mat"]]
        sec = s["id"] + 1
        ops.section("ElasticMembranePlateSection", sec, m["E"], m["nu"], s["t"], 0.0)
        ops.element("ShellMITC4", s["id"] + 1, *s["n"], sec)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for l in model.get("nloads", []):
        ops.load(l["node"], *[float(x) for x in l["comp"]])
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1.0e-12, 100)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    ok = ops.analyze(1)
    if ok != 0:
        return {"ok": ok, "disp": {}}
    return {"ok": 0, "disp": {n["id"]: [ops.nodeDisp(n["id"], k) for k in range(1, 7)] for n in model["nodes"]}}


def scordelis_model(n: int) -> tuple[dict, int]:
    R, L, phi0, t, g = 25.0, 50.0, 40.0 * math.pi / 180.0, 0.25, 90.0
    E, nu = 4.32e8, 0.0
    base = 3000

    def nid(i, j):
        return base + j * (n + 1) + i

    nodes = []
    for j in range(n + 1):
        for i in range(n + 1):
            phi = phi0 * i / n
            y = (L * 0.5) * j / n
            fix = [0, 0, 0, 0, 0, 0]
            if j == 0:
                fix[UX] = fix[UZ] = 1
            if j == n:
                fix[UY] = fix[RX] = fix[RZ] = 1
            if i == 0:
                fix[UX] = fix[RY] = fix[RZ] = 1
            nodes.append({"id": nid(i, j), "x": R * math.sin(phi), "y": y, "z": R * math.cos(phi), "fix": fix})
    nodes_by_id = {nd["id"]: nd for nd in nodes}
    shells = []
    fz: dict[int, float] = {}
    sid = 0
    for j in range(n):
        for i in range(n):
            conn = [nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)]
            shells.append({"id": 9000 + sid, "n": conn, "mat": 0, "t": t})
            fn = -g * quad_area(nodes_by_id, conn) * 0.25
            for node_id in conn:
                fz[node_id] = fz.get(node_id, 0.0) + fn
            sid += 1
    nloads = [{"node": k, "comp": [0.0, 0.0, v, 0.0, 0.0, 0.0]} for k, v in sorted(fz.items()) if v]
    return {"mats": [mat(E, nu)], "nodes": nodes, "shells": shells, "spress": [], "nloads": nloads}, nid(n, n)


def pinched_cylinder_model(n: int) -> tuple[dict, int]:
    R, L, t, P = 300.0, 600.0, 3.0, 1.0
    E, nu = 3.0e6, 0.3
    base = 4000

    def nid(i, j):
        return base + j * (n + 1) + i

    nodes = []
    for j in range(n + 1):
        for i in range(n + 1):
            th = 0.5 * math.pi * i / n
            z = (L * 0.5) * j / n
            fix = [0, 0, 0, 0, 0, 0]
            if i == 0:
                fix[UY] = fix[RX] = fix[RZ] = 1
            if i == n:
                fix[UX] = fix[RY] = fix[RZ] = 1
            if j == n:
                fix[UZ] = fix[RX] = fix[RY] = 1
            if j == 0:
                fix[UX] = fix[UY] = 1
            nodes.append({"id": nid(i, j), "x": R * math.cos(th), "y": R * math.sin(th), "z": z, "fix": fix})
    shells = []
    sid = 0
    for j in range(n):
        for i in range(n):
            shells.append({"id": 10000 + sid, "n": [nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)], "mat": 0, "t": t})
            sid += 1
    load_node = nid(0, n)
    nloads = [{"node": load_node, "comp": [-0.25 * P, 0.0, 0.0, 0.0, 0.0, 0.0]}]
    return {"mats": [mat(E, nu)], "nodes": nodes, "shells": shells, "spress": [], "nloads": nloads}, load_node


def bad_model(kind: str) -> dict:
    nodes = [
        {"id": 10, "x": 0.0, "y": 0.0, "z": 0.0, "fix": [1, 1, 1, 1, 1, 1]},
        {"id": 20, "x": 100.0, "y": 0.0, "z": 0.0, "fix": [1, 1, 1, 1, 1, 1]},
        {"id": 30, "x": 100.0, "y": 100.0, "z": 0.0, "fix": [0, 0, 0, 0, 0, 0]},
        {"id": 40, "x": 0.0, "y": 100.0, "z": 0.0, "fix": [1, 1, 1, 1, 1, 1]},
    ]
    if kind == "warped":
        nodes[2]["z"] = 35.0
        conn = [10, 20, 30, 40]
    elif kind == "concave":
        nodes[2]["x"] = 50.0
        nodes[2]["y"] = 30.0
        conn = [10, 20, 30, 40]
    elif kind == "bowtie":
        conn = [10, 30, 20, 40]
    elif kind == "duplicate":
        conn = [10, 20, 20, 40]
    else:
        conn = [10, 20, 30, 40]
    return {
        "mats": [mat()],
        "nodes": nodes,
        "shells": [{"id": 1, "n": conn, "mat": 0, "t": 10.0}],
        "spress": [{"shell": 1, "p": -0.01}],
        "nloads": [{"node": 30, "comp": [0.0, 0.0, -1.0, 0.0, 0.0, 0.0]}],
    }


def duplicate_shell_id_model() -> dict:
    model, _ = square_plate(2, "simple", pressure=False, id_base=5000, shell_base=1)
    for s in model["shells"]:
        s["id"] = 77
    model["spress"] = [{"shell": 77, "p": -0.01}]
    return model


def main() -> int:
    audit = Audit()
    ensure_cli()
    print("=" * 88)
    print("MITC4 shell deep audit (black-box frame_cli + independent Python model generation)")
    print("=" * 88)

    D = plate_D(E0, NU0, 10.0)
    w_ss_ref = 0.00406 * 0.01 * 1000.0 ** 4 / D
    ss = []
    for n in (4, 8, 16, 24):
        model, c = square_plate(n, "simple")
        r = run_cli(model)
        w = abs(r["disp"][c][UZ]) if not r["singular"] else float("inf")
        ss.append((n, w, abs(w - w_ss_ref) / w_ss_ref))
    print("\n[flat square plate: simply supported, uniform pressure]")
    for n, w, e in ss:
        print(f"  N={n:<2d} w={w:.8g} mm  err={100.0 * e:.3f}%")
    audit.check("SS plate N=16 within 1% of Kirchhoff", ss[2][2] < 0.01, f"err={ss[2][2]:.3e}")
    audit.check("SS plate mesh convergence improves from N=4 to N=24", ss[-1][2] < ss[0][2], f"e4={ss[0][2]:.3e} e24={ss[-1][2]:.3e}")

    w_cl_ref = 0.00126 * 0.01 * 1000.0 ** 4 / D
    cl = []
    for n in (8, 16, 24):
        model, c = square_plate(n, "clamped")
        r = run_cli(model)
        w = abs(r["disp"][c][UZ]) if not r["singular"] else float("inf")
        cl.append((n, w, abs(w - w_cl_ref) / w_cl_ref))
    print("\n[flat square plate: clamped, uniform pressure]")
    for n, w, e in cl:
        print(f"  N={n:<2d} w={w:.8g} mm  err={100.0 * e:.3f}%")
    audit.check("clamped plate N=16 within 2% of theory", cl[1][2] < 0.02, f"err={cl[1][2]:.3e}")
    audit.check("clamped plate mesh convergence improves from N=8 to N=24", cl[-1][2] < cl[0][2], f"e8={cl[0][2]:.3e} e24={cl[-1][2]:.3e}")

    thin_t = 1.0
    thin_D = plate_D(E0, NU0, thin_t)
    thin_ref = 0.00406 * 0.01 * 1000.0 ** 4 / thin_D
    model, c = square_plate(16, "simple", t=thin_t)
    r = run_cli(model)
    thin_w = abs(r["disp"][c][UZ])
    thin_err = abs(thin_w - thin_ref) / thin_ref
    print("\n[thin plate locking check]")
    print(f"  t/a=0.001  w={thin_w:.8g} mm  ref={thin_ref:.8g} mm  err={100.0 * thin_err:.3f}%")
    audit.check("thin plate does not shear-lock", thin_err < 0.02, f"err={thin_err:.3e}")

    print("\n[load path, reaction equilibrium, and orientation]")
    mp, _ = square_plate(10, "clamped", pressure=True, lumped=False)
    ml, _ = square_plate(10, "clamped", pressure=False, lumped=True)
    rp = run_cli(mp)
    rl = run_cli(ml)
    rel, absd, scale = max_disp_rel(rp["disp"], rl["disp"])
    audit.check("SPRESS equals independently lumped nodal pressure", rel < 1e-9, f"rel={rel:.3e} abs={absd:.3e}")
    rfz = sum(v[UZ] for v in rp["rf"].values())
    total_load = 0.01 * 1000.0 * 1000.0
    audit.check("reaction sum balances shell pressure", abs(rfz - total_load) / total_load < 1e-9, f"sumRz={rfz:.8g} load={total_load:.8g}")
    mr, _ = square_plate(8, "simple", reverse=True)
    mb, _ = square_plate(8, "simple", reverse=False)
    rr = run_cli(mr)
    rb = run_cli(mb)
    rel, absd, _ = max_disp_rel(rr["disp"], rb["disp"], dofs=(UX, UY, UZ))
    audit.check("reversed shell orientation with adjusted pressure gives same translations", rel < 1e-9, f"rel={rel:.3e} abs={absd:.3e}")

    print("\n[3D rotation equivariance]")
    base, _ = square_plate(8, "clamped")
    R = rodrigues((1.0, 2.0, 4.0), 0.73)
    rot = rotate_model(base, R)
    rb = run_cli(base)
    rr = run_cli(rot)
    max_abs = 0.0
    scale = 0.0
    for node_id, u0 in rb["disp"].items():
        u0t = rotate_xyz(R, u0[:3])
        u0r = rotate_xyz(R, u0[3:6])
        u1 = rr["disp"][node_id]
        max_abs = max(max_abs, max(abs(u1[i] - u0t[i]) for i in range(3)))
        max_abs = max(max_abs, max(abs(u1[i + 3] - u0r[i]) for i in range(3)))
        scale = max(scale, max(abs(x) for x in u0t + u0r + u1))
    rel = max_abs / max(scale, 1e-30)
    audit.check("full displacement/rotation vector transforms under rigid 3D rotation", rel < 1e-8, f"rel={rel:.3e} abs={max_abs:.3e}")

    print("\n[distorted quadrilateral mesh]")
    dist = []
    for n in (8, 16):
        model, c = square_plate(n, "simple", perturb=0.18)
        r = run_cli(model)
        w = abs(r["disp"][c][UZ])
        e = abs(w - w_ss_ref) / w_ss_ref
        dist.append((n, w, e))
        print(f"  N={n:<2d} w={w:.8g} mm  err={100.0 * e:.3f}%")
    audit.check("distorted mesh remains accurate at N=16", dist[1][2] < 0.025, f"err={dist[1][2]:.3e}")
    audit.check("distorted mesh improves with refinement", dist[1][2] < dist[0][2], f"e8={dist[0][2]:.3e} e16={dist[1][2]:.3e}")

    print("\n[OpenSees ShellMITC4 extra cross-check]")
    m_os = cantilever_shell_model()
    r_cli = run_cli(m_os)
    r_os = run_opensees_shell(m_os)
    if r_os is None:
        audit.warn("OpenSeesPy not available", "skipping extra shell cross-check")
    elif r_os["ok"] != 0:
        audit.check("OpenSees solved skew/tilted cantilever shell", False, f"ok={r_os['ok']}")
    else:
        rel, absd, scale = max_disp_rel(r_cli["disp"], r_os["disp"])
        print(f"  max displacement diff vs OpenSees = {rel:.3e} rel  ({absd:.3e} abs / {scale:.3e} scale)")
        audit.check("skew/tilted cantilever shell matches OpenSees", rel < 1e-7, f"rel={rel:.3e}")

    print("\n[curved-shell benchmark: Scordelis-Lo roof]")
    roof_ref = 0.3024
    roof = []
    for n in (8, 16, 24):
        model, c = scordelis_model(n)
        r = run_cli(model)
        w = abs(r["disp"][c][UZ])
        e = abs(w - roof_ref) / roof_ref
        roof.append((n, w, e))
        print(f"  N={n:<2d} w={w:.8g}  err={100.0 * e:.3f}%")
    audit.check("Scordelis-Lo N=24 within 3% of reference", roof[-1][2] < 0.03, f"err={roof[-1][2]:.3e}")
    audit.check("Scordelis-Lo converges toward reference", roof[-1][2] < roof[0][2], f"e8={roof[0][2]:.3e} e24={roof[-1][2]:.3e}")

    print("\n[curved-shell benchmark: pinched cylinder]")
    cyl_ref = 1.8248e-5
    cyl = []
    for n in (8, 16, 24, 32):
        model, c = pinched_cylinder_model(n)
        r = run_cli(model)
        w = abs(r["disp"][c][UX])
        cyl.append((n, w, w / cyl_ref))
        print(f"  N={n:<2d} w={w:.8e}  ratio={100.0 * w / cyl_ref:.2f}%")
    audit.check("pinched cylinder converges upward", cyl[3][1] > cyl[2][1] > cyl[1][1] > cyl[0][1], "")
    audit.check("pinched cylinder N=32 in [90%,105%] of reference", 0.90 <= cyl[3][2] <= 1.05, f"ratio={cyl[3][2]:.3f}")

    print("\n[invalid/ambiguous input handling]")
    miss = square_plate(2, "simple", pressure=False)[0]
    miss["spress"] = [{"shell": 999999, "p": -0.01}]
    audit.check("missing shell pressure reference is rejected", run_cli(miss)["singular"] == 1, "")
    audit.check("duplicate shell corner node is rejected", run_cli(bad_model("duplicate"))["singular"] == 1, "")
    audit.check("concave shell quad is rejected", run_cli(bad_model("concave"))["singular"] == 1, "")
    audit.check("bow-tie/self-intersecting quad is rejected or singular", run_cli(bad_model("bowtie"))["singular"] == 1, "")
    dup = run_cli(duplicate_shell_id_model())
    if dup["singular"] == 0:
        audit.warn("duplicate ShellQuad ids are accepted", "ShellPressure by id becomes ambiguous")
    else:
        audit.check("duplicate ShellQuad ids rejected", True, "")
    warped = run_cli(bad_model("warped"))
    if warped["singular"] == 0:
        audit.warn("non-coplanar four-node shell is accepted", "flat-shell facet silently projects warped geometry")
    else:
        audit.check("non-coplanar shell rejected", True, "")
    bad_nu = bad_model("normal")
    bad_nu["mats"] = [mat(E0, 0.60)]
    audit.check("invalid shell Poisson ratio is rejected", run_cli(bad_nu)["singular"] == 1, "")

    print("\n[CLI end-to-end timing: clamped pressure plate]")
    for n in (8, 16, 24, 32):
        model, _ = square_plate(n, "clamped")
        t0 = time.perf_counter()
        r = run_cli(model)
        ms = (time.perf_counter() - t0) * 1000.0
        ndof = 6 * len(model["nodes"])
        audit.report(f"N={n} mesh timing", f"nodes={len(model['nodes'])} shells={len(model['shells'])} dof={ndof} singular={r['singular']} wall={ms:.1f} ms")

    print("-" * 88)
    print(f"checks failed={audit.failures} warnings={audit.warnings}")
    return 1 if audit.failures else 0


if __name__ == "__main__":
    sys.exit(main())
