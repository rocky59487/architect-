#!/usr/bin/env python
"""
Independent precision audit for FrameCore's line-elastic beam solver.

This is intentionally separate from the UE automation fixtures. It treats
frame_cli.exe as a black-box executable, feeds it hand-built and generated models,
then compares the output against:
  1. closed-form beam formulas where available;
  2. an independent dense direct-stiffness implementation in Python/numpy;
  3. rigid-body rotation equivariance.

Units follow FrameCore: N, mm, MPa.
"""
from __future__ import annotations

import math
import os
import subprocess
import copy
from dataclasses import dataclass
from typing import Iterable

import numpy as np


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CLI = os.path.join(ROOT, "Plugins", "FrameSolver", "Standalone", "frame_cli.exe")

DOF = 6
Ux, Uy, Uz, Rx, Ry, Rz = range(6)


@dataclass
class Check:
    name: str
    value: float
    tol: float
    ok: bool
    note: str = ""


def norm(v: np.ndarray) -> float:
    return float(np.linalg.norm(v))


def rect_section(b: float, d: float) -> dict:
    # b is local-z width, d is local-y depth.
    a = b * d
    iz = b * d**3 / 12.0
    iy = d * b**3 / 12.0
    bb, tt = max(b, d), min(b, d)
    j = bb * tt**3 * (1.0 / 3.0 - 0.21 * (tt / bb) * (1.0 - tt**4 / (12.0 * bb**4)))
    return dict(A=a, Iy=iy, Iz=iz, J=j, cy=d / 2.0, cz=b / 2.0, Asy=0.0, Asz=0.0)


def square_section(side: float) -> dict:
    s = rect_section(side, side)
    s["J"] = 0.1406 * side**4
    return s


def refvec_for(pi: Iterable[float], pj: Iterable[float]) -> tuple[float, float, float]:
    pi = np.array(tuple(pi), dtype=float)
    pj = np.array(tuple(pj), dtype=float)
    x = pj - pi
    x /= norm(x)
    # Prefer global Z, but avoid degeneracy.
    if abs(float(np.dot(x, np.array([0.0, 0.0, 1.0])))) < 0.92:
        return (0.0, 0.0, 1.0)
    return (1.0, 0.0, 0.0)


def model_to_cli(model: dict) -> str:
    lines: list[str] = []
    for m in model["materials"]:
        lines.append(f"MAT {m['E']:.17g} {m['G']:.17g} {m.get('rho', 0.0):.17g}")
    for s in model["sections"]:
        vals = (s["A"], s["Iy"], s["Iz"], s["J"], s["cy"], s["cz"], s.get("Asy", 0.0), s.get("Asz", 0.0))
        lines.append("SEC " + " ".join(f"{v:.17g}" for v in vals))
    for n in model["nodes"]:
        f = n["fix"]
        lines.append(
            f"NODE {n['id']} {n['x']:.17g} {n['y']:.17g} {n['z']:.17g} "
            + " ".join(str(int(v)) for v in f)
        )
    for e in model["members"]:
        rv = e["refvec"]
        lines.append(
            f"MEMBER {e['id']} {e['i']} {e['j']} {e['mat']} {e['sec']} "
            + " ".join(f"{v:.17g}" for v in rv)
        )
    for l in model.get("nloads", []):
        lines.append(f"NLOAD {l['node']} " + " ".join(f"{v:.17g}" for v in l["comp"]))
    for u in model.get("udls", []):
        lines.append(f"UDL {u['member']} " + " ".join(f"{v:.17g}" for v in u["w_local"]))
    opt = model.get("opt")
    if opt:
        lines.append(f"OPT {int(opt.get('enableReleases', False))} {int(opt.get('useTimoshenko', False))} {opt.get('pivotTol', 1e-12):.17g}")
    lines.append("END")
    return "\n".join(lines) + "\n"


def run_cli(model: dict) -> dict:
    if not os.path.exists(CLI):
        raise FileNotFoundError(CLI)
    p = subprocess.run([CLI], input=model_to_cli(model), capture_output=True, text=True, check=False)
    if p.returncode != 0:
        raise RuntimeError(f"frame_cli failed rc={p.returncode}\nstdout={p.stdout}\nstderr={p.stderr}")
    out = {"singular": None, "disp": {}, "mf": {}}
    for line in p.stdout.splitlines():
        toks = line.split()
        if not toks:
            continue
        if toks[0] == "SINGULAR":
            out["singular"] = int(toks[1])
        elif toks[0] == "DISP":
            out["disp"][int(toks[1])] = np.array([float(v) for v in toks[2:]], dtype=float)
        elif toks[0] == "MF":
            out["mf"][int(toks[1])] = np.array([float(v) for v in toks[2:]], dtype=float)
    return out


def local_axes(pi: np.ndarray, pj: np.ndarray, refvec: Iterable[float]) -> np.ndarray:
    x = pj - pi
    x /= norm(x)
    ref = np.array(tuple(refvec), dtype=float)
    z = np.cross(x, ref)
    if norm(z) < 1.0e-8:
        ref = np.array([0.0, 1.0, 0.0])
        z = np.cross(x, ref)
        if norm(z) < 1.0e-8:
            ref = np.array([1.0, 0.0, 0.0])
            z = np.cross(x, ref)
    z /= norm(z)
    y = np.cross(z, x)
    y /= norm(y)
    return np.vstack([x, y, z])


def local_stiffness(E: float, G: float, A: float, Iy: float, Iz: float, J: float, L: float) -> np.ndarray:
    k = np.zeros((12, 12), dtype=float)
    EA = E * A / L
    GJ = G * J / L
    L2 = L * L
    L3 = L2 * L
    az, bz, cz, dz = 12.0 * E * Iz / L3, 6.0 * E * Iz / L2, 4.0 * E * Iz / L, 2.0 * E * Iz / L
    ay, by, cy, dy = 12.0 * E * Iy / L3, 6.0 * E * Iy / L2, 4.0 * E * Iy / L, 2.0 * E * Iy / L

    k[0, 0] = EA
    k[0, 6] = -EA
    k[6, 0] = -EA
    k[6, 6] = EA
    k[3, 3] = GJ
    k[3, 9] = -GJ
    k[9, 3] = -GJ
    k[9, 9] = GJ

    k[1, 1] = az
    k[1, 5] = bz
    k[1, 7] = -az
    k[1, 11] = bz
    k[5, 1] = bz
    k[5, 5] = cz
    k[5, 7] = -bz
    k[5, 11] = dz
    k[7, 1] = -az
    k[7, 5] = -bz
    k[7, 7] = az
    k[7, 11] = -bz
    k[11, 1] = bz
    k[11, 5] = dz
    k[11, 7] = -bz
    k[11, 11] = cz

    k[2, 2] = ay
    k[2, 4] = -by
    k[2, 8] = -ay
    k[2, 10] = -by
    k[4, 2] = -by
    k[4, 4] = cy
    k[4, 8] = by
    k[4, 10] = dy
    k[8, 2] = -ay
    k[8, 4] = by
    k[8, 8] = ay
    k[8, 10] = by
    k[10, 2] = -by
    k[10, 4] = dy
    k[10, 8] = by
    k[10, 10] = cy
    return k


def local_stiffness_timoshenko(E: float, G: float, A: float, Iy: float, Iz: float,
                               J: float, L: float, Asy: float, Asz: float) -> np.ndarray:
    k = np.zeros((12, 12), dtype=float)
    EA = E * A / L
    GJ = G * J / L
    L2 = L * L

    phiz = 12.0 * E * Iz / (G * Asy * L2) if Asy > 0.0 else 0.0
    bz0 = E * Iz / ((1.0 + phiz) * L * L2)
    az, bz = 12.0 * bz0, 6.0 * bz0 * L
    cz, dz = (4.0 + phiz) * bz0 * L2, (2.0 - phiz) * bz0 * L2

    phiy = 12.0 * E * Iy / (G * Asz * L2) if Asz > 0.0 else 0.0
    by0 = E * Iy / ((1.0 + phiy) * L * L2)
    ay, by = 12.0 * by0, 6.0 * by0 * L
    cy, dy = (4.0 + phiy) * by0 * L2, (2.0 - phiy) * by0 * L2

    k[0, 0] = EA
    k[0, 6] = -EA
    k[6, 0] = -EA
    k[6, 6] = EA
    k[3, 3] = GJ
    k[3, 9] = -GJ
    k[9, 3] = -GJ
    k[9, 9] = GJ

    k[1, 1] = az
    k[1, 5] = bz
    k[1, 7] = -az
    k[1, 11] = bz
    k[5, 1] = bz
    k[5, 5] = cz
    k[5, 7] = -bz
    k[5, 11] = dz
    k[7, 1] = -az
    k[7, 5] = -bz
    k[7, 7] = az
    k[7, 11] = -bz
    k[11, 1] = bz
    k[11, 5] = dz
    k[11, 7] = -bz
    k[11, 11] = cz

    k[2, 2] = ay
    k[2, 4] = -by
    k[2, 8] = -ay
    k[2, 10] = -by
    k[4, 2] = -by
    k[4, 4] = cy
    k[4, 8] = by
    k[4, 10] = dy
    k[8, 2] = -ay
    k[8, 4] = by
    k[8, 8] = ay
    k[8, 10] = by
    k[10, 2] = -by
    k[10, 4] = dy
    k[10, 8] = by
    k[10, 10] = cy
    return k


def transform12(R: np.ndarray) -> np.ndarray:
    T = np.zeros((12, 12), dtype=float)
    for b in range(4):
        T[3 * b : 3 * b + 3, 3 * b : 3 * b + 3] = R
    return T


def fixed_end_udl(w_local: Iterable[float], L: float) -> np.ndarray:
    wx, wy, wz = (-float(v) for v in w_local)
    q = np.zeros(12, dtype=float)
    q[0] += wx * L / 2.0
    q[6] += wx * L / 2.0
    q[1] += wy * L / 2.0
    q[7] += wy * L / 2.0
    q[5] += wy * L * L / 12.0
    q[11] += -wy * L * L / 12.0
    q[2] += wz * L / 2.0
    q[8] += wz * L / 2.0
    q[4] += -wz * L * L / 12.0
    q[10] += wz * L * L / 12.0
    return q


def py_solve(model: dict) -> dict:
    nodes = model["nodes"]
    members = model["members"]
    node_index = {n["id"]: i for i, n in enumerate(nodes)}
    ndof = len(nodes) * DOF
    K = np.zeros((ndof, ndof), dtype=float)
    F = np.zeros(ndof, dtype=float)
    member_cache = []

    for load in model.get("nloads", []):
        ni = node_index[load["node"]]
        F[ni * DOF : ni * DOF + DOF] += np.array(load["comp"], dtype=float)

    for e in members:
        ni = node_index[e["i"]]
        nj = node_index[e["j"]]
        pi = np.array([nodes[ni]["x"], nodes[ni]["y"], nodes[ni]["z"]], dtype=float)
        pj = np.array([nodes[nj]["x"], nodes[nj]["y"], nodes[nj]["z"]], dtype=float)
        L = norm(pj - pi)
        mat = model["materials"][e["mat"]]
        sec = model["sections"][e["sec"]]
        use_timoshenko = bool(model.get("opt", {}).get("useTimoshenko", False))
        if use_timoshenko and sec.get("Asy", 0.0) > 0.0 and sec.get("Asz", 0.0) > 0.0:
            kl = local_stiffness_timoshenko(mat["E"], mat["G"], sec["A"], sec["Iy"], sec["Iz"],
                                            sec["J"], L, sec["Asy"], sec["Asz"])
        else:
            kl = local_stiffness(mat["E"], mat["G"], sec["A"], sec["Iy"], sec["Iz"], sec["J"], L)
        T = transform12(local_axes(pi, pj, e["refvec"]))
        qf = np.zeros(12, dtype=float)
        for udl in model.get("udls", []):
            if udl["member"] == e["id"]:
                qf += fixed_end_udl(udl["w_local"], L)
        kg = T.T @ kl @ T
        dofs = list(range(ni * DOF, ni * DOF + DOF)) + list(range(nj * DOF, nj * DOF + DOF))
        for a, ga in enumerate(dofs):
            for b, gb in enumerate(dofs):
                K[ga, gb] += kg[a, b]
        if np.any(qf):
            peq = -(T.T @ qf)
            for a, ga in enumerate(dofs):
                F[ga] += peq[a]
        member_cache.append((e, kl, T, qf, dofs))

    fixed = []
    free = []
    for i, n in enumerate(nodes):
        for d, is_fixed in enumerate(n["fix"]):
            (fixed if is_fixed else free).append(i * DOF + d)
    if not free:
        return {"singular": 1, "u": np.zeros(ndof), "reactions": np.zeros(ndof), "mf": {}, "cond": math.inf}

    Kff = K[np.ix_(free, free)]
    Ff = F[free]
    try:
        cond = float(np.linalg.cond(Kff))
    except np.linalg.LinAlgError:
        cond = math.inf
    try:
        uf = np.linalg.solve(Kff, Ff)
        singular = 0
    except np.linalg.LinAlgError:
        uf = np.zeros(len(free), dtype=float)
        singular = 1
    u = np.zeros(ndof, dtype=float)
    if not singular:
        u[free] = uf
    reactions = K @ u - F

    mf = {}
    if not singular:
        for idx, (e, kl, T, qf, dofs) in enumerate(member_cache):
            ue = u[dofs]
            Q = kl @ (T @ ue) + qf
            mf[idx] = np.array([Q[0], Q[1], Q[2], Q[3], Q[4], Q[5], -Q[6], Q[7], Q[8], Q[9], Q[10], Q[11]])
    return {"singular": singular, "u": u, "reactions": reactions, "mf": mf, "cond": cond}


def disp_array_cli(out: dict, nnode: int) -> np.ndarray:
    arr = np.zeros(nnode * DOF, dtype=float)
    for i in range(nnode):
        arr[i * DOF : i * DOF + DOF] = out["disp"][i]
    return arr


def max_scaled_diff(a: np.ndarray, b: np.ndarray) -> float:
    scale = max(float(np.max(np.abs(a))), float(np.max(np.abs(b))), 1.0e-30)
    return float(np.max(np.abs(a - b)) / scale)


def add_check(checks: list[Check], name: str, value: float, tol: float, note: str = "") -> None:
    checks.append(Check(name=name, value=value, tol=tol, ok=bool(value <= tol), note=note))


def clone_model(model: dict) -> dict:
    return copy.deepcopy(model)


def scale_loads(model: dict, factor: float) -> dict:
    m = clone_model(model)
    for load in m.get("nloads", []):
        load["comp"] = [factor * float(v) for v in load["comp"]]
    for udl in m.get("udls", []):
        udl["w_local"] = tuple(factor * float(v) for v in udl["w_local"])
    return m


def disp_by_node_id(model: dict, out: dict) -> dict[int, np.ndarray]:
    return {int(node["id"]): out["disp"][int(node["id"])] for node in model["nodes"]}


def mf_by_member_id(model: dict, out: dict) -> dict[int, np.ndarray]:
    return {int(member["id"]): out["mf"][int(member["id"])] for member in model["members"] if int(member["id"]) in out["mf"]}


def stacked_disp_by_id(model: dict, out: dict, ids: Iterable[int] | None = None) -> np.ndarray:
    d = disp_by_node_id(model, out)
    use_ids = sorted(d) if ids is None else list(ids)
    return np.concatenate([d[i] for i in use_ids])


def stacked_mf_by_id(model: dict, out: dict, ids: Iterable[int] | None = None) -> np.ndarray:
    mf = mf_by_member_id(model, out)
    use_ids = sorted(mf) if ids is None else list(ids)
    return np.concatenate([mf[i] for i in use_ids])


def load_work(model: dict, out: dict, loads: list[dict]) -> float:
    d = disp_by_node_id(model, out)
    total = 0.0
    for load in loads:
        total += float(np.dot(np.array(load["comp"], dtype=float), d[int(load["node"])]))
    return total


def compare_cli_to_python(checks: list[Check], name: str, model: dict, tol_disp: float, tol_force: float) -> None:
    cli = run_cli(model)
    py = py_solve(model)
    add_check(checks, f"{name}: singular flag", abs(float(cli["singular"]) - float(py["singular"])), 0.0)
    if cli["singular"] or py["singular"]:
        return
    u_cli = disp_array_cli(cli, len(model["nodes"]))
    add_check(checks, f"{name}: displacement vs independent DSM", max_scaled_diff(u_cli, py["u"]), tol_disp)
    all_cli = np.concatenate([cli["mf"][i] for i in sorted(cli["mf"])])
    all_py = np.concatenate([py["mf"][i] for i in sorted(py["mf"])])
    add_check(checks, f"{name}: member forces vs independent DSM", max_scaled_diff(all_cli, all_py), tol_force)

    # Reactions are not printed by CLI. Check global force and moment equilibrium
    # via the independent solve.
    total_force = np.zeros(3, dtype=float)
    total_moment = np.zeros(3, dtype=float)
    for l in model.get("nloads", []):
        ni = next(i for i, n in enumerate(model["nodes"]) if n["id"] == l["node"])
        r = np.array([model["nodes"][ni]["x"], model["nodes"][ni]["y"], model["nodes"][ni]["z"]], dtype=float)
        c = np.array(l["comp"], dtype=float)
        total_force += c[:3]
        total_moment += c[3:] + np.cross(r, c[:3])
    for e in model["members"]:
        ni = next(i for i, n in enumerate(model["nodes"]) if n["id"] == e["i"])
        nj = next(i for i, n in enumerate(model["nodes"]) if n["id"] == e["j"])
        pi = np.array([model["nodes"][ni]["x"], model["nodes"][ni]["y"], model["nodes"][ni]["z"]], dtype=float)
        pj = np.array([model["nodes"][nj]["x"], model["nodes"][nj]["y"], model["nodes"][nj]["z"]], dtype=float)
        L = norm(pj - pi)
        R = local_axes(pi, pj, e["refvec"])
        for udl in model.get("udls", []):
            if udl["member"] == e["id"]:
                w_global = R.T @ np.array(udl["w_local"], dtype=float)
                f = w_global * L
                total_force += f
                total_moment += np.cross((pi + pj) / 2.0, f)
    reaction_force = np.zeros(3, dtype=float)
    reaction_moment = np.zeros(3, dtype=float)
    for i, n in enumerate(model["nodes"]):
        r = np.array([n["x"], n["y"], n["z"]], dtype=float)
        rr = py["reactions"][i * DOF : i * DOF + DOF]
        reaction_force += rr[:3]
        reaction_moment += rr[3:] + np.cross(r, rr[:3])
    force_scale = max(float(np.max(np.abs(total_force))), 1.0)
    moment_scale = max(float(np.max(np.abs(total_moment))), 1.0)
    add_check(checks, f"{name}: global force equilibrium",
              float(np.max(np.abs(reaction_force + total_force)) / force_scale), 1.0e-9)
    add_check(checks, f"{name}: global moment equilibrium",
              float(np.max(np.abs(reaction_moment + total_moment)) / moment_scale), 1.0e-9)


def cantilever_closed_form_case() -> dict:
    E, G = 199500.0, 76800.0
    L = 2173.25
    sec = rect_section(137.0, 271.0)
    model = {
        "materials": [dict(E=E, G=G, rho=7850.0)],
        "sections": [sec],
        "nodes": [
            dict(id=0, x=0.0, y=0.0, z=0.0, fix=[1, 1, 1, 1, 1, 1]),
            dict(id=1, x=L, y=0.0, z=0.0, fix=[0, 0, 0, 0, 0, 0]),
        ],
        "members": [dict(id=0, i=0, j=1, mat=0, sec=0, refvec=(0.0, 0.0, 1.0))],
        "nloads": [dict(node=1, comp=[1234.0, -777.0, 991.0, 50000.0, 0.0, 0.0])],
    }
    return model


def check_cantilever_closed_form(checks: list[Check]) -> None:
    model = cantilever_closed_form_case()
    cli = run_cli(model)
    sec = model["sections"][0]
    mat = model["materials"][0]
    L = model["nodes"][1]["x"]
    Fx, Fy, Fz, Mx, _, _ = model["nloads"][0]["comp"]
    d = cli["disp"][1]

    # With refVec = global Z, local y is global Z and local z is -global Y.
    ux = Fx * L / (mat["E"] * sec["A"])
    uy = Fy * L**3 / (3.0 * mat["E"] * sec["Iy"])
    uz = Fz * L**3 / (3.0 * mat["E"] * sec["Iz"])
    rx = Mx * L / (mat["G"] * sec["J"])
    ry = -Fz * L * L / (2.0 * mat["E"] * sec["Iz"])
    rz = Fy * L * L / (2.0 * mat["E"] * sec["Iy"])
    exact = np.array([ux, uy, uz, rx, ry, rz], dtype=float)
    add_check(checks, "closed form cantilever: tip DOF vector", max_scaled_diff(d, exact), 2.0e-10,
              "limited by frame_cli %.12g output")

    mf = cli["mf"][0]
    root = np.array([-Fx, -Fz, Fy, -Mx, -Fy * L, -Fz * L], dtype=float)
    add_check(checks, "closed form cantilever: root force vector", max_scaled_diff(mf[:6], root), 2.0e-10,
              "sign convention checked in local axes")


def simply_supported_udl_case() -> dict:
    E, G = 210000.0, 80769.0
    L = 4321.0
    w = 3.7
    sec = rect_section(120.0, 240.0)
    return {
        "materials": [dict(E=E, G=G, rho=7850.0)],
        "sections": [sec],
        "nodes": [
            dict(id=0, x=0.0, y=0.0, z=0.0, fix=[1, 1, 1, 1, 0, 0]),
            dict(id=1, x=L / 2.0, y=0.0, z=0.0, fix=[0, 0, 0, 0, 0, 0]),
            dict(id=2, x=L, y=0.0, z=0.0, fix=[0, 1, 1, 0, 0, 0]),
        ],
        "members": [
            dict(id=0, i=0, j=1, mat=0, sec=0, refvec=(0.0, 0.0, 1.0)),
            dict(id=1, i=1, j=2, mat=0, sec=0, refvec=(0.0, 0.0, 1.0)),
        ],
        "udls": [
            dict(member=0, w_local=(0.0, -w, 0.0)),
            dict(member=1, w_local=(0.0, -w, 0.0)),
        ],
    }


def check_simply_supported_udl(checks: list[Check]) -> None:
    model = simply_supported_udl_case()
    cli = run_cli(model)
    sec = model["sections"][0]
    mat = model["materials"][0]
    L = model["nodes"][2]["x"]
    w = -model["udls"][0]["w_local"][1]
    exact_defl = 5.0 * w * L**4 / (384.0 * mat["E"] * sec["Iz"])
    got = abs(float(cli["disp"][1][Uz]))
    add_check(checks, "closed form simply-supported UDL: midspan deflection",
              abs(got - exact_defl) / exact_defl, 2.0e-10, "limited by frame_cli %.12g output")
    exact_m = w * L**2 / 8.0
    got_m = abs(float(cli["mf"][0][11]))
    add_check(checks, "closed form simply-supported UDL: midspan moment",
              abs(got_m - exact_m) / exact_m, 2.0e-10, "limited by frame_cli %.12g output")


def skew_frame_case() -> dict:
    steel = dict(E=205000.0, G=79000.0, rho=7850.0)
    rc = dict(E=24870.0, G=10362.5, rho=2400.0)
    nodes = [
        dict(id=0, x=0.0, y=0.0, z=0.0, fix=[1, 1, 1, 1, 1, 1]),
        dict(id=1, x=4400.0, y=-250.0, z=0.0, fix=[1, 1, 1, 1, 1, 1]),
        dict(id=2, x=-300.0, y=3350.0, z=0.0, fix=[1, 1, 1, 1, 1, 1]),
        dict(id=3, x=120.0, y=80.0, z=3150.0, fix=[0, 0, 0, 0, 0, 0]),
        dict(id=4, x=4560.0, y=-130.0, z=3020.0, fix=[0, 0, 0, 0, 0, 0]),
        dict(id=5, x=-160.0, y=3480.0, z=3310.0, fix=[0, 0, 0, 0, 0, 0]),
        dict(id=6, x=2100.0, y=1710.0, z=5280.0, fix=[0, 0, 0, 0, 0, 0]),
    ]
    sections = [rect_section(220.0, 390.0), rect_section(160.0, 260.0), square_section(310.0)]
    raw_members = [
        (0, 0, 3, 1, 2),
        (1, 1, 4, 1, 2),
        (2, 2, 5, 1, 2),
        (3, 3, 4, 0, 0),
        (4, 4, 5, 0, 1),
        (5, 5, 3, 0, 1),
        (6, 3, 6, 0, 0),
        (7, 4, 6, 0, 1),
        (8, 5, 6, 0, 0),
    ]
    by_id = {n["id"]: n for n in nodes}
    members = []
    for mid, i, j, mat, sec in raw_members:
        pi = (by_id[i]["x"], by_id[i]["y"], by_id[i]["z"])
        pj = (by_id[j]["x"], by_id[j]["y"], by_id[j]["z"])
        members.append(dict(id=mid, i=i, j=j, mat=mat, sec=sec, refvec=refvec_for(pi, pj)))
    return {
        "materials": [steel, rc],
        "sections": sections,
        "nodes": nodes,
        "members": members,
        "nloads": [
            dict(node=3, comp=[18000.0, -7000.0, -62000.0, 0.0, 1.2e6, -0.4e6]),
            dict(node=4, comp=[-12000.0, 3000.0, -54000.0, -0.8e6, 0.0, 0.2e6]),
            dict(node=5, comp=[6000.0, 9000.0, -48000.0, 0.3e6, -0.1e6, 0.0]),
            dict(node=6, comp=[4000.0, -5000.0, -35000.0, 0.0, 0.0, 0.6e6]),
        ],
        "udls": [
            dict(member=3, w_local=(0.0, -1.25, 0.35)),
            dict(member=4, w_local=(0.0, -0.95, -0.20)),
            dict(member=7, w_local=(0.0, -0.55, 0.10)),
        ],
    }


def rotate_matrix() -> np.ndarray:
    axis = np.array([1.0, 2.0, 3.0], dtype=float)
    axis /= norm(axis)
    a = 0.73
    c, s = math.cos(a), math.sin(a)
    x, y, z = axis
    K = np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]], dtype=float)
    return c * np.eye(3) + s * K + (1.0 - c) * np.outer(axis, axis)


def rotated_cantilever_model(R: np.ndarray | None = None) -> dict:
    E, G = 205000.0, 79000.0
    L = 1987.0
    p1 = np.array([L, 0.0, 0.0])
    load = np.array([0.0, -850.0, 1300.0])
    ref = np.array([0.0, 0.0, 1.0])
    if R is not None:
        p1 = R @ p1
        load = R @ load
        ref = R @ ref
    return {
        "materials": [dict(E=E, G=G, rho=7850.0)],
        "sections": [rect_section(90.0, 210.0)],
        "nodes": [
            dict(id=0, x=0.0, y=0.0, z=0.0, fix=[1, 1, 1, 1, 1, 1]),
            dict(id=1, x=float(p1[0]), y=float(p1[1]), z=float(p1[2]), fix=[0, 0, 0, 0, 0, 0]),
        ],
        "members": [dict(id=0, i=0, j=1, mat=0, sec=0, refvec=tuple(float(v) for v in ref))],
        "nloads": [dict(node=1, comp=[float(load[0]), float(load[1]), float(load[2]), 0.0, 0.0, 0.0])],
    }


def check_rotation_equivariance(checks: list[Check]) -> None:
    R = rotate_matrix()
    base = run_cli(rotated_cantilever_model(None))
    rot = run_cli(rotated_cantilever_model(R))
    u0 = base["disp"][1][0:3]
    r0 = base["disp"][1][3:6]
    u1 = rot["disp"][1][0:3]
    r1 = rot["disp"][1][3:6]
    add_check(checks, "rotation equivariance: translations", max_scaled_diff(u1, R @ u0), 2.0e-10,
              "limited by frame_cli %.12g output")
    add_check(checks, "rotation equivariance: rotations", max_scaled_diff(r1, R @ r0), 2.0e-10,
              "limited by frame_cli %.12g output")


def timoshenko_deep_cantilever_case() -> dict:
    E, G = 210000.0, 80769.0
    L, P = 320.0, 1000.0
    sec = rect_section(100.0, 100.0)
    sec["Asy"] = sec["A"] * 5.0 / 6.0
    sec["Asz"] = sec["A"] * 5.0 / 6.0
    return {
        "materials": [dict(E=E, G=G, rho=7850.0)],
        "sections": [sec],
        "nodes": [
            dict(id=0, x=0.0, y=0.0, z=0.0, fix=[1, 1, 1, 1, 1, 1]),
            dict(id=1, x=L, y=0.0, z=0.0, fix=[0, 0, 0, 0, 0, 0]),
        ],
        "members": [dict(id=0, i=0, j=1, mat=0, sec=0, refvec=(0.0, 0.0, 1.0))],
        "nloads": [dict(node=1, comp=[0.0, 0.0, P, 0.0, 0.0, 0.0])],
        "opt": dict(useTimoshenko=True),
    }


def check_timoshenko_deep_beam(checks: list[Check]) -> None:
    model = timoshenko_deep_cantilever_case()
    cli = run_cli(model)
    py = py_solve(model)
    sec, mat = model["sections"][0], model["materials"][0]
    L = model["nodes"][1]["x"]
    P = model["nloads"][0]["comp"][2]
    exact = P * L**3 / (3.0 * mat["E"] * sec["Iz"]) + P * L / (mat["G"] * sec["Asy"])
    got = abs(float(cli["disp"][1][Uz]))
    add_check(checks, "Timoshenko deep cantilever: closed form deflection",
              abs(got - exact) / exact, 2.0e-10, "EB bending + shear term")
    add_check(checks, "Timoshenko deep cantilever: independent DSM displacement",
              max_scaled_diff(disp_array_cli(cli, 2), py["u"]), 2.0e-10)


def vertical_refvec_fallback_case() -> dict:
    E, G = 200000.0, 77000.0
    H = 2750.0
    return {
        "materials": [dict(E=E, G=G, rho=7850.0)],
        "sections": [rect_section(180.0, 260.0)],
        "nodes": [
            dict(id=0, x=0.0, y=0.0, z=0.0, fix=[1, 1, 1, 1, 1, 1]),
            dict(id=1, x=0.0, y=0.0, z=H, fix=[0, 0, 0, 0, 0, 0]),
        ],
        "members": [dict(id=0, i=0, j=1, mat=0, sec=0, refvec=(0.0, 0.0, 1.0))],
        "nloads": [dict(node=1, comp=[2400.0, -900.0, -12000.0, 0.0, 0.0, 150000.0])],
    }


def check_vertical_refvec_fallback(checks: list[Check]) -> None:
    model = vertical_refvec_fallback_case()
    compare_cli_to_python(checks, "vertical member with parallel refVec fallback", model, 2.0e-10, 2.0e-10)


def mechanism_case() -> dict:
    return {
        "materials": [dict(E=210000.0, G=80769.0, rho=7850.0)],
        "sections": [rect_section(100.0, 100.0)],
        "nodes": [
            dict(id=0, x=0.0, y=0.0, z=0.0, fix=[1, 1, 1, 0, 0, 0]),
            dict(id=1, x=1600.0, y=0.0, z=0.0, fix=[0, 0, 0, 0, 0, 0]),
        ],
        "members": [dict(id=0, i=0, j=1, mat=0, sec=0, refvec=(0.0, 0.0, 1.0))],
        "nloads": [dict(node=1, comp=[0.0, 0.0, 1000.0, 0.0, 0.0, 0.0])],
    }


def check_mechanism_detection(checks: list[Check]) -> None:
    model = mechanism_case()
    cli = run_cli(model)
    py = py_solve(model)
    add_check(checks, "mechanism detection: CLI singular flag", abs(float(cli["singular"] - 1)), 0.0)
    add_check(checks, "mechanism detection: independent DSM singular flag", abs(float(py["singular"] - 1)), 0.0)


def split_load_cases(model: dict) -> tuple[dict, dict]:
    a = clone_model(model)
    b = clone_model(model)
    a["nloads"] = [copy.deepcopy(model["nloads"][0]), copy.deepcopy(model["nloads"][2])]
    b["nloads"] = [copy.deepcopy(model["nloads"][1]), copy.deepcopy(model["nloads"][3])]
    a["udls"] = [copy.deepcopy(model["udls"][0]), copy.deepcopy(model["udls"][2])]
    b["udls"] = [copy.deepcopy(model["udls"][1])]
    return a, b


def check_linearity_and_scaling(checks: list[Check]) -> None:
    base_model = skew_frame_case()
    base = run_cli(base_model)
    a_model, b_model = split_load_cases(base_model)
    a = run_cli(a_model)
    b = run_cli(b_model)

    ids_n = sorted(int(n["id"]) for n in base_model["nodes"])
    ids_m = sorted(int(m["id"]) for m in base_model["members"])
    u_base = stacked_disp_by_id(base_model, base, ids_n)
    u_sum = stacked_disp_by_id(a_model, a, ids_n) + stacked_disp_by_id(b_model, b, ids_n)
    mf_base = stacked_mf_by_id(base_model, base, ids_m)
    mf_sum = stacked_mf_by_id(a_model, a, ids_m) + stacked_mf_by_id(b_model, b, ids_m)
    add_check(checks, "linearity: split load cases displacement superposition",
              max_scaled_diff(u_base, u_sum), 5.0e-10)
    add_check(checks, "linearity: split load cases member-force superposition",
              max_scaled_diff(mf_base, mf_sum), 5.0e-10)

    factor = -2.75
    scaled_model = scale_loads(base_model, factor)
    scaled = run_cli(scaled_model)
    add_check(checks, "linearity: negative load scaling displacement",
              max_scaled_diff(stacked_disp_by_id(scaled_model, scaled, ids_n), factor * u_base), 5.0e-10)
    add_check(checks, "linearity: negative load scaling member force",
              max_scaled_diff(stacked_mf_by_id(scaled_model, scaled, ids_m), factor * mf_base), 5.0e-10)


def check_maxwell_betti(checks: list[Check]) -> None:
    base = skew_frame_case()
    base["udls"] = []
    loads_a = [
        dict(node=6, comp=[12345.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
        dict(node=4, comp=[0.0, 0.0, 0.0, 0.0, 0.0, 2.5e6]),
    ]
    loads_b = [
        dict(node=5, comp=[0.0, -8000.0, 0.0, 0.0, 0.0, 0.0]),
        dict(node=3, comp=[0.0, 0.0, 0.0, 0.0, 1.1e6, 0.0]),
    ]
    case_a = clone_model(base)
    case_b = clone_model(base)
    case_a["nloads"] = loads_a
    case_b["nloads"] = loads_b
    out_a = run_cli(case_a)
    out_b = run_cli(case_b)
    wab = load_work(case_a, out_b, loads_a)
    wba = load_work(case_b, out_a, loads_b)
    scale = max(abs(wab), abs(wba), 1.0)
    add_check(checks, "Maxwell-Betti reciprocity: Fa^T ub == Fb^T ua",
              abs(wab - wba) / scale, 5.0e-10)


def check_order_invariance(checks: list[Check]) -> None:
    model = skew_frame_case()
    perm = clone_model(model)
    perm["nodes"] = [perm["nodes"][i] for i in [2, 0, 1, 5, 3, 6, 4]]
    perm["members"] = [perm["members"][i] for i in [4, 0, 8, 2, 6, 1, 5, 7, 3]]
    out0 = run_cli(model)
    out1 = run_cli(perm)
    ids_n = sorted(int(n["id"]) for n in model["nodes"])
    ids_m = sorted(int(m["id"]) for m in model["members"])
    add_check(checks, "input order invariance: node displacement by node id",
              max_scaled_diff(stacked_disp_by_id(model, out0, ids_n), stacked_disp_by_id(perm, out1, ids_n)), 5.0e-10)
    add_check(checks, "input order invariance: member forces by member id",
              max_scaled_diff(stacked_mf_by_id(model, out0, ids_m), stacked_mf_by_id(perm, out1, ids_m)), 5.0e-10)


def ill_conditioned_cantilever_case() -> dict:
    sec = rect_section(22.0, 44.0)
    return {
        "materials": [dict(E=205000.0, G=79000.0, rho=7850.0)],
        "sections": [sec],
        "nodes": [
            dict(id=0, x=0.0, y=0.0, z=0.0, fix=[1, 1, 1, 1, 1, 1]),
            dict(id=1, x=85000.0, y=0.0, z=0.0, fix=[0, 0, 0, 0, 0, 0]),
        ],
        "members": [dict(id=0, i=0, j=1, mat=0, sec=0, refvec=(0.0, 0.0, 1.0))],
        "nloads": [dict(node=1, comp=[120.0, -45.0, 60.0, 1000.0, 0.0, 0.0])],
    }


def check_ill_conditioned_case(checks: list[Check]) -> None:
    model = ill_conditioned_cantilever_case()
    py = py_solve(model)
    cli = run_cli(model)
    add_check(checks, "ill-conditioned slender cantilever: condition number lower bound",
              1.0e7 / max(py["cond"], 1.0), 1.0, f"cond(Kff)={py['cond']:.3e}")
    if not cli["singular"] and not py["singular"]:
        add_check(checks, "ill-conditioned slender cantilever: displacement vs DSM",
                  max_scaled_diff(disp_array_cli(cli, 2), py["u"]), 1.0e-8, f"cond(Kff)={py['cond']:.3e}")
        add_check(checks, "ill-conditioned slender cantilever: member force vs DSM",
                  max_scaled_diff(stacked_mf_by_id(model, cli), np.concatenate([py["mf"][0]])), 1.0e-8)


def print_report(checks: list[Check]) -> int:
    print("=" * 88)
    print("FrameCore independent precision audit")
    print("=" * 88)
    failures = 0
    for c in checks:
        status = "PASS" if c.ok else "FAIL"
        print(f"{status:4}  {c.name:<64} {c.value:.3e}  tol={c.tol:.1e}")
        if c.note:
            print(f"      note: {c.note}")
        if not c.ok:
            failures += 1
    print("-" * 88)
    print(f"checks={len(checks)} failures={failures}")
    print("CLI precision note: frame_cli prints %.12g, so black-box comparisons cannot prove")
    print("better than roughly 1e-11 to 1e-12 even if the internal solve is more accurate.")
    return 1 if failures else 0


def main() -> int:
    checks: list[Check] = []
    check_cantilever_closed_form(checks)
    check_simply_supported_udl(checks)
    compare_cli_to_python(checks, "closed-form cantilever model", cantilever_closed_form_case(), 2.0e-10, 2.0e-10)
    compare_cli_to_python(checks, "simply-supported UDL model", simply_supported_udl_case(), 2.0e-10, 2.0e-10)
    compare_cli_to_python(checks, "skew 3D frame with nodal loads and UDLs", skew_frame_case(), 5.0e-10, 5.0e-10)
    check_rotation_equivariance(checks)
    check_timoshenko_deep_beam(checks)
    check_vertical_refvec_fallback(checks)
    check_mechanism_detection(checks)
    check_linearity_and_scaling(checks)
    check_maxwell_betti(checks)
    check_order_invariance(checks)
    check_ill_conditioned_case(checks)
    return print_report(checks)


if __name__ == "__main__":
    raise SystemExit(main())
