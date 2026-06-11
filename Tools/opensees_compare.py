#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
#14 OpenSees offline cross-validation harness (NOT shipped in the game; licensing).

Solves the SAME reference structures three ways and reports the numeric differences:
  (1) our engine  -> frame_cli.exe (stdin/stdout driver around frame::solve)
  (2) OpenSees    -> openseespy (the independent reference solver)
  (3) analytic    -> closed-form oracle, where one exists

Both solvers are linear-elastic Euler-Bernoulli, so they must agree to a tight
tolerance. We compare GLOBAL node displacements (frame-invariant, so independent of
each solver's local-axis / Iy-Iz convention) plus the axial force and resultant end
moment (also convention-free). Reference models use SQUARE sections (Iy == Iz) so the
3D geomTransf vecxz / Iy-Iz ambiguity cannot bias the comparison.

Reference set includes the Taiwan RC code modulus E_c = 4700*sqrt(f'c) (MPa), matching
the 112-year concrete design code (ACI 318-19 basis).

Usage:  python Tools/opensees_compare.py        (exit 0 iff every diff is within tol)
"""
import os, sys, math, subprocess

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CLI  = os.path.join(ROOT, "Plugins", "FrameSolver", "Standalone", "frame_cli.exe")
BUILD_CLI = os.path.join(ROOT, "Plugins", "FrameSolver", "Standalone", "build_cli.bat")

try:
    import openseespy.opensees as ops
except Exception as e:  # pragma: no cover
    print("[FAIL] openseespy not importable:", e)
    print("       install with:  python -m pip install openseespy")
    sys.exit(2)


# ----------------------------------------------------------------- section helpers
def square_section(side):
    """A, Iy, Iz, J, cy, cz, Asy, Asz for a solid square (Iy == Iz; St-Venant J)."""
    A = side * side
    I = side ** 4 / 12.0
    J = 0.1406 * side ** 4            # St-Venant torsion constant of a square
    return dict(A=A, Iy=I, Iz=I, J=J, cy=side / 2, cz=side / 2, Asy=0.0, Asz=0.0)


def rect_section(b, d):
    """Rectangular section, b = width (local z), d = depth (local y): Iz = b d^3/12
    (x-y bending), Iy = d b^3/12 (x-z bending). Iy != Iz so it actually exercises the
    local-axis convention."""
    A = b * d
    Iz = b * d ** 3 / 12.0
    Iy = d * b ** 3 / 12.0
    bb, tt = max(b, d), min(b, d)            # St-Venant J of a rectangle (approx series)
    J = bb * tt ** 3 * (1.0 / 3.0 - 0.21 * (tt / bb) * (1.0 - (tt ** 4) / (12.0 * bb ** 4)))
    return dict(A=A, Iy=Iy, Iz=Iz, J=J, cy=d / 2, cz=b / 2, Asy=0.0, Asz=0.0)


def refvec_for(pi, pj):
    """(1,0,0) for a vertical member else (0,0,1) — dodges the local-axis degeneracy."""
    vertical = abs(pi[0] - pj[0]) < 1e-9 and abs(pi[1] - pj[1]) < 1e-9
    return (1.0, 0.0, 0.0) if vertical else (0.0, 0.0, 1.0)


# ----------------------------------------------------------------- our engine (frame_cli)
def run_frame_cli(model):
    lines = []
    for m in model["materials"]:
        lines.append(f"MAT {m['E']} {m['G']} {m.get('rho', 7850.0)}")
    for s in model["sections"]:
        lines.append("SEC {A} {Iy} {Iz} {J} {cy} {cz} {Asy} {Asz}".format(**s))
    for n in model["nodes"]:
        fp = list(n["fix"]) + list(n.get("presc", [0, 0, 0, 0, 0, 0]))
        lines.append("NODE {} {} {} {} {}".format(n["id"], n["x"], n["y"], n["z"],
                                                  " ".join(str(v) for v in fp)))
    for e in model["members"]:
        rv = e["refvec"]
        lines.append(f"MEMBER {e['id']} {e['i']} {e['j']} {e['mat']} {e['sec']} {rv[0]} {rv[1]} {rv[2]} {e.get('active', 1)}")
    for l in model.get("nloads", []):
        c = l["comp"]
        lines.append(f"NLOAD {l['node']} {c[0]} {c[1]} {c[2]} {c[3]} {c[4]} {c[5]}")
    for h in model.get("hinges", []):
        lines.append(f"HINGE {h['member']} {h['dof']} {h['Mp']}")
    lines.append("END")
    text = "\n".join(lines) + "\n"

    p = subprocess.run([CLI], input=text, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError("frame_cli failed: " + p.stderr)
    disp, mf, singular = {}, {}, None
    for ln in p.stdout.splitlines():
        t = ln.split()
        if not t:
            continue
        if t[0] == "SINGULAR":
            singular = int(t[1])
        elif t[0] == "DISP":
            disp[int(t[1])] = [float(x) for x in t[2:8]]
        elif t[0] == "MF":
            mf[int(t[1])] = [float(x) for x in t[2:14]]
    return singular, disp, mf


# ----------------------------------------------------------------- OpenSees
def run_opensees(model):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    free_ry = set(model.get("os_free_ry", []))   # hinge-state twin: support Ry left free
    for n in model["nodes"]:
        ops.node(n["id"], float(n["x"]), float(n["y"]), float(n["z"]))
        f = n["fix"]; p = n.get("presc", [0, 0, 0, 0, 0, 0])
        # A prescribed (nonzero) DOF is imposed via sp() below, NOT fixed to 0 here.
        mask = [int(f[d]) if (f[d] and p[d] == 0.0) else 0 for d in range(6)]
        if n["id"] in free_ry:
            mask[4] = 0                          # global Ry freed (the hinge rotation)
        ops.fix(n["id"], *mask)

    for e in model["members"]:
        if not e.get("active", 1):
            continue   # our side deactivates the member; OpenSees simply never builds it
        pi = next(nn for nn in model["nodes"] if nn["id"] == e["i"])
        pj = next(nn for nn in model["nodes"] if nn["id"] == e["j"])
        axis = (pj["x"] - pi["x"], pj["y"] - pi["y"], pj["z"] - pi["z"])
        # vecxz: any vector not parallel to the member axis (mirror our refVec rule).
        vert = abs(axis[0]) < 1e-9 and abs(axis[1]) < 1e-9
        vecxz = (1.0, 0.0, 0.0) if vert else (0.0, 0.0, 1.0)
        ops.geomTransf("Linear", e["id"] + 1, *vecxz)
        m = model["materials"][e["mat"]]
        s = model["sections"][e["sec"]]
        # CONVENTION BRIDGE: our refVec defines the local x-y plane (local y toward refVec),
        # whereas OpenSees geomTransf vecxz gives local y = vecxz x axis (perpendicular to
        # the x-vecxz plane). Using the SAME vector for both therefore swaps the roles of
        # the two bending axes, so our Iz (x-y bending) is OpenSees' Iy and vice versa.
        # We pass (Iy_os, Iz_os) = (our Iz, our Iy). For square sections this is a no-op;
        # the rectangular cantilever below positively verifies the mapping is correct.
        # element elasticBeamColumn $tag $i $j $A $E $G $Jx $Iy $Iz $transfTag
        ops.element("elasticBeamColumn", e["id"] + 1, e["i"], e["j"],
                    s["A"], m["E"], m["G"], s["J"], s["Iz"], s["Iy"], e["id"] + 1)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for l in model.get("nloads", []):
        c = l["comp"]
        ops.load(l["node"], c[0], c[1], c[2], c[3], c[4], c[5])
    # hinge-state twin: the constant plastic moment acting at the freed support rotation
    for hm in model.get("os_hinge_moments", []):
        ops.load(hm["node"], 0.0, 0.0, 0.0, 0.0, float(hm["My"]), 0.0)
    # prescribed (imposed) support displacements -> single-point constraints in the pattern
    for n in model["nodes"]:
        f = n["fix"]; p = n.get("presc", [0, 0, 0, 0, 0, 0])
        for d in range(6):
            if f[d] and p[d] != 0.0:
                ops.sp(n["id"], d + 1, float(p[d]))

    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1.0e-10, 20)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    ok = ops.analyze(1)

    disp = {}
    for n in model["nodes"]:
        disp[n["id"]] = [ops.nodeDisp(n["id"], k) for k in range(1, 7)]
    mf = {}
    for e in model["members"]:
        # localForce (NOT eleForce, which is GLOBAL) -> axial along local x + local moments,
        # directly comparable to our compression-positive local end forces (abs / resultant).
        # An inactive member has no OpenSees element; it reads zero, which is exactly our
        # removed-member force convention -- so the comparison also gates that convention.
        mf[e["id"]] = ops.eleResponse(e["id"] + 1, "localForce") if e.get("active", 1) else [0.0] * 12
    return ok, disp, mf


# ----------------------------------------------------------------- comparison
def max_disp_norm(disp):
    return max((max(abs(v) for v in d) for d in disp.values()), default=1.0) or 1.0


def compare_disp(d_mine, d_os):
    scale = max(max_disp_norm(d_mine), max_disp_norm(d_os), 1e-30)
    worst = 0.0
    for nid in d_mine:
        for k in range(6):
            worst = max(worst, abs(d_mine[nid][k] - d_os[nid][k]) / scale)
    return worst


def resultant_moment(mf12, end):  # end 0 -> i (My,Mz at idx 4,5), end 1 -> j (idx 10,11)
    base = 0 if end == 0 else 6
    return math.hypot(mf12[base + 4], mf12[base + 5])


def compare_forces(mf_mine, mf_os):
    worstN, worstM = 0.0, 0.0
    Nscale = max((max(abs(mf_mine[e][0]), abs(mf_mine[e][6])) for e in mf_mine), default=1.0) or 1.0
    Mscale = max((max(resultant_moment(mf_mine[e], 0), resultant_moment(mf_mine[e], 1)) for e in mf_mine), default=1.0) or 1.0
    for e in mf_mine:
        worstN = max(worstN, abs(abs(mf_mine[e][0]) - abs(mf_os[e][0])) / Nscale)
        for end in (0, 1):
            worstM = max(worstM, abs(resultant_moment(mf_mine[e], end) - resultant_moment(mf_os[e], end)) / Mscale)
    return worstN, worstM


# ----------------------------------------------------------------- reference models
def model_cantilever3d():
    sec = square_section(100.0)
    nodes = [
        dict(id=0, x=0.0, y=0.0, z=0.0, fix=[1, 1, 1, 1, 1, 1]),
        dict(id=1, x=2000.0, y=0.0, z=0.0, fix=[0, 0, 0, 0, 0, 0]),
    ]
    members = [dict(id=0, i=0, j=1, mat=0, sec=0, refvec=refvec_for((0, 0, 0), (2000, 0, 0)))]
    P = 1000.0
    nloads = [dict(node=1, comp=[0, 0, -P, 0, 0, 0])]
    mats = [dict(E=210000.0, G=80769.0, rho=7850.0)]
    analytic = dict(node=1, dof=2, value=-P * 2000.0 ** 3 / (3.0 * 210000.0 * sec["Iz"]))
    return dict(name="cantilever3d (steel, tip load)", materials=mats, sections=[sec],
                nodes=nodes, members=members, nloads=nloads, analytic=analytic)


def model_cantilever_rect():
    # rectangular section (Iy != Iz) loaded transversely in BOTH global Y and Z, so the
    # two distinct inertias are both exercised -> positively verifies the refVec / Iy-Iz
    # convention bridge (not just masked by a square section).
    sec = rect_section(100.0, 200.0)          # b=100 (z), d=200 (y)
    L, E = 2000.0, 210000.0
    nodes = [
        dict(id=0, x=0.0, y=0.0, z=0.0, fix=[1, 1, 1, 1, 1, 1]),
        dict(id=1, x=L,   y=0.0, z=0.0, fix=[0, 0, 0, 0, 0, 0]),
    ]
    members = [dict(id=0, i=0, j=1, mat=0, sec=0, refvec=refvec_for((0, 0, 0), (L, 0, 0)))]
    Py, Pz = -800.0, -1000.0
    nloads = [dict(node=1, comp=[0, Py, Pz, 0, 0, 0])]
    mats = [dict(E=E, G=80769.0, rho=7850.0)]
    # global-Z deflection uses our Iz (x-y bending); checked analytically.
    analytic = dict(node=1, dof=2, value=Pz * L ** 3 / (3.0 * E * sec["Iz"]))
    return dict(name="cantilever rect (Iy!=Iz, biaxial tip load)", materials=mats, sections=[sec],
                nodes=nodes, members=members, nloads=nloads, analytic=analytic)


def model_portal_rc():
    fc = 28.0
    Ec = 4700.0 * math.sqrt(fc)              # Taiwan 112 / ACI 318-19, MPa  (~24870)
    Gc = Ec / (2.0 * (1.0 + 0.2))            # concrete nu ~ 0.2
    sec = square_section(500.0)
    H, L = 3500.0, 6000.0
    nodes = [
        dict(id=0, x=0.0, y=0.0, z=0.0, fix=[1, 1, 1, 1, 1, 1]),
        dict(id=1, x=L,   y=0.0, z=0.0, fix=[1, 1, 1, 1, 1, 1]),
        dict(id=2, x=0.0, y=0.0, z=H,   fix=[0, 0, 0, 0, 0, 0]),
        dict(id=3, x=L,   y=0.0, z=H,   fix=[0, 0, 0, 0, 0, 0]),
    ]
    members = [
        dict(id=0, i=0, j=2, mat=0, sec=0, refvec=refvec_for((0, 0, 0), (0, 0, H))),
        dict(id=1, i=1, j=3, mat=0, sec=0, refvec=refvec_for((L, 0, 0), (L, 0, H))),
        dict(id=2, i=2, j=3, mat=0, sec=0, refvec=refvec_for((0, 0, H), (L, 0, H))),
    ]
    nloads = [
        dict(node=2, comp=[50000.0, 0, -100000.0, 0, 0, 0]),   # lateral + gravity
        dict(node=3, comp=[0, 0, -100000.0, 0, 0, 0]),
    ]
    mats = [dict(E=Ec, G=Gc, rho=2400.0)]
    return dict(name=f"portal RC (f'c={fc}, Ec={Ec:.0f} MPa)", materials=mats, sections=[sec],
                nodes=nodes, members=members, nloads=nloads, analytic=None)


def model_settlement():
    # Fixed-fixed beam whose far end settles by delta (Uz); free midspan node so the reduced
    # system is non-empty. Cross-checks our prescribed-displacement path against OpenSees sp().
    # Analytic: end moment 6*E*I*delta/L^2, reaction 12*E*I*delta/L^3.
    sec = square_section(100.0)
    L, delta = 2000.0, 1.0
    nodes = [
        dict(id=0, x=0.0,     y=0.0, z=0.0, fix=[1, 1, 1, 1, 1, 1]),
        dict(id=1, x=L / 2.0, y=0.0, z=0.0, fix=[0, 0, 0, 0, 0, 0]),
        dict(id=2, x=L,       y=0.0, z=0.0, fix=[1, 1, 1, 1, 1, 1], presc=[0, 0, -delta, 0, 0, 0]),
    ]
    members = [dict(id=0, i=0, j=1, mat=0, sec=0, refvec=(0, 0, 1)),
               dict(id=1, i=1, j=2, mat=0, sec=0, refvec=(0, 0, 1))]
    mats = [dict(E=210000.0, G=80769.0, rho=7850.0)]
    return dict(name="prescribed settlement (fixed-fixed, end settles)", materials=mats,
                sections=[sec], nodes=nodes, members=members, nloads=[], analytic=None)


def model_hinge_states():
    # Collapse stage 4a cross-check: a formed-hinge STATE solved by an independent solver.
    # OUR side: HINGE token (release + residual Mp baked into the element condensation).
    # OPENSEES side: there is no member-end release on elasticBeamColumn, so the same physics
    # is modelled the textbook way -- free the support's Ry (the hinge rotation about global
    # -y is our local z; the beam runs along +x with refVec (0,0,1)) and apply the constant
    # hinge moment externally at that node. The interior node response must then match.
    # Mp is read from OUR elastic run (50% of the current end moment -> a valid mid-history
    # state); the support-node Ry differs BY CONSTRUCTION (ours: condensed inside the element,
    # support stays fixed; OpenSees: a free dof reading the hinge rotation), so only the
    # interior node is compared, and member forces are skipped (our hinged end recovers 0 by
    # the condensation contract while OpenSees' element carries Mp there).
    sec = square_section(100.0)
    L, P = 4000.0, 40000.0
    nodes = [
        dict(id=0, x=0.0,     y=0.0, z=0.0, fix=[1, 1, 1, 1, 1, 1]),
        dict(id=1, x=L / 2.0, y=0.0, z=0.0, fix=[0, 0, 0, 0, 0, 0]),
        dict(id=2, x=L,       y=0.0, z=0.0, fix=[1, 1, 1, 1, 1, 1]),
    ]
    members = [dict(id=0, i=0, j=1, mat=0, sec=0, refvec=(0, 0, 1)),
               dict(id=1, i=1, j=2, mat=0, sec=0, refvec=(0, 0, 1))]
    nloads = [dict(node=1, comp=[0, 0, -P, 0, 0, 0])]
    mats = [dict(E=210000.0, G=80769.0, rho=7850.0)]
    elastic = dict(name="hinge probe (elastic)", materials=mats, sections=[sec],
                   nodes=nodes, members=members, nloads=nloads)
    sing, _d, mf = run_frame_cli(elastic)
    if sing:
        raise RuntimeError("hinge probe solve unexpectedly singular")
    mp0 = 0.5 * mf[0][5]            # 50% of the elastic endI.Mz of member 0 (signed, local)
    # Our local z = (0,-1,0) here, so an internal local moment Mz = Mp is a global moment
    # -Mp about +y; the freed OpenSees support balances it with an external My = -Mp.
    hinged = dict(name="hinge state: support hinge with residual Mp (vs OpenSees freed Ry + moment)",
                  materials=mats, sections=[sec], nodes=nodes, members=members, nloads=nloads,
                  hinges=[dict(member=0, dof=5, Mp=mp0)],
                  os_free_ry=[0], os_hinge_moments=[dict(node=0, My=-mp0)],
                  compare_nodes=[1], skip_forces=True, analytic=None)
    return [hinged]


def model_collapse_states():
    # Collapse stage 3c cross-check: the driver's per-step states are ordinary linear solves of
    # a partially-removed topology, so each state must match OpenSees built WITHOUT the removed
    # members. State A = full propped cantilever (indeterminate); state B = the same model with
    # the prop DEACTIVATED on our side (active=0 token) vs simply omitted on the OpenSees side
    # (its foot node remains, fully fixed and element-less). The removal SEQUENCE itself is
    # gated by the F30 closed-form oracles; this leg gates the per-step linear states.
    sec = square_section(100.0)
    L, h, P = 3000.0, 1000.0, 40000.0
    nodes = [
        dict(id=0, x=0.0,     y=0.0, z=0.0, fix=[1, 1, 1, 1, 1, 1]),
        dict(id=1, x=L / 2.0, y=0.0, z=0.0, fix=[0, 0, 0, 0, 0, 0]),
        dict(id=2, x=L,       y=0.0, z=0.0, fix=[0, 0, 0, 0, 0, 0]),
        dict(id=3, x=L,       y=0.0, z=-h,  fix=[1, 1, 1, 1, 1, 1]),
    ]
    def members(prop_active):
        return [dict(id=0, i=0, j=1, mat=0, sec=0, refvec=(0, 0, 1)),
                dict(id=1, i=1, j=2, mat=0, sec=0, refvec=(0, 0, 1)),
                dict(id=2, i=3, j=2, mat=0, sec=0, refvec=refvec_for((L, 0, -h), (L, 0, 0)),
                     active=prop_active)]
    nloads = [dict(node=1, comp=[0, 0, -P, 0, 0, 0])]
    mats = [dict(E=210000.0, G=80769.0, rho=7850.0)]
    full = dict(name="collapse state A: propped cantilever (full)", materials=mats,
                sections=[sec], nodes=nodes, members=members(1), nloads=nloads, analytic=None)
    cut = dict(name="collapse state B: prop removed by flag (vs OpenSees omission)", materials=mats,
               sections=[sec], nodes=nodes, members=members(0), nloads=nloads, analytic=None)
    return [full, cut]


# ----------------------------------------------------------------- shells (MITC4)
def run_frame_cli_shell(model, dkq=False):
    lines = []
    if dkq:
        lines.append("OPT 0 0 1e-12 0 1")   # S8-8b: useDKQPlate = 1 (releases/timoshenko/pivot default)
    for m in model["smats"]:
        lines.append(f"SMAT {m['E']} {m['nu']} {m['G']}")
    for n in model["nodes"]:
        f = n["fix"]
        lines.append("NODE {id} {x} {y} {z} {0} {1} {2} {3} {4} {5}".format(*f, **n))
    for s in model["shells"]:
        nn = s["n"]
        lines.append(f"SHELL {s['id']} {nn[0]} {nn[1]} {nn[2]} {nn[3]} {s['mat']} {s['t']} {s.get('active', 1)}")
    for l in model.get("nloads", []):
        c = l["comp"]
        lines.append(f"NLOAD {l['node']} {c[0]} {c[1]} {c[2]} {c[3]} {c[4]} {c[5]}")
    lines.append("END")
    p = subprocess.run([CLI], input="\n".join(lines) + "\n", capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError("frame_cli failed: " + p.stderr)
    disp, singular = {}, None
    for ln in p.stdout.splitlines():
        t = ln.split()
        if not t:
            continue
        if t[0] == "SINGULAR":
            singular = int(t[1])
        elif t[0] == "DISP":
            disp[int(t[1])] = [float(x) for x in t[2:8]]
    return singular, disp


def run_opensees_shell(model, dkq=False):
    # ShellMITC4 (default) / ShellDKGQ (dkq=True) + ElasticMembranePlateSection. ShellMITC4 IS the
    # Reissner-Mindlin MITC4 we implement; ShellDKGQ shares the Batoz DKQ thin-plate bending with
    # our opt-in DKQ path. ShellDKGQ's GQ12 membrane differs from our Q4, so DKQ legs use FLAT
    # plates (pure bending) where the in-plane membrane is inactive.
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    for n in model["nodes"]:
        ops.node(n["id"], float(n["x"]), float(n["y"]), float(n["z"]))
        ops.fix(n["id"], *[int(x) for x in n["fix"]])
    for s in model["shells"]:
        if not s.get("active", 1):
            continue   # our side deactivates the facet; OpenSees simply never builds it
        m = model["smats"][s["mat"]]
        secTag = s["id"] + 1
        ops.section("ElasticMembranePlateSection", secTag, m["E"], m["nu"], s["t"], 0.0)
        nn = s["n"]
        ops.element("ShellDKGQ" if dkq else "ShellMITC4", s["id"] + 1, nn[0], nn[1], nn[2], nn[3], secTag)
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
    disp = {n["id"]: [ops.nodeDisp(n["id"], k) for k in range(1, 7)] for n in model["nodes"]}
    return ok, disp


def shell_plate_model(name, nx, ny, Lx, Ly, t, E, nu, tilt_deg, Fz_total):
    """Cantilever plate clamped along the x=0 edge, point loads (global Fz_total split
    over the free-edge nodes). `tilt_deg` rigidly rotates the whole plate about global Y,
    so a non-zero tilt exercises the 3D facet rotation AND couples membrane+bending — the
    same global load is fed to both solvers, so the comparison stays valid."""
    G = E / (2.0 * (1.0 + nu))
    th = math.radians(tilt_deg)
    c, s = math.cos(th), math.sin(th)
    idx = lambda i, j: j * (nx + 1) + i
    nodes = []
    for j in range(ny + 1):
        for i in range(nx + 1):
            x0, y0, z0 = Lx * i / nx, Ly * j / ny, 0.0
            x = c * x0 + s * z0
            z = -s * x0 + c * z0
            fix = [1, 1, 1, 1, 1, 1] if i == 0 else [0, 0, 0, 0, 0, 0]
            nodes.append(dict(id=idx(i, j), x=x, y=y0, z=z, fix=fix))
    shells = []
    sid = 0
    for j in range(ny):
        for i in range(nx):
            shells.append(dict(id=sid, n=[idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)],
                               mat=0, t=t))
            sid += 1
    nedge = ny + 1
    nloads = [dict(node=idx(nx, j), comp=[0, 0, Fz_total / nedge, 0, 0, 0]) for j in range(ny + 1)]
    return dict(name=name, smats=[dict(E=E, nu=nu, G=G)], nodes=nodes, shells=shells, nloads=nloads)


def shell_models():
    # Element-removal cross-check (collapse stage 3a): OUR side keeps the facet in the model but
    # deactivated (ShellQuad::active=false); the OpenSees side simply never builds that element.
    # The two must agree like any other plate — removal-by-flag == removal-by-omission. The
    # deactivated facet is INTERIOR, so every node keeps at least one active facet (no mechanism).
    holed = shell_plate_model("shell cantilever plate (interior facet deactivated)", 8, 4,
                              1000.0, 500.0, 10.0, 30000.0, 0.3, 0.0, -100.0)
    holed["shells"][11]["active"] = 0   # quad (i=3, j=1) — interior, all 4 nodes shared
    return [
        shell_plate_model("shell cantilever plate (flat)", 8, 4, 1000.0, 500.0, 10.0,
                          30000.0, 0.3, 0.0, -100.0),
        shell_plate_model("shell cantilever plate (tilted 30deg, 3D + membrane)", 8, 4,
                          1000.0, 500.0, 10.0, 30000.0, 0.3, 30.0, -100.0),
        holed,
    ]


# ----------------------------------------------------------------- modal (eigen)
def run_frame_cli_modal(model, nModes):
    lines = []
    for m in model["materials"]:
        lines.append(f"MAT {m['E']} {m['G']} {m.get('rho', 7850.0)}")
    for s in model["sections"]:
        lines.append("SEC {A} {Iy} {Iz} {J} {cy} {cz} {Asy} {Asz}".format(**s))
    for n in model["nodes"]:
        fp = list(n["fix"]) + list(n.get("presc", [0, 0, 0, 0, 0, 0]))
        lines.append("NODE {} {} {} {} {}".format(n["id"], n["x"], n["y"], n["z"], " ".join(str(v) for v in fp)))
    for e in model["members"]:
        rv = e["refvec"]
        lines.append(f"MEMBER {e['id']} {e['i']} {e['j']} {e['mat']} {e['sec']} {rv[0]} {rv[1]} {rv[2]}")
    lines.append(f"EIGEN {nModes}")
    lines.append("END")
    p = subprocess.run([CLI], input="\n".join(lines) + "\n", capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError("frame_cli modal failed: " + p.stderr)
    for ln in p.stdout.splitlines():
        t = ln.split()
        if t and t[0] == "FREQ":
            return sorted(float(x) for x in t[2:2 + int(t[1])])
    return []


def run_opensees_modal(model, nModes):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    for n in model["nodes"]:
        ops.node(n["id"], float(n["x"]), float(n["y"]), float(n["z"]))
        ops.fix(n["id"], *[int(x) for x in n["fix"]])
    for e in model["members"]:
        pi = next(nn for nn in model["nodes"] if nn["id"] == e["i"])
        pj = next(nn for nn in model["nodes"] if nn["id"] == e["j"])
        vert = abs(pj["x"] - pi["x"]) < 1e-9 and abs(pj["y"] - pi["y"]) < 1e-9
        vecxz = (1.0, 0.0, 0.0) if vert else (0.0, 0.0, 1.0)
        ops.geomTransf("Linear", e["id"] + 1, *vecxz)
        m = model["materials"][e["mat"]]; s = model["sections"][e["sec"]]
        massPerLen = m.get("rho", 0.0) * 1e-12 * s["A"]      # consistent units (tonne/mm^3 * mm^2)
        # consistent (-cMass) mass to match our localMass12; Iy/Iz swapped per the convention bridge.
        ops.element("elasticBeamColumn", e["id"] + 1, e["i"], e["j"],
                    s["A"], m["E"], m["G"], s["J"], s["Iz"], s["Iy"], e["id"] + 1,
                    "-mass", massPerLen, "-cMass")
    eigs = ops.eigen(nModes)
    return sorted(math.sqrt(abs(l)) for l in eigs)


def modal_models():
    def beam(name, n, L, supports):
        sec = square_section(100.0)
        nodes = []
        for i in range(n + 1):
            fix = supports(i, n)
            nodes.append(dict(id=i, x=L * i / n, y=0.0, z=0.0, fix=fix))
        members = [dict(id=i, i=i, j=i + 1, mat=0, sec=0, refvec=(0, 0, 1)) for i in range(n)]
        mats = [dict(E=210000.0, G=80769.0, rho=7850.0)]
        return dict(name=name, materials=mats, sections=[sec], nodes=nodes, members=members, nModes=4)
    cant = beam("modal cantilever", 12, 3000.0,
                lambda i, n: [1, 1, 1, 1, 1, 1] if i == 0 else [0, 0, 0, 0, 0, 0])
    ss = beam("modal simply-supported", 12, 4000.0,
              lambda i, n: [1, 1, 1, 1, 0, 0] if i == 0 else ([0, 1, 1, 0, 0, 0] if i == n else [0, 0, 0, 0, 0, 0]))
    return [cant, ss]


# ----------------------------------------------------------------- P-Delta (second order)
def pdelta_column_model(frac):
    """Cantilever column along +Z (8 elements), base encastre, tip carries axial compression P
    (global -Z) + lateral H (global +X). Mirrors the standalone F40 fixture. The beam-column
    stability-function closed form gives the exact second-order tip sway."""
    side, E = 200.0, 200000.0
    G = E / (2.0 * (1.0 + 0.3))
    L, H = 6000.0, 1000.0
    I = side ** 4 / 12.0
    Pcr = math.pi ** 2 * E * I / (4.0 * L * L)
    P = frac * Pcr
    nE = 8
    sec = square_section(side)
    nodes = [dict(id=k, x=0.0, y=0.0, z=L * k / nE,
                  fix=[1, 1, 1, 1, 1, 1] if k == 0 else [0, 0, 0, 0, 0, 0]) for k in range(nE + 1)]
    members = [dict(id=k, i=k, j=k + 1, mat=0, sec=0, refvec=(1.0, 0.0, 0.0)) for k in range(nE)]
    nloads = [dict(node=nE, comp=[H, 0.0, -P, 0, 0, 0])]
    mats = [dict(E=E, G=G, rho=7850.0)]
    k = math.sqrt(P / (E * I))
    exact = H * (math.tan(k * L) - k * L) / (P * k)
    return dict(name=f"P-Delta column (P/Pcr={frac})", materials=mats, sections=[sec],
                nodes=nodes, members=members, nloads=nloads, tip_node=nE, exact=exact)


def run_frame_cli_pdelta(model, path):
    """Drive frame_cli in P-Delta mode (path 0=frozen reuse, 1=K_T reference). Returns the
    (converged, diverged, iters) status tuple and the second-order node displacements."""
    lines = []
    for m in model["materials"]:
        lines.append(f"MAT {m['E']} {m['G']} {m.get('rho', 7850.0)}")
    for s in model["sections"]:
        lines.append("SEC {A} {Iy} {Iz} {J} {cy} {cz} {Asy} {Asz}".format(**s))
    for n in model["nodes"]:
        fp = list(n["fix"]) + list(n.get("presc", [0, 0, 0, 0, 0, 0]))
        lines.append("NODE {} {} {} {} {}".format(n["id"], n["x"], n["y"], n["z"],
                                                  " ".join(str(v) for v in fp)))
    for e in model["members"]:
        rv = e["refvec"]
        lines.append(f"MEMBER {e['id']} {e['i']} {e['j']} {e['mat']} {e['sec']} {rv[0]} {rv[1]} {rv[2]}")
    for l in model.get("nloads", []):
        c = l["comp"]
        lines.append(f"NLOAD {l['node']} {c[0]} {c[1]} {c[2]} {c[3]} {c[4]} {c[5]}")
    lines.append(f"PDELTA {path}")
    lines.append("END")
    p = subprocess.run([CLI], input="\n".join(lines) + "\n", capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError("frame_cli pdelta failed: " + p.stderr)
    disp, status = {}, None
    for ln in p.stdout.splitlines():
        t = ln.split()
        if not t:
            continue
        if t[0] == "PDSTATUS":
            status = (int(t[1]), int(t[2]), int(t[3]))
        elif t[0] == "DISP":
            disp[int(t[1])] = [float(x) for x in t[2:8]]
    return status, disp


def run_opensees_pdelta(model):
    """OpenSees second-order solve: PDelta geomTransf + Newton (axial force updated to the
    converged state). This is the P-large-Delta linearization — its closed-form target is the
    SAME beam-column exact value, reached from a slightly different formulation than ours."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    for n in model["nodes"]:
        ops.node(n["id"], float(n["x"]), float(n["y"]), float(n["z"]))
        ops.fix(n["id"], *[int(x) for x in n["fix"]])
    for e in model["members"]:
        pi = next(nn for nn in model["nodes"] if nn["id"] == e["i"])
        pj = next(nn for nn in model["nodes"] if nn["id"] == e["j"])
        vert = abs(pj["x"] - pi["x"]) < 1e-9 and abs(pj["y"] - pi["y"]) < 1e-9
        vecxz = (1.0, 0.0, 0.0) if vert else (0.0, 0.0, 1.0)
        ops.geomTransf("PDelta", e["id"] + 1, *vecxz)
        m = model["materials"][e["mat"]]
        s = model["sections"][e["sec"]]
        ops.element("elasticBeamColumn", e["id"] + 1, e["i"], e["j"],
                    s["A"], m["E"], m["G"], s["J"], s["Iz"], s["Iy"], e["id"] + 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for l in model.get("nloads", []):
        ops.load(l["node"], *[float(x) for x in l["comp"]])
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1.0e-12, 100)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    ok = ops.analyze(1)
    disp = {n["id"]: [ops.nodeDisp(n["id"], k) for k in range(1, 7)] for n in model["nodes"]}
    return ok, disp


# ----------------------------------------------------------------- driver
def main():
    if not os.path.exists(CLI):
        print("[build] frame_cli.exe missing -> building...")
        rc = subprocess.run(["cmd", "/c", BUILD_CLI]).returncode
        if rc != 0 or not os.path.exists(CLI):
            print("[FAIL] could not build frame_cli.exe")
            return 2

    # Gate tolerances are STRICT by default (the beam models match OpenSees to ~1e-12, and the
    # MITC4 shell to ~1e-10 on the flat/tilted plates exercised here — see TOL_SHELL_VS_OS below;
    # skewed/warped meshes are ~1e-7-1e-8, cross-checked in shell_mitc4_deep_audit.py — so a
    # loose gate would not catch a precision regression). Pass --relaxed for a cross-platform
    # report run where floating-point/BLAS differences make sub-1e-8 agreement unrealistic.
    relaxed = ("--relaxed" in sys.argv)
    if relaxed:
        TOL_DISP_VS_OS, TOL_FORCE_VS_OS, TOL_VS_ANALYTIC = 1e-4, 1e-3, 2e-3
        TOL_SHELL_VS_OS = 1e-3
        TOL_MODAL = 1e-2
        TOL_PDELTA = 2e-2
    else:
        TOL_DISP_VS_OS, TOL_FORCE_VS_OS, TOL_VS_ANALYTIC = 1e-8, 1e-8, 1e-6
        # Our flat-shell and OpenSees ShellMITC4 are the same element and agree on node
        # displacements to ~1e-10 (membrane + MITC4 bending + drilling + 3D rotation all
        # match). The 1e-7 gate leaves headroom for FP/OpenSees-version variation while
        # still catching any real regression.
        TOL_SHELL_VS_OS = 1e-7
        # Our consistent mass vs OpenSees elasticBeamColumn -cMass: same EB consistent mass,
        # natural frequencies agree closely.
        TOL_MODAL = 1e-4
        # P-Delta: OpenSees PDelta geomTransf (P-large-Delta linearization, axial updated by Newton)
        # vs our consistent geometric stiffness (axial frozen at first order). Same closed-form
        # target but different formulation -> agreement is discretization-level, not a precision gate.
        # The TIGHT P-Delta oracle is the analytic beam-column F40 (frozen==reference 1e-12,
        # reference==exact ~1e-5); this leg is a cross-tool consistency check only.
        TOL_PDELTA = 1e-2
    print(f"  tolerances: {'RELAXED (report)' if relaxed else 'STRICT (gate)'}  "
          f"disp={TOL_DISP_VS_OS:.0e} force={TOL_FORCE_VS_OS:.0e} analytic={TOL_VS_ANALYTIC:.0e}")

    models = [model_cantilever3d(), model_cantilever_rect(), model_portal_rc(), model_settlement()] \
             + model_collapse_states() + model_hinge_states()
    failures = 0
    print("=" * 64)
    print(" #14 OpenSees offline cross-validation")
    print("=" * 64)
    for M in models:
        sing, d_mine, mf_mine = run_frame_cli(M)
        ok_os, d_os, mf_os = run_opensees(M)
        print(f"\n[{M['name']}]")
        if sing or ok_os != 0:
            print(f"  [FAIL] solve flagged singular (ours={sing}, openseesAnalyze={ok_os})")
            failures += 1
            continue
        cn = M.get("compare_nodes")
        if cn:   # hinge-state models: the support's hinge rotation lives in different places
            d_mine = {k: v for k, v in d_mine.items() if k in cn}
            d_os   = {k: v for k, v in d_os.items()   if k in cn}
        wd = compare_disp(d_mine, d_os)
        okd = wd < TOL_DISP_VS_OS
        print(f"  disp diff (ours vs OpenSees)     = {wd:.3e}   {'PASS' if okd else 'FAIL'} (tol {TOL_DISP_VS_OS:.0e})")
        failures += (0 if okd else 1)
        if M.get("skip_forces"):   # hinged end recovers 0 by contract vs OpenSees carrying Mp
            print("  member forces: skipped by design (hinged-end recovery convention differs)")
        else:
            wN, wM = compare_forces(mf_mine, mf_os)
            okf = wN < TOL_FORCE_VS_OS and wM < TOL_FORCE_VS_OS
            print(f"  axial diff (ours vs OpenSees)    = {wN:.3e}   {'PASS' if wN < TOL_FORCE_VS_OS else 'FAIL'}")
            print(f"  moment diff (ours vs OpenSees)   = {wM:.3e}   {'PASS' if wM < TOL_FORCE_VS_OS else 'FAIL'}")
            failures += (0 if okf else 1)

        a = M.get("analytic")
        if a:
            mine = d_mine[a["node"]][a["dof"]]
            osv  = d_os[a["node"]][a["dof"]]
            em = abs(mine - a["value"]) / max(abs(a["value"]), 1e-30)
            eo = abs(osv - a["value"]) / max(abs(a["value"]), 1e-30)
            oka = em < TOL_VS_ANALYTIC and eo < TOL_VS_ANALYTIC
            print(f"  vs analytic: ours={mine:.6g} openSees={osv:.6g} exact={a['value']:.6g}")
            print(f"               relerr ours={em:.3e} openSees={eo:.3e}   {'PASS' if oka else 'FAIL'} (tol {TOL_VS_ANALYTIC:.0e})")
            failures += (0 if oka else 1)

    # ---- MITC4 shell cross-validation vs OpenSees ShellMITC4 ----
    print("\n" + "-" * 64)
    print(" MITC4 shell vs OpenSees ShellMITC4")
    print("-" * 64)
    for M in shell_models():
        sing, d_mine = run_frame_cli_shell(M)
        ok_os, d_os = run_opensees_shell(M)
        print(f"\n[{M['name']}]")
        if sing or ok_os != 0:
            print(f"  [FAIL] solve flagged singular (ours={sing}, openseesAnalyze={ok_os})")
            failures += 1
            continue
        wd = compare_disp(d_mine, d_os)
        okd = wd < TOL_SHELL_VS_OS
        print(f"  disp diff (ours vs OpenSees ShellMITC4) = {wd:.3e}   "
              f"{'PASS' if okd else 'FAIL'} (tol {TOL_SHELL_VS_OS:.0e})")
        failures += (0 if okd else 1)

    # ---- DKQ thin-plate cross-validation vs OpenSees ShellDKGQ (S8-8b) ----
    print("\n" + "-" * 64)
    print(" DKQ thin plate vs OpenSees ShellDKGQ")
    print("-" * 64)
    # DKQ is thin-plate-only; use a FLAT plate (pure bending -> in-plane membrane DOF ~ 0, so the
    # GQ12-vs-Q4 membrane difference is inactive and only the SHARED Batoz DKQ bending is compared).
    dkq_models = [
        shell_plate_model("DKQ thin cantilever plate (flat, pure bending)", 8, 4,
                          1000.0, 500.0, 5.0, 30000.0, 0.3, 0.0, -20.0),
    ]
    for M in dkq_models:
        sing, d_mine = run_frame_cli_shell(M, dkq=True)
        ok_os, d_os = run_opensees_shell(M, dkq=True)
        print(f"\n[{M['name']}]")
        if sing or ok_os != 0:
            print(f"  [FAIL] solve flagged singular (ours={sing}, openseesAnalyze={ok_os})")
            failures += 1
            continue
        wd = compare_disp(d_mine, d_os)
        okd = wd < TOL_SHELL_VS_OS
        print(f"  disp diff (ours vs OpenSees ShellDKGQ) = {wd:.3e}   "
              f"{'PASS' if okd else 'FAIL'} (tol {TOL_SHELL_VS_OS:.0e})")
        failures += (0 if okd else 1)

    # ---- modal cross-validation vs OpenSees eigen ----
    print("\n" + "-" * 64)
    print(" Modal (natural frequencies) vs OpenSees eigen")
    print("-" * 64)
    for M in modal_models():
        f_mine = run_frame_cli_modal(M, M["nModes"])
        f_os = run_opensees_modal(M, M["nModes"])
        print(f"\n[{M['name']}]")
        k = min(len(f_mine), len(f_os), 3)
        worst = max((abs(f_mine[i] - f_os[i]) / max(abs(f_os[i]), 1e-30) for i in range(k)), default=1.0)
        okm = k > 0 and worst < TOL_MODAL
        print(f"  omega(rad/s) ours={[round(x, 4) for x in f_mine[:k]]}")
        print(f"               OS  ={[round(x, 4) for x in f_os[:k]]}")
        print(f"  max rel diff = {worst:.3e}   {'PASS' if okm else 'FAIL'} (tol {TOL_MODAL:.0e})")
        failures += (0 if okm else 1)

    # ---- P-Delta (second order) cross-validation vs OpenSees PDelta geomTransf ----
    print("\n" + "-" * 64)
    print(" P-Delta (second order) vs OpenSees PDelta geomTransf + Newton")
    print("-" * 64)
    print("  scope: OpenSees PDelta transf (axial UPDATED by Newton) vs our consistent geometric")
    print("         stiffness (axial FROZEN at first order); same beam-column target, loose tol.")
    for M in [pdelta_column_model(0.3)]:
        st_f, d_f = run_frame_cli_pdelta(M, 0)        # frozen reuse path
        st_r, d_r = run_frame_cli_pdelta(M, 1)        # K_T reference path
        ok_os, d_os = run_opensees_pdelta(M)
        print(f"\n[{M['name']}]")
        conv = st_f is not None and st_f[0] == 1 and st_r is not None and st_r[0] == 1
        if not conv or ok_os != 0:
            print(f"  [FAIL] pdelta solve issue (frozen={st_f}, reference={st_r}, openseesAnalyze={ok_os})")
            failures += 1
            continue
        tn = M["tip_node"]
        tip_frozen, tip_ref, tip_os, tip_exact = d_f[tn][0], d_r[tn][0], d_os[tn][0], M["exact"]
        e_fr = abs(tip_frozen - tip_ref) / max(abs(tip_ref), 1e-30)        # internal two-path lock
        e_fx = abs(tip_frozen - tip_exact) / max(abs(tip_exact), 1e-30)    # ours vs analytic
        e_ox = abs(tip_os - tip_exact) / max(abs(tip_exact), 1e-30)        # OpenSees vs analytic
        e_mo = abs(tip_frozen - tip_os) / max(abs(tip_os), 1e-30)          # ours vs OpenSees
        okp = e_fr < 1e-9 and e_fx < TOL_PDELTA and e_ox < TOL_PDELTA and e_mo < TOL_PDELTA
        print(f"  tip sway  ours(frozen)={tip_frozen:.6g}  ours(ref)={tip_ref:.6g}  "
              f"OpenSees={tip_os:.6g}  exact={tip_exact:.6g}")
        print(f"  frozen-vs-ref={e_fr:.3e}  ours-vs-exact={e_fx:.3e}  "
              f"OS-vs-exact={e_ox:.3e}  ours-vs-OS={e_mo:.3e}   "
              f"{'PASS' if okp else 'FAIL'} (tol {TOL_PDELTA:.0e})")
        failures += (0 if okp else 1)

    print("\n" + "=" * 64)
    if failures == 0:
        print(" OPENSEES GATE: PASS  (engine matches OpenSees + analytic)")
        return 0
    print(f" OPENSEES GATE: FAIL  ({failures} check(s) out of tolerance)")
    return 1


if __name__ == "__main__":
    sys.exit(main())
