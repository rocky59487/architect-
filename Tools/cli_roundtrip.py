#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
S6 CLI round-trip integration test (NOT shipped in the game). Builds frame_cli and drives the
J1 / J1b / J1.5 text bridge end-to-end, checking invariants against analytic / engine oracles.
Pure Python (no openseespy). This is gate leg [5/5].

Covers the S6 protocol: VERSION handshake, MEMBER ... tonly token, TONLY / SIZEOPT / DYNC commands,
the MAT cap token (real allowables -> the 10-bar weight now matches standalone F44), DYNC per-frame
DFRAME streaming, and the daemon block loop (many models through one process, EOR-delimited, ==
independent cli per block). The rigorous mechanics live in the F-tests / UE / deep audit; this leg
proves the text bridge round-trips so a Grasshopper / external client can drive the engine.

Usage:  python Tools/cli_roundtrip.py        (exit 0 iff every check passes)
"""
import os, sys, math, subprocess, ctypes

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STAND = os.path.join(ROOT, "Plugins", "FrameSolver", "Standalone")
CLI   = os.path.join(STAND, "frame_cli.exe")
DLL   = os.path.join(STAND, "frame_capi.dll")
BUILD_CLI  = os.path.join(STAND, "build_cli.bat")
BUILD_CAPI = os.path.join(STAND, "build_capi.bat")
IN = 25.4   # mm per inch

_fails = 0
def check(tag, ok, detail=""):
    global _fails
    if not ok:
        _fails += 1
    print(("  [PASS] " if ok else "  [FAIL] ") + tag + (("  " + detail) if detail else ""))


# ------------------------------------------------------------------ section helpers
def sq(side):
    A = side * side
    I = side ** 4 / 12.0
    J = 0.1406 * side ** 4
    return "{} {} {} {} {} {} {} {}".format(A, I, I, J, side / 2, side / 2, 5.0 / 6 * A, 5.0 / 6 * A)

def rect(b, d):
    A = b * d
    Iz = b * d ** 3 / 12.0
    Iy = d * b ** 3 / 12.0
    bb, tt = max(b, d), min(b, d)
    r = tt / bb
    J = bb * tt ** 3 * (1.0 / 3.0 - 0.21 * r * (1.0 - r ** 4 / 12.0))
    return "{} {} {} {} {} {} {} {}".format(A, Iy, Iz, J, b / 2, d / 2, 5.0 / 6 * A, 5.0 / 6 * A)


# ------------------------------------------------------------------ driver (EOR-delimited blocks)
def _parse(rows):
    d = {"VERSION": None, "SINGULAR": None, "DISP": {}, "MF": {},
         "TONLY": None, "SLACK": [], "SIZEOPT": None, "AREA": {}, "WEIGHTVOL": None,
         "DYNC": None, "DEVENT": [], "DFRAME": [], "COROT": None, "_first": None}
    for t in rows:
        if d["_first"] is None:
            d["_first"] = t[0]
        tag = t[0]
        if tag == "VERSION":     d["VERSION"] = t[1] if len(t) > 1 else ""
        elif tag == "SINGULAR":  d["SINGULAR"] = int(t[1])
        elif tag == "DISP":      d["DISP"][int(t[1])] = [float(x) for x in t[2:8]]
        elif tag == "MF":        d["MF"][int(t[1])] = [float(x) for x in t[2:14]]
        elif tag == "TONLY":     d["TONLY"] = (int(t[1]), int(t[2]), int(t[3]))
        elif tag == "SLACK":     d["SLACK"] = [int(x) for x in t[1:]]
        elif tag == "SIZEOPT":   d["SIZEOPT"] = (int(t[1]), int(t[2]), int(t[3]))
        elif tag == "AREA":      d["AREA"][int(t[1])] = (float(t[2]), float(t[3]))
        elif tag == "WEIGHTVOL": d["WEIGHTVOL"] = float(t[1])
        elif tag == "DYNC":      d["DYNC"] = (int(t[1]), int(t[2]), int(t[3]), float(t[4]))
        elif tag == "DEVENT":    d["DEVENT"].append([float(t[1]), int(t[2]), int(t[3]), int(t[4])])
        elif tag == "DFRAME":    d["DFRAME"].append([float(t[1]), float(t[2])])
        elif tag == "COROT":     d["COROT"] = (int(t[1]), int(t[2]), int(t[3]), int(t[4]))
    return d


def run_blocks(blocks):
    """blocks = list of (list-of-lines). One process; returns one parsed dict per EOR-delimited block."""
    text = "".join("\n".join(blk) + "\nEND\n" for blk in blocks)
    p = subprocess.run([CLI], input=text, capture_output=True, text=True, errors="replace")
    if p.returncode != 0:
        raise RuntimeError("frame_cli rc=%d: %s" % (p.returncode, p.stderr))
    results, cur = [], []
    for ln in p.stdout.splitlines():
        t = ln.split()
        if not t:
            continue
        if t[0] == "EOR":
            results.append(_parse(cur)); cur = []
        else:
            cur.append(t)
    if cur:
        results.append(_parse(cur))
    return results


def run(lines):
    return run_blocks([lines])[0]


def ensure_built():
    for bat, art in ((BUILD_CLI, CLI), (BUILD_CAPI, DLL)):
        p = subprocess.run(["cmd", "/c", bat], capture_output=True, text=True, errors="replace")
        if p.returncode != 0 or not os.path.exists(art):
            print("[FAIL] %s failed (rc=%d)" % (os.path.basename(bat), p.returncode))
            print(p.stdout[-800:]); print(p.stderr[-400:])
            sys.exit(1)


def capi_solve(text):
    """Drive the C ABI DLL (frame_capi_solve_text) -- two-call size-then-fill, like a real client."""
    lib = ctypes.CDLL(DLL)
    lib.frame_capi_solve_text.restype = ctypes.c_int
    lib.frame_capi_solve_text.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_int]
    b = (text if text.rstrip().endswith("END") else text.rstrip() + "\nEND\n").encode()
    n = lib.frame_capi_solve_text(b, None, 0)           # query length
    buf = ctypes.create_string_buffer(n + 1)
    lib.frame_capi_solve_text(b, buf, n + 1)
    return buf.value.decode(errors="replace")


# ------------------------------------------------------------------ models
def cantilever(P=1000.0, L=2000.0, side=100.0, E=210000.0, G=80769.0):
    return ["MAT {} {} 7850".format(E, G), "SEC " + sq(side),
            "NODE 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0",
            "NODE 1 {} 0 0 0 0 0 0 0 0 0 0 0 0 0 0".format(L),
            "MEMBER 0 0 1 0 0 0 0 1",
            "NLOAD 1 0 0 {} 0 0 0".format(P)]

def x_portal():
    return ["MAT 200000 76923.07692307692 7850",
            "SEC " + rect(250, 250), "SEC " + rect(60, 60),
            "NODE 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0",
            "NODE 1 6000 0 0 1 1 1 1 1 1 0 0 0 0 0 0",
            "NODE 2 0 0 3000 0 1 0 0 0 0 0 0 0 0 0 0",
            "NODE 3 6000 0 3000 0 1 0 0 0 0 0 0 0 0 0 0",
            "MEMBER 0 0 2 0 0 0 1 0", "MEMBER 1 1 3 0 0 0 1 0", "MEMBER 2 2 3 0 0 0 1 0",
            "MEMBER 3 0 3 0 1 0 1 0 1 1", "MEMBER 4 1 2 0 1 0 1 0 1 1",
            "NLOAD 2 50000 0 0 0 0 0"]

def ten_bar(cap=None):
    bay = 360.0 * IN
    side10 = math.sqrt(10.0 * IN * IN)
    capstr = "" if cap is None else " {} {}".format(cap, cap)
    bars = [(5, 3), (3, 1), (6, 4), (4, 2), (3, 4), (1, 2), (5, 4), (6, 3), (3, 2), (4, 1)]
    Pn = 1.0e5 * 4.4482216152605
    L = ["MAT 68947.572931684 26518.0 7850" + capstr, "SEC " + sq(side10)]
    def nd(i, x, z, sup):
        return "NODE {} {} 0 {} {} 0 0 0 0 0 0".format(i, x, z, "1 1 1" if sup else "0 1 0")
    L += [nd(1, 2 * bay, bay, False), nd(2, 2 * bay, 0, False), nd(3, bay, bay, False),
          nd(4, bay, 0, False), nd(5, 0, bay, True), nd(6, 0, 0, True)]
    for k, (i, j) in enumerate(bars):
        L.append("MEMBER {} {} {} 0 0 0 1 0".format(k, i, j))
    L += ["NLOAD 2 0 0 {} 0 0 0".format(-Pn), "NLOAD 4 0 0 {} 0 0 0".format(-Pn)]
    return L


def main():
    ensure_built()
    print("CLI round-trip integration (frame_cli J1 / J1b / J1.5 bridge)")

    # ---- 1: VERSION handshake + baseline cantilever DISP == PL^3/3EI ----
    side, E, P, L = 100.0, 210000.0, 1000.0, 2000.0
    d = run(cantilever(P, L, side, E))
    check("VERSION handshake is the first stdout line", d["_first"] == "VERSION" and bool(d["VERSION"]),
          "first=%s sha=%s" % (d["_first"], d["VERSION"]))
    exact = P * L ** 3 / (3.0 * E * (side ** 4 / 12.0))
    got = d["DISP"][1][2]
    rel = abs(got - exact) / abs(exact)
    check("baseline DISP cantilever == PL^3/3EI (rel<1e-6)", rel < 1e-6,
          "got=%.6g exact=%.6g rel=%.2e" % (got, exact, rel))

    # ---- 2: TONLY -- X-braced portal, one diagonal goes slack (F42 plumbing) ----
    d = run(x_portal() + ["TONLY"])
    conv = d["TONLY"][0] if d["TONLY"] else -1
    check("TONLY converged + exactly one brace slack",
          d["TONLY"] is not None and conv == 1 and len(d["SLACK"]) == 1,
          "TONLY=%s slack=%s" % (d["TONLY"], d["SLACK"]))

    # ---- 3: SIZEOPT plumbing (no cap token -> hardcoded allowable; converges, 10 areas, vol>0) ----
    d = run(ten_bar(cap=None) + ["SIZEOPT {} 100 1e-8".format(0.1 * IN * IN)])
    so = d["SIZEOPT"]
    check("SIZEOPT converged + 10 sized areas + positive volume",
          so is not None and so[0] == 1 and so[2] == 0 and len(d["AREA"]) == 10
          and d["WEIGHTVOL"] and d["WEIGHTVOL"] > 0,
          "SIZEOPT=%s nAreas=%d vol=%s" % (so, len(d["AREA"]), d["WEIGHTVOL"]))

    # ---- 4: J1b MAT cap token -> real allowables; literature-INDEPENDENT FSD invariants ----
    # (avoid a self-referential "==engine value" oracle: assert the pin-jointed literature lower
    #  bound 1593.2 lb, the fully-stressed invariant on the sized bars, AND consistency with F44.)
    sigA = 25000.0 * 0.0068947572931684          # 25000 psi -> MPa
    d = run(ten_bar(cap=sigA) + ["SIZEOPT {} 100 1e-10".format(0.1 * IN * IN)])
    weight_lb = 0.1 * d["WEIGHTVOL"] / (IN ** 3)  # 0.1 lb/in^3 * sum(A_in2*L_in) = 0.1*WEIGHTVOL/IN^3
    fs_bars = sum(1 for (a, dc) in d["AREA"].values() if 0.999 <= dc <= 1.001)  # fully-stressed bars
    check("MAT cap token -> 10-bar FSD: weight >= lit 1593.2, sized bars D/C=1, == F44",
          d["SIZEOPT"][0] == 1 and weight_lb >= 1593.0 and fs_bars >= 5 and abs(weight_lb - 1608.49) < 2.0,
          "weight_lb=%.4f fs_bars=%d" % (weight_lb, fs_bars))

    # ---- 5: DYNC per-frame DFRAME streaming (loaded cantilever vibrates; frames stored) ----
    dynbeam = ["MAT 210000 80769 7850", "SEC " + sq(100.0),
               "NODE 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0",
               "NODE 1 1000 0 0 0 0 0 0 0 0 0 0 0 0 0 0",
               "NODE 2 2000 0 0 0 0 0 0 0 0 0 0 0 0 0 0",
               "MEMBER 0 0 1 0 0 0 0 1", "MEMBER 1 1 2 0 0 0 0 1",
               "NLOAD 2 0 0 -5000 0 0 0", "DYNC 0.001 0.03"]
    d = run(dynbeam)
    dy = d["DYNC"]
    check("DYNC streams per-frame DFRAME (nFrames>0, DFRAME count matches)",
          dy is not None and dy[0] in (0, 1, 2) and dy[2] > 0 and len(d["DFRAME"]) == dy[2]
          and all(d["DFRAME"][k][0] <= d["DFRAME"][k + 1][0] for k in range(len(d["DFRAME"]) - 1)),
          "DYNC=%s nDFRAME=%d" % (dy, len(d["DFRAME"])))

    # ---- 6: J1.5 daemon -- two models through ONE process == two independent cli runs ----
    multi = run_blocks([cantilever(P, L, side, E), x_portal() + ["TONLY"]])
    indep0 = run(cantilever(P, L, side, E))
    indep1 = run(x_portal() + ["TONLY"])
    ok_count = len(multi) == 2
    same0 = ok_count and multi[0]["DISP"].get(1) == indep0["DISP"].get(1)
    same1 = ok_count and multi[1]["TONLY"] == indep1["TONLY"] and multi[1]["SLACK"] == indep1["SLACK"]
    check("daemon: 2 blocks/1 process == 2 independent cli runs (byte-identical)",
          ok_count and same0 and same1,
          "nBlocks=%d same0=%s same1=%s" % (len(multi), same0, same1))

    # ---- 7: J2 C API DLL -- frame_capi_solve_text == frame_cli.exe (byte-identical protocol) ----
    model_lines = cantilever(P, L, side, E)
    capi_out = capi_solve("\n".join(model_lines))
    # the CLI stdout for the same single block (processAll backs both, so output is byte-identical)
    p = subprocess.run([CLI], input="\n".join(model_lines) + "\nEND\n",
                        capture_output=True, text=True, errors="replace")
    cli_out = p.stdout
    same = capi_out.strip() == cli_out.strip()
    capi_ver = capi_out.split()[1] if capi_out.split() and capi_out.split()[0] == "VERSION" else ""
    check("C API DLL frame_capi_solve_text == frame_cli.exe (byte-identical)",
          same and bool(capi_ver),
          "identical=%s capi_sha=%s len=%d" % (same, capi_ver, len(capi_out)))

    # ---- 8: COROT planar cantilever elastica (S9) -- large-displacement tip vs Mattiasson table ----
    def planar_cantilever(n, Lp, Pp, E=1.0, A=100.0, Iz=1e-3):
        lines = ["MAT {} 0.4 0".format(E), "SEC {} {} {} {} 1 1 0 0".format(A, Iz, Iz, Iz)]
        for k in range(n + 1):
            x = Lp * k / n
            # node 0 encastre; others: in-plane (Ux,Uy,Rz) free, out-of-plane (Uz,Rx,Ry) fixed.
            flags = "1 1 1 1 1 1" if k == 0 else "0 0 1 1 1 0"
            lines.append("NODE {} {} 0 0 {} 0 0 0 0 0 0".format(k, x, flags))
        for k in range(n):
            lines.append("MEMBER {} {} {} 0 0 0 0 1".format(k, k, k + 1))
        lines.append("NLOAD {} 0 {} 0 0 0 0".format(n, Pp))   # transverse +Y tip load
        return lines

    alpha, Lp, Ee, Iz, n = 5.0, 1.0, 1.0, 1e-3, 16
    Pp = alpha * Ee * Iz / Lp ** 2
    d = run(planar_cantilever(n, Lp, Pp) + ["COROT 20 80 1e-9"])
    co = d["COROT"]
    dv = d["DISP"][n][1] / Lp     # tip Uy / L (transverse)
    dh = -d["DISP"][n][0] / Lp    # -tip Ux / L (axial shortening)
    dv_exact, dh_exact = 0.7137915236, 0.3876283607   # WS_F F-3 elastica table, alpha=5
    check("COROT planar elastica alpha=5 (tip dv/L,dh/L vs Mattiasson; converged, not diverged)",
          co is not None and co[0] == 1 and co[1] == 0
          and abs(dv - dv_exact) < 2e-3 and abs(dh - dh_exact) < 5e-3,
          "COROT=%s dv=%.6f(exp %.6f) dh=%.6f(exp %.6f)" % (co, dv, dv_exact, dh, dh_exact))

    print("\n%s  (failures=%d)" % ("ALL PASS" if _fails == 0 else "FAILURES", _fails))
    return 0 if _fails == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
