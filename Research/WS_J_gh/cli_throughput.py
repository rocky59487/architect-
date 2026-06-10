"""cli_throughput.py — WS-J end-to-end frame_cli throughput measurement (scratch, NOT engine).

Measures what a Grasshopper CLI-bridge component would actually pay per solve:
  spawn frame_cli.exe + write model text to stdin + engine solve + parse stdout.
Sizes: trivial (process overhead), ~1.6k DOF, ~19.5k DOF, ~107k DOF towers
(same topology family as Standalone/frame_perf.cpp's makeTower, ported here).
"""

import os
import subprocess
import time

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
CLI = os.path.join(ROOT, "Plugins", "FrameSolver", "Standalone", "frame_cli.exe")


def square_sec(side):
    A = side * side
    I = side ** 4 / 12.0
    J = 0.1406 * side ** 4
    return f"SEC {A} {I} {I} {J} {side/2} {side/2} 0 0"


def tower_text(nx, ny, st):
    sx, sy, h = 6000.0, 5000.0, 3300.0
    nid = lambda i, j, k: k * (nx + 1) * (ny + 1) + j * (nx + 1) + i
    out = ["MAT 200000 76923.07692307692 7850"]
    out += [square_sec(520.0), square_sec(360.0), square_sec(220.0)]
    pos = {}
    for k in range(st + 1):
        dx, dy = 35.0 * k, -18.0 * k
        for j in range(ny + 1):
            for i in range(nx + 1):
                x, y, z = i * sx + dx, j * sy + dy, k * h
                pos[nid(i, j, k)] = (x, y, z)
                fix = "1 1 1 1 1 1" if k == 0 else "0 0 0 0 0 0"
                out.append(f"NODE {nid(i,j,k)} {x} {y} {z} {fix}")
    mid = 0

    def member(a, b, sec):
        nonlocal mid
        (xa, ya, za), (xb, yb, zb) = pos[a], pos[b]
        dx, dy, dz = xb - xa, yb - ya, zb - za
        n = (dx * dx + dy * dy + dz * dz) ** 0.5
        ref = "0 0 1" if abs(dz / n) < 0.92 else "1 0 0"
        out.append(f"MEMBER {mid} {a} {b} 0 {sec} {ref}")
        mid += 1

    for k in range(st):
        for j in range(ny + 1):
            for i in range(nx + 1):
                member(nid(i, j, k), nid(i, j, k + 1), 0)
    for k in range(1, st + 1):
        for j in range(ny + 1):
            for i in range(nx):
                member(nid(i, j, k), nid(i + 1, j, k), 1)
        for j in range(ny):
            for i in range(nx + 1):
                member(nid(i, j, k), nid(i, j + 1, k), 1)
    for k in range(st):
        for i in range(nx):
            for face in (0, 1):
                jf = ny if face else 0
                a = nid(i, jf, k) if (i + k) % 2 == 0 else nid(i + 1, jf, k)
                b = nid(i + 1, jf, k + 1) if (i + k) % 2 == 0 else nid(i, jf, k + 1)
                member(a, b, 2)
        for j in range(ny):
            for face in (0, 1):
                ifc = nx if face else 0
                a = nid(ifc, j, k) if (j + k) % 2 == 0 else nid(ifc, j + 1, k)
                b = nid(ifc, j + 1, k + 1) if (j + k) % 2 == 0 else nid(ifc, j, k + 1)
                member(a, b, 2)
    for k in range(1, st + 1):
        sf = k / st
        for j in range(ny + 1):
            for i in range(nx + 1):
                edge = (1 if i in (0, nx) else 0) + (1 if j in (0, ny) else 0)
                trib = 1.0 if edge == 0 else (0.65 if edge == 1 else 0.42)
                fx = 900.0 * sf * (1.0 + 0.08 * (j - ny / 2.0))
                fy = -350.0 * sf * (1.0 + 0.05 * (i - nx / 2.0))
                out.append(f"NLOAD {nid(i,j,k)} {fx} {fy} {-38000.0*trib} 0 0 0")
    out.append("END")
    return "\n".join(out) + "\n", (nx + 1) * (ny + 1) * (st + 1) * 6


def trivial_text():
    return ("MAT 200000 76923 7850\n" + square_sec(200.0) +
            "\nNODE 0 0 0 0 1 1 1 1 1 1\nNODE 1 0 0 3000 0 0 0 0 0 0\n"
            "MEMBER 0 0 1 0 0 1 0 0\nNLOAD 1 1000 0 0 0 0 0\nEND\n"), 12


def bench(name, text, dof, runs=3):
    times = []
    out_bytes = 0
    for _ in range(runs):
        t0 = time.perf_counter()
        p = subprocess.run([CLI], input=text, capture_output=True, text=True)
        t1 = time.perf_counter()
        if p.returncode not in (0,):
            print(f"[cli] {name} FAILED rc={p.returncode} stderr={p.stderr[:200]}")
            return
        out_bytes = len(p.stdout)
        times.append((t1 - t0) * 1000.0)
    times.sort()
    print(f"[cli] {name} dof={dof} inKB={len(text)//1024} outKB={out_bytes//1024} "
          f"runs={runs} medianMs={times[len(times)//2]:.1f} minMs={times[0]:.1f}")


def main():
    if not os.path.exists(CLI):
        print("[cli] frame_cli.exe missing — build with Standalone\\build_cli.bat first")
        return
    t, d = trivial_text()
    bench("trivial(overhead)", t, d, runs=5)
    for (nx, ny, st), label in (((5, 4, 8), "tower-1.6kDOF"),
                                ((12, 9, 24), "tower-19.5kDOF"),
                                ((22, 18, 40), "tower-107kDOF")):
        t0 = time.perf_counter()
        text, dof = tower_text(nx, ny, st)
        genMs = (time.perf_counter() - t0) * 1000.0
        print(f"[gen] {label} genMs={genMs:.1f}")
        bench(label, text, dof, runs=3)


if __name__ == "__main__":
    main()
