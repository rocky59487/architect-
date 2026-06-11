#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
P-Delta (second-order) detailed cross-validation sweep (S3). A standalone REPORT tool — it is NOT
part of the four-leg gate (the gate's OpenSees leg already runs one P-Delta case inside
opensees_compare.py). This sweeps P/Pcr and tabulates, for the cantilever column:

    ours(frozen)   ours(reference)   OpenSees(PDelta+Newton)   beam-column exact

plus the frozen-path iteration count (which grows ~ log(tol)/log(P/Pcr) toward the Euler load).

Honest scope: our consistent element geometric stiffness (axial FROZEN at first order) reaches the
analytic beam-column closed form to ~1e-5..1e-7 on an 8-element mesh AT EVERY P/Pcr below 1 — that
is the engine's hard oracle and the gate condition here. OpenSees PDelta geomTransf is the
P-large-Delta linearization (axial UPDATED by Newton), which omits the in-element P-small-delta and
so drifts from the closed form as P/Pcr grows (~1e-3 at 0.3, larger near the load) — it is reported
as a cross-tool reference column, NOT gated.

Usage:  python Tools/pdelta_compare.py     (exit 0 iff every case matches the analytic oracle)
"""
import os, sys, subprocess

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import opensees_compare as oc   # reuses pdelta_column_model / run_frame_cli_pdelta / run_opensees_pdelta


def main():
    if not os.path.exists(oc.CLI):
        print("[build] frame_cli.exe missing -> building...")
        if subprocess.run(["cmd", "/c", oc.BUILD_CLI]).returncode != 0 or not os.path.exists(oc.CLI):
            print("[FAIL] could not build frame_cli.exe")
            return 2

    TOL_LOCK   = 1e-9    # frozen vs reference: the two engine paths share one fixed point
    TOL_EXACT  = 2e-3    # ours vs the analytic beam-column closed form (the hard oracle)
    print("=" * 86)
    print(" P-Delta second-order sweep: cantilever column, 8 elements (gate = ours-vs-exact)")
    print("=" * 86)
    print(f"  {'P/Pcr':>6} | {'frozen':>11} {'reference':>11} {'OpenSees':>11} {'exact':>11} "
          f"| {'frz=ref':>8} {'o-vs-ex':>8} {'OS-vs-ex':>8} {'iters':>5}  status")
    failures = 0
    for frac in [0.1, 0.3, 0.5, 0.7, 0.9]:
        M = oc.pdelta_column_model(frac)
        st_f, d_f = oc.run_frame_cli_pdelta(M, 0)
        st_r, d_r = oc.run_frame_cli_pdelta(M, 1)
        ok_os, d_os = oc.run_opensees_pdelta(M)
        tn = M["tip_node"]
        tf, tr, tos, tex = d_f[tn][0], d_r[tn][0], d_os[tn][0], M["exact"]
        e_fr = abs(tf - tr) / max(abs(tr), 1e-30)
        e_fx = abs(tf - tex) / max(abs(tex), 1e-30)
        e_ox = abs(tos - tex) / max(abs(tex), 1e-30)
        iters = st_f[2] if st_f else -1
        ok = (st_f and st_f[0] == 1 and st_r and st_r[0] == 1 and ok_os == 0
              and e_fr < TOL_LOCK and e_fx < TOL_EXACT)
        print(f"  {frac:6.2f} | {tf:11.6g} {tr:11.6g} {tos:11.6g} {tex:11.6g} "
              f"| {e_fr:8.1e} {e_fx:8.1e} {e_ox:8.1e} {iters:5d}  {'PASS' if ok else 'FAIL'}")
        failures += (0 if ok else 1)
    print("=" * 86)
    if failures == 0:
        print(" PDELTA SWEEP: PASS  (every case matches the analytic beam-column oracle)")
        return 0
    print(f" PDELTA SWEEP: FAIL  ({failures} case(s) off the analytic oracle)")
    return 1


if __name__ == "__main__":
    sys.exit(main())
