#!/usr/bin/env python3
r"""FrameCore high-performance benchmark/tuning harness prototype.

This script intentionally treats FrameCore as a black box. It runs the existing
standalone/research executables, parses their stable text output, writes JSONL,
and summarizes repeated runs with median/geomean statistics.

Examples:
  python Research\WS_B_solver\run_hp_tuning.py --suite smoke --runs 2
  python Research\WS_B_solver\run_hp_tuning.py --suite frame_perf --grid "3x2x4,5x4x8" --runs 3
  python Research\WS_B_solver\run_hp_tuning.py --suite ws_b --include-large --timeout-sec 1800
"""

from __future__ import annotations

import argparse
import csv
import datetime as _dt
import itertools
import json
import math
import os
import platform
import re
import statistics
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_FRAME_PERF = ROOT / "Plugins" / "FrameSolver" / "Standalone" / "frame_perf.exe"
DEFAULT_RESEARCH_BIN = ROOT / "Research" / "bin"
DEFAULT_OUT_ROOT = ROOT / "Research" / "out" / "hp_tuning"

SCHEMA_VERSION = 1
KEY_VALUE_RE = re.compile(r'([A-Za-z_][A-Za-z0-9_]*)=("[^"]*"|[^\s|]+)')
BRACKET_RE = re.compile(r"^\[([^\]]+)\]\s*(.*)$")

PRIMARY_METRICS = {
    "frame_perf.solve_median_ms",
    "hpfem.mf_apply_ms",
    "hpfem.pcg_iters",
    "hpfem.pcg_ms",
    "hpfem.framecore_bsr6_apply_ms",
    "hpfem.threaded_apply_speedup",
    "hpfem.parallel_pcg_speedup_ms",
    "hpfem.combined_hpfem_speedup_ms",
    "hpfem.combined_hpfem_apply_ms_avg",
    "hpfem.combined_hpfem_precond_ms_avg",
    "hpfem.combined_hpfem_other_ms_avg",
    "hpfem.combined_hpfem_initial_rel",
    "hpfem.combined_hpfem_seed_ms",
    "hpfem.factor_bypass_setup_speedup",
    "hpfem.factor_bypass_first_solve_speedup",
    "hpfem.factor_bypass_batch_speedup",
    "hpfem.recycle_speedup_ms",
    "hpfem.recycle_rec_iters_avg",
    "hpfem.line_schwarz_speedup_ms",
    "hpfem.line_schwarz_iters",
    "hpfem.full_apply_passed",
    "hpfem.full_apply_max_reaction_abs",
    "hpfem.mechanism_guard_passed",
    "scale.factor_plus_solve_ms",
    "solver_compare.total_ms",
    "sparse_buckling.sparse_ms",
}


@dataclass(frozen=True)
class Job:
    suite: str
    label: str
    exe: Path
    args: tuple[str, ...]
    params: dict[str, Any]


def now_stamp() -> str:
    return _dt.datetime.now().strftime("%Y%m%d_%H%M%S")


def iso_now() -> str:
    return _dt.datetime.now().astimezone().isoformat(timespec="seconds")


def command_text(cmd: Iterable[str]) -> str:
    if os.name == "nt":
        return subprocess.list2cmdline(list(cmd))
    import shlex

    return shlex.join(list(cmd))


def parse_scalar(raw: str) -> Any:
    value = raw.strip().rstrip(",")
    if len(value) >= 2 and value[0] == '"' and value[-1] == '"':
        return value[1:-1]
    if value.lower() == "nan":
        return float("nan")
    if re.fullmatch(r"[+-]?\d+", value):
        try:
            return int(value)
        except ValueError:
            pass
    try:
        return float(value)
    except ValueError:
        return value


def extract_kv(text: str) -> dict[str, Any]:
    return {m.group(1): parse_scalar(m.group(2)) for m in KEY_VALUE_RE.finditer(text)}


def finite_number(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value))


def context_key(context: dict[str, Any]) -> str:
    if not context:
        return ""
    return "|".join(f"{key}={context[key]}" for key in sorted(context))


def emit_metric(
    out: list[dict[str, Any]],
    metric: str,
    value: Any,
    unit: str,
    line: str,
    context: dict[str, Any] | None = None,
    fields: dict[str, Any] | None = None,
) -> None:
    if not finite_number(value):
        return
    out.append(
        {
            "metric": metric,
            "value": float(value),
            "unit": unit,
            "context": context or {},
            "line": line,
            "fields": fields or {},
        }
    )


def parse_frame_perf(stdout: str) -> list[dict[str, Any]]:
    metrics: list[dict[str, Any]] = []
    case_context: dict[str, Any] = {}
    for raw_line in stdout.splitlines():
        line = raw_line.strip()
        if line.startswith("case="):
            fields = extract_kv(line)
            case_context = {"case": fields.get("case", "frame_perf")}
            continue
        if line.startswith("nodes="):
            fields = extract_kv(line)
            raw_mib = re.search(r"\(~([0-9.]+)\s+MiB raw triplets\)", line)
            if raw_mib:
                fields["tripletRawMiB"] = parse_scalar(raw_mib.group(1))
            emit_metric(metrics, "frame_perf.build_ms", fields.get("buildMs"), "ms", line, case_context, fields)
            emit_metric(metrics, "frame_perf.free_dof", fields.get("freeDof"), "count", line, case_context, fields)
            emit_metric(metrics, "frame_perf.nodes", fields.get("nodes"), "count", line, case_context, fields)
            emit_metric(metrics, "frame_perf.members", fields.get("members"), "count", line, case_context, fields)
            emit_metric(metrics, "frame_perf.triplet_raw_mib", fields.get("tripletRawMiB"), "MiB", line, case_context, fields)
            continue
        if line.startswith("solveMs "):
            fields = extract_kv(line)
            emit_metric(metrics, "frame_perf.solve_median_ms", fields.get("median"), "ms", line, case_context, fields)
            emit_metric(metrics, "frame_perf.solve_min_ms", fields.get("min"), "ms", line, case_context, fields)
            emit_metric(metrics, "frame_perf.solve_max_ms", fields.get("max"), "ms", line, case_context, fields)
            emit_metric(metrics, "frame_perf.solve_mean_ms", fields.get("mean"), "ms", line, case_context, fields)
            emit_metric(metrics, "frame_perf.checksum", fields.get("checksum"), "checksum", line, case_context, fields)
    return metrics


def parse_scale(stdout: str) -> list[dict[str, Any]]:
    metrics: list[dict[str, Any]] = []
    for raw_line in stdout.splitlines():
        line = raw_line.strip()
        match = BRACKET_RE.match(line)
        if not match:
            continue
        tag, rest = match.group(1), match.group(2)
        fields = extract_kv(rest)
        context = {"case": f"{fields.get('nx')}x{fields.get('ny')}x{fields.get('st')}"}
        if tag == "scale-begin":
            emit_metric(metrics, "scale.build_ms", fields.get("buildMs"), "ms", line, context, fields)
            emit_metric(metrics, "scale.mem_mib", fields.get("memMiB"), "MiB", line, context, fields)
            emit_metric(metrics, "scale.nf", fields.get("nf"), "count", line, context, fields)
        elif tag == "scale":
            factor = fields.get("factorMs")
            solve = fields.get("solveMs")
            emit_metric(metrics, "scale.factor_ms", factor, "ms", line, context, fields)
            emit_metric(metrics, "scale.solve_ms", solve, "ms", line, context, fields)
            if finite_number(factor) and finite_number(solve):
                emit_metric(
                    metrics,
                    "scale.factor_plus_solve_ms",
                    float(factor) + float(solve),
                    "ms",
                    line,
                    context,
                    fields,
                )
            emit_metric(metrics, "scale.mem_factor_mib", fields.get("memFactorMiB"), "MiB", line, context, fields)
            emit_metric(metrics, "scale.peak_mib", fields.get("peakMiB"), "MiB", line, context, fields)
            emit_metric(metrics, "scale.nf", fields.get("nf"), "count", line, context, fields)
            emit_metric(metrics, "scale.nnz_k", fields.get("nnzK"), "count", line, context, fields)
            emit_metric(metrics, "scale.umax", fields.get("umax"), "model-unit", line, context, fields)
        elif tag == "oom":
            context["status"] = "oom"
            emit_metric(metrics, "scale.peak_mib", fields.get("peakMiB"), "MiB", line, context, fields)
            emit_metric(metrics, "scale.nf", fields.get("nf"), "count", line, context, fields)
    return metrics


def parse_solver_compare(stdout: str) -> list[dict[str, Any]]:
    metrics: list[dict[str, Any]] = []
    for raw_line in stdout.splitlines():
        line = raw_line.strip()
        match = BRACKET_RE.match(line)
        if not match:
            continue
        tag, rest = match.group(1), match.group(2)
        fields = extract_kv(rest)
        parts = rest.split()
        if tag == "case" and parts:
            context = {"case": parts[0]}
            emit_metric(metrics, "solver_compare.nf", fields.get("nf"), "count", line, context, fields)
            emit_metric(metrics, "solver_compare.nnz", fields.get("nnz"), "count", line, context, fields)
        elif tag == "solver" and parts:
            case_name = parts[0]
            solver_parts: list[str] = []
            for token in parts[1:]:
                if "=" in token:
                    break
                solver_parts.append(token)
            solver_name = " ".join(solver_parts) if solver_parts else "unknown"
            context = {"case": case_name, "solver": solver_name}
            factor = fields.get("factorMs")
            setup = fields.get("setupMs")
            solve = fields.get("solveMs")
            emit_metric(metrics, "solver_compare.factor_ms", factor, "ms", line, context, fields)
            emit_metric(metrics, "solver_compare.setup_ms", setup, "ms", line, context, fields)
            emit_metric(metrics, "solver_compare.solve_ms", solve, "ms", line, context, fields)
            base = factor if finite_number(factor) else setup
            if finite_number(base) and finite_number(solve):
                emit_metric(
                    metrics,
                    "solver_compare.total_ms",
                    float(base) + float(solve),
                    "ms",
                    line,
                    context,
                    fields,
                )
            emit_metric(metrics, "solver_compare.iters", fields.get("iters"), "count", line, context, fields)
            emit_metric(metrics, "solver_compare.cg_err", fields.get("cgErr"), "ratio", line, context, fields)
            emit_metric(metrics, "solver_compare.residual", fields.get("res"), "ratio", line, context, fields)
    return metrics


def parse_sparse_buckling(stdout: str) -> list[dict[str, Any]]:
    metrics: list[dict[str, Any]] = []
    for raw_line in stdout.splitlines():
        line = raw_line.strip()
        match = BRACKET_RE.match(line)
        if not match:
            continue
        tag, rest = match.group(1), match.group(2)
        if tag in {"done", "scale", "solver", "case"}:
            continue
        fields = extract_kv(rest)
        context = {"case": tag}
        emit_metric(metrics, "sparse_buckling.nf", fields.get("nf"), "count", line, context, fields)
        emit_metric(metrics, "sparse_buckling.dense_run", fields.get("denseRun"), "bool", line, context, fields)
        emit_metric(metrics, "sparse_buckling.lam_dense", fields.get("lamDense"), "factor", line, context, fields)
        emit_metric(metrics, "sparse_buckling.lam_sparse", fields.get("lamSparse"), "factor", line, context, fields)
        emit_metric(metrics, "sparse_buckling.dense_ms", fields.get("tDenseMs"), "ms", line, context, fields)
        emit_metric(metrics, "sparse_buckling.sparse_ms", fields.get("tSparseMs"), "ms", line, context, fields)
        emit_metric(metrics, "sparse_buckling.residual", fields.get("residual"), "ratio", line, context, fields)
        emit_metric(
            metrics,
            "sparse_buckling.rel_err_sparse_vs_dense",
            fields.get("relErrSparseVsDense"),
            "ratio",
            line,
            context,
            fields,
        )
        emit_metric(
            metrics,
            "sparse_buckling.rel_err_vs_analytic",
            fields.get("relErrVsAnalytic"),
            "ratio",
            line,
            context,
            fields,
        )
        emit_metric(metrics, "sparse_buckling.dense_mem_est_mib", fields.get("denseMemEstMiB"), "MiB", line, context, fields)
        emit_metric(metrics, "sparse_buckling.kg_empty", fields.get("KgEmpty"), "bool", line, context, fields)
        emit_metric(metrics, "sparse_buckling.engine_singular", fields.get("engineSingular"), "bool", line, context, fields)
    return metrics


def parse_hpfem(stdout: str) -> list[dict[str, Any]]:
    metrics: list[dict[str, Any]] = []
    for raw_line in stdout.splitlines():
        line = raw_line.strip()
        match = BRACKET_RE.match(line)
        if not match:
            continue
        tag, rest = match.group(1), match.group(2)
        if tag not in {"mfop", "framecore_bsr6", "threaded_apply", "parallel_pcg", "combined_hpfem", "recycle", "line_schwarz", "full_apply", "full_apply_summary", "mechanism_guard", "mechanism_guard_summary"}:
            continue
        fields = extract_kv(rest)
        if tag == "mechanism_guard_summary":
            context = {"case": "all", "experiment": "mechanism_guard_oracle"}
            emit_metric(metrics, "hpfem.mechanism_guard_passed", fields.get("passed"), "count", line, context, fields)
            emit_metric(metrics, "hpfem.mechanism_guard_total", fields.get("total"), "count", line, context, fields)
            continue
        if tag == "full_apply_summary":
            context = {"case": "all", "experiment": "full_apply_oracle"}
            emit_metric(metrics, "hpfem.full_apply_passed", fields.get("passed"), "count", line, context, fields)
            emit_metric(metrics, "hpfem.full_apply_total", fields.get("total"), "count", line, context, fields)
            emit_metric(metrics, "hpfem.full_apply_max_apply_rel", fields.get("maxApplyRel"), "ratio", line, context, fields)
            emit_metric(metrics, "hpfem.full_apply_max_rhs_rel", fields.get("maxRhsRel"), "ratio", line, context, fields)
            emit_metric(metrics, "hpfem.full_apply_max_u_rel", fields.get("maxURel"), "ratio", line, context, fields)
            emit_metric(metrics, "hpfem.full_apply_max_reaction_rel", fields.get("maxReactionRel"), "ratio", line, context, fields)
            emit_metric(metrics, "hpfem.full_apply_max_reaction_abs", fields.get("maxReactionAbs"), "force", line, context, fields)
            continue
        experiment = {
            "mfop": "matrix_free_operator",
            "framecore_bsr6": "framecore_bsr6",
            "threaded_apply": "threaded_apply",
            "parallel_pcg": "parallel_pcg",
            "combined_hpfem": "combined_hpfem",
            "recycle": "pcg_recycling",
            "line_schwarz": "line_schwarz",
            "full_apply": "full_apply_oracle",
            "mechanism_guard": "mechanism_guard_oracle",
        }[tag]
        context = {
            "case": fields.get("preset", "unknown"),
            "experiment": experiment,
        }
        if tag == "mfop":
            context["precond"] = fields.get("precond", "unknown")
            context["coarseBins"] = fields.get("coarseBins", "1x1")
        if tag in {"threaded_apply", "parallel_pcg", "combined_hpfem"}:
            context["threads"] = fields.get("threads", "unknown")
        if tag in {"parallel_pcg", "combined_hpfem"}:
            context["precond"] = fields.get("precond", "unknown")
            context["coarseBins"] = fields.get("coarseBins", "1x1")
        if tag == "combined_hpfem":
            context["basisMax"] = fields.get("basisMax", "unknown")
            context["rhs"] = fields.get("rhs", "unknown")
            context["seedLoadBasis"] = fields.get("seedLoadBasis", 0)
        if tag == "recycle":
            context["basisMax"] = fields.get("basisMax", "unknown")
        if tag == "line_schwarz":
            context["mode"] = "floor_lines" if fields.get("floorLines", 1) else "columns_only"
        emit_metric(metrics, "hpfem.nf", fields.get("nf"), "count", line, context, fields)
        emit_metric(metrics, "hpfem.apply_rel", fields.get("applyRel"), "ratio", line, context, fields)
        if "sparseApplyMs" in fields:
            emit_metric(metrics, "hpfem.sparse_apply_ms", fields.get("sparseApplyMs"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.speedup_vs_sparse_apply", fields.get("speedup"), "ratio", line, context, fields)
        if tag == "mfop":
            emit_metric(metrics, "hpfem.coarse_dofs", fields.get("coarseDofs"), "count", line, context, fields)
            emit_metric(metrics, "hpfem.mf_apply_ms", fields.get("mfApplyMs"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.pcg_ms", fields.get("pcgMs"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.pcg_iters", fields.get("pcgIters"), "count", line, context, fields)
            emit_metric(metrics, "hpfem.pcg_rel", fields.get("pcgRel"), "ratio", line, context, fields)
            emit_metric(metrics, "hpfem.pcg_vs_ldlt", fields.get("pcgVsLdlt"), "ratio", line, context, fields)
            emit_metric(metrics, "hpfem.ldlt_solve_ms", fields.get("ldltSolveMs"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.operator_build_ms", fields.get("opBuildMs"), "ms", line, context, fields)
        elif tag == "framecore_bsr6":
            emit_metric(metrics, "hpfem.framecore_bsr6_apply_ms", fields.get("bsrApplyMs"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.framecore_bsr6_build_ms", fields.get("buildMs"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.block_nnz", fields.get("blockNnz"), "count", line, context, fields)
        elif tag == "threaded_apply":
            emit_metric(metrics, "hpfem.threaded_apply_seq_ms", fields.get("seqMs"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.threaded_apply_par_ms", fields.get("parMs"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.threaded_apply_speedup", fields.get("speedup"), "ratio", line, context, fields)
            emit_metric(metrics, "hpfem.threaded_apply_spawn_par_ms", fields.get("spawnParMs"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.threaded_apply_spawn_speedup", fields.get("spawnSpeedup"), "ratio", line, context, fields)
            emit_metric(metrics, "hpfem.blocks", fields.get("blocks"), "count", line, context, fields)
        elif tag == "parallel_pcg":
            emit_metric(metrics, "hpfem.parallel_pcg_coarse_dofs", fields.get("coarseDofs"), "count", line, context, fields)
            emit_metric(metrics, "hpfem.hp_setup_ms", fields.get("hpSetupMs"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.preconditioner_build_ms", fields.get("precondBuildMs"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.factor_bypass_setup_speedup", fields.get("factorBypassSetupSpeedup"), "ratio", line, context, fields)
            emit_metric(metrics, "hpfem.factor_bypass_first_solve_speedup", fields.get("factorBypassFirstSolveSpeedup"), "ratio", line, context, fields)
            emit_metric(metrics, "hpfem.direct_first_ms", fields.get("directFirstMs"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.hp_first_ms", fields.get("hpFirstMs"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.ldlt_solve_ms", fields.get("ldltSolveMs"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.parallel_pcg_serial_iters", fields.get("serialIters"), "count", line, context, fields)
            emit_metric(metrics, "hpfem.parallel_pcg_iters", fields.get("parallelIters"), "count", line, context, fields)
            emit_metric(metrics, "hpfem.parallel_pcg_serial_ms", fields.get("serialMs"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.parallel_pcg_serial_apply_ms", fields.get("serialApplyMs"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.parallel_pcg_serial_precond_ms", fields.get("serialPrecondMs"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.parallel_pcg_serial_other_ms", fields.get("serialOtherMs"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.parallel_pcg_ms", fields.get("parallelMs"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.parallel_pcg_apply_ms", fields.get("parallelApplyMs"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.parallel_pcg_precond_ms", fields.get("parallelPrecondMs"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.parallel_pcg_other_ms", fields.get("parallelOtherMs"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.parallel_pcg_speedup_ms", fields.get("speedupMs"), "ratio", line, context, fields)
            emit_metric(metrics, "hpfem.parallel_pcg_true_rel", fields.get("parallelTrueRel"), "ratio", line, context, fields)
            emit_metric(metrics, "hpfem.parallel_pcg_vs_ldlt", fields.get("pcgVsLdlt"), "ratio", line, context, fields)
        elif tag == "combined_hpfem":
            emit_metric(metrics, "hpfem.combined_hpfem_coarse_dofs", fields.get("coarseDofs"), "count", line, context, fields)
            emit_metric(metrics, "hpfem.combined_hpfem_seed_count", fields.get("seedCount"), "count", line, context, fields)
            emit_metric(metrics, "hpfem.combined_hpfem_seed_iters", fields.get("seedIters"), "count", line, context, fields)
            emit_metric(metrics, "hpfem.combined_hpfem_seed_ms", fields.get("seedMs"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.combined_hpfem_seed_apply_ms", fields.get("seedApplyMs"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.combined_hpfem_seed_precond_ms", fields.get("seedPrecondMs"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.hp_setup_ms", fields.get("hpSetupMs"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.preconditioner_build_ms", fields.get("precondBuildMs"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.ldlt_solve_ms_avg", fields.get("ldltSolveMsAvg"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.factor_bypass_setup_speedup", fields.get("factorBypassSetupSpeedup"), "ratio", line, context, fields)
            emit_metric(metrics, "hpfem.factor_bypass_first_solve_speedup", fields.get("factorBypassFirstSolveSpeedup"), "ratio", line, context, fields)
            emit_metric(metrics, "hpfem.factor_bypass_batch_speedup", fields.get("factorBypassBatchSpeedup"), "ratio", line, context, fields)
            emit_metric(metrics, "hpfem.direct_batch_ms", fields.get("directBatchMs"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.hp_batch_ms", fields.get("hpBatchMs"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.combined_hpfem_serial_iters_avg", fields.get("serialItersAvg"), "count", line, context, fields)
            emit_metric(metrics, "hpfem.combined_hpfem_iters_avg", fields.get("combinedItersAvg"), "count", line, context, fields)
            emit_metric(metrics, "hpfem.combined_hpfem_iters_skip1_avg", fields.get("combinedItersSkip1Avg"), "count", line, context, fields)
            emit_metric(metrics, "hpfem.combined_hpfem_speedup_iters", fields.get("speedupIters"), "ratio", line, context, fields)
            emit_metric(metrics, "hpfem.combined_hpfem_serial_ms_avg", fields.get("serialMsAvg"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.combined_hpfem_serial_apply_ms_avg", fields.get("serialApplyMsAvg"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.combined_hpfem_serial_precond_ms_avg", fields.get("serialPrecondMsAvg"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.combined_hpfem_serial_other_ms_avg", fields.get("serialOtherMsAvg"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.combined_hpfem_ms_avg", fields.get("combinedMsAvg"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.combined_hpfem_apply_ms_avg", fields.get("combinedApplyMsAvg"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.combined_hpfem_precond_ms_avg", fields.get("combinedPrecondMsAvg"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.combined_hpfem_other_ms_avg", fields.get("combinedOtherMsAvg"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.combined_hpfem_speedup_ms", fields.get("speedupMs"), "ratio", line, context, fields)
            emit_metric(metrics, "hpfem.combined_hpfem_project_ms_avg", fields.get("projectMsAvg"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.combined_hpfem_basis_add_ms", fields.get("basisAddMs"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.combined_hpfem_initial_rel", fields.get("maxCombinedInitialRel"), "ratio", line, context, fields)
            emit_metric(metrics, "hpfem.combined_hpfem_seed_true_rel", fields.get("maxSeedTrueRel"), "ratio", line, context, fields)
            emit_metric(metrics, "hpfem.combined_hpfem_true_rel", fields.get("maxCombinedTrueRel"), "ratio", line, context, fields)
            emit_metric(metrics, "hpfem.combined_hpfem_vs_ldlt", fields.get("maxCombinedVsLdlt"), "ratio", line, context, fields)
        elif tag == "recycle":
            emit_metric(metrics, "hpfem.recycle_base_iters_avg", fields.get("baseItersAvg"), "count", line, context, fields)
            emit_metric(metrics, "hpfem.recycle_rec_iters_avg", fields.get("recItersAvg"), "count", line, context, fields)
            emit_metric(metrics, "hpfem.recycle_rec_iters_skip1_avg", fields.get("recItersSkip1Avg"), "count", line, context, fields)
            emit_metric(metrics, "hpfem.recycle_speedup_iters", fields.get("speedupIters"), "ratio", line, context, fields)
            emit_metric(metrics, "hpfem.recycle_base_ms_avg", fields.get("baseMsAvg"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.recycle_rec_ms_avg", fields.get("recMsAvg"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.recycle_speedup_ms", fields.get("speedupMs"), "ratio", line, context, fields)
            emit_metric(metrics, "hpfem.recycle_project_ms_avg", fields.get("projectMsAvg"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.recycle_basis_add_ms", fields.get("basisAddMs"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.recycle_max_rel_x", fields.get("maxRelXRec"), "ratio", line, context, fields)
            emit_metric(metrics, "hpfem.recycle_max_true_rel", fields.get("maxTrueRec"), "ratio", line, context, fields)
        elif tag == "line_schwarz":
            emit_metric(metrics, "hpfem.line_schwarz_base_iters", fields.get("baseIters"), "count", line, context, fields)
            emit_metric(metrics, "hpfem.line_schwarz_iters", fields.get("schwarzIters"), "count", line, context, fields)
            emit_metric(metrics, "hpfem.line_schwarz_base_ms", fields.get("baseMs"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.line_schwarz_ms", fields.get("schwarzMs"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.line_schwarz_speedup_iters", fields.get("speedupIters"), "ratio", line, context, fields)
            emit_metric(metrics, "hpfem.line_schwarz_speedup_ms", fields.get("speedupMs"), "ratio", line, context, fields)
            emit_metric(metrics, "hpfem.line_schwarz_setup_ms", fields.get("setupMs"), "ms", line, context, fields)
            emit_metric(metrics, "hpfem.line_schwarz_subdomains", fields.get("subdomains"), "count", line, context, fields)
            emit_metric(metrics, "hpfem.line_schwarz_skipped", fields.get("skipped"), "count", line, context, fields)
            emit_metric(metrics, "hpfem.line_schwarz_pcg_vs_ldlt", fields.get("pcgVsLdlt"), "ratio", line, context, fields)
            emit_metric(metrics, "hpfem.line_schwarz_true_rel", fields.get("maxTrueRel"), "ratio", line, context, fields)
        elif tag == "full_apply":
            emit_metric(metrics, "hpfem.full_apply_case_ok", 1 if fields.get("status") == "ok" else 0, "bool", line, context, fields)
            emit_metric(metrics, "hpfem.full_apply_apply_rel", fields.get("applyRel"), "ratio", line, context, fields)
            emit_metric(metrics, "hpfem.full_apply_rhs_rel", fields.get("rhsRel"), "ratio", line, context, fields)
            emit_metric(metrics, "hpfem.full_apply_u_rel", fields.get("uRel"), "ratio", line, context, fields)
            emit_metric(metrics, "hpfem.full_apply_reaction_rel", fields.get("reactionRel"), "ratio", line, context, fields)
            emit_metric(metrics, "hpfem.full_apply_reaction_abs", fields.get("reactionMaxAbs"), "force", line, context, fields)
        elif tag == "mechanism_guard":
            emit_metric(metrics, "hpfem.mechanism_guard_case_ok", 1 if fields.get("status") == "ok" else 0, "bool", line, context, fields)
            emit_metric(metrics, "hpfem.mechanism_guard_expected_singular", fields.get("expectSingular"), "bool", line, context, fields)
            emit_metric(metrics, "hpfem.mechanism_guard_prepared_singular", fields.get("preparedSingular"), "bool", line, context, fields)
            emit_metric(metrics, "hpfem.mechanism_guard_solve_singular", fields.get("solveSingular"), "bool", line, context, fields)
            emit_metric(metrics, "hpfem.mechanism_guard_pivot_margin", fields.get("pivotMargin"), "ratio", line, context, fields)
    return metrics


def parse_metrics(suite: str, stdout: str) -> list[dict[str, Any]]:
    if suite == "frame_perf":
        return parse_frame_perf(stdout)
    if suite == "hpfem":
        return parse_hpfem(stdout)
    if suite == "scale":
        return parse_scale(stdout)
    if suite == "solver_compare":
        return parse_solver_compare(stdout)
    if suite == "sparse_buckling":
        return parse_sparse_buckling(stdout)
    return []


def split_csv(value: str) -> list[str]:
    return [part.strip() for part in value.split(",") if part.strip()]


def split_coarse_bins(value: str) -> list[tuple[int, int, str]]:
    out: list[tuple[int, int, str]] = []
    for part in split_csv(value):
        match = re.fullmatch(r"(\d+)\s*x\s*(\d+)", part, flags=re.IGNORECASE)
        if not match:
            raise ValueError(f"coarse bin item must look like NxM, got: {part}")
        bx, by = int(match.group(1)), int(match.group(2))
        if bx <= 0 or by <= 0:
            raise ValueError(f"coarse bins must be positive, got: {part}")
        out.append((bx, by, f"{bx}x{by}"))
    return out or [(1, 1, "1x1")]


def canonical_size_key(key: str) -> str:
    lowered = key.strip().lower()
    if lowered in {"st", "story", "stories"}:
        return "stories"
    if lowered in {"x", "nx"}:
        return "nx"
    if lowered in {"y", "ny"}:
        return "ny"
    raise ValueError(f"unknown grid key: {key}")


def parse_size_grid(spec: str | None, default_spec: str) -> list[dict[str, int]]:
    text = (spec or default_spec).strip()
    if not text:
        return []

    if "=" not in text:
        out: list[dict[str, int]] = []
        for item in re.split(r"[;,]\s*", text):
            item = item.strip()
            if not item:
                continue
            match = re.fullmatch(r"(\d+)\s*x\s*(\d+)\s*x\s*(\d+)", item, flags=re.IGNORECASE)
            if not match:
                raise ValueError(f"grid item must look like NxNyxStories, got: {item}")
            nx, ny, stories = (int(match.group(i)) for i in range(1, 4))
            out.append({"nx": nx, "ny": ny, "stories": stories})
        return out

    base = {"nx": 3, "ny": 2, "stories": 4}
    values: dict[str, list[int]] = {}
    for chunk in text.split(";"):
        if not chunk.strip():
            continue
        if "=" not in chunk:
            raise ValueError(f"grid chunk must be key=value[,value], got: {chunk}")
        key, raw_values = chunk.split("=", 1)
        canon = canonical_size_key(key)
        parsed_values = [int(v.strip()) for v in raw_values.split(",") if v.strip()]
        if not parsed_values:
            raise ValueError(f"grid key has no values: {key}")
        values[canon] = parsed_values

    keys = list(values)
    out = []
    for combo in itertools.product(*(values[key] for key in keys)):
        item = dict(base)
        item.update({key: value for key, value in zip(keys, combo)})
        out.append(item)
    return out


def parse_suites(raw: list[str] | None) -> list[str]:
    requested = raw or ["smoke"]
    expanded: list[str] = []
    aliases = {
        "smoke": ["frame_perf", "scale"],
        "hpfem_smoke": ["hpfem"],
        "ws_b": ["scale", "solver_compare", "sparse_buckling"],
        "all": ["frame_perf", "hpfem", "scale", "solver_compare", "sparse_buckling"],
    }
    valid = {"frame_perf", "hpfem", "scale", "solver_compare", "sparse_buckling", *aliases}
    for group in requested:
        for item in split_csv(group):
            key = item.strip().lower().replace("-", "_")
            if key not in valid:
                raise ValueError(f"unknown suite '{item}'; expected one of {sorted(valid)}")
            expanded.extend(aliases.get(key, [key]))

    seen: set[str] = set()
    out: list[str] = []
    for suite in expanded:
        if suite not in seen:
            out.append(suite)
            seen.add(suite)
    return out


def size_label(prefix: str, params: dict[str, int]) -> str:
    return f"{prefix}_nx{params['nx']}_ny{params['ny']}_st{params['stories']}"


def ensure_exe(path: Path, strict: bool, warnings: list[str]) -> bool:
    if path.exists():
        return True
    msg = f"missing executable: {path}"
    if strict:
        raise FileNotFoundError(msg)
    warnings.append(msg)
    return False


def build_jobs(args: argparse.Namespace, warnings: list[str]) -> list[Job]:
    suites = parse_suites(args.suite)
    jobs: list[Job] = []

    frame_perf = args.frame_perf_exe.resolve()
    research_bin = args.research_bin.resolve()

    if "frame_perf" in suites and ensure_exe(frame_perf, args.strict_missing, warnings):
        presets = split_csv(args.frame_presets)
        if presets:
            for preset in presets:
                exe_args = ["--preset", preset, "--repeat", str(args.repeat), "--warmup", str(args.warmup)]
                if args.frame_dry:
                    exe_args.append("--dry")
                jobs.append(
                    Job(
                        suite="frame_perf",
                        label=f"frame_perf_preset_{preset}",
                        exe=frame_perf,
                        args=tuple(exe_args),
                        params={"preset": preset, "repeat": args.repeat, "warmup": args.warmup, "dry": args.frame_dry},
                    )
                )
        else:
            for params in parse_size_grid(args.frame_grid or args.grid, "3x2x4"):
                exe_args = [
                    "--nx",
                    str(params["nx"]),
                    "--ny",
                    str(params["ny"]),
                    "--stories",
                    str(params["stories"]),
                    "--repeat",
                    str(args.repeat),
                    "--warmup",
                    str(args.warmup),
                ]
                if args.frame_dry:
                    exe_args.append("--dry")
                job_params = dict(params)
                job_params.update({"repeat": args.repeat, "warmup": args.warmup, "dry": args.frame_dry})
                jobs.append(
                    Job(
                        suite="frame_perf",
                        label=size_label("frame_perf", params),
                        exe=frame_perf,
                        args=tuple(exe_args),
                        params=job_params,
                    )
                )

    if "scale" in suites:
        exe = research_bin / "exp_million_dof.exe"
        if ensure_exe(exe, args.strict_missing, warnings):
            for params in parse_size_grid(args.scale_grid or args.grid, "3x2x4"):
                jobs.append(
                    Job(
                        suite="scale",
                        label=size_label("scale", params),
                        exe=exe,
                        args=(
                            "--nx",
                            str(params["nx"]),
                            "--ny",
                            str(params["ny"]),
                            "--stories",
                            str(params["stories"]),
                        ),
                        params=dict(params),
                    )
                )

    if "hpfem" in suites:
        cases = split_csv(args.hpfem_cases)
        coarse_bins = split_coarse_bins(args.hpfem_coarse_bins)
        for case in cases:
            exe = research_bin / "exp_framecore_bsr6_matvec.exe"
            if ensure_exe(exe, args.strict_missing, warnings):
                jobs.append(
                    Job(
                        suite="hpfem",
                        label=f"framecore_bsr6_{case}",
                        exe=exe,
                        args=("--preset", case, "--repeat", str(args.hpfem_repeat)),
                        params={"experiment": "framecore_bsr6", "preset": case, "repeat": args.hpfem_repeat},
                    )
                )
            exe = research_bin / "exp_matrix_free_operator.exe"
            if ensure_exe(exe, args.strict_missing, warnings):
                pcg_max = args.hpfem_pcg_max_iter
                if pcg_max <= 0:
                    pcg_max = 5000 if case == "small" else 1500
                for precond in split_csv(args.hpfem_preconds):
                    bins_to_run = coarse_bins if "coarse" in precond else [(1, 1, "1x1")]
                    for bx, by, bins_label in bins_to_run:
                        jobs.append(
                            Job(
                                suite="hpfem",
                                label=f"matrix_free_operator_{case}_{precond}_{bins_label}",
                                exe=exe,
                                args=(
                                    "--preset",
                                    case,
                                    "--repeat",
                                    str(args.hpfem_repeat),
                                    "--pcgMaxIter",
                                    str(pcg_max),
                                    "--pcgTol",
                                    str(args.hpfem_pcg_tol),
                                    "--precond",
                                    precond,
                                    "--coarseBinsX",
                                    str(bx),
                                    "--coarseBinsY",
                                    str(by),
                                    "--maxCoarseDofs",
                                    str(args.hpfem_max_coarse_dofs),
                                ),
                                params={
                                    "experiment": "matrix_free_operator",
                                    "preset": case,
                                    "precond": precond,
                                    "coarseBins": bins_label,
                                    "repeat": args.hpfem_repeat,
                                    "pcgMaxIter": pcg_max,
                                    "pcgTol": args.hpfem_pcg_tol,
                                    "maxCoarseDofs": args.hpfem_max_coarse_dofs,
                                },
                        )
                    )
            exe = research_bin / "exp_threaded_element_apply.exe"
            if ensure_exe(exe, args.strict_missing, warnings):
                for threads in split_csv(args.hpfem_thread_counts):
                    jobs.append(
                        Job(
                            suite="hpfem",
                            label=f"threaded_apply_{case}_t{threads}",
                            exe=exe,
                            args=("--preset", case, "--repeat", str(args.hpfem_repeat), "--threads", threads),
                            params={
                                "experiment": "threaded_apply",
                                "preset": case,
                                "threads": int(threads),
                                "repeat": args.hpfem_repeat,
                            },
                        )
                    )
            exe = research_bin / "exp_parallel_pcg.exe"
            if ensure_exe(exe, args.strict_missing, warnings):
                pcg_max = args.hpfem_pcg_max_iter
                if pcg_max <= 0:
                    pcg_max = 5000 if case == "small" else 2000
                for threads in split_csv(args.hpfem_thread_counts):
                    for precond in split_csv(args.hpfem_parallel_pcg_preconds):
                        bins_to_run = coarse_bins if "coarse" in precond else [(1, 1, "1x1")]
                        for bx, by, bins_label in bins_to_run:
                            jobs.append(
                                Job(
                                    suite="hpfem",
                                    label=f"parallel_pcg_{case}_t{threads}_{precond}_{bins_label}",
                                    exe=exe,
                                    args=(
                                        "--preset",
                                        case,
                                        "--threads",
                                        threads,
                                        "--pcgMaxIter",
                                        str(pcg_max),
                                        "--pcgTol",
                                        str(args.hpfem_pcg_tol),
                                        "--precond",
                                        precond,
                                        "--coarseBinsX",
                                        str(bx),
                                        "--coarseBinsY",
                                        str(by),
                                        "--maxCoarseDofs",
                                        str(args.hpfem_max_coarse_dofs),
                                    ),
                                    params={
                                        "experiment": "parallel_pcg",
                                        "preset": case,
                                        "threads": int(threads),
                                        "precond": precond,
                                        "coarseBins": bins_label,
                                        "pcgMaxIter": pcg_max,
                                        "pcgTol": args.hpfem_pcg_tol,
                                    },
                                )
                            )
            exe = research_bin / "exp_pcg_recycling.exe"
            if ensure_exe(exe, args.strict_missing, warnings):
                pcg_max = args.hpfem_pcg_max_iter
                if pcg_max <= 0:
                    pcg_max = 5000 if case == "small" else 2000
                for basis_max in split_csv(args.hpfem_recycle_basis):
                    jobs.append(
                        Job(
                            suite="hpfem",
                            label=f"pcg_recycling_{case}_b{basis_max}",
                            exe=exe,
                            args=(
                                "--preset",
                                case,
                                "--rhs",
                                str(args.hpfem_recycle_rhs),
                                "--basisMax",
                                basis_max,
                                "--pcgMaxIter",
                                str(pcg_max),
                                "--pcgTol",
                                str(args.hpfem_pcg_tol),
                            ),
                            params={
                                "experiment": "pcg_recycling",
                                "preset": case,
                                "rhs": args.hpfem_recycle_rhs,
                                "basisMax": int(basis_max),
                                "pcgMaxIter": pcg_max,
                                "pcgTol": args.hpfem_pcg_tol,
                            },
                        )
                    )
            if args.hpfem_combined_rhs > 0:
                exe = research_bin / "exp_parallel_pcg.exe"
                if ensure_exe(exe, args.strict_missing, warnings):
                    pcg_max = args.hpfem_pcg_max_iter
                    if pcg_max <= 0:
                        pcg_max = 5000 if case == "small" else 2000
                    for threads in split_csv(args.hpfem_thread_counts):
                        for precond in split_csv(args.hpfem_combined_preconds):
                            bins_to_run = coarse_bins if "coarse" in precond else [(1, 1, "1x1")]
                            for bx, by, bins_label in bins_to_run:
                                for basis_max in split_csv(args.hpfem_combined_basis):
                                    exe_args = [
                                        "--preset",
                                        case,
                                        "--threads",
                                        threads,
                                        "--pcgMaxIter",
                                        str(pcg_max),
                                        "--pcgTol",
                                        str(args.hpfem_pcg_tol),
                                        "--precond",
                                        precond,
                                        "--coarseBinsX",
                                        str(bx),
                                        "--coarseBinsY",
                                        str(by),
                                        "--maxCoarseDofs",
                                        str(args.hpfem_max_coarse_dofs),
                                        "--rhs",
                                        str(args.hpfem_combined_rhs),
                                        "--basisMax",
                                        basis_max,
                                    ]
                                    if args.hpfem_combined_seed_load_basis:
                                        exe_args.append("--seedLoadBasis")
                                    seed_suffix = "_seed" if args.hpfem_combined_seed_load_basis else ""
                                    jobs.append(
                                        Job(
                                            suite="hpfem",
                                            label=f"combined_hpfem_{case}_t{threads}_{precond}_{bins_label}_b{basis_max}_r{args.hpfem_combined_rhs}{seed_suffix}",
                                            exe=exe,
                                            args=tuple(exe_args),
                                            params={
                                                "experiment": "combined_hpfem",
                                                "preset": case,
                                                "threads": int(threads),
                                                "precond": precond,
                                                "coarseBins": bins_label,
                                                "rhs": args.hpfem_combined_rhs,
                                                "basisMax": int(basis_max),
                                                "seedLoadBasis": args.hpfem_combined_seed_load_basis,
                                                "pcgMaxIter": pcg_max,
                                                "pcgTol": args.hpfem_pcg_tol,
                                            },
                                        )
                                    )
            exe = research_bin / "exp_line_schwarz_pcg.exe"
            if ensure_exe(exe, args.strict_missing, warnings):
                pcg_max = args.hpfem_pcg_max_iter
                if pcg_max <= 0:
                    pcg_max = 5000 if case == "small" else 3000
                for mode in split_csv(args.hpfem_line_schwarz_modes):
                    exe_args = [
                        "--preset",
                        case,
                        "--pcgMaxIter",
                        str(pcg_max),
                        "--pcgTol",
                        str(args.hpfem_pcg_tol),
                    ]
                    if mode == "columns":
                        exe_args.append("--columnsOnly")
                    elif mode == "floor":
                        exe_args.append("--floorLines")
                    else:
                        raise ValueError(f"unknown line schwarz mode: {mode}")
                    jobs.append(
                        Job(
                            suite="hpfem",
                            label=f"line_schwarz_{case}_{mode}",
                            exe=exe,
                            args=tuple(exe_args),
                            params={
                                "experiment": "line_schwarz",
                                "preset": case,
                                "mode": mode,
                                "pcgMaxIter": pcg_max,
                                "pcgTol": args.hpfem_pcg_tol,
                            },
                        )
                    )
        exe = research_bin / "exp_full_apply_oracle.exe"
        if ensure_exe(exe, args.strict_missing, warnings):
            jobs.append(
                Job(
                    suite="hpfem",
                    label="full_apply_oracle",
                    exe=exe,
                    args=(),
                    params={"experiment": "full_apply_oracle"},
                )
            )
        exe = research_bin / "exp_mechanism_guard_oracle.exe"
        if ensure_exe(exe, args.strict_missing, warnings):
            jobs.append(
                Job(
                    suite="hpfem",
                    label="mechanism_guard_oracle",
                    exe=exe,
                    args=(),
                    params={"experiment": "mechanism_guard_oracle"},
                )
            )

    if "solver_compare" in suites:
        exe = research_bin / "exp_solver_compare.exe"
        if ensure_exe(exe, args.strict_missing, warnings):
            exe_args = () if args.include_large else ("--noLarge",)
            jobs.append(
                Job(
                    suite="solver_compare",
                    label="solver_compare_large" if args.include_large else "solver_compare_small",
                    exe=exe,
                    args=exe_args,
                    params={"include_large": args.include_large},
                )
            )

    if "sparse_buckling" in suites:
        exe = research_bin / "exp_sparse_buckling.exe"
        if ensure_exe(exe, args.strict_missing, warnings):
            exe_args = () if args.include_large else ("--noLarge",)
            jobs.append(
                Job(
                    suite="sparse_buckling",
                    label="sparse_buckling_large" if args.include_large else "sparse_buckling_small",
                    exe=exe,
                    args=exe_args,
                    params={"include_large": args.include_large},
                )
            )

    return jobs


def git_snapshot() -> dict[str, Any]:
    def run_git(*git_args: str) -> str:
        try:
            proc = subprocess.run(
                ["git", "-C", str(ROOT), *git_args],
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=5,
            )
        except Exception:
            return ""
        return proc.stdout.strip() if proc.returncode == 0 else ""

    status = run_git("status", "--short")
    return {
        "commit": run_git("rev-parse", "--short", "HEAD"),
        "branch": run_git("rev-parse", "--abbrev-ref", "HEAD"),
        "dirty_count": len([line for line in status.splitlines() if line.strip()]),
    }


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8", errors="replace")


def run_one(
    job: Job,
    iteration: int,
    run_index: int,
    args: argparse.Namespace,
    out_dir: Path,
    git_info: dict[str, Any],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    run_id = f"run_{run_index:04d}_{job.label}_i{iteration + 1}"
    cmd = [str(job.exe), *job.args]
    started_at = iso_now()
    t0 = time.perf_counter()
    timeout_hit = False
    stdout = ""
    stderr = ""
    returncode: int | None

    try:
        proc = subprocess.run(
            cmd,
            cwd=str(ROOT),
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=args.timeout_sec,
        )
        returncode = proc.returncode
        stdout = proc.stdout
        stderr = proc.stderr
    except subprocess.TimeoutExpired as exc:
        timeout_hit = True
        returncode = None
        stdout = exc.stdout or ""
        stderr = exc.stderr or ""

    elapsed_ms = (time.perf_counter() - t0) * 1000.0

    stdout_log = ""
    stderr_log = ""
    if not args.no_raw_logs:
        log_dir = out_dir / "logs"
        stdout_path = log_dir / f"{run_id}.stdout.txt"
        stderr_path = log_dir / f"{run_id}.stderr.txt"
        write_text(stdout_path, stdout)
        write_text(stderr_path, stderr)
        stdout_log = str(stdout_path)
        stderr_log = str(stderr_path)

    parsed = parse_metrics(job.suite, stdout)
    ok = (returncode == 0) and not timeout_hit
    run_record = {
        "type": "run",
        "schema": SCHEMA_VERSION,
        "run_id": run_id,
        "suite": job.suite,
        "label": job.label,
        "params": job.params,
        "iteration": iteration + 1,
        "cmd": cmd,
        "cmd_text": command_text(cmd),
        "cwd": str(ROOT),
        "started_at": started_at,
        "elapsed_ms": elapsed_ms,
        "returncode": returncode,
        "ok": ok,
        "timeout": timeout_hit,
        "stdout_log": stdout_log,
        "stderr_log": stderr_log,
        "host": {
            "node": platform.node(),
            "platform": platform.platform(),
            "python": sys.version.split()[0],
        },
        "git": git_info,
    }
    if args.embed_raw:
        run_record["stdout"] = stdout
        run_record["stderr"] = stderr

    metric_records: list[dict[str, Any]] = []
    for metric in parsed:
        metric_records.append(
            {
                "type": "metric",
                "schema": SCHEMA_VERSION,
                "run_id": run_id,
                "suite": job.suite,
                "label": job.label,
                "params": job.params,
                "iteration": iteration + 1,
                **metric,
            }
        )

    return run_record, metric_records


def geomean(values: list[float]) -> float | None:
    positives = [value for value in values if value > 0 and math.isfinite(value)]
    if not positives:
        return None
    return math.exp(sum(math.log(value) for value in positives) / len(positives))


def aggregate(metric_records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[tuple[str, str, str, str, str], list[float]] = {}
    for record in metric_records:
        value = record.get("value")
        if not finite_number(value):
            continue
        key = (
            record.get("suite", ""),
            record.get("label", ""),
            context_key(record.get("context", {})),
            record.get("metric", ""),
            record.get("unit", ""),
        )
        groups.setdefault(key, []).append(float(value))

    rows: list[dict[str, Any]] = []
    for (suite, label, ctx, metric, unit), values in sorted(groups.items()):
        values_sorted = sorted(values)
        gm = geomean(values_sorted)
        rows.append(
            {
                "suite": suite,
                "label": label,
                "context": ctx,
                "metric": metric,
                "unit": unit,
                "count": len(values_sorted),
                "median": statistics.median(values_sorted),
                "geomean": gm,
                "min": values_sorted[0],
                "max": values_sorted[-1],
            }
        )
    return rows


def write_summary_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["suite", "label", "context", "metric", "unit", "count", "median", "geomean", "min", "max"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def fmt_num(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        if abs(value) >= 1000:
            return f"{value:.1f}"
        if abs(value) >= 10:
            return f"{value:.2f}"
        if abs(value) >= 0.01:
            return f"{value:.3f}"
        return f"{value:.3e}"
    return str(value)


def primary_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    selected = [row for row in rows if row["metric"] in PRIMARY_METRICS]
    return sorted(selected, key=lambda row: (row["suite"], row["metric"], float(row["median"])))


def write_summary_md(path: Path, rows: list[dict[str, Any]], jsonl_path: Path, summary_csv: Path) -> None:
    selected = primary_rows(rows)
    lines = [
        "# HP Tuning Summary",
        "",
        f"- generated_at: {iso_now()}",
        f"- jsonl: {jsonl_path}",
        f"- csv: {summary_csv}",
        "",
        "## Primary Metrics",
        "",
        "| suite | label | context | metric | n | median | geomean | unit |",
        "|---|---|---|---|---:|---:|---:|---|",
    ]
    for row in selected:
        lines.append(
            "| {suite} | {label} | {context} | {metric} | {count} | {median} | {geomean} | {unit} |".format(
                suite=row["suite"],
                label=row["label"],
                context=row["context"],
                metric=row["metric"],
                count=row["count"],
                median=fmt_num(row["median"]),
                geomean=fmt_num(row["geomean"]),
                unit=row["unit"],
            )
        )
    write_text(path, "\n".join(lines) + "\n")


def print_primary_summary(rows: list[dict[str, Any]]) -> None:
    selected = primary_rows(rows)
    if not selected:
        print("[summary] no primary metrics parsed")
        return
    print("[summary] primary metrics:")
    for row in selected[:24]:
        ctx = f" {row['context']}" if row["context"] else ""
        print(
            f"  {row['suite']} {row['label']}{ctx} "
            f"{row['metric']} n={row['count']} median={fmt_num(row['median'])} "
            f"geomean={fmt_num(row['geomean'])} {row['unit']}"
        )


def make_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run existing FrameCore benchmark/research executables and collect JSONL tuning data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--suite",
        action="append",
        help="Suite(s): smoke, frame_perf, scale, solver_compare, sparse_buckling, ws_b, all. Comma-separated or repeatable.",
    )
    parser.add_argument("--grid", help="Size grid for frame_perf/scale, e.g. '3x2x4,5x4x8' or 'nx=3,5;ny=2;stories=4,8'.")
    parser.add_argument("--frame-grid", help="Override size grid for frame_perf.")
    parser.add_argument("--scale-grid", help="Override size grid for exp_million_dof.")
    parser.add_argument("--frame-presets", default="", help="Comma-separated frame_perf presets, e.g. xxl,mega. Overrides frame grid.")
    parser.add_argument("--runs", type=int, default=1, help="Process-level repeats per job.")
    parser.add_argument("--repeat", type=int, default=3, help="frame_perf internal repeat count.")
    parser.add_argument("--warmup", type=int, default=1, help="frame_perf internal warmup count.")
    parser.add_argument("--frame-dry", action="store_true", help="Pass --dry to frame_perf.")
    parser.add_argument("--hpfem-cases", default="small,xxl", help="Comma-separated cases for hpfem suite: small,xxl,mega.")
    parser.add_argument("--hpfem-repeat", type=int, default=20, help="Internal apply repeat count for hpfem experiments.")
    parser.add_argument("--hpfem-preconds", default="diag,block6,diag_coarse,block6_coarse", help="Comma-separated matrix-free preconditioner modes.")
    parser.add_argument("--hpfem-coarse-bins", default="1x1,2x2", help="Comma-separated coarse aggregate grids for *_coarse preconditioners, e.g. 1x1,2x2,3x3.")
    parser.add_argument("--hpfem-max-coarse-dofs", type=int, default=4096, help="Guardrail for dense coarse inverse size in matrix-free experiments.")
    parser.add_argument("--hpfem-thread-counts", default="4", help="Comma-separated thread counts for exp_threaded_element_apply.")
    parser.add_argument("--hpfem-parallel-pcg-preconds", default="diag", help="Comma-separated preconditioners for exp_parallel_pcg: diag,block6,block6_coarse.")
    parser.add_argument("--hpfem-recycle-basis", default="0,4,8", help="Comma-separated basisMax values for exp_pcg_recycling.")
    parser.add_argument("--hpfem-recycle-rhs", type=int, default=12, help="Number of related RHS vectors for exp_pcg_recycling.")
    parser.add_argument("--hpfem-combined-rhs", type=int, default=0, help="Enable combined_hpfem jobs with this many RHS; 0 disables.")
    parser.add_argument("--hpfem-combined-basis", default="8", help="Comma-separated basisMax values for combined_hpfem jobs.")
    parser.add_argument("--hpfem-combined-preconds", default="block6_coarse", help="Comma-separated preconditioners for combined_hpfem jobs.")
    parser.add_argument("--hpfem-combined-seed-load-basis", action="store_true", help="Seed combined_hpfem with solved parametric load modes before the RHS sequence.")
    parser.add_argument("--hpfem-line-schwarz-modes", default="columns,floor", help="Comma-separated line Schwarz modes: columns,floor.")
    parser.add_argument("--hpfem-pcg-max-iter", type=int, default=0, help="PCG max iterations for hpfem matrix-free experiment; 0 selects per-case defaults.")
    parser.add_argument("--hpfem-pcg-tol", type=float, default=1e-10, help="PCG tolerance for hpfem matrix-free experiment.")
    parser.add_argument("--include-large", action="store_true", help="Allow large research cases; default passes --noLarge where supported.")
    parser.add_argument("--timeout-sec", type=float, default=600.0, help="Per-process timeout.")
    parser.add_argument("--out-dir", type=Path, help="Output directory. Defaults to Research/out/hp_tuning/<timestamp>.")
    parser.add_argument("--frame-perf-exe", type=Path, default=DEFAULT_FRAME_PERF, help="Path to frame_perf.exe.")
    parser.add_argument("--research-bin", type=Path, default=DEFAULT_RESEARCH_BIN, help="Directory containing research exe files.")
    parser.add_argument("--strict-missing", action="store_true", help="Fail if an expected executable is missing instead of skipping it.")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without executing them.")
    parser.add_argument("--no-raw-logs", action="store_true", help="Do not write per-run stdout/stderr log files.")
    parser.add_argument("--embed-raw", action="store_true", help="Embed stdout/stderr text in the JSONL run records.")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = make_arg_parser()
    args = parser.parse_args(argv)
    if args.runs < 1:
        parser.error("--runs must be >= 1")
    if args.repeat < 1:
        parser.error("--repeat must be >= 1")
    if args.warmup < 0:
        parser.error("--warmup must be >= 0")

    warnings: list[str] = []
    try:
        jobs = build_jobs(args, warnings)
    except (ValueError, FileNotFoundError) as exc:
        parser.error(str(exc))

    for warning in warnings:
        print(f"[warn] {warning}")

    if not jobs:
        print("[error] no jobs to run")
        return 2

    for job in jobs:
        print(f"[job] {job.suite} {job.label}: {command_text([str(job.exe), *job.args])}")

    if args.dry_run:
        return 0

    out_dir = args.out_dir.resolve() if args.out_dir else (DEFAULT_OUT_ROOT / now_stamp()).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    jsonl_path = out_dir / "results.jsonl"
    summary_csv = out_dir / "summary.csv"
    summary_md = out_dir / "summary.md"
    git_info = git_snapshot()

    all_metric_records: list[dict[str, Any]] = []
    run_index = 1
    with jsonl_path.open("w", encoding="utf-8") as jsonl:
        for job in jobs:
            for iteration in range(args.runs):
                print(f"[run] {job.label} iteration={iteration + 1}/{args.runs}")
                run_record, metric_records = run_one(job, iteration, run_index, args, out_dir, git_info)
                jsonl.write(json.dumps(run_record, ensure_ascii=False, sort_keys=True) + "\n")
                for metric_record in metric_records:
                    jsonl.write(json.dumps(metric_record, ensure_ascii=False, sort_keys=True) + "\n")
                all_metric_records.extend(metric_records)
                run_index += 1
                status = "ok" if run_record["ok"] else f"rc={run_record['returncode']}"
                print(f"[done] {run_record['run_id']} {status} elapsedMs={run_record['elapsed_ms']:.1f} metrics={len(metric_records)}")

    rows = aggregate(all_metric_records)
    write_summary_csv(summary_csv, rows)
    write_summary_md(summary_md, rows, jsonl_path, summary_csv)

    print(f"[out] jsonl={jsonl_path}")
    print(f"[out] summary_csv={summary_csv}")
    print(f"[out] summary_md={summary_md}")
    print_primary_summary(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
