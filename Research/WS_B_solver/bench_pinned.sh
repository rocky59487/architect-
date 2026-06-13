#!/usr/bin/env bash
# bench_pinned.sh -- stable repeated/interleaved benchmark for the HP solver.
#
# Single runs of exp_parallel_pcg are noisy (the one-time LDLT factor and machine
# drift dominate short batches). This harness runs each labelled config REPS times,
# INTERLEAVED round-robin so machine drift hits every config equally, then reports
# median / min / max of the chosen metric. Use it to compare OLD vs NEW honestly.
#
# Usage:  bench_pinned.sh [REPS] [METRIC]
#   REPS    repeats per config (default 7)
#   METRIC  printf token to extract, e.g. factorBypassBatchSpeedup (default) or
#           combinedMsAvg, deflItersAvg, speedupMs, ...
#
# Edit the CONFIGS array below to choose what to compare.
set -u
cd "$(dirname "$0")/../.." || exit 1
EXE=./Research/bin/exp_parallel_pcg.exe
REPS="${1:-7}"
METRIC="${2:-factorBypassBatchSpeedup}"

COMMON="--preset xxl --threads 16 --precond block6_coarse --coarseBinsX 2 --coarseBinsY 2 --pcgMaxIter 2000 --pcgTol 1e-10"

# label | extra args   (one per line; '|' separated)
CONFIGS=(
  "nonseeded_OLD|--rhs 8 --basisMax 8"
  "nonseeded_NEW|--rhs 8 --basisMax 8 --parallelPrecond --coarseSolve banded"
  "seeded32_OLD|--rhs 32 --basisMax 8 --seedLoadBasis"
  "seeded32_NEW|--rhs 32 --basisMax 8 --seedLoadBasis --parallelPrecond --coarseSolve banded"
)

declare -A SAMPLES
for c in "${CONFIGS[@]}"; do SAMPLES["${c%%|*}"]=""; done

echo "# bench_pinned: REPS=$REPS METRIC=$METRIC  ($(date '+%H:%M:%S'))"
for ((r=1; r<=REPS; r++)); do
  for c in "${CONFIGS[@]}"; do
    label="${c%%|*}"; extra="${c#*|}"
    val=$($EXE $COMMON $extra 2>/dev/null | tr ' ' '\n' | grep -oE "${METRIC}=[0-9.eE+-]+" | head -1 | cut -d= -f2)
    SAMPLES["$label"]="${SAMPLES["$label"]} $val"
    printf "  rep%-2d %-16s %s=%s\n" "$r" "$label" "$METRIC" "${val:-NA}"
  done
done

echo "# summary (median / min / max of $METRIC over $REPS reps)"
for c in "${CONFIGS[@]}"; do
  label="${c%%|*}"
  vals=$(echo "${SAMPLES["$label"]}" | tr ' ' '\n' | grep -E '^[0-9]' | sort -g)
  n=$(echo "$vals" | wc -l)
  if [ "$n" -gt 0 ]; then
    med=$(echo "$vals" | awk -v n="$n" 'NR==int((n+1)/2){print}')
    mn=$(echo "$vals" | head -1); mx=$(echo "$vals" | tail -1)
    printf "  %-16s median=%-10s min=%-10s max=%-10s n=%s\n" "$label" "$med" "$mn" "$mx" "$n"
  fi
done
