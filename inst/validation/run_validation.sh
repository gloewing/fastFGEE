#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
OUT=${FGEE_VALIDATION_OUT:-"$HERE/results"}
mkdir -p "$OUT/logs" "$OUT/resources"
export FGEE_VALIDATION_OUT="$OUT"

run_one() {
  local script=$1
  local stem
  stem=$(basename "$script" .R)
  echo "== $stem =="
  if command -v /usr/bin/time >/dev/null 2>&1; then
    /usr/bin/time -v -o "$OUT/resources/${stem}.time.txt" \
      Rscript "$HERE/$script" >"$OUT/logs/${stem}.log" 2>&1
  else
    Rscript "$HERE/$script" >"$OUT/logs/${stem}.log" 2>&1
  fi
  tail -20 "$OUT/logs/${stem}.log"
}

run_one run_golden_checks.R
run_one stress_working_stats.R
run_one benchmark_correlation_operators.R
run_one benchmark_working_stats.R
run_one benchmark_tuning.R

if [[ "${FGEE_RUN_ENGINE_BENCHMARK:-1}" == "1" ]]; then
  run_one benchmark_fgee_engines.R
fi

echo "Validation outputs: $OUT"
