#!/bin/bash
# Read-only accounting, progress and live sampler diagnostics.
set -euo pipefail
cd "${DCM_ROOT:-/ibex/user/meleg/DCM_MODEL}"
experiment="${1:-$(cat logs/latest_k150_experiment.txt)}"
date -Is
prepare_job=$(cat "$experiment/preparation_job_id.txt")
sacct -X -j "$prepare_job" --format=JobID,JobName%24,State,ExitCode,Elapsed
ids=$(awk 'NR>1 {print $3}' "$experiment/manifest.tsv" | paste -sd, -)
if [[ -n "$ids" ]]; then
  sacct -X -j "$ids" --format=JobID,JobName%24,State,ExitCode,Elapsed,Timelimit
fi
if [[ ! -f "$experiment/PREPARATION_COMPLETE.txt" ]]; then
  echo "Preparation not yet complete. Its log is logs/dcm-k150-prepare-${prepare_job}.out/.err"
fi
shopt -s nullglob
printf '%4s %5s %-8s %6s %3s %3s %6s %7s %8s %5s %10s\n' \
  case chain phase draws cap max hits 'hits%' 'recent%' div mean_leaps
while IFS=$'\t' read -r case_id label fit_job check_job cpus run_dir; do
  [[ "$case_id" =~ ^[0-9]+$ ]] || continue
  files=("$run_dir"/stan_csv/*.csv)
  for path in "${files[@]}"; do
    awk -F, -v case_id="$case_id" -f "$experiment/live_diagnostics.awk" "$path"
  done
  if [[ -f "$run_dir/FIT_FAILED.txt" ]]; then
    echo "Case $case_id fitting error: $(cat "$run_dir/FIT_FAILED.txt")"
  fi
done < "$experiment/manifest.tsv"
