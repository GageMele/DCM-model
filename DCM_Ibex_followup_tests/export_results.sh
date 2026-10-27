#!/bin/bash
# Snapshot compact results from every campaign; leave original files in place.
set -euo pipefail
cd -- "${DCM_ROOT:-/ibex/user/meleg/DCM_MODEL}"
mkdir -p logs
stamp=$(date -u +%Y%m%dT%H%M%SZ)
archive="$PWD/logs/dcm_followup_comparison_${stamp}_$$.tar.gz"
file_list=$(mktemp)
trap 'rm -f -- "$file_list"' EXIT
roots=()
for directory in experiments benchmarks runs DCM_Ibex_targeted_tests/original_review; do
  [[ -d "$directory" ]] && roots+=("$directory")
done
(( ${#roots[@]} > 0 )) || { echo "No results directories found."; exit 1; }
find "${roots[@]}" \
  \( -type d \( -name stan_csv -o -name .git \) -prune \) -o \
  \( -type f \( -name '*.csv' -o -name '*.tsv' -o -name '*.txt' -o -name '*.md' -o -name '*.png' -o -name '*.pdf' -o -name '*.stan' -o -name '*.R' -o -name 'run_config.rds' \) \
    ! -name '*_observations.csv' ! -name '*predictive_exceedance.csv' \
    ! -name 'sampler_iterations_including_warmup.csv' ! -name 'sampler_diagnostics.csv' \
    -size -20971521c -print0 \) > "$file_list"
[[ -s "$file_list" ]] || { echo "No compact results found."; exit 1; }
status=0
tar --null --files-from="$file_list" -czf "$archive" || status=$?
(( status <= 1 )) || exit "$status"
if (( status == 1 )); then echo "Some active files changed while collecting; completed files remain included."; fi
printf '\nComparison archive:\n%s\n' "$archive"
ls -lh -- "$archive"
