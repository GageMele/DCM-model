#!/bin/bash
# Export all comparison tables while leaving plots and large raw files on Ibex.
set -euo pipefail
cd "${DCM_ROOT:-/ibex/user/meleg/DCM_MODEL}"
experiment="${ONSET_EXPERIMENT:-$(cat logs/latest_onset150_experiment.txt)}"
archive="$PWD/logs/dcm_onset_compact_$(date -u +%Y%m%dT%H%M%SZ)_$$.tar.gz"
list=$(mktemp)
trap 'rm -f -- "$list"' EXIT
find "$experiment/results" -maxdepth 1 -type f -name '*.csv' -print0 > "$list"
find experiments benchmarks \
  -type d \( -name stan_csv -o -name repairs -o -name 'scheduler_backup_*' \) -prune -o \
  -type f \( -name cases.tsv -o -name manifest.tsv -o -name source_controls.tsv \
    -o -name existing_evaluations.tsv -o -name test_summary.csv \
    -o -name diagnostic_summary.csv -o -name run_config.txt \
    -o -name CHECKS_FAILED.txt -o -name FIT_FAILED.txt -o -name CHECKS_COMPLETE.txt \) -print0 >> "$list"
find "$experiment/repairs/onset_support_r1" -maxdepth 1 -type f \
  \( -name actual_cohort_support_audit.csv -o -name repair_manifest.tsv -o -name README.md \) -print0 >> "$list"
find "$experiment/runs" -maxdepth 3 -type f -name parameter_diagnostics.csv -print0 >> "$list"
for name in campaign.R base_models.R data_identity_audit.csv data_checksums.csv; do
  [[ ! -f "$experiment/$name" ]] || printf '%s\0' "$experiment/$name" >> "$list"
done
status=0
tar --null --files-from="$list" -czf "$archive" || status=$?
(( status<=1 )) || exit "$status"
printf '\nCompact comparison archive:\n%s\n' "$archive"
ls -lh -- "$archive"
