#!/bin/bash
# Export compact comparisons; leave all original fitted models in place.
set -euo pipefail
root="${DCM_ROOT:-/ibex/user/meleg/DCM_MODEL}"
experiment=$(cat "$root/logs/latest_shape150_experiment.txt")
stamp=$(date -u +%Y%m%dT%H%M%SZ)
archive="$root/logs/dcm_shape_comparison_${stamp}_$$.tar.gz"
list=$(mktemp)
trap 'rm -f -- "$list"' EXIT
cd -- "$root"
relative=${experiment#"$root/"}
find "$relative" \
  \( -type d \( -name stan_csv -o -name data -o -name .git \) -prune \) -o \
  \( -type f \( -name '*.csv' -o -name '*.tsv' -o -name '*.txt' -o -name '*.json' \
       -o -name '*.R' -o -name '*.md' -o -name '*.png' -o -name '*.stan' \) \
       ! -name 'sampler_iterations_including_warmup.csv' -size -20971521c -print0 \) > "$list"
code=0
tar --null --files-from="$list" -czf "$archive" || code=$?
(( code <= 1 )) || exit "$code"
printf '\nComparison archive:\n%s\n' "$archive"
ls -lh -- "$archive"
