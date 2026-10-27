#!/bin/bash
# Collect all saved campaigns; never submits or reruns a model.
set -euo pipefail
root="${DCM_ROOT:-/ibex/user/meleg/DCM_MODEL}"
scripts=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd -- "$root"
if ! command -v Rscript >/dev/null 2>&1; then
    module load "${R_MODULE:-R/4.5.0/gnu-12.2.0}"
fi
stamp=$(date -u +%Y%m%dT%H%M%SZ)
name="dcm_all_run_review_${stamp}_$$"
Rscript --vanilla "$scripts/collect_all_runs.R" "$root" "$root/logs/$name"
tar -czf "$root/logs/$name.tar.gz" -C "$root/logs" "$name"
printf '\nAll-run review archive: %s\n' "$root/logs/$name.tar.gz"
printf 'Upload this archive for the full cross-campaign comparison.\n'
