#!/bin/bash
# Recover checking only; repeating resumes the same accepted registrations.
set -euo pipefail
export DCM_ROOT="${DCM_ROOT:-/ibex/user/meleg/DCM_MODEL}"
cd "$DCM_ROOT"
package=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
export EXPERIMENT_DIR="${ONSET_EXPERIMENT:-$(cat logs/latest_onset150_experiment.txt)}"
export ONSET_REPAIR_DIR="$EXPERIMENT_DIR/repairs/onset_support_r1"
mkdir -p "$ONSET_REPAIR_DIR/scheduler"
exec 9>"$ONSET_REPAIR_DIR/submission.lock"
flock -n 9 || { echo 'Another recovery submission is active.' >&2; exit 1; }
source "$EXPERIMENT_DIR/ibex_environment.sh"
(cd "$package"; sha256sum --check --quiet SHA256SUMS.txt)
for file in "$package/"*; do
  [[ -f "$file" ]] || continue
  name=$(basename "$file")
  if [[ -f "$ONSET_REPAIR_DIR/$name" ]]; then
    cmp --silent "$file" "$ONSET_REPAIR_DIR/$name" || { echo "Different recovery version: $name" >&2; exit 1; }
  else
    cp -- "$file" "$ONSET_REPAIR_DIR/$name"
  fi
done
for file in "$ONSET_REPAIR_DIR/"*.sh "$ONSET_REPAIR_DIR/"*.slurm; do bash -n "$file"; done
# Test new scoring on fixtures and every real held-out cast before patching.
Rscript --vanilla "$ONSET_REPAIR_DIR/test_and_audit.R" "$ONSET_REPAIR_DIR" "$EXPERIMENT_DIR"
bash "$ONSET_REPAIR_DIR/patch_files.sh"
Rscript --vanilla "$ONSET_REPAIR_DIR/submit_repair.R" "$DCM_ROOT" "$ONSET_REPAIR_DIR" "$EXPERIMENT_DIR"
printf '\nRecovery directory: %s\nSupport audit: %s/actual_cohort_support_audit.csv\n' "$ONSET_REPAIR_DIR" "$ONSET_REPAIR_DIR"
