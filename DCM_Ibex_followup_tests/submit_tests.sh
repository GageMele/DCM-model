#!/bin/bash
# Resume the same experiment on repeat; explicitly set FOLLOWUP_NEW=1 for a new one.
set -euo pipefail
export LC_ALL=C
export DCM_ROOT="${DCM_ROOT:-/ibex/user/meleg/DCM_MODEL}"
package=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd -- "$DCM_ROOT"
mkdir -p logs
command -v flock >/dev/null || { echo 'Need flock to prevent concurrent submissions.' >&2; exit 1; }
exec 9>logs/followup150_submission.lock
flock -n 9 || { echo 'Another follow-up submission is active; no jobs submitted.' >&2; exit 1; }
source_experiment="${FOLLOWUP_SOURCE:-$(cat logs/latest_targeted150_experiment.txt)}"
source "$source_experiment/ibex_environment.sh"
for file in "$package/"*.sh "$package/"*.slurm; do bash -n "$file"; done
(cd -- "$package"; sha256sum --check --quiet SHA256SUMS.txt)
Rscript --vanilla "$package/test_offline.R" "$package"
Rscript --vanilla "$package/submit.R" "$DCM_ROOT" "$package"
