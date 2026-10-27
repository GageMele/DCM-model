#!/bin/bash
# Queue the validated matrix and legacy prediction checks in a fresh experiment.
set -euo pipefail
export LC_ALL=C
export DCM_ROOT="${DCM_ROOT:-/ibex/user/meleg/DCM_MODEL}"
package=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd -- "$DCM_ROOT"
source_experiment="${SOURCE_EXPERIMENT:-$(cat logs/latest_extensive150_experiment.txt)}"
source "$source_experiment/ibex_environment.sh"
for file in "$package/"*.sh "$package/"*.slurm; do bash -n "$file"; done
Rscript --vanilla "$package/test_submission.R" "$package"
Rscript --vanilla "$package/test_reporting.R" "$package"
Rscript --vanilla "$package/test_existing.R" "$package"
Rscript --vanilla "$package/test_dependency_recovery.R" "$package"
Rscript --vanilla "$package/submit.R" "$DCM_ROOT" "$package"
