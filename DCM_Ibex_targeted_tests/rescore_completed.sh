#!/bin/bash
# Add prediction checks for earlier tests that finished since the campaign launch.
set -euo pipefail
export LC_ALL=C
root="${DCM_ROOT:-/ibex/user/meleg/DCM_MODEL}"
cd -- "$root"
experiment="${EXPERIMENT_DIR:-$(cat logs/latest_targeted150_experiment.txt)}"
source "$experiment/ibex_environment.sh"
Rscript --vanilla "$experiment/submit_existing.R" "$experiment" "$root"
