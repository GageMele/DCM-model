#!/bin/bash
# Refresh completed results; this never starts a model fit.
set -euo pipefail
cd "${DCM_ROOT:-/ibex/user/meleg/DCM_MODEL}"
export EXPERIMENT_DIR="$(cat logs/latest_onset150_experiment.txt)"
source "$EXPERIMENT_DIR/ibex_environment.sh"
Rscript --vanilla "$EXPERIMENT_DIR/04_collect_results.R" "$EXPERIMENT_DIR"
