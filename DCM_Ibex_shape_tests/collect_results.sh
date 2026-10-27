#!/bin/bash
# Update summaries without refitting any model.
set -euo pipefail
root="${DCM_ROOT:-/ibex/user/meleg/DCM_MODEL}"
experiment=$(cat "$root/logs/latest_shape150_experiment.txt")
source "$experiment/ibex_environment.sh"
Rscript --vanilla "$experiment/04_collect_results.R" "$experiment"
