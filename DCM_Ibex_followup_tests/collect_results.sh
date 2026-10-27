#!/bin/bash
# Refresh compact all-campaign comparisons without fitting or prediction.
set -euo pipefail
export LC_ALL=C
cd -- "${DCM_ROOT:-/ibex/user/meleg/DCM_MODEL}"
experiment=$(cat logs/latest_followup150_experiment.txt)
source "$experiment/ibex_environment.sh"
Rscript --vanilla "$experiment/04_collect_results.R" "$experiment"
