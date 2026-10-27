#!/bin/bash
# Gather all cases, including failed or timed-out fits, after checks have ended.
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
set -euo pipefail
cd "${DCM_ROOT:?Missing project root}"
source "${EXPERIMENT_DIR:?Missing experiment directory}/ibex_environment.sh"
Rscript --vanilla "$EXPERIMENT_DIR/04_collect_results.R" "$EXPERIMENT_DIR"
