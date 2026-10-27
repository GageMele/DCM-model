#!/bin/bash
# This allocation prepares the subsets; no data processing runs on the login node.
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
set -euo pipefail
: "${DCM_ROOT:?Missing project root}"
: "${EXPERIMENT_DIR:?Missing experiment directory}"
cd "$DCM_ROOT"
source "$EXPERIMENT_DIR/ibex_environment.sh"
export LC_ALL=C
srun --ntasks=1 --cpus-per-task=1 --cpu-bind=cores \
  Rscript --vanilla "$EXPERIMENT_DIR/01_prepare_tests.R" "$EXPERIMENT_DIR"
