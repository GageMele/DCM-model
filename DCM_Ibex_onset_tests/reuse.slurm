#!/bin/bash
# Diagnostics and figures run as a separate job after its fit has ended.
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=16:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
set -euo pipefail
cd "${DCM_ROOT:?Missing project root}"
source "${EXPERIMENT_DIR:?Missing experiment directory}/ibex_environment.sh"
export LC_ALL=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
srun --ntasks=1 --cpus-per-task=1 --cpu-bind=cores \
  Rscript --vanilla "$EXPERIMENT_DIR/reuse_check.R"
