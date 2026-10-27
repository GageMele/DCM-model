#!/bin/bash
# Prepare nested subsets and validate both prior choices on every dataset.
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
set -euo pipefail
cd "${DCM_ROOT:?Missing project root}"
source "${EXPERIMENT_DIR:?Missing experiment directory}/ibex_environment.sh"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
srun --ntasks=1 --cpus-per-task=1 --cpu-bind=cores \
  Rscript --vanilla "$EXPERIMENT_DIR/01_prepare.R"
