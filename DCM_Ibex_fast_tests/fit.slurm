#!/bin/bash
# One independent model fit; CPUs are overridden for its actual case.
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
set -euo pipefail
: "${DCM_ROOT:?Missing project root}"
: "${EXPERIMENT_DIR:?Missing experiment directory}"
: "${CASE_ID:?Missing case}"
cd "$DCM_ROOT"
source "$EXPERIMENT_DIR/ibex_environment.sh"
export LC_ALL=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
date -u
hostname
echo "Exploratory case $CASE_ID; CPUs $SLURM_CPUS_PER_TASK"
srun --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
  Rscript --vanilla "$EXPERIMENT_DIR/02_fit_case.R"
