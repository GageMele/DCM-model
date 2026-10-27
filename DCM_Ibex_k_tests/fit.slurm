#!/bin/bash
# Fit one frozen case; CPU, memory and walltime are set by the submit helper.
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
set -euo pipefail
cd "${DCM_ROOT:?Missing project root}"
source "${EXPERIMENT_DIR:?Missing experiment directory}/ibex_environment.sh"
export LC_ALL=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
date -u
hostname
echo "Basis-size case ${CASE_ID:?Missing case}; CPUs $SLURM_CPUS_PER_TASK"
srun --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
  Rscript --vanilla "$EXPERIMENT_DIR/02_fit_case.R"
