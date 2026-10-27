#!/bin/bash
# Stage-specific CPUs, memory and time are supplied by the dispatcher.
set -euo pipefail
cd -- "${DCM_ROOT:?Missing project directory}"
source "${EXPERIMENT_DIR:?Missing experiment directory}/ibex_environment.sh"
export LC_ALL=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
case "${STAGE:?Missing stage}" in
  prepare) Rscript --vanilla "$EXPERIMENT_DIR/test_offline.R" "$EXPERIMENT_DIR"
           srun --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
             Rscript --vanilla "$EXPERIMENT_DIR/01_prepare.R" ;;
  controller) python3 "$EXPERIMENT_DIR/scheduler.py" control "$EXPERIMENT_DIR" ;;
  fit) srun --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
         Rscript --vanilla "$EXPERIMENT_DIR/02_fit_case.R" ;;
  check) srun --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=cores \
           Rscript --vanilla "$EXPERIMENT_DIR/03_check_case.R" ;;
  *) echo "Unknown stage: $STAGE" >&2; exit 2 ;;
esac
