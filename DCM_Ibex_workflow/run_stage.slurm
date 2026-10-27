#!/bin/bash -l
# Resource defaults for a pilot FIT; override for preparation/checking/production.
#SBATCH --job-name=dcm-model
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail
export DCM_ROOT="${DCM_ROOT:-/ibex/user/$USER/DCM_MODEL}"
cd "$DCM_ROOT"
source ./ibex_environment.sh
stage="${1:-fit}"
export MODEL_DATA="${MODEL_DATA:-$DCM_ROOT/data/processed/model_data.rds}"
export MODEL_THREADS="${MODEL_THREADS:-2}"
export STAN_NUM_THREADS="$MODEL_THREADS"
# Temporary compilation files can use /tmp; input/CSV/fit outputs stay persistent.
export TMPDIR="/tmp/dcm-${USER}-${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"
echo "Stage: $stage; job: $SLURM_JOB_ID; root: $DCM_ROOT"
df -h "$DCM_ROOT"
module list
case "$stage" in
  prepare) script=01_prepare_data.R ;;
  fit) script=02_fit_model.R ;;
  check)
    : "${RUN_DIR:?Export RUN_DIR to the completed run directory before submitting checks.}"
    script=03_check_and_plot.R ;;
  *) echo "Stage must be prepare, fit, or check" >&2; exit 2 ;;
esac
srun --ntasks=1 --cpus-per-task="$SLURM_CPUS_PER_TASK" Rscript --vanilla "$script"
