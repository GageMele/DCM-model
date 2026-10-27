#!/bin/bash
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
set -euo pipefail
cd "${DCM_ROOT:?}"
source "${EXPERIMENT_DIR:?}/ibex_environment.sh"
export LC_ALL=C
srun --ntasks=1 --cpus-per-task=1 --cpu-bind=cores Rscript --vanilla "$EXPERIMENT_DIR/score_existing.R"
