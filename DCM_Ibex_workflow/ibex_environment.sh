#!/bin/bash
# Source this file; use the same R module/library for setup and every job.
export DCM_ROOT="${DCM_ROOT:-/ibex/user/$USER/DCM_MODEL}"
# Confirmed in your live module list and successful installation log.
export R_MODULE="${R_MODULE:-R/4.5.0/gnu-12.2.0}"
module purge
module load "$R_MODULE"
DCM_R_SERIES="$(Rscript --vanilla -e 'cat(paste(R.version$major, strsplit(R.version$minor,".",fixed=TRUE)[[1]][1],sep="."))')"
export R_LIBS_USER="$DCM_ROOT/software/R-$DCM_R_SERIES"
mkdir -p "$R_LIBS_USER"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
# Setup records its exact installation path here (independent of login node).
if [[ -s "$DCM_ROOT/software/cmdstan_path.txt" ]]; then
  export CMDSTAN="$(cat "$DCM_ROOT/software/cmdstan_path.txt")"
fi
