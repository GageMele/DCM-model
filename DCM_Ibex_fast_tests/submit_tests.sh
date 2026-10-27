#!/bin/bash
# Freeze an experiment and submit all cases as independently schedulable jobs.
set -euo pipefail
export DCM_ROOT="${DCM_ROOT:-/ibex/user/meleg/DCM_MODEL}"
cd "$DCM_ROOT"
package="$DCM_ROOT/DCM_Ibex_fast_tests"
export FULL_MODEL_DATA="${FULL_MODEL_DATA:-$DCM_ROOT/data/processed/full_unsmoothed_v2.rds}"
for file in "$FULL_MODEL_DATA" "$DCM_ROOT/02_fit_model.R" "$DCM_ROOT/ibex_environment.sh"; do
  [[ -f "$file" ]] || { echo "Required file absent: $file" >&2; exit 1; }
done
mkdir -p logs experiments
export EXPERIMENT_DIR="$DCM_ROOT/experiments/fast150_$(date -u +%Y%m%dT%H%M%SZ)_$$"
mkdir "$EXPERIMENT_DIR"
cp "$package/"*.R "$package/"*.slurm "$package/cases.tsv" "$EXPERIMENT_DIR/"
cp "$DCM_ROOT/02_fit_model.R" "$DCM_ROOT/ibex_environment.sh" "$EXPERIMENT_DIR/"
# Freeze the full input too, so other preparation work cannot change these tests.
printf '%s\n' "$FULL_MODEL_DATA" > "$EXPERIMENT_DIR/original_full_bundle_path.txt"
cp "$FULL_MODEL_DATA" "$EXPERIMENT_DIR/full_input_250m.rds"
export FULL_MODEL_DATA="$EXPERIMENT_DIR/full_input_250m.rds"
source "$EXPERIMENT_DIR/ibex_environment.sh"
# Check syntax and the case manifest before submitting jobs.
for file in "$EXPERIMENT_DIR/"*.slurm; do bash -n "$file"; done
Rscript --vanilla -e 'for (f in commandArgs(TRUE)) parse(f); cat("R syntax checks passed\n")' \
  "$EXPERIMENT_DIR/common.R" "$EXPERIMENT_DIR/01_prepare_tests.R" \
  "$EXPERIMENT_DIR/02_fit_case.R" "$EXPERIMENT_DIR/03_check_case.R" \
  "$EXPERIMENT_DIR/04_collect_results.R" "$EXPERIMENT_DIR/02_fit_model.R"
Rscript --vanilla -e 'd<-read.delim(commandArgs(TRUE)[1]); stopifnot(nrow(d)==24L,!anyDuplicated(d$case_id),all(d$iter>d$warmup),all(d$chains>=2),all(d$cpus==d$chains*pmax(1,d$threads))); cat("Case manifest passed\n")' \
  "$EXPERIMENT_DIR/cases.tsv"
submission=$(sbatch --parsable --job-name=dcm-f150-prepare --chdir="$DCM_ROOT" \
  --export=ALL "$EXPERIMENT_DIR/prepare.slurm")
prepare_job="${submission%%;*}"
printf '%s\n' "$prepare_job" > "$EXPERIMENT_DIR/preparation_job_id.txt"
printf 'case_id\tlabel\tfit_job\tcheck_job\tcpus\trun_dir\n' > "$EXPERIMENT_DIR/manifest.tsv"
printf 'case_id\tstage\tjob_id\n' > "$EXPERIMENT_DIR/job_events.tsv"
# Save the pointer before submission, so a partial submission remains discoverable.
printf '%s\n' "$EXPERIMENT_DIR" > logs/latest_fast150_experiment.txt
checks=()
while IFS=$'\t' read -r case_id label subset prior mode chains iterations warmup delta tree threads cpus; do
  [[ "$case_id" == case_id ]] && continue
  [[ -n "$case_id" ]] || continue
  export CASE_ID="$case_id"
  submission=$(sbatch --parsable --job-name="dcm-f150-c${case_id}" \
    --chdir="$DCM_ROOT" --dependency="afterok:$prepare_job" \
    --cpus-per-task="$cpus" --constraint="${TEST_CONSTRAINT:-intel&cascadelake}" \
    --hint=nomultithread --mem="${TEST_MEMORY:-32G}" --time="${TEST_TIME:-12:00:00}" \
    --export=ALL "$EXPERIMENT_DIR/fit.slurm")
  fit_job="${submission%%;*}"
  printf '%s\tfit\t%s\n' "$case_id" "$fit_job" >> "$EXPERIMENT_DIR/job_events.tsv"
  submission=$(sbatch --parsable --job-name="dcm-f150-check${case_id}" \
    --chdir="$DCM_ROOT" --dependency="afterany:$fit_job" \
    --export=ALL "$EXPERIMENT_DIR/check.slurm")
  check_job="${submission%%;*}"
  printf '%s\tchecks\t%s\n' "$case_id" "$check_job" >> "$EXPERIMENT_DIR/job_events.tsv"
  checks+=("$check_job")
  printf -v run_name 'case_%02d' "$case_id"
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$case_id" "$label" "$fit_job" "$check_job" "$cpus" \
    "$EXPERIMENT_DIR/runs/$run_name" >> "$EXPERIMENT_DIR/manifest.tsv"
  echo "Case $case_id: $label; fit $fit_job; checks $check_job; $cpus CPUs"
done < "$EXPERIMENT_DIR/cases.tsv"
dependencies=$(IFS=:; echo "${checks[*]}")
submission=$(sbatch --parsable --job-name=dcm-f150-results --chdir="$DCM_ROOT" \
  --dependency="afterany:$dependencies" --export=ALL "$EXPERIMENT_DIR/collect.slurm")
printf '%s\n' "${submission%%;*}" > "$EXPERIMENT_DIR/collection_job_id.txt"
echo "Preparation job: $prepare_job"
echo "Experiment: $EXPERIMENT_DIR"
echo "All 24 fits submitted independently; Slurm controls simultaneous starts."
