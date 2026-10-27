#!/bin/bash
# Freeze a separate experiment from the existing 150 m campaign, then submit.
set -euo pipefail
export DCM_ROOT="${DCM_ROOT:-/ibex/user/meleg/DCM_MODEL}"
cd "$DCM_ROOT"
package="$DCM_ROOT/DCM_Ibex_k_tests"
source_experiment="${SOURCE_EXPERIMENT:-$(cat logs/latest_fast150_experiment.txt)}"
for file in "$source_experiment/PREPARATION_COMPLETE.txt" \
    "$source_experiment/data/p10.rds" "$source_experiment/02_fit_model.R" \
    "$source_experiment/ibex_environment.sh" \
    "$source_experiment/full_response_statistics.rds"; do
  [[ -f "$file" ]] || { echo "Required frozen source absent: $file" >&2; exit 1; }
done
mkdir -p logs experiments
export EXPERIMENT_DIR="$DCM_ROOT/experiments/k150_$(date -u +%Y%m%dT%H%M%SZ)_$$"
mkdir "$EXPERIMENT_DIR"
mkdir "$EXPERIMENT_DIR/data"
cp "$package/"*.R "$package/"*.slurm "$package/"*.awk \
   "$package/cases.tsv" "$package/monitor.sh" "$EXPERIMENT_DIR/"
cp "$source_experiment/02_fit_model.R" "$source_experiment/ibex_environment.sh" \
   "$source_experiment/full_response_statistics.rds" "$EXPERIMENT_DIR/"
cp "$source_experiment/data/p10.rds" "$EXPERIMENT_DIR/data/p10.rds"
# Verify the input copy byte-for-byte; do not reselect or renormalize observations.
cmp -s "$source_experiment/data/p10.rds" "$EXPERIMENT_DIR/data/p10.rds"
printf '%s\n' "$source_experiment" > "$EXPERIMENT_DIR/source_experiment.txt"
md5sum "$EXPERIMENT_DIR/data/p10.rds" "$EXPERIMENT_DIR/02_fit_model.R" \
  > "$EXPERIMENT_DIR/frozen_inputs.md5"
source "$EXPERIMENT_DIR/ibex_environment.sh"
for file in "$EXPERIMENT_DIR/"*.slurm "$EXPERIMENT_DIR/monitor.sh"; do bash -n "$file"; done
Rscript --vanilla -e 'for (f in commandArgs(TRUE)) parse(f); cat("R syntax checks passed\n")' \
  "$EXPERIMENT_DIR/"*.R
# A local dry run can stop here; no Slurm commands have run at this point.
if [[ "${KTEST_VALIDATE_ONLY:-0}" == 1 ]]; then
  echo "Frozen experiment ready for validation: $EXPERIMENT_DIR"
  exit 0
fi
export CASE_ID=1
fit_time="${KTEST_TIME:-7-00:00:00}"
fit_memory="${KTEST_MEMORY:-32G}"
constraint="${KTEST_CONSTRAINT:-intel&cascadelake}"
# Check that Slurm accepts the new resource request before submitting any jobs.
sbatch --test-only --job-name=dcm-k150-resource-check --chdir="$DCM_ROOT" \
  --cpus-per-task=8 --constraint="$constraint" --hint=nomultithread \
  --mem="$fit_memory" --time="$fit_time" --export=ALL "$EXPERIMENT_DIR/fit.slurm"
printf 'case_id\tlabel\tfit_job\tcheck_job\tcpus\trun_dir\n' > "$EXPERIMENT_DIR/manifest.tsv"
printf 'case_id\tstage\tjob_id\n' > "$EXPERIMENT_DIR/job_events.tsv"
# Keep a partial campaign discoverable if a later submission fails.
printf '%s\n' "$EXPERIMENT_DIR" > logs/latest_k150_experiment.txt
submission=$(sbatch --parsable --job-name=dcm-k150-prepare --chdir="$DCM_ROOT" \
  --export=ALL "$EXPERIMENT_DIR/prepare.slurm")
prepare_job="${submission%%;*}"
[[ "$prepare_job" =~ ^[0-9]+$ ]] || { echo "Invalid preparation job ID" >&2; exit 1; }
printf '%s\n' "$prepare_job" > "$EXPERIMENT_DIR/preparation_job_id.txt"
checks=()
while IFS=$'\t' read -r case_id label subset prior mode chains iterations warmup delta tree threads cpus kd kdist kmld; do
  [[ "$case_id" == case_id ]] && continue
  [[ -n "$case_id" ]] || continue
  export CASE_ID="$case_id"
  submission=$(sbatch --parsable --job-name="dcm-k150-c${case_id}" \
    --chdir="$DCM_ROOT" --dependency="afterok:$prepare_job" \
    --cpus-per-task="$cpus" --constraint="$constraint" --hint=nomultithread \
    --mem="$fit_memory" --time="$fit_time" --export=ALL "$EXPERIMENT_DIR/fit.slurm")
  fit_job="${submission%%;*}"
  [[ "$fit_job" =~ ^[0-9]+$ ]] || { echo "Invalid fit job ID" >&2; exit 1; }
  printf '%s\tfit\t%s\n' "$case_id" "$fit_job" >> "$EXPERIMENT_DIR/job_events.tsv"
  submission=$(sbatch --parsable --job-name="dcm-k150-check${case_id}" \
    --chdir="$DCM_ROOT" --dependency="afterany:$fit_job" \
    --export=ALL "$EXPERIMENT_DIR/check.slurm")
  check_job="${submission%%;*}"
  [[ "$check_job" =~ ^[0-9]+$ ]] || { echo "Invalid check job ID" >&2; exit 1; }
  checks+=("$check_job")
  printf '%s\tchecks\t%s\n' "$case_id" "$check_job" >> "$EXPERIMENT_DIR/job_events.tsv"
  printf -v run_name 'case_%02d' "$case_id"
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$case_id" "$label" "$fit_job" "$check_job" "$cpus" \
    "$EXPERIMENT_DIR/runs/$run_name" >> "$EXPERIMENT_DIR/manifest.tsv"
  echo "Case $case_id: $label; k=($kd,$kdist,$kmld); tree=$tree; fit=$fit_job; checks=$check_job"
done < "$EXPERIMENT_DIR/cases.tsv"
dependencies=$(IFS=:; echo "${checks[*]}")
submission=$(sbatch --parsable --job-name=dcm-k150-results --chdir="$DCM_ROOT" \
  --dependency="afterany:$dependencies" --export=ALL "$EXPERIMENT_DIR/collect.slurm")
printf '%s\n' "${submission%%;*}" > "$EXPERIMENT_DIR/collection_job_id.txt"
echo "Experiment: $EXPERIMENT_DIR"
echo "Preparation job: $prepare_job"
echo "12 independent fits submitted; each requests 8 CPUs, $fit_memory, and $fit_time."
echo "Current jobs are unchanged. Run bash DCM_Ibex_k_tests/monitor.sh for progress."
