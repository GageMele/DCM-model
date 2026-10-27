#!/bin/bash
# Freeze the existing data and submit 12 separately schedulable baseline fits.
set -euo pipefail
export LC_ALL=C
export DCM_ROOT="${DCM_ROOT:-/ibex/user/meleg/DCM_MODEL}"
cd "$DCM_ROOT"
package="$DCM_ROOT/DCM_Ibex_baseline_tests"
if [[ "${BASE_SELF_TEST_RUNNING:-0}" != 1 ]]; then
  bash "$package/self_test_submission.sh"
fi
source_experiment="${SOURCE_EXPERIMENT:-$(cat logs/latest_fast150_experiment.txt)}"
for file in "$source_experiment/PREPARATION_COMPLETE.txt" \
    "$source_experiment/full_input_250m.rds" "$source_experiment/data/p10.rds" \
    "$source_experiment/data/p20.rds" "$source_experiment/ibex_environment.sh"; do
  [[ -f "$file" ]] || { echo "Required frozen source absent: $file" >&2; exit 1; }
done
mkdir -p logs experiments
export EXPERIMENT_DIR="$DCM_ROOT/experiments/baseline150_$(date -u +%Y%m%dT%H%M%SZ)_$$"
mkdir "$EXPERIMENT_DIR"
mkdir "$EXPERIMENT_DIR/data"
cp "$package/"*.R "$package/"*.slurm "$package/"*.awk \
   "$package/cases.tsv" "$package/monitor_all_dcm_v3.sh" "$EXPERIMENT_DIR/"
cp "$source_experiment/ibex_environment.sh" "$EXPERIMENT_DIR/"
cp "$source_experiment/full_input_250m.rds" "$EXPERIMENT_DIR/full_input_250m.rds"
for subset in p10 p20; do
  cp "$source_experiment/data/$subset.rds" "$EXPERIMENT_DIR/data/$subset.rds"
  cmp -s "$source_experiment/data/$subset.rds" "$EXPERIMENT_DIR/data/$subset.rds"
done
cmp -s "$source_experiment/full_input_250m.rds" "$EXPERIMENT_DIR/full_input_250m.rds"
printf '%s\n' "$source_experiment" > "$EXPERIMENT_DIR/source_experiment.txt"
md5sum "$EXPERIMENT_DIR/full_input_250m.rds" "$EXPERIMENT_DIR/data/"*.rds \
  "$EXPERIMENT_DIR/baseline.R" "$EXPERIMENT_DIR/cases.tsv" > "$EXPERIMENT_DIR/frozen_inputs.md5"
source "$EXPERIMENT_DIR/ibex_environment.sh"
for file in "$package/"*.sh "$EXPERIMENT_DIR/"*.slurm; do bash -n "$file"; done
Rscript --vanilla -e 'for (f in commandArgs(TRUE)) parse(f); cat("R syntax checks passed\n")' \
  "$EXPERIMENT_DIR/"*.R
Rscript --vanilla -e 'a<-commandArgs(TRUE); source(file.path(a[1],"baseline.R")); validate_baseline_cases(read.delim(file.path(a[1],"cases.tsv"))); cat("12 baseline cases passed\n")' "$EXPERIMENT_DIR"
# Useful for checking the frozen files without submitting any Slurm jobs.
if [[ "${BASE_VALIDATE_ONLY:-0}" == 1 ]]; then
  echo "Frozen experiment ready: $EXPERIMENT_DIR"
  exit 0
fi
export CASE_ID=1
fit_time="${BASE_TIME:-7-00:00:00}"
constraint="${BASE_CONSTRAINT:-intel&cascadelake}"
# Check the largest CPU/memory request before creating any actual jobs.
sbatch --test-only --job-name=dcm-base150-resource-check --chdir="$DCM_ROOT" \
  --cpus-per-task=8 --constraint="$constraint" --hint=nomultithread \
  --mem=64G --time="$fit_time" --export=ALL "$EXPERIMENT_DIR/fit.slurm"
printf 'case_id\tlabel\tfit_job\tcheck_job\tcpus\trun_dir\n' > "$EXPERIMENT_DIR/manifest.tsv"
printf 'case_id\tstage\tjob_id\n' > "$EXPERIMENT_DIR/job_events.tsv"
printf '%s\n' "$EXPERIMENT_DIR" > logs/latest_baseline150_experiment.txt
submission=$(sbatch --parsable --job-name=dcm-base150-prepare --chdir="$DCM_ROOT" \
  --export=ALL "$EXPERIMENT_DIR/prepare.slurm")
prepare_job="${submission%%;*}"
[[ "$prepare_job" =~ ^[0-9]+$ ]] || { echo "Invalid preparation job ID" >&2; exit 1; }
printf '%s\n' "$prepare_job" > "$EXPERIMENT_DIR/preparation_job_id.txt"
checks=()
while IFS=$'\t' read -r case_id label subset prior mode chains iterations warmup delta tree threads cpus kd equation cores refresh variant memory; do
  [[ "$case_id" == case_id ]] && continue
  [[ -n "$case_id" ]] || continue
  export CASE_ID="$case_id"
  submission=$(sbatch --parsable --job-name="dcm-base150-c${case_id}" \
    --chdir="$DCM_ROOT" --dependency="afterok:$prepare_job" \
    --cpus-per-task="$cpus" --constraint="$constraint" --hint=nomultithread \
    --mem="${memory}G" --time="$fit_time" --export=ALL "$EXPERIMENT_DIR/fit.slurm")
  fit_job="${submission%%;*}"
  [[ "$fit_job" =~ ^[0-9]+$ ]] || { echo "Invalid fit job ID" >&2; exit 1; }
  printf '%s\tfit\t%s\n' "$case_id" "$fit_job" >> "$EXPERIMENT_DIR/job_events.tsv"
  # Record the fit immediately so a failed check submission cannot hide it.
  printf -v run_name 'case_%02d' "$case_id"
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$case_id" "$label" "$fit_job" "NA" "$cpus" \
    "$EXPERIMENT_DIR/runs/$run_name" >> "$EXPERIMENT_DIR/manifest.tsv"
  submission=$(sbatch --parsable --job-name="dcm-base150-check${case_id}" \
    --chdir="$DCM_ROOT" --dependency="afterany:$fit_job" \
    --export=ALL "$EXPERIMENT_DIR/check.slurm")
  check_job="${submission%%;*}"
  [[ "$check_job" =~ ^[0-9]+$ ]] || { echo "Invalid check job ID" >&2; exit 1; }
  checks+=("$check_job")
  printf '%s\tchecks\t%s\n' "$case_id" "$check_job" >> "$EXPERIMENT_DIR/job_events.tsv"
  # Replace only this manifest row, through a new file in the same directory.
  awk -F '\t' -v OFS='\t' -v id="$case_id" -v check="$check_job" \
    'NR>1 && $1==id {$4=check} {print}' "$EXPERIMENT_DIR/manifest.tsv" \
    > "$EXPERIMENT_DIR/manifest.next.tsv"
  mv "$EXPERIMENT_DIR/manifest.next.tsv" "$EXPERIMENT_DIR/manifest.tsv"
  echo "Case $case_id: $label; fit=$fit_job; checks=$check_job; CPUs=$cpus; memory=${memory}G"
done < "$EXPERIMENT_DIR/cases.tsv"
dependencies=$(IFS=:; echo "${checks[*]}")
submission=$(sbatch --parsable --job-name=dcm-base150-results --chdir="$DCM_ROOT" \
  --dependency="afterany:$dependencies" --export=ALL "$EXPERIMENT_DIR/collect.slurm")
collection_job="${submission%%;*}"
[[ "$collection_job" =~ ^[0-9]+$ ]] || { echo "Invalid collection job ID" >&2; exit 1; }
printf '%s\n' "$collection_job" > "$EXPERIMENT_DIR/collection_job_id.txt"
echo "Experiment: $EXPERIMENT_DIR"
echo "Preparation job: $prepare_job; collector: $collection_job"
echo "12 fits submitted: 8 CPUs each, 32G (10%/20%) or 64G (50%), time limit $fit_time."
echo "Monitor: bash DCM_Ibex_baseline_tests/monitor_all_dcm_v3.sh --diagnostics"
