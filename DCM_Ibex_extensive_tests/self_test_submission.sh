#!/bin/bash
# Exercise submission wiring with fake Slurm/R tools in an isolated directory.
# This test never calls the real sbatch or fits a model.
set -euo pipefail
package=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
sandbox=$(mktemp -d "${TMPDIR:-/tmp}/dcm-extensive-selftest.XXXXXXXX")
trap 'rm -rf -- "$sandbox"' EXIT
mkdir -p "$sandbox/project/DCM_Ibex_extensive_tests" "$sandbox/project/source/data" \
  "$sandbox/project/logs" "$sandbox/bin"
cp "$package/"* "$sandbox/project/DCM_Ibex_extensive_tests/"
for name in PREPARATION_COMPLETE.txt full_input_250m.rds ibex_environment.sh; do
  : > "$sandbox/project/source/$name"
done
for subset in p10 p20 p50; do : > "$sandbox/project/source/data/$subset.rds"; done
cat > "$sandbox/bin/Rscript" <<'MOCK_R'
#!/bin/bash
# R syntax and scientific construction are checked separately with real R.
exit 0
MOCK_R
cat > "$sandbox/bin/sbatch" <<'MOCK_SLURM'
#!/bin/bash
set -euo pipefail
job_name= dependency= cpus= memory= time_limit= test_only=0
for arg in "$@"; do
  case "$arg" in
    --test-only) test_only=1 ;;
    --job-name=*) job_name=${arg#*=} ;;
    --dependency=*) dependency=${arg#*=} ;;
    --cpus-per-task=*) cpus=${arg#*=} ;;
    --mem=*) memory=${arg#*=} ;;
    --time=*) time_limit=${arg#*=} ;;
  esac
done
if (( test_only )); then
  [[ "$cpus" == 16 && "$memory" == 64G ]] || exit 1
  : > "$MOCK_STATE/resource_check"
  exit 0
fi
counter=$(cat "$MOCK_STATE/counter")
counter=$((counter + 1))
printf '%s\n' "$counter" > "$MOCK_STATE/counter"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$counter" "$job_name" "$dependency" \
  "$cpus" "$memory" "$time_limit" "${CASE_ID:-}" >> "$MOCK_STATE/submissions.tsv"
printf '%s\n' "$counter"
MOCK_SLURM
chmod +x "$sandbox/bin/Rscript" "$sandbox/bin/sbatch"
printf '1000\n' > "$sandbox/counter"
export MOCK_STATE="$sandbox"
export PATH="$sandbox/bin:$PATH"
export DCM_ROOT="$sandbox/project"
export SOURCE_EXPERIMENT="$sandbox/project/source"
export EXT_SELF_TEST_RUNNING=1 EXT_VALIDATE_ONLY=0 EXT_TIME=7-00:00:00
bash "$sandbox/project/DCM_Ibex_extensive_tests/submit_tests.sh" > "$sandbox/submission.out"
experiment=$(cat "$sandbox/project/logs/latest_extensive150_experiment.txt")
[[ -f "$sandbox/resource_check" && $(cat "$sandbox/counter") == 1206 ]]
[[ $(cat "$experiment/preparation_job_id.txt") == 1001 ]]
[[ $(cat "$experiment/collection_job_id.txt") == 1206 ]]
[[ $(wc -l < "$experiment/manifest.tsv") == 103 ]]
[[ $(wc -l < "$experiment/job_events.tsv") == 205 ]]
while IFS=$'\t' read -r id label fit check cpus run; do
  [[ "$id" == case_id ]] && continue
  [[ "$fit" == $((1000 + 2 * id)) && "$check" == $((1001 + 2 * id)) && "$cpus" -ge 4 ]]
  awk -F '\t' -v fit="$fit" -v check="$check" -v id="$id" -v cpus="$cpus" '
    $1==fit {if($3!="afterok:1001" || $4!=cpus || $7!=id) exit 1; f=1}
    $1==check {if($3!="afterany:"fit || $7!=id) exit 1; c=1}
    END {if(!f || !c) exit 1}' "$sandbox/submissions.tsv"
done < "$experiment/manifest.tsv"
awk -F '\t' '
  $2~/^dcm-x150-c[0-9]+$/ {
    if($6!="7-00:00:00") exit 1
    expected=($7==85 || $7==87 || $7==89 || $7==91)?"64G":"32G"
    if($5!=expected) exit 1
    fits++
  }
  $2=="dcm-x150-results" {if($3!~/^afterany:1003:1005:1007:/) exit 1; collector=1}
  END {if(fits!=102 || !collector) exit 1}' "$sandbox/submissions.tsv"
echo "Submission self-test passed: 102 fits, correct resources/dependencies, 206 jobs total."
