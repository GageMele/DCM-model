#!/bin/bash
# Exercise Slurm submission wiring and live diagnostics without submitting jobs.
set -euo pipefail
export PATH="/c/rtools45/usr/bin:$PATH"
workspace="$PWD"
package="$workspace/DCM_Ibex_k_tests"
for script in "$package/"*.sh "$package/"*.slurm; do bash -n "$script"; done
test_root=$(mktemp -d "$workspace/validation_k150_shell_XXXXXX")
mkdir -p "$test_root/root/logs" "$test_root/root/experiments/source/data" "$test_root/bin"
cp -r "$package" "$test_root/root/DCM_Ibex_k_tests"
parent="$test_root/root/experiments/source"
printf 'frozen test input\n' > "$parent/data/p10.rds"
printf 'completed\n' > "$parent/PREPARATION_COMPLETE.txt"
printf '# frozen model stub\n' > "$parent/02_fit_model.R"
printf '# frozen environment stub\n' > "$parent/ibex_environment.sh"
printf 'stats stub\n' > "$parent/full_response_statistics.rds"
printf '%s\n' "$parent" > "$test_root/root/logs/latest_fast150_experiment.txt"
cat > "$test_root/bin/sbatch" <<'MOCK'
#!/bin/bash
set -euo pipefail
printf '%s\n' "$*" >> "$MOCK_LOG"
if [[ " $* " == *' --test-only '* ]]; then
    [[ "${FAIL_RESOURCE:-0}" == 0 ]] || exit 1
    exit 0
fi
id=$(cat "$MOCK_COUNTER")
id=$((id+1))
printf '%s\n' "$id" > "$MOCK_COUNTER"
printf '%s;mock-cluster\n' "$id"
MOCK
cat > "$test_root/bin/Rscript" <<'MOCK'
#!/bin/bash
# Real R parsing and basis construction are checked by validate_k150.R.
exit 0
MOCK
chmod +x "$test_root/bin/sbatch" "$test_root/bin/Rscript"
export PATH="$test_root/bin:$PATH"
export MOCK_LOG="$test_root/sbatch.log" MOCK_COUNTER="$test_root/counter"
printf '9000\n' > "$MOCK_COUNTER"
export DCM_ROOT="$test_root/root"

# Resource rejection must occur before any real submission or campaign pointer.
if FAIL_RESOURCE=1 bash "$DCM_ROOT/DCM_Ibex_k_tests/submit_tests.sh"; then
    echo 'Expected resource rejection did not occur.' >&2; exit 1
fi
[[ "$(cat "$MOCK_COUNTER")" == 9000 ]]
[[ ! -f "$DCM_ROOT/logs/latest_k150_experiment.txt" ]]

# Successful campaign: 1 prepare + 12 fits + 12 checks + 1 collector.
bash "$DCM_ROOT/DCM_Ibex_k_tests/submit_tests.sh"
[[ "$(cat "$MOCK_COUNTER")" == 9026 ]]
experiment=$(cat "$DCM_ROOT/logs/latest_k150_experiment.txt")
[[ "$(cat "$experiment/preparation_job_id.txt")" == 9001 ]]
[[ "$(cat "$experiment/collection_job_id.txt")" == 9026 ]]
cmp -s "$parent/data/p10.rds" "$experiment/data/p10.rds"
awk -F '\t' 'NR>1 {n++; if(NF!=6 || $5!=8 || $3!~/^[0-9]+$/ || $4!~/^[0-9]+$/) exit 1}
    END {if(n!=12) exit 1}' "$experiment/manifest.tsv"
[[ "$(grep -c -- '--dependency=afterok:9001' "$MOCK_LOG")" == 12 ]]
[[ "$(grep -c -- '--dependency=afterany:' "$MOCK_LOG")" == 13 ]]

# Fixture includes warmup, post-warmup, and a partial final CSV row.
fixture="$test_root/chain.csv"
{
  printf '# id = 2\n# max_depth = 12\n'
  printf 'lp__,accept_stat__,stepsize__,treedepth__,n_leapfrog__,divergent__,energy__,theta\n'
  for ((i=1;i<=60;i++)); do
    if (( i<=10 )); then depth=1; else depth=12; fi
    printf '%s\n' "-10,0.95,0.01,$depth,15,0,12,0.1"
  done
  printf '# Adaptation terminated\n'
  printf '%s\n' '-10,0.95,0.01,12,20,0,12,0.1' \
    '-10,0.95,0.01,11,10,1,12,0.1' '-10,0.95,0.01,12,30,0,12,0.1'
  printf '%s' '-10,0.95,0.01,12'
} > "$fixture"
awk -F, -v case_id=9 -f "$package/live_diagnostics.awk" "$fixture" > "$test_root/live.txt"
cat "$test_root/live.txt"
awk 'NR==1 {if($1!=9 || $2!=2 || $3!="warmup" || $4!=60 || $7!=50 || $8!=83.3 || $9!=100) exit 1}
     NR==2 {if($3!="sampling" || $4!=3 || $7!=2 || $10!=1 || $11!=20) exit 1}
     END {if(NR!=2) exit 1}' "$test_root/live.txt"
echo 'Shell syntax, resource rejection, 26-job wiring, input isolation and live phase diagnostics passed.'
