#!/bin/bash
# Validate Slurm wiring and the reporter fix using fake submission commands.
set -euo pipefail
export PATH="/c/rtools45/usr/bin:$PATH"
workspace="$PWD"
package="$workspace/DCM_Ibex_equation_tests"
for script in "$package/"*.sh "$package/"*.slurm; do bash -n "$script"; done
test_root=$(mktemp -d "$workspace/validation_equation150_shell_XXXXXX")
mkdir -p "$test_root/root/logs" "$test_root/root/experiments/source/data" "$test_root/bin"
cp -r "$package" "$test_root/root/DCM_Ibex_equation_tests"
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
# Actual R parsing and model construction are tested separately.
exit 0
MOCK
chmod +x "$test_root/bin/sbatch" "$test_root/bin/Rscript"
export PATH="$test_root/bin:$PATH"
export MOCK_LOG="$test_root/sbatch.log" MOCK_COUNTER="$test_root/counter"
printf '9000\n' > "$MOCK_COUNTER"
export DCM_ROOT="$test_root/root"
if FAIL_RESOURCE=1 bash "$DCM_ROOT/DCM_Ibex_equation_tests/submit_tests.sh"; then
    echo 'Expected rejection did not occur.' >&2; exit 1
fi
[[ "$(cat "$MOCK_COUNTER")" == 9000 ]]
[[ ! -f "$DCM_ROOT/logs/latest_equation150_experiment.txt" ]]
bash "$DCM_ROOT/DCM_Ibex_equation_tests/submit_tests.sh"
# One prepare + 18 fits + 18 checks + one collector.
[[ "$(cat "$MOCK_COUNTER")" == 9038 ]]
experiment=$(cat "$DCM_ROOT/logs/latest_equation150_experiment.txt")
[[ "$(cat "$experiment/preparation_job_id.txt")" == 9001 ]]
[[ "$(cat "$experiment/collection_job_id.txt")" == 9038 ]]
cmp -s "$parent/data/p10.rds" "$experiment/data/p10.rds"
awk -F '\t' 'NR>1 {n++; if(NF!=6 || ($5!=8 && $5!=16)) exit 1; cpus+=$5}
    END {if(n!=18 || cpus!=168) exit 1}' "$experiment/manifest.tsv"
[[ "$(grep -c -- '--dependency=afterok:9001' "$MOCK_LOG")" == 18 ]]
[[ "$(grep -c -- '--dependency=afterany:' "$MOCK_LOG")" == 19 ]]
echo 'Shell syntax, resource rejection, 38-job wiring and byte-identical input copies passed.'
