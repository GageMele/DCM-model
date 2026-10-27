#!/bin/bash
# Read-only monitor: original pilot, threading benchmarks, fast150, k150 and equation150.
# Optional --diagnostics scans saved CSVs in all 150 m campaigns.
set -uo pipefail
cd "${DCM_ROOT:-/ibex/user/meleg/DCM_MODEL}" || exit 1
shopt -s nullglob
declare -A known_jobs

# Include the original pilot and benchmark jobs even if manifests are absent.
for job in 53749517 53759873 53759874 53759876 53759877 53759878; do
    known_jobs["$job"]=1
done

# Discover all recorded campaigns, including completed fits and helper jobs.
manifests=(benchmarks/threading_*/manifest.tsv
           experiments/fast150_*/manifest.tsv
           experiments/k150_*/manifest.tsv experiments/equation150_*/manifest.tsv)
for manifest in "${manifests[@]}"; do
    while IFS= read -r job; do
        [[ "$job" =~ ^[0-9]+$ ]] && known_jobs["$job"]=1
    done < <(awk -F '\t' '
        NR==1 {for(i=1;i<=NF;i++)
            if($i=="job_id" || $i=="fit_job" || $i=="check_job") keep[i]=1; next}
        {for(i in keep) print $i}' "$manifest")
    campaign=$(dirname "$manifest")
    for name in preparation_job_id.txt collection_job_id.txt comparison_job_id.txt; do
        [[ -f "$campaign/$name" ]] || continue
        job=$(head -n 1 "$campaign/$name")
        job="${job%%;*}"
        [[ "$job" =~ ^[0-9]+$ ]] && known_jobs["$job"]=1
    done
done
ids=$(printf '%s\n' "${!known_jobs[@]}" | sort -n | paste -sd, -)

printf '\nDCM snapshot: %s\n' "$(date -Is)"
printf '\nACTIVE DCM JOBS: scheduling reasons and time limits\n'
# Filter by the complete job name, excluding CROCO and other workloads.
squeue --me -o '%.12i %.30j %.10M %.12l %.2t %.40R' |
    awk 'NR==1 || $2 ~ /^dcm-/'

printf '\nALL RECORDED DCM JOBS: includes completed/failed fits and checks\n'
# Explicit IDs retain jobs that have disappeared from squeue.
sacct -X --starttime=2026-10-13 -j "$ids" \
    --format=JobID,JobName%26,State%16,ExitCode,Elapsed,Timelimit,AllocCPUS

printf '\nLATEST PRINTED CHAIN PROGRESS\n'
logs=(logs/dcm-pilot-*.out logs/dcm-bC*.out
      logs/dcm-f150-c[0-9]*.out logs/dcm-k150-c[0-9]*.out logs/dcm-e150-c[0-9]*.out)
for path in "${logs[@]}"; do
    printf '\n%s\n' "$path"
    # Read only a short tail; retain the latest checkpoint for each chain.
    tail -n 200 "$path" | awk '
        /^Chain [0-9]+ Iteration:/ {
            latest[$2]=$0
            if($2>maximum) maximum=$2
        }
        END {
            found=0
            for(i=1;i<=maximum;i++) if(i in latest) {print latest[i]; found=1}
            if(!found) print "No iteration checkpoint in the log tail yet."
        }'
done

printf '\nRECORDED FITTING OR CHECKING ERRORS\n'
errors=0
for manifest in "${manifests[@]}"; do
    while IFS= read -r run; do
        [[ -d "$run" ]] || continue
        for name in FIT_FAILED.txt CHECKS_FAILED.txt; do
            if [[ -f "$run/$name" ]]; then
                errors=$((errors+1))
                printf '\n%s\n' "$run/$name"
                cat "$run/$name"
            fi
        done
    done < <(awk -F '\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="run_dir") col=i; next}
        col {print $col}' "$manifest")
done
(( errors > 0 )) || echo "No FIT_FAILED/CHECKS_FAILED markers found; check Slurm states above too."

if [[ "${1:-}" == --diagnostics ]]; then
    diagnostics_awk="DCM_Ibex_equation_tests/live_diagnostics.awk"
    if [[ ! -f "$diagnostics_awk" ]]; then
        echo "CSV diagnostics unavailable: $diagnostics_awk is missing."
        exit 0
    fi
    printf '\nLIVE CSV DIAGNOSTICS: saved warmup and sampling, fast150, k150 and equation150\n'
    printf 'These campaigns save warmup; original pilot/benchmarks use the log progress above.\n'
    for manifest in experiments/fast150_*/manifest.tsv experiments/k150_*/manifest.tsv experiments/equation150_*/manifest.tsv; do
        printf '\nCampaign: %s\n' "$(dirname "$manifest")"
        while IFS=$'\t' read -r case_id label fit_job check_job cpus run; do
            [[ "$case_id" =~ ^[0-9]+$ ]] || continue
            files=("$run"/stan_csv/*.csv)
            (( ${#files[@]} > 0 )) || continue
            printf '\nCase %s: %s\n' "$case_id" "$label"
            printf '%4s %5s %-8s %6s %3s %3s %6s %7s %8s %5s %10s\n' \
                case chain phase draws cap max hits 'hits%' 'recent%' div mean_leaps
            for csv in "${files[@]}"; do
                awk -F, -v case_id="$case_id" -f "$diagnostics_awk" "$csv"
            done
        done < "$manifest"
    done
fi
printf '\nSnapshot complete. This script does not modify or submit jobs.\n'
