#!/bin/bash
# Leave jobs running. Inspect limits; optionally attempt an authorized extension.
# Slurm normally requires an administrator to increase a RUNNING job's limit.
set -u
jobs=(53749517 53759873 53759874 53759876 53759877)
ids=$(IFS=,; echo "${jobs[*]}")
squeue -j "$ids" -o '%.12i %.18j %.10M %.12l %.2t'
if [[ "${1:-}" != --attempt ]]; then
  echo "An administrator is normally required to extend these running jobs."
  echo "Requested target: 14 days TOTAL walltime per job, not 14 extra days."
  echo "The authorized extension attempt is available with --attempt."
  exit 0
fi
for job in "${jobs[@]}"; do
  echo "Requesting 14-day total walltime for job $job"
  if ! scontrol update "JobId=$job" "TimeLimit=14-00:00:00"; then
    echo "Job $job was not extended. Request an extension from Ibex Support."
  fi
done
squeue -j "$ids" -o '%.12i %.18j %.10M %.12l %.2t'
