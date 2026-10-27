# Verified Ibex queue snapshot: 27 October 2026

Captured for `meleg` at **10:47:40 UTC / 13:47:40 Saudi time**: **335 account entries, 20 RUNNING, 315 PENDING**. Exact DCM manifest matching identifies **19 running shape fits (304 CPUs)**. The other running job is `53515610`, `C_D1_202202_run` (96 CPUs). No pending DCM match, DCM controller, or DCM check appears in the capture.

Read the [complete annotated inventory](../../../docs/CURRENT_JOBS.md) for all running jobs, planned model settings, pending reasons, operational interpretation, and the earlier historical snapshot. The October 25 scientific result tables remain unchanged.

## Original evidence

The six files in `raw/` are the original archive members, without byte changes:

- [all_queued_jobs.csv](raw/all_queued_jobs.csv): every account queue entry, all 13 Slurm fields.
- [running_jobs.csv](raw/running_jobs.csv): the exact RUNNING subset, all 20 entries.
- [dcm_manifest_matches.csv](raw/dcm_manifest_matches.csv): exact DCM job-ID associations with campaign, stage, case, label, manifest and run directory.
- [squeue_raw.txt](raw/squeue_raw.txt): original pipe-delimited output.
- [metadata.json](raw/metadata.json): query timestamps, command, account, totals and scope.
- [README.md](raw/README.md): original readable table containing all 335 entries.

The [previously reported metadata](reported_capture_metadata.json) is preserved and equals the archive metadata. The original remote archive path is:

~~~text
/ibex/user/meleg/DCM_MODEL/logs/ibex_jobs_20261027T104740Z_1793098060392524450.tar.gz
~~~

Archive bytes: **10,014**. SHA256:

~~~text
d49aae3377f68af41d2fa664872ca1c4061a570a6c7d2ec46aa05cbdc010ddcc
~~~

The source archive was uploaded by the user. Its original member contents, per-member hashes and archive hash are preserved here; the compressed container remains external.

## Derived views

- [classified_queue.csv](classified_queue.csv): all 335 entries with explicit match evidence.
- [dcm_running_with_settings.csv](dcm_running_with_settings.csv): all 19 DCM fits joined to every field of the frozen shape case matrix, with `plan_` prefixes.
- [pending_jobs.csv](pending_jobs.csv): all 315 pending entries.
- [running_without_dcm_match.csv](running_without_dcm_match.csv): the other running job.
- [summary.json](summary.json): counts, resource totals, name groups, and interpretation limits.
- [validation.json](validation.json): checks of the supplied evidence.
- [provenance.json](provenance.json): file hashes, source paths, join method and execution boundary.
- [build_inventory.py](build_inventory.py): exact offline source used to validate, preserve and derive the structured snapshot files.

To reproduce the structured import, keep this repository at `DCM-model/` and the source archive at `github_handover_inputs/jobs_20261027/` under the same parent, then run `python DCM-model/operations/jobs/20261027T104740Z/build_inventory.py` from that parent. The code uses only Python's standard library and requires Python 3.9 or later. It verifies the expected archive hash and fails on inconsistent evidence. It writes only the dated snapshot's structured evidence files and a copy of itself, and prints a summary. It does not regenerate hand-edited documentation, repository-wide checksums, or publish to GitHub. It makes no SSH or Slurm calls.

## Interpretation limits

Exact manifest matching identifies the DCM fits; unmatched names are retained without asserting another project's ownership. All visible DCM names match the manifest. The missing controller needs accounting/log inspection; this capture does not establish a failure or require an automatic restart. The queue does not show completed-job history, chain progress, model validity or new predictive scores. Slurm start times have no recorded offset; query start/end are explicitly UTC. CPU and memory fields are allocations/requests rather than utilization.

Keep the [October 25 registrations](../20261025/shape_all_282_job_registrations.csv) and all scientific history. Never infer completed case counts by subtracting live jobs from the 282-case plan.
