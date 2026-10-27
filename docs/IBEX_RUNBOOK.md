# Ibex operational runbook

## Execution locations

Project root: /ibex/user/meleg/DCM_MODEL

Windows connects with:

~~~powershell
ssh meleg@ilogin.ibex.kaust.edu.sa
~~~

After authentication, run the Bash commands below on Ibex. Local Downloads paths and PowerShell environment syntax do not apply inside that Bash session. Login nodes in the supplied logs include login509-02-l and login509-02-r.

This repository is documentation and preserved deployment source. It is not a live mount of Ibex. Do not assume a GitHub update changes cluster code or jobs.

## Current frozen experiment and pointers

Latest verified shape experiment:

~~~text
/ibex/user/meleg/DCM_MODEL/experiments/shape150_20261022T102851Z_3426481
~~~

The logs/latest_shape150_experiment.txt pointer should resolve to the intended campaign. Other pointers include latest_fast150_experiment.txt, latest_equation150_experiment.txt, latest_followup150_experiment.txt, and latest_onset150_experiment.txt. Read their content rather than guessing timestamps.

A frozen experiment carries cases.tsv, copied source, configuration, prepared data/memberships, Slurm ledger/manifest, run directories, logs, result tables, and completion/failure markers. Shape uses three-digit run names (case_034); older campaigns often use two digits. Inspect the actual manifest.

## Environment

Recorded cluster stack:

| Component | Recorded value |
|---|---|
| R module | R/4.5.0/gnu-12.2.0 |
| Dependencies loaded by module | JDK 19.0.1, GNU 12.2.0 |
| R user library | $DCM_ROOT/software/R-4.5 |
| CmdStan path pointer | software/cmdstan_path.txt |
| CmdStan | 2.37.0 |
| Representative shape OS | Rocky Linux 9.4 |
| brms | 2.23.0 |
| cmdstanr | 0.9.0 |
| mgcv | 1.9-1 |
| posterior | 1.7.1 |

These are observed versions, not a complete package lockfile. Source the existing experiment's ibex_environment.sh instead of installing/upgrading packages during an active campaign.

BLAS/OpenMP helper environment variables OMP_NUM_THREADS, OPENBLAS_NUM_THREADS, MKL_NUM_THREADS, and VECLIB_MAXIMUM_THREADS are capped at one. Stan's requested within-chain threading is separately configured. Do not silently change the compiler, CmdStan, package versions, factor levels, or data ordering when resuming.

## Job inventory and read-only monitoring

See [CURRENT_JOBS.md](CURRENT_JOBS.md) for the successful user-reported October 27 capture: 20 RUNNING jobs and 335 total entries at 10:47 UTC / 13:47 Saudi time, across the entire meleg account. The archive with individual jobs is awaiting transfer. Direct SSH from the authoring workspace was blocked earlier; the October 25 archived classifications remain historical.

~~~bash
cd /ibex/user/meleg/DCM_MODEL
experiment=$(cat logs/latest_shape150_experiment.txt)

python3 "$experiment/scheduler.py" status "$experiment"
bash DCM_Ibex_shape_tests/monitor_all_dcm_v8.sh --diagnostics

squeue -u meleg -o "%.18i %.45j %.10T %.12M %.6C %R"
sacct -X -j 54250664,54250665 \
  --format=JobID,JobName%45,State,Elapsed,ExitCode
~~~

The scheduler status command prints its last saved dispatch_status.json. Check its UTC timestamp; it is not itself a fresh scheduler pass. squeue gives live queue state, sacct includes recorded completed jobs. The two IDs above are historical preparation/initial-controller IDs; use the current ledger for successors.

The v8 monitor scans manifests under experiments/benchmarks and queries Slurm. Its “all” label does not guarantee every standalone pilot directory or existing_evaluations.tsv prediction-only row is enumerated. The all-campaign collector is the broader scientific registry. Earlier v4 was in the extensive package, not the baseline package; the missing-path console error did not stop model fits.

Read logs without resubmitting:

~~~bash
for log_file in "$experiment"/logs/*-54250665.out \
                "$experiment"/logs/*-54250665.err; do
    [[ -f "$log_file" ]] || continue
    printf '\n%s\n' "$log_file"
    tail -n 60 "$log_file"
done
~~~

Replace that job ID with the relevant current ID. Inspect FIT_COMPLETE.txt/FIT_FAILED.txt and CHECKS_COMPLETE.txt/CHECKS_FAILED.txt separately in each run.

## Refresh and export completed results

~~~bash
cd /ibex/user/meleg/DCM_MODEL
bash DCM_Ibex_shape_tests/collect_results.sh
bash DCM_Ibex_shape_tests/export_results.sh
~~~

Collection reads existing outputs. It does not refit models or magically generate missing predictions. Export prints a timestamped archive path. Preserve previous exports.

The compact export omits raw chains, data directories, .git, large sampler tables, and individual files over its size limit. It does not constitute a complete backup of every shell/Python script or RDS object. Active files may change during tar; inspect reported status. A tar status of one can mean files changed during collection, not necessarily a missing archive.

The known October 25 export was:

~~~text
/ibex/user/meleg/DCM_MODEL/logs/dcm_shape_comparison_20261025T111507Z_160796.tar.gz
~~~

Download in a separate local PowerShell session:

~~~powershell
scp meleg@ilogin.ibex.kaust.edu.sa:/ibex/user/meleg/DCM_MODEL/logs/dcm_shape_comparison_20261025T111507Z_160796.tar.gz "$env:USERPROFILE\Downloads\"
~~~

Use the newly printed archive name for a fresh snapshot. An SSH connection reset does not establish that Slurm jobs stopped; reconnect and query Slurm.

## Shape dispatcher configuration

The verified campaign configuration used partition batch and fit constraint intel&cascadelake. Its package identity record was 0859014f1537ded5a7298d47e124ff97af5a271bb89b54314fcf51e170c90db6; the historical self-extracting launcher's ZIP payload hash is a different object and must not be substituted for that identifier.

Defaults are 48 active fitting/checking jobs and 1024 requested CPUs. The dispatcher counts accepted pending/running model stages against limits. allocated_model_cpus is accounted requested capacity, not measured CPU utilization. Slurm determines actual start times and resource allocation.

Fits request 16 CPUs through 50% and 32 CPUs at 75/100%, memory 64–256 GiB. Fit walltime requests are seven days, prior checks four hours, preparation eight hours, evaluation two days, controller seven days; actual walltime is capped to partition MaxTime. Preparation requests 1 CPU/96 GiB, controller 1 CPU/4 GiB, checking 4 CPUs/64 GiB.

The scheduler:
- Freezes the campaign and preflights formulas/data on a compute node.
- Prioritizes 20% anchors and quick prior checks.
- Records accepted job IDs immediately and uses locks/unique names.
- Recovers uncertain submission acknowledgments by exact job identity or stops for inspection.
- Retries explicit submission-capacity rejections without changing scientific settings.
- Checks completed fits and refreshes summaries.
- Arranges a successor controller near walltime; if blocked by quota, may requeue only the lightweight controller.
- Does not automatically repeat accepted model fits or cancel unrelated jobs.

Array throttling alone would not bypass a MaxSubmitJobs limit. The dispatcher therefore limits actual submissions.

## Resume versus create a new campaign

The historical queue_shape_dcm.sh installs a versioned package and normally resumes its matching frozen campaign. SHAPE_NEW=1 explicitly requests a new campaign. Do not use it to recover an existing fit.

Relevant controls are DCM_ROOT, SHAPE_SOURCE, SHAPE_MAX_ACTIVE, SHAPE_MAX_CPUS, SHAPE_PARTITION, and SHAPE_CONSTRAINT. Verify actual accepted names in scheduler.py before changing them. Default source resolves through latest_followup150_experiment.txt.

Before any restart:
1. Confirm the latest pointer, existing controller state, ledger, and locks.
2. Inspect active fit/check jobs. Never clear a lock simply because it exists.
3. Compare the frozen campaign source with its recorded identity and any applied repairs.
4. Preserve the collector identifier fix described below.
5. Use the existing documented resume path only after establishing that another controller is not active.

The root launchers can submit expensive work. They are not monitoring commands. Historical installers may refuse to overwrite an existing package; that protects the frozen release.

## Collector identifier repair

The preserved original shape collector contains:

~~~r
key<-paste0(basename(campaign),":",sprintf("%03d",as.integer(id)))
~~~

Threading IDs are strings, so the live collector was repaired to preserve them:

~~~r
id_text <- as.character(id)
id_key <- if (grepl("^[0-9]+$", id_text)) {
  sprintf("%03d", as.integer(id_text))
} else {
  id_text
}
key <- paste0(basename(campaign), ":", id_key)
~~~

This restores separate C4_T0/C4_T2/C4_T4/C4_T8 rows and the 586 count. The original package in this repository remains unchanged for checksum fidelity. Do not overwrite a patched live collector with it. The [repair utility](../scripts/repair_shape_collector_ids.py) inspects by default; explicit --apply creates a backup and changes only this exact known expression. It submits no jobs. Do not rewrite frozen package hashes to hide a change.

## Older evaluations and repairs

DCM_Ibex_targeted_tests/rescore_completed.sh reuses eligible completed fits to submit predictions and a collector. It is not purely read-only, even though it performs zero refits. Eligibility audits raw data/membership and excludes incomplete/ineligible fits.

repair_targeted_dcm_dependencies.sh repairs the documented interrupted scheduling state and preserves registered fits. repair_onset_checks.sh installs a support repair and schedules evaluations of existing fits. Review their target paths, ledgers, and current completion state before reusing them; they were designed for known historical failures.

## Recovery evidence to collect

For a failed stage save job ID/state/exit code, .out/.err tail, its failure marker, run config, exact campaign, data identity, and current ledger. For a fit retain complete/partial sampler CSV metadata. For a check retain the fit and repair only the scoring cause where possible.

Never infer a scientific failure from a zero-row report crash, unsupported profile onset, missing monitor path, or integer-coercion warning. Never relabel an unconverged fit as valid merely because Slurm exited zero.

## Backup and handover

Preserve raw inputs, prepared RDS files, scaling/memberships, hydrated fits, chain CSVs, frozen source, ledgers, and exports outside Git. Git holds compact source and documentation. A compact comparison archive alone cannot reproduce posterior predictions without saved fits and data. Record UTC (Ibex display is commonly UTC+3) to avoid confusing local log timestamps.
