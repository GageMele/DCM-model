# DCM vertical-profile campaign: 102 fits

This package submits **96 posterior fits and 6 prior-predictive checks**, followed
by a separate checking job for each fit and a final collector: **206 Slurm jobs**.
It creates a new `extensive150_*` experiment and preserves existing campaigns.
All fits become eligible together after one preparation job passes. Slurm controls
how many actually run concurrently, according to available resources and account limits.

## Submit on Ibex

Upload `DCM_Ibex_extensive_tests_150m.zip` into `/ibex/user/meleg/DCM_MODEL`, then:

```bash
cd /ibex/user/meleg/DCM_MODEL
unzip DCM_Ibex_extensive_tests_150m.zip
bash DCM_Ibex_extensive_tests/submit_tests.sh
```

The submitter uses `logs/latest_baseline150_experiment.txt` to locate the frozen
full input and the already prepared 10%, 20%, and 50% bundles. It reuses the
environment script from that experiment. The baseline fits do not need to finish:
only baseline preparation must have completed, as it had in the latest log.

The submitter first runs a mock submission test, shell/R syntax checks, verifies
the 102-row matrix, and asks Slurm to check the largest resource request. It then
submits preparation, all fits, their checking jobs, and the collector. Preparation
builds Stan data/code and resolves priors using the installed Ibex brms version;
no fit starts unless that succeeds. All submitted job IDs are recorded immediately.

Fit resources: 4, 8, or 16 CPUs, 32 GiB at 10%/20%, 64 GiB at 50%, and a seven-day
limit per fit. The sum across all 102 fit requests is **808 CPUs and 3,392 GiB**;
these are requested resources if all fits overlapped, not a promise of simultaneous
allocation or an estimate of consumed memory. Most fits request 8 CPUs. Checking
jobs request 1 CPU/32 GiB/4 hours. Compilation is included in reported fit time.
These limits can be changed explicitly with `EXT_TIME` and `EXT_CONSTRAINT`.
Default CPU constraint: `intel&cascadelake`.

Run the submit command once. If submission stops partway, the existing
`manifest.tsv` and `job_events.tsv` identify exactly which jobs were accepted;
starting another experiment would duplicate those accepted cases.

## Monitor and collect

```bash
cd /ibex/user/meleg/DCM_MODEL
bash DCM_Ibex_extensive_tests/monitor_all_dcm_v4.sh
# Detailed CSV scan, including warmup and post-warmup tree hits separately:
bash DCM_Ibex_extensive_tests/monitor_all_dcm_v4.sh --diagnostics
```

The monitor also includes the pilot, threading benchmarks, fast150, k150,
equation150 and baseline150 campaigns. The normal snapshot is cheaper than the
full CSV scan. For a periodically refreshed terminal view:

```bash
watch -n 120 bash DCM_Ibex_extensive_tests/monitor_all_dcm_v4.sh
```

The final collector runs automatically. To refresh partial results without
refitting any model:

```bash
cd /ibex/user/meleg/DCM_MODEL
export DCM_ROOT="$PWD"
export EXPERIMENT_DIR="$(cat logs/latest_extensive150_experiment.txt)"
source "$EXPERIMENT_DIR/ibex_environment.sh"
Rscript --vanilla "$EXPERIMENT_DIR/04_collect_results.R"
```

See `PLAN.md` for equations, design rationale, interpretation and next decisions.
`cases.tsv` is the exact executable matrix. `basis_counts.csv` and
`split_counts.csv` contain locally verified actual-data construction counts.
`VALIDATION.md` states what was tested and what requires execution on Ibex.

## Main outputs

Under the new experiment:

- `preflight/`: generated Stan programs and resolved priors for distinct designs.
- `scaling.rds`, `validation_casts.csv`, `data_checksums.csv`: frozen predictor
  transformations, validation identities and source hashes.
- `runs/case_XX/fit.rds`: portable, hydrated brms draws; raw CmdStan CSVs retained.
- `runs/case_XX/checks/parameter_diagnostics.csv`: R-hat, ESS and mean MCSE.
- `sampler_summary_by_phase.csv`: per-chain warmup/sampling divergences,
  tree hits and E-BFMI; original CSVs support the live recent-50-draw monitor.
- `prediction_scores.csv`: cast-balanced training/validation scores, both
  conditional and population predictions.
- `*_cast_metrics.csv`, `*_by_region.csv`, `*_strata.csv`: errors, shape errors,
  intervals and residual correlation, stratified by Region, depth, MLD, distance,
  and whether the date/ENSO group was present in the training subset.
- `*_observations.csv`, `profile_checks.pdf`: predictions at observed depths and
  observed-versus-predicted profiles with 90% predictive intervals.
- `results/predictive_ranking.csv`: completed posterior candidates, convergence
  screens first, then validation error. Tree hits remain a separate column.
- `results/matched_contrasts.csv`: sensitivity comparisons against each model's
  k=7 reference, including hardware matching and CPU-cost ratios.
- `results/paired_vs_original_k15.csv`: changes in errors on identical casts
  relative to the new original-equation k=15 reference.
- `results/case_status.csv`: all cases, including incomplete and failed checks.

Prior-only runs have separate predictive quantile outputs and are excluded from
the predictive ranking. A successful diagnostic screen alone does not establish
scientific adequacy.
