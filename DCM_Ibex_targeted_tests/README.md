# Targeted DCM tests - 16 October 2026

This package prepares **63 new fits: 57 posterior fits and 6 prior-predictive fits**. It also queues prediction checks for every currently completed, compatible posterior fit from the fast150, k150, equation150, baseline150 and extensive150 campaigns. Poor diagnostics do not exclude a completed fit from these checks. Earlier fits are read without refitting or modifying them.

## Queue the batch

Place `queue_targeted_dcm.sh` in `/ibex/user/meleg/DCM_MODEL`, then run on Ibex:

```bash
cd /ibex/user/meleg/DCM_MODEL
bash queue_targeted_dcm.sh
```

The installer extracts this package, verifies its checksums and invokes `submit_tests.sh`. It needs the existing `logs/latest_extensive150_experiment.txt` pointer and that campaign's frozen full input, p10/p20/p50 data bundles, validation data and environment script. It creates a fresh `experiments/targeted150_<UTC>_<pid>` directory and records every submitted job ID. It leaves existing fits and queued jobs in place. An existing package directory stops installation; do not repeat the submission command unless you intend to start another 63-fit experiment.

Every fit waits for a preparation job that verifies the data split, priors and Stan data/code construction for **all 63 cases**. All fit jobs can then run independently, subject to Ibex scheduling and account limits. Fit requests use four parallel chains, 4/8/16 CPUs, 32 GB or 64 GB, and a seven-day walltime. The large AR and 50% cases request 64 GB. Prediction checks use one CPU and 32 GB. Hardware is recorded, and fit jobs request the same CPU architecture (`intel&cascadelake`).

## Data policy

Percentages refer to the full set of 1,080 casts, keeping every eligible observation in each selected cast. No new model trains on 5%, and no 10% subset is split into two 5% pieces.

| Split | Training casts | Independent validation casts | Training observations | Validation observations |
| --- | ---: | ---: | ---: | ---: |
| Primary 10% | 108 | 108 | 15,228 | 15,226 |
| Primary 20% | 216 | Same 108 | 30,455 | Same 15,226 |
| Primary 50% | 540 | Same 108 | 76,139 | Same 15,226 |
| Date-blocked 10% | 108 | 117 | 15,228 | 16,497 |

These are counts verified locally against the frozen source; the Ibex preparation verifies them again. The primary training subsets are nested, and validation lies outside the 50% training set. Main fits and all common validation use 10-150 m. The historical depth comparison additionally fits 10-250 m on the same 108 training casts (26,028 observations), then scores only 10-150 m.

The primary holdout tests new casts and may share dates and ENSO states with training. The separate date-blocked tests hold out complete `Ocean:Region:date` groups and still use at least 10% of the full casts for both training and validation. They test new dates, rather than new cruises or chronological forecasting. Their scores use a different cohort and are not paired with the primary scores.

MLD and log-distance scaling uses training casts only. Primary 10% scaling stays fixed for the 20% and 50% comparisons; date-blocked scaling uses its own training data. Response values, depth observations, cast membership and raw covariates are preserved. No response smoothing or depth binning is added.

## What to read

- `PLAN.md`: questions, equations, priors, exact case ranges and interpretation.
- `cases.tsv`: every parameter for every new fit; cases 2 and 8 are the two main controls.
- `ORIGINAL_REVIEW.md` and `original_review/`: diagnostics and plots from your original saved fit, including comparison curves against current completed fits.
- `VALIDATION.md`: completed local checks and the required full Ibex preparation gate.

## Monitor and collect

```bash
cd /ibex/user/meleg/DCM_MODEL
bash DCM_Ibex_targeted_tests/monitor_all_dcm_v5.sh --diagnostics
```

The v5 monitor includes all previous recorded campaigns, the original pilot and threading benchmarks, this new batch, preparation/collection jobs and old-fit prediction checks. It reports chain progress and separates warmup from sampling tree-depth hits and divergences. A terminal disconnection does not stop submitted Slurm jobs.

To collect a partial report without rerunning any model:

```bash
experiment=$(cat logs/latest_targeted150_experiment.txt)
source "$experiment/ibex_environment.sh"
Rscript --vanilla "$experiment/04_collect_results.R" "$experiment"
```

To add common validation checks for older models that finish after this batch was submitted:

```bash
bash DCM_Ibex_targeted_tests/rescore_completed.sh
```

This command queues predictions only. It skips previously registered evaluations and never submits a fit. It should be run after preparation succeeds. Earlier 5% fits remain in the previous all-run comparison but are not used for these new >=10% training comparisons.

## Outputs

`results/all_results.csv` retains every available result, including diagnostic failures. `case_status.csv` retains every new case and registered old-fit evaluation, including incomplete or failed runs. `all_posterior_comparisons.csv`, `prior_checks.csv`, `all_paired_cast_comparisons.csv` and `planned_controlled_contrasts.csv` provide the corresponding comparisons. There is no shortlist based only on error.

Pairwise results require an identical validation-data checksum and identical cast membership. Casts receive equal weight in summary errors. Validation response values are removed from the data passed into prediction functions and used only for scoring. All pairwise differences, including those involving poorly converged fits, carry convergence flags. The paired standard error is descriptive: casts can share dates, so it is not an independent-sample significance test.

Each run's `checks/` directory contains regional/stratified scores, every held-out cast's errors, observed/predicted example profiles, 90% predictive intervals, expected population profiles and uncertainty, expected-peak summaries, ocean profiles and enhancement probability plots. Scoring uses a fixed selection of up to 1,000 posterior draws. Training reconstruction scores use the same core 108 training casts across the primary 10/20/50% fits, so their reconstruction comparison has a common cohort; they are not scores over all 20/50% training rows.

Conditional predictions retain learned group effects and integrate unseen group effects by simulation. Population predictions set all group effects to zero, including any group effects on distributional parameters; they are not a marginal average over new groups. Validation uncertainty and expected-mean uncertainty are labelled separately. Do not compare prior-only simulations to posterior predictive accuracy.

To collect compact reports and plots for comparison:

```bash
bash DCM_Ibex_targeted_tests/export_results.sh
```

This creates a new timestamped archive under `logs/` and prints its full path. Raw chains, saved fit objects and large observation-level tables remain on Ibex.

Optional resource overrides, set before launching, are `TARGET_TIME`, `TARGET_CONSTRAINT`, `DCM_ROOT` and `SOURCE_EXPERIMENT`. `TARGET_FREEZE_ONLY=1` copies and verifies inputs without submitting jobs; it is not a full Stan preflight.
