# Compare every DCM run

This collects every registered run, every available completed result, and small supporting tables and plots across all experiment generations. There is no error threshold, diagnostic exclusion, top-N selection, or score ordering. It includes the fast150, k150, equation150, baseline150 and extensive150 campaigns, registered threading benchmarks, and original/unregistered runs with a saved configuration under runs/.

Copy this folder to /ibex/user/meleg/DCM_MODEL, then run:

```bash
cd /ibex/user/meleg/DCM_MODEL
bash DCM_Ibex_all_run_review/run_review.sh
```

Upload the resulting `logs/dcm_all_run_review_*.tar.gz` archive itself. Pasted terminal output only shows that an archive exists; it does not provide its contents. This script reads small existing outputs and does not load saved fits, refit models, submit jobs, or recalculate predictions. It can be run again into a new timestamped snapshot as checks finish.

## Tables

- `all_registered_runs.csv`: every case, including incomplete runs and failed checks. All original summary fields are retained.
- `all_completed_posterior_runs.csv`: every posterior fit with its completion marker, including poor diagnostics and partial checks.
- `all_prior_predictive_runs.csv`: separate scale/plausibility checks.
- `model_and_compute.csv`: equation, data fraction, priors, spline k, sampler settings, hardware and compute cost.
- `sampler_and_precision.csv`: R-hat, ESS, MCSE when available, divergences, depth limits, E-BFMI and recorded screens.
- `prediction_and_profiles.csv`: every available training and validation metric, including conditional and population predictions, bias, profile shape, peak depth/height, interval coverage and width, and residual dependence.
- `all_*by_region.csv`, `all_*strata.csv`, `all_*cast_metrics.csv`: every available regional, subgroup and individual-cast result, with the run identifier attached.
- `coverage_counts.csv`, `missing_or_incomplete_metrics.csv`, `collection_warnings.txt`: coverage and gaps, rather than silently dropping runs.
- `evidence/`: compact original summaries, resolved priors, formulas/Stan code, diagnostics, logs saved inside run directories, and plots. `evidence_index.csv` maps these to their original paths. Large raw draws and models are excluded.

## How the full comparison should work

1. Review all completed fits in campaign/case order, keeping diagnostics alongside scientific metrics. The tables retain runs that need additional sampling or have failed checks.
2. Compare model equations and k settings for profile shape, DCM location/height, bias and calibration by Region, depth, MLD, and distance. An overall RMSE can hide systematic differences.
3. Compare priors, warmup/retained draws, adapt_delta, tree limits, threads and data fractions while checking which other settings changed. Use paired, identical-cast evaluations for equation comparisons; do not infer isolated effects from runs with multiple differences.
4. Compare predictive benefits with wall time, allocated core-hours, ESS and Monte Carlo precision. More data does not replace sufficient posterior sampling.
5. Report every run and every available metric. A shortlist may be useful later, but is a separate view of this complete inventory.

## Important limits

Older campaigns mostly have bounded in-sample predictive checks, while extensive150 has common training and held-out cast evaluations. These scores cannot be combined into one error ranking. Held-out comparisons require the same validation data, prediction mode, response preparation and depth range. Saved-file hashes help identify shared files; a differing hash does not by itself prove differing observations.

This package exposes missing old-run scores. Backfilling them requires reading the saved fits and first verifying that evaluation casts were excluded from each training set. It is additional prediction work, not a model refit. The collector does not pretend those results already exist.

The two historical diagnostic screens differ: older diagnostic_screen_pass included zero tree-depth hits; newer convergence_screen_pass separates tree efficiency. Review their component metrics. Never treat a missing flag as a pass. Prior-predictive runs are retained separately and have no accuracy score.

Completion status comes from saved marker files. It is a snapshot, not Slurm accounting; incomplete runs can include timed-out or cancelled runs. Registered runs outside the three discovery roots or benchmark manifests without resolvable paths are explicitly outside automatic coverage and need manual reconciliation.
