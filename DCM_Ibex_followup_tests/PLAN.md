# Focused DCM follow-up: 18 October 2026

This campaign has 30 posterior fits and four prior-predictive checks. It follows the comparison of all 187 completed posterior fits, rather than selecting runs by RMSE alone. Older campaigns and their pending evaluations continue to be included in monitoring and result collection.

## Equations and default settings

The anchors are the original mean equation without MLD (`original_z`), shared depth-varying effects of MLD and distance (`both_z`), and those effects plus a shared depth-varying product interaction (`interaction_z`). All retain Region, the regional depth smooth and the nested date/ENSO random intercepts. Predictors use the existing training-derived standardization. The interaction multiplies standardized MLD and distance; it is not a general three-dimensional tensor surface.

Defaults are four chains, 1,000 warmup plus 2,000 retained iterations per chain (iter=3,000; 8,000 retained draws total), adapt_delta=.95, max_treedepth=14, two threads per chain, k_depth=7 and numeric-by k=5. The regional residual model is `sigma ~ 0 + Region`. Proper priors match the targeted campaign: mean Intercept student_t(3,.73,5.448), slopes normal(0,2), mean group/smooth SD student_t(3,0,1.5), and log-sigma coefficients normal(0,1).

| Cases | Runs | Question and direct controls |
|---|---:|---|
| 1-6 | 6 | Each anchor with seeds 1234 and 20261018. Compare equations with equal settings and inspect seed sensitivity. |
| 7-9 | 3 | Change only adapt_delta to .90 relative to cases 1/3/5. |
| 10-12 | 3 | Change only threads per chain to four relative to cases 1/3/5. |
| 13-15 | 3 | Change training from 108 to 216 casts (10% to20%); same frozen scaling and validation. |
| 16-18 | 3 | Both-varying equation: regional sigma control, regional depth-varying sigma, then stronger sigma-smooth SD prior. All use 2,000 warmup,2,000 retained and delta=.99. Compare17 versus16 for structure,18 versus17 for the prior. |
| 19-21 | 3 | The same variance comparisons for the interaction equation. |
| 22-24 | 3 | All anchors on the existing disjoint nested-date split. Compare equations within this cohort. |
| 25-27 | 3 | Alternative108 training casts: primary20% minus primary10%. Same size and predictor scaling, different composition. |
| 28-30 | 3 | Interaction equation: depth k=10; numeric-by k=3; numeric-by k=7. Each compares with case5. |
| 31-34 | 4 | Prior-predictive checks for both/interaction equations and sigma-smooth SD half-normal scales .5/.25. No accuracy ranking. |

The sigma smooth is `s(depth, by=Region, bs='cr', k=5)` in log-sigma. Its positive SD prior is half-normal(0,.5) or half-normal(0,.25). No new50% fits, large tensor expansions, cast curves, Student likelihoods or AR fits are added: the current50% jobs and their baseline evaluations must inform those decisions first.

## Data and comparisons

The frozen source is the successfully prepared `targeted150` campaign. Main training sets have108/216 casts of1,080; the original primary holdout has108 and the date-blocked holdout117. Training and validation contain entire casts. Responses remain unsmoothed and unbinned. All new fitting and scoring is confined to10-150m.

A new108-cast reserve holdout is selected by fixed seed and predictor/Region metadata from casts outside primary50% training, date-blocked training, and both previous holdouts. It is shared by all new fits; each split uses its existing scaler. It is independent of these new fits' training, but reserve results will be used for development. Historical full-data fits cannot be assigned independent reserve-test accuracy. Preserve a further untouched test or use outer grouped cross-validation for final claims.

Every posterior fit is scored in conditional and population modes on training, existing validation, and reserve data. The20% training check uses its original108-cast training core so training scores remain comparable; alternative10% uses its own108 training casts. These training checks are labelled separately from validation.

Scores use up to2,000 sampled posterior draws. Full scores and two1,000-draw batch summaries reveal prediction Monte Carlo sensitivity without extra posterior-predict calls. Batch differences are diagnostics, not formal MCSEs. No held-out chlorophyll values are supplied to prediction functions. New random-effect levels are simulated; population curves set group effects to zero.

Keep every case, including failed convergence screens. Compare cast-balanced RMSE, bias, profile shape, observed-peak errors, positive-enhancement area, regional and depth calibration, proper90% interval scores, and runtime/core-hours. Primary and reserve comparisons require identical validation hashes and cast membership. Date-blocked versus primary comparisons change training and validation cohorts and cannot isolate a single date effect.

The collector reports all matching conditional-cast pairs, planned controls, regional contrasts and a500-resample Region-stratified nested-date cluster bootstrap for planned contrasts. This is exploratory uncertainty on a development set, not a selection-adjusted final interval. Very small differences need stable sampling and predictive batches before interpretation.

## Execution and recovery

Submit all34 fits concurrently subject to normal Slurm limits. Physical CPU requests are8 per default fit and16 for cases10-12. Memory is32G or64G; posterior walltime is7days (prior checks2hours). Checks request1CPU/32G/8hours. These are resource limits, not runtime forecasts.

An Ibex preparation job validates all34 formula/prior/data combinations before fitting is released. Every case saves warmup and has a separate diagnostics/prediction job. Final collection waits for active recorded fits, checks and older evaluations across all campaigns. It does not wait for other collectors, avoiding a dependency cycle.

Repeat the same installer or submit helper to resume registration. A lock prevents concurrent submitters. Accepted job IDs are recorded immediately and recovered by unique job names if needed. Failed fits are not automatically rerun. An uncertain submission outcome stops safely instead of risking duplicate work. Completed preparation does not keep stale Slurm dependencies. Expired fit dependencies are omitted only after Slurm proves the fit terminal. Explicitly setting FOLLOWUP_NEW=1 starts a genuinely new experiment; normal repetition does not.

Useful commands, from `/ibex/user/meleg/DCM_MODEL`:

```bash
# Include all older and new runs; inspect warmup and sampling separately.
bash DCM_Ibex_followup_tests/monitor_all_dcm_v6.sh --diagnostics

# Refresh existing compact summaries without fitting or prediction.
bash DCM_Ibex_followup_tests/collect_results.sh

# Package compact results from every campaign for comparison.
bash DCM_Ibex_followup_tests/export_results.sh
```

Scientific interpretation follows [Stan diagnostics](https://mc-stan.org/learn-stan/diagnostics-warnings.html) and the [brms distributional-model documentation](https://paulbuerkner.com/brms/articles/brms_distreg.html). Fit settings are experimental configurations, not established optima.
