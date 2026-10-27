# Evaluation and fair comparison

## Keep every model in the register

Preserve complete, incomplete, failed, prior-only, reused, and diagnostically invalid entries. A completed posterior fit may be awaiting checks; a checks_failed marker can reflect an evaluation bug rather than failed MCMC. Exclude inappropriate rows from a particular numerical comparison, while retaining their status and exclusion reason in the registry.

A common case number does not identify a model across campaigns. Use campaign-qualified strings. Do not coerce threading identifiers C4_T0/C4_T2/C4_T4/C4_T8 to integers.

## Cohorts and prediction targets

Training samples contain whole casts. Validation must exclude those casts; date-blocked validation additionally excludes its nested date groups from training. Do not randomly split rows from one profile across training and validation.

For below-100% shape runs, the primary 10% validation and 10% reserve cohorts are excluded from training. At 100%, those casts are included, and outputs are training_full, reconstruction_primary, and reconstruction_reserve. They are reconstruction, never held-out validation.

Conditional predictions retain group effects; new levels use the specified Gaussian group-effect sampling. Population predictions with re_formula = NA set group effects to zero. They are not marginal predictions integrating over unknown date/ENSO effects.

The current scoring batches observations by Ocean–Region–date so a newly drawn date effect is shared inside that batch. These marginal score calculations are not a full joint predictive distribution over all unseen groups. Do not infer complete joint uncertainty from them.

Construct prediction frames from predictors, excluding held-out chl_norm. Use posterior_epred for expected means and posterior_predict for observation-level intervals. The uncertainty of an expected profile is not the uncertainty of a new chlorophyll observation.

## Basic error and calibration

For each cast c with n_c observations, MSE_c is the mean squared residual. The primary aggregate is:

RMSE_cast = sqrt(mean_c(MSE_c)).

This gives each cast equal weight. It is not pooled row RMSE or the arithmetic mean of per-cast RMSE. Save row-weighted and cast-weighted quantities under distinct labels.

MAE uses analogous cast-level summaries. Shape RMSE first subtracts each cast's mean from its observed and predicted profile, isolating deviations around the cast offset. A mean-offset correction can improve whole-profile RMSE without appreciably improving shape RMSE.

For 90% predictive intervals [L,U], interval score is:

(U−L) + 20 max(L−y,0) + 20 max(y−U,0).

Lower is better. Always inspect interval width and empirical coverage together. Narrow intervals can look attractive while undercovering.

## Profile-level outcomes

Save peak depth and height, depth-integrated positive enhancement, conditional/population scores, regional summaries, and individual-cast errors.

The original profile helper defines a DCM class by peak depth >15 m. It does not automatically impose an enhancement-amplitude threshold. Later onset/width rules use separate thresholds. Boundary maxima are flagged, because a curve still rising at 150 m does not establish a true interior maximum.

Positive area uses trapezoidal integration and can bridge observed gaps. Interpret it with vertical support. Adjacent-row residual lag-one correlation is an observational dependence diagnostic, not MCMC chain autocorrelation.

## Onset definition and missing vertical support

Onset thresholds are +.25, +.5, and +1 enhancement, corresponding to 1.25×, 1.5×, and 2× surface chlorophyll. The detector searches at and below 15 m depth on a supported 1 m grid, using a 5 m running median and five consecutive points satisfying the threshold.

It does not extrapolate outside observed support or bridge gaps greater than 3 m. It distinguishes detection, already-above-threshold/left censoring, gap censoring, not reached, and insufficient support. A missing onset is not assigned 150 m.

The reserve singleton cast has only one retained depth (10 m). It remains in ordinary point errors, while onset/width is unidentifiable. This was the cause of the earlier common checking failure and must not be “fixed” by dropping the cast from all evaluation.

Five-metre scoring medians do not alter unsmoothed model-fitting observations. Keep observed and predicted diagnostic processing consistent.

## Sharpness and transition width

The shape campaign measures maximum positive 4 m centered slope and the depth of strongest rise. The rapid-rise rule uses the recorded .02 slope and .5 enhancement criteria.

Transition width measures the rising segment from 25% to 75% of a profile's own positive enhancement peak, with three-point persistence. Low enhancement, insufficient support, boundary peaks, and censoring are explicit. This is an operational diagnostic, not a universal biological DCM definition.

Expected model curves are smoother than individual noisy profiles. Inspect posterior predictive profiles as well as expected curves when deciding whether a slope discrepancy indicates missing mean structure, observation variability, or both.

## Draw selection and reproducibility

Current point/interval scoring uses up to 2000 posterior draws; earlier checks often used 1000. Onset uses 400 expected-curve draws; MLD slices use up to 1000. Exact choices and seeds live in the preserved code.

Two 1000-draw evaluation batches in a follow-up study are a sensitivity check, not a formal Monte Carlo standard-error analysis. Keep posterior-draw subsampling distinct from MCMC chain length and effective sample size.

## Paired contrasts

Compare the same response definition, depth interval, raw/prepared data version, training membership, evaluation cohort, prediction mode, and outcome definition. Identify every exception explicitly.

For ordinary point errors retain all eligible casts. For onset/width compare common finite detections in both models and report the denominator, non-detections, and censoring separately. A favorable average on a selectively detected subset is insufficient by itself.

The October 25 analysis used 5000 paired cast-bootstrap replicates. It did not cluster those replicates by date. Earlier onset analyses also used date-clustered procedures; distinguish them when comparing uncertainty intervals.

Inspect influential-cast contributions and Region-level changes. Do not remove difficult Red Sea casts to improve a primary score. The reserve has repeatedly informed model development and cannot serve as an untouched final holdout.

## MLD profile slices

Plot expected profiles at regional MLD 10th/50th/90th percentiles, holding regional distance fixed and specifying group-effect handling. Use saved training scaling and factor levels.

Pairwise training-hull flags indicate some unsupported combinations but do not prove support in the full joint predictor space. A slice is an associative model summary, not evidence that changing MLD experimentally causes the predicted deformation.

## Reporting standard

For every completed fit report configuration, training/evaluation identity, diagnostics, fit cost, overall/region/cast errors, shape/onset/peak/width/slope, uncertainty, and support flags. Show all cases, with filters for valid comparisons, rather than publishing only an error-ranked shortlist.
