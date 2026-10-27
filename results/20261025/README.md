# DCM model comparison — partial campaign review, 25 October 2026

## Main conclusion

Continue the regional MLD-by-depth model, the half-shift model and the baseline controls. The regional MLD model is the strongest completed **20% primary-training** candidate for overall conditional RMSE, especially in Central/South Red Sea. It is not a universal winner for onset, peak depth or sharpness. Half-shift remains useful for onset. Shared MLD effects and the prescribed sharpening warps are less convincing. The current regional tensor fits and all new 40–100% training cases are unfinished in this snapshot.

No model has been fitted, re-fitted, submitted, cancelled or changed for this review.

## Evidence and scope

Source: `dcm_shape_comparison_20261025T111507Z_160796.tar.gz`, especially its frozen `cases.tsv`, `models.R`, per-run checks and `results/` tables. This is an archive snapshot, not a live Ibex query.

- 586 registered entries across campaigns; 356 completed entries. Registrations include reused controls and prior-only cases; these are not 356 independent posterior fits.
- Current campaign: 282 cases; 67 completed posterior fits, 12 completed prior checks, 48 accepted fitting stages and 155 unsubmitted cases at the reported dispatcher timestamp.
- The eligible source contains 1,080 casts across eight regions, with existing preprocessing and the 10–150 m fitting interval. Primary 20% training uses 216 casts / 30,455 observations. Primary validation and reserve each contain 108 whole casts.
- All 67 completed new posterior fits are retained in `all_67_completed_posterior_runs.csv`; all current cases and every older registration remain in separate complete tables, including incomplete cases.
- Older runs are reviewed from the combined aggregate summaries. The archive does not contain complete raw prediction evidence for all older runs, so a new cast-by-cast audit of every old model is not claimed.
- Recalculated 2,412 aggregate values; maximum discrepancy 5.51e-14. Recalculated all 344 planned matched-cast comparisons. Verified identical primary validation and reserve cast sets across the four 20% cast selections.
- These repeatedly inspected cohorts are development evaluation sets. They are not an untouched final test.

## Like-for-like primary 20% comparison

Each row below uses the same training data and validation/reserve casts. Residual SD varies by Region in every row. Thus “baseline” means the original depth **mean equation**, not every setting of the original historical fit. RMSE is the square root of the mean of within-cast mean squared errors: every cast gets equal weight.

| Case | Mean equation family | Validation RMSE | Reserve RMSE | Fit wall hours | Divergences | Tree hits |
|---|---|---:|---:|---:|---:|---:|
| 19 | Baseline regional depth + linear distance | 1.540 | 2.739 | 19.7 | 0 | 0 |
| 22 | Shared MLD/distance depth interactions | 1.534 | 2.736 | 35.7 | 0 | 8 |
| 25 | Half MLD shift | 1.489 | 2.648 | 27.2 | 0 | 2 |
| 28 | Full MLD shift | 1.505 | 2.624 | 25.7 | 1 | 0 |
| 31 | Shared MLD-by-depth effect | 1.537 | 2.738 | 19.9 | 0 | 0 |
| 34 | Region-specific MLD-by-depth effect | 1.409 | 2.559 | 33.2 | 0 | 0 |
| 129 | Regional MLD-by-depth; k=15 Central/South RS, k=7 elsewhere | 1.401 | 2.506 | 33.9 | 0 | 0 |

Case 34 reduces validation RMSE by 8.5% and reserve RMSE by 6.6% relative to case 19. It improves within-cast MSE for 54% of validation casts and 62% of reserve casts. Its overall improvement is influenced by a few high-error casts, rather than uniformly distributed: three casts contribute about 56% / 67% of the total absolute paired MSE changes in validation / reserve. The exploratory paired cast bootstrap interval for its RMSE difference includes zero in both cohorts; this is not conclusive evidence of superiority. The bootstrap does not account for shared sampling dates or model-selection uncertainty.

The four 20% selections retain the same validation and reserve casts and MCMC seed. Across these selections, average validation / reserve RMSE is 1.613 / 2.846 for the baseline and 1.532 / 2.686 for regional MLD. Regional MLD improves reserve RMSE in all four selections but validation RMSE in only three. Two of its repeat-selection fits have 1 and 4 divergences. Do not average these scores into an ensemble prediction: they are means of separate model scores.

## Shape, onset and sharpness are different objectives

All comparisons below use casts eligible in BOTH compared models. Missing/censored onsets and widths remain explicit; they are not scored as zero error or discarded silently.

- **Half shift, 10% (case 7 vs 1):** onset absolute error decreases by 4.77 m on 60 common validation casts and 3.32 m on 59 common reserve casts. Width error decreases by 1.85 m / 1.49 m on 46 / 45 common casts. Peak-depth error increases by 2.69 m / 2.35 m. Maximum-rise slope error changes little and slightly worsens.
- **Half shift, 20% (25 vs 19):** onset error decreases by 2.03 m / 0.76 m on 64 / 59 common casts. Width error decreases by 1.80 m / 1.46 m on 54 / 50 common casts. Peak-depth error increases by 3.23 m / 2.60 m.
- **Regional MLD, 20% (34 vs 19):** onset error decreases by 1.62 m on 60 common validation casts but increases by 1.05 m on 56 common reserve casts. Width error increases by 0.49 m / 0.38 m on 53 / 50 common casts. Slope error changes by +0.0123 / -0.0065 chl_norm units per metre. Overall RMSE gains therefore do not establish that the onset or steepness problem is solved.
- **De-meaned profile shape RMSE:** case 19 → 34 changes 1.279 → 1.271 on validation, versus 2.230 → 1.970 on reserve. The validation gain in total RMSE is largely a reduction in profile-level offset error, not a large reduction in this de-meaned shape error.
- **Maximum-rise underprediction:** on the same 70 validation and 71 reserve eligible casts, median per-cast predicted/observed maximum-slope ratios are 0.20 / 0.24 for baseline, 0.20 / 0.22 for half shift, and 0.21 / 0.29 for regional MLD. With higher Red Sea k the ratios are 0.24 / 0.29. These operational maximum-slope metrics use a 1 m grid, 5 m running median and 4 m centred differences. They are not necessarily the slope immediately before the main DCM peak. Expected curves naturally smooth noisy individual profiles; the result still warrants checking transition sharpness directly.
- Width uses the distance between 25% and 75% of each profile's own positive enhancement peak. Cases with left-censoring, insufficient enhancement or boundary peaks have fewer comparable widths. A lower width MAE on a different set of casts is not directly a better model.

## Regional findings and k

The regional MLD effect helps chiefly in Red Sea casts. For case 19 → 34, Central RS RMSE changes from 3.720 → 3.083 on validation and 7.925 → 7.248 on reserve; South RS changes 2.497 → 2.377 and 3.721 → 3.225. Some Pacific regions show little gain or modest worsening. There is no justification for claiming one universal MLD response across all regions.

Case 129 adds k=15 in Central/South RS while keeping other regions at 7. Relative to case 34, South RS RMSE improves 2.377 → 2.195 on validation and 3.225 → 2.958 on reserve. Central RS changes 3.083 → 3.134 (slightly worse) and 7.248 → 7.082 (better). Its small overall advantage does not establish that higher k is always needed in Central RS. The basis dimension is a flexibility limit, not the fitted effective degrees of freedom.

At 10% training, increasing the half-shift basis from 7 to 11 improves validation RMSE 1.400 → 1.379 and reserve RMSE 2.978 → 2.959; k=15 gives 1.381 / 2.972. That supports retaining k=11 but does not demonstrate a universal optimum. The older `k150` screen has R-hat approximately 1.84–3.03 and severe tree-depth limitations; it cannot establish an optimal k.

`mld_shape_comparison.png` shows three models fitted to the same primary 20% subset. Solid lines are inside the saved pairwise depth–MLD and MLD–distance training hulls; dashed lines are outside at least one hull. Shading is a 90% expected-mean interval, excluding residual noise and date/ENSO effects. In Central RS, some extreme MLD slices at fixed median distance are unsupported: their dramatic peaks must not be interpreted as established biological responses. Pairwise hull inclusion is itself weaker than full multivariate support. The curves describe conditional associations, not causal MLD effects.

## Fixed sharpening versus learned shape change

The prescribed warp models multiply half-shift depth by exp(lambda * standardized MLD). Positive lambda narrows profiles as MLD grows; this is imposed rather than learned freely.

Relative to half shift case 25, lambda=+0.25 (case 221) changes validation RMSE 1.489 → 1.476 but reserve RMSE 2.648 → 2.745. Lambda=+0.50 (case 223) gives 1.580 / 2.846 and worse matched-cast slope errors. Lambda=-0.25 (225) gives 1.556 / 2.636. The present data do not support imposing a universal deeper-MLD-means-sharper-profile rule. Regional varying-coefficient models and pending tensors are more appropriate ways to learn the association.

## Sampler settings, priors and computational cost

All 67 completed new posterior fits have R-hat <1.006, minimum bulk ESS >1,300, minimum tail ESS >1,800 and minimum E-BFMI >0.62. Of these, 56 pass the convergence screen, while 27 pass the stricter screen that also requires zero maximum-depth hits. Eleven fits have divergences. A small number of tree hits is an efficiency warning distinct from the inferential concern raised by divergences. Neither good R-hat nor a zero-hit result proves the scientific mean or likelihood is adequate.

For half/full shifts, changing delta .95/.97/.99, doubling warmup, or doubling retained draws changes RMSE very little. For example half-shift validation RMSE stays approximately 1.488; delta=.99 has 136 tree hits versus 2 at the baseline. Higher delta is not automatically faster or less tree-limited.

A defensible provisional setting is four chains, 1,000 warmup plus 2,000 retained draws per chain, delta=.97, tree depth 14. Use diagnostics to decide extensions; do not shorten MCMC merely because there are more observations. Existing .95 runs may offer a useful speed tradeoff when diagnostics pass. There is no current evidence that all runs need higher tree depth or longer warmup.

For the regional MLD equation, case 241 at delta=.95 fits in 28.1 h with zero divergences/tree hits and essentially the same RMSE as case 34's 33.2 h at .97. Case 249 with tighter normal(0,1) fixed/smooth-population coefficients fits in 24.1 h and gives RMSE 1.406 / 2.553; zero divergences/tree hits. This small predictive change supports further prior regularization comparisons, not a definitive prior winner. All three run on the same reported Xeon Gold 6248 CPU model and 16 allocated CPUs, but on different nodes; these are observed walltimes, not replicated speed benchmarks.

Prior-only simulations are very broad: baseline central 90% predictive range approximately -19 to +20 chl_norm units, with regional MLD simulations broader. Twelve prior checks are not posterior accuracy results. Gaussian posterior prediction also puts approximately 6–7% of simulated responses below -1 in several anchors; the mean and observation-level predictive distribution need separate scrutiny if -1 is the physical lower bound of the enhancement response. Preserve the original response definition when evaluating any alternative likelihood.

## Generalization and residual limitations

Conditional prediction includes estimated date/ENSO offsets where levels were seen; unseen levels are integrated using simulated group effects. Population prediction sets these group effects to zero. For the primary 20% baseline → regional MLD comparison, validation population RMSE changes only 1.606 → 1.586, much less than the conditional 1.540 → 1.409. Reserve population RMSE improves 3.087 → 2.811. The model preferred for reconstructing known-date casts need not be the best for a new-date population prediction.

On the two date-blocked selections with completed comparisons, half shift beats baseline validation RMSE: 2.002 → 1.963 and 2.127 → 2.046. Full shift is slightly worse on these validation cohorts (2.033 and 2.155) despite better reserve scores. Regional MLD date-blocked fits are not complete. Do not compare these cohort RMSE levels directly with primary-cohort scores as if only the split method changed.

Residual lag-one correlations remain around 0.97 for the primary anchors. Much within-profile structure remains; both mean misspecification and residual dependence should be investigated before treating 1 m observations as independent information. Good sampler diagnostics do not resolve this. Do not claim calibrated onset uncertainty solely from aggregate observation-level 90% coverage.

## Provisional choice and next decisions

For reconstruction of held-out casts in represented dates, prioritize case 34 / 129 regional MLD-by-depth, with case 249 as a promising more regularized version. Keep half shift as a strong onset comparator and the baseline for peak location. The learned equation is:

`mean = Region intercept + distance_z coefficient + regional depth smooth + mld_z * regional depth-effect smooth + nested date intercept + nested ENSO intercept`, with `sigma ~ 0 + Region`.

At any fixed depth the added MLD effect is linear in standardized MLD, while its coefficient varies smoothly with depth and independently by Region. It can therefore change shape as well as level. It does not impose a universal sharpening response. The exact implemented formula construction is preserved in the archive's `models.R`.

Keep the already-submitted tensors, higher-k and larger-data tests running. Do not select a final model before the date-blocked regional MLD and tensor results, the 40–100% fits, and explicit onset/sharpness comparisons are available. No completed new case yet answers whether 40%, 50%, 75% or 100% is preferable. Once the mean structure is clearer, the strong residual dependence and Gaussian predictive support deserve targeted investigation.

## Files

- `all_586_registrations.csv`: every registered old/new row, unfiltered source fields.
- `all_282_new_cases.csv`: complete current campaign including unfinished cases.
- `all_67_completed_posterior_runs.csv`: settings, diagnostics and multiple performance measures for EVERY completed new posterior fit.
- `12_prior_predictive_checks.csv`: prior-only cases, explicitly separate.
- `paired_comparisons_common_casts.csv`: all 344 recomputed planned contrasts; onset/slope/width common-pair counts included.
- `four_20percent_cast_selections.csv`: all four selections of six core equations.
- `sharpness_common_casts.csv`: matched-cast maximum-slope comparison.
- `exploratory_cast_bootstrap.csv`: descriptive paired-cast resampling; not an untouched final-test analysis.
- `cast_selection_comparison.png`, `mld_shape_comparison.png`: scientific comparison figures.
- `audit.json`: scope, counts and numerical reconciliation.
- `analyze_shape_review.py`: annotated reproduction script; requires an extracted experiment directory, numpy, pandas and matplotlib. It never fits or submits models.
