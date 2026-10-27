# What we have learned

**Evidence cutoff: 25 October 2026.** These findings synthesize the supplied comparison archives and the original saved-model review. They are not a live query of Ibex. The full compact review, all 586 registrations, all 282 new cases, all 67 completed posterior rows, planned comparisons, analysis script, and figures are included in [results/20261025](../results/20261025/README.md). The original external source archive is needed to rerun its complete cast-level audit.

## Overall interpretation

Region-specific vertical structure is essential. MLD can help, but how it enters matters. The latest completed results favor region-specific depth-dependent MLD effects for overall conditional prediction error. Half-shift remains useful for onset. Shared linear MLD effects offer smaller and less consistent improvements. A single universally imposed sharpening warp has not worked reliably.

This does not yet select a final model. Sharp transitions, peak depth, residual dependence, negative predictive support, and generalization to unseen dates remain concerns. All new 40–100% cases and regional tensor comparisons were unfinished in the latest archive.

## Original saved final fit

The saved original model used 255,369 observations at 10–250 m, four chains, iter 5000/warmup 1000, delta .97, and tree limit 14. It had no divergences or tree hits; its largest observed post-warmup tree was 13. Maximum R-hat was 1.001655, minimum bulk ESS 3680.014, minimum tail ESS 6854.821, and minimum chain E-BFMI .701435.

It sampled well. That does not prove optimal prediction or identify which notebook preprocessing blocks preceded it. Its data lack MLD and cast IDs, and date is a scaled numeric matrix without the original calendar origin. Historical reconstructions cannot be declared held-out validation.

The original mean log-distance coefficient was about −.21009; residual sigma 1.48825, date-group SD .33434, ENSO-group SD .68584. These belong to that dataset and equation. Its bounded 640-observation training reconstruction RMSE 2.407 and 90% coverage .930 are not comparable to new held-out cohort scores. See the [full original review](../DCM_Ibex_targeted_tests/ORIGINAL_REVIEW.md), including original/current curves and enhancement probabilities.

## Lessons from the early campaigns

| Question | Observed evidence | Consequence |
|---|---|---|
| Were tree limits hit? | Early warmup saturated; later equation B cap 10 hit 997/1000 posterior draws, C 419/1000 | Separate warmup from sampling; raise limits selectively |
| Are short runs enough? | Many two-chain equation fits had R-hat around 1.01–1.03; four-chain longer tests reached about 1.006 | Do not choose equations from short unconverged runs |
| Can MLD replace regional curves? | MLD-only error about 1.925 versus regional structure about 1.463 in the earlier common comparison | Retain regional shape |
| Does regional sigma help? | Additive RMSE about 1.465 → 1.4301; interval score 7.340 → 5.492 | Regional residual SD is a strong baseline improvement |
| Does sigma need a depth smooth? | Interaction interval score 5.440 → 4.913, but RMSE 1.4161 → 1.4442, coverage 90.4% → 85.8%, fit about 25 → 51 h | Calibration/complexity tradeoff, not an automatic improvement |
| Is a large regional 3D tensor best? | Early RMSE roughly 4.69–8.88 and no convergence-screen passes | Those fits are not credible winners; smaller tensors remain open |
| Do cast effects solve generalization? | Cast intercepts did not give a clear held-out gain; cast curves worsened it | In-sample flexibility does not guarantee new-cast prediction |
| Does Student likelihood help? | Estimated nu near 1 caused concerns; fixed nu 4 had worse error and about 82% interval coverage | Retain as negative findings, not default replacements |
| Does AR(1) help? | Matched scalar-sigma RMSE 1.4642 → 1.4430, interval score 7.346 → 7.447 | Modest mean gain, no broad calibration resolution |
| Does depth domain matter? | Matched 250 m versus 150 m bridge RMSE about 1.4408 versus 1.4628 | Deeper observations can alter fitted structure in the shallower domain |
| Did early low-k fits resolve k? | Some k150 R-hats were about 1.84–3.03 | Those runs cannot establish the best k |

These rounded earlier contrasts belong to their original cohorts/configurations; do not compare their absolute values directly to every later table. Prior-only cases 23/24 in the fast campaign were not posterior model successes.

## Shifts and the onset campaign

Full-shift at 10% improved one common primary/reserve comparison from about 1.4272/3.0406 to 1.4010/2.9365. Half-shift showed useful onset gains. Ratio and hybrid coordinates did not give consistent improvements. Shift plus distance had a favorable point error around 1.3904 but six divergences, so it was not a valid winner.

The initial all-model “Insufficient vertical support” checking failure was a scoring bug triggered by a singleton-depth cast. Repaired checks reused existing fits. The failure did not imply that every model had failed to sample. Missing/censored onset remains explicitly recorded after repair.

## Latest primary 20% comparison

All rows below use the same 216 training casts, fixed 108-cast primary validation cohort and 108-cast reserve cohort, regional sigma, and the latest scoring convention. Errors are conditional cast-balanced RMSE. Case numbers refer to shape150_20261022T102851Z_3426481.

| Case | Equation/basis | Validation | Reserve | Fit hours | Divergences | Tree hits |
|---:|---|---:|---:|---:|---:|---:|
| 19 | No MLD, k7 | 1.540173 | 2.738907 | 19.714 | 0 | 0 |
| 22 | Shared MLD/distance interaction, k7 baseline | 1.534243 | 2.736312 | 35.744 | 0 | 8 |
| 25 | Half-shift, k7 | 1.488747 | 2.647956 | 27.153 | 0 | 2 |
| 28 | Full-shift, k7 | 1.505464 | 2.623645 | 25.657 | 1 | 0 |
| 31 | Shared MLD varying, k7 | 1.537314 | 2.737896 | 19.885 | 0 | 0 |
| 34 | Regional MLD varying, k7 | 1.408783 | 2.558899 | 33.226 | 0 | 0 |
| 129 | Regional MLD varying, Central/South RS k15; others k7 | 1.400621 | 2.506452 | 33.928 | 0 | 0 |

Case 34 improves RMSE by about 8.5% on validation and 6.6% on reserve relative to case 19. It is a promising equation, especially compared with the weak shared MLD version. Full-shift case 28 still needs its divergence addressed.

The gain is not universal across casts. Case 34 improved cast MSE for about 54% of validation and 62% of reserve casts. The three largest cast contributions represented roughly 56%/67% of the absolute MSE change. Paired cast-bootstrap intervals included zero. Do not delete those casts; examine their data quality, Region, environmental support, and profile errors.

## Offset, shape, and regional gains

Demeaned shape RMSE for case 19 → 34 changed only 1.279265 → 1.271310 on primary validation, versus about 2.2295 → 1.9695 on reserve. Some whole-profile improvement therefore comes from correcting offsets rather than resolving all shape details.

| Region/cohort | No MLD case 19 | Regional MLD case 34 | Regional k case 129 |
|---|---:|---:|---:|
| Central RS validation | 3.720 | 3.083 | 3.134 |
| Central RS reserve | 7.925 | 7.248 | 7.082 |
| South RS validation | 2.497 | 2.377 | 2.195 |
| South RS reserve | 3.721 | 3.225 | 2.958 |

Higher regional k helped South RS in both cohorts; Central RS moved in opposite directions across cohorts. This is evidence for targeted flexibility, not a blanket requirement for k15 everywhere.

Conditional versus zero-group-effect predictions also matters. For case 19 → 34, validation population RMSE changed only about 1.606 → 1.586, compared with conditional 1.540 → 1.409. Reserve population RMSE improved about 3.087 → 2.811. Known/shared date and ENSO effects account for some of the conditional advantage.

## Sensitivity to 20% cast selection

These are averages over four separate training selections, not an ensemble:

| Equation | Mean validation RMSE | Mean reserve RMSE |
|---|---:|---:|
| No MLD | 1.6132 | 2.8458 |
| Shared interaction | 1.6051 | 2.8425 |
| Shared MLD varying | 1.6099 | 2.8441 |
| Half-shift | 1.5805 | 2.7605 |
| Full-shift | 1.5941 | 2.7323 |
| Regional MLD varying | 1.5319 | 2.6856 |

Regional MLD improved reserve in all four selections and validation in three. Repeated regional cases 210 and 218 had one and four divergences, respectively. This supports robustness of the scientific direction while preventing a clean claim that every repeat is computationally valid.

## Onset, peak, and sharpness are separate outcomes

For jointly comparable detected casts, half-shift case 7 versus no-MLD case 1 at 10% reduced onset MAE by about 4.77 m in validation and 3.32 m in reserve, but increased peak-depth error by about 2.69/2.35 m. Width error improved about 1.85/1.49 m. Comparable onset counts were 60/59 and width counts 46/45.

At 20%, half-shift 25 versus 19 reduced onset error about 2.03/.76 m and width error 1.80/1.46 m, while worsening peak-depth error about 3.23/2.60 m. Onset counts were 64/59 and width counts 54/50.

Regional MLD 34 versus 19 improved validation onset about 1.62 m but worsened reserve onset about 1.05 m; width error worsened about .49/.38 m. Joint onset counts were 60/56, width counts 53/50. Whole-profile RMSE improvements therefore do not prove onset superiority.

On common comparable casts, median predicted/observed strongest positive slopes were approximately:

| Case | Validation predicted/observed | Reserve predicted/observed |
|---:|---|---|
| 19 | .1997 | .2388 |
| 25 | .2037 | .2236 |
| 34 | .2073 | .2896 |
| 129 | .2421 | .2919 |

These are ratios, not raw slopes; samples were about 70 validation/71 reserve casts. Expected model curves remain substantially smoother than the observed diagnostic curves. Noisy individual observations and expected means are different targets, but the sharpness shortfall remains a useful model-development signal.

## k, warps, priors, and sampler settings

For half-shift at 10%, k7/k11/k15 gave validation errors 1.40023/1.37893/1.38052 and reserve 2.97771/2.95924/2.97250. More k was not monotonically better.

Relative to half-shift case 25, fixed positive warps gave:
- lambda +.25, case 221: 1.47624 validation / 2.74478 reserve.
- lambda +.50, case 223: 1.57953 / 2.84601, one divergence.
- lambda −.25, case 225: 1.55558 / 2.63620.

No warp consistently improved both cohorts. This argues against adopting a universal imposed sharpening coefficient now. It does not reject learned regional shape changes.

Regional MLD delta .95, case 241, had 1.408155/2.558737 at about 28.09 h with clean diagnostics. Tight mean b priors, case 249, had 1.406492/2.552677 at about 24.107 h, also clean. Both deserve follow-through at larger sizes. Half-shift delta .99 increased tree hits to 136 without useful error gains.

## Remaining scientific concerns

- Adjacent-depth residual correlation is around .97 in leading models; random intercepts and regional sigma do not remove vertical dependence.
- Gaussian predictive draws below chl_norm = −1 occur around 6–7%, implying negative chlorophyll on the response scale.
- Broad prior predictive 90% ranges around −19 to 20 for baseline models, and broader regional effects, motivate further calibration.
- Earlier date-blocked half-shift comparisons improved roughly 2.002080 → 1.962904 and 2.127208 → 2.045935, while full shift was slightly worse at about 2.0333/2.1546. Current regional MLD date-blocked cases were pending.
- Both validation and reserve have been repeatedly inspected. They are development cohorts, not an untouched final test.
- No new 40–100% or regional tensor conclusions can be drawn from this snapshot.

## Current shortlist and what would change it

Keep regional MLD varying with k7, the Central/South RS k15 variant, no-MLD controls, half-shift, and full-shift. Evaluate completed regional tensor results as they arrive. Compare regional error, onset/width/slope, predictive intervals, and new-date behavior alongside RMSE.

The next decision should follow collection of the already queued larger-data and unseen-date results. It should not be another unstructured batch of tests or a declaration that case 129 is the final model.
