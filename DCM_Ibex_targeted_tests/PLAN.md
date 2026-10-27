# Questions and controlled comparisons

The completed results do not identify one optimal equation yet. Many mean-equation changes gave small overall error differences, while regional calibration, cast structure and training-set sensitivity raised larger questions. This follow-up separates those explanations before expanding large tensor models further.

## Evidence behind the choices

On the earlier common 108-cast holdout, the simple additive model at depth k=7 had conditional RMSE about 1.465. Its 90% interval coverage was only about 63% in Central RS, 53% in North RS and 65% in South RS, despite overall coverage around 90%. Several other regions had nearly complete coverage. That suggests testing regional residual scales instead of making every region share one residual SD. These are exploratory findings from a holdout already used for model development.

For the same additive equation, increasing training from 10% to 20% and 50% changed RMSE to about 1.563 and 1.602, but improved peak-depth error and individual-cast errors for many casts. Residual SD and predictive widths increased. This is not evidence that more data is intrinsically worse: it raises questions about heterogeneity, extreme profiles, coverage and learned random effects. The new regional-variance learning curves test these directly on the same holdout.

Depth-adjacent residual correlation was near 0.97 in representative summaries. Cast intercepts, cast curves and an AR residual model test different explanations. They are separate comparisons: mean profile variation and residual dependence are not interchangeable.

One earlier Student-t model estimated degrees of freedom near 1, with a very small scale and poor held-out coverage. The follow-up fixes nu=4 rather than allowing nearly Cauchy tails to dominate that comparison.

The longer four-chain equation runs had much better R-hat and ESS than the short two-chain runs. Many completed four-chain models did not hit the tree limit. This motivates controlled warmup, retained-draw and delta tests rather than automatically increasing every limit.

## Main controls

All new equations retain `Region`, `(1|Ocean:Region:date)` and `(1|Ocean:Region:enso)`.

Control **case 2**:

```r
bf(
  chl_norm ~ 1 + Region + distance_z + mld_z +
    s(depth, by = Region, bs = "cr", k = 7) +
    (1 | Ocean:Region:date) + (1 | Ocean:Region:enso),
  sigma ~ 0 + Region,
  family = gaussian()
)
```

Control **case 8** replaces `distance_z + mld_z` with:

```r
s(depth, by = mld_z, bs = "cr", k = 5) +
s(depth, by = distance_z, bs = "cr", k = 5)
```

This permits linear effects of MLD and distance whose size changes smoothly with depth. It is not a general nonlinear response to MLD or distance. The numeric-by smooths contain their corresponding constant-in-depth effect; an additional duplicate main term is omitted. Both controls have four chains, 1,000 warmup and 1,000 retained iterations per chain, adapt_delta=0.95, max_treedepth=12, two threads per chain, seed 1234, 10% training and the same 10% whole-cast holdout.

## New case matrix

| Cases | Question | Direct comparison |
| --- | --- | --- |
| 1-12 | Are regional SD, depth-varying regional SD, cast intercepts/curves or fixed-nu Student-t errors useful? | Six configurations for each control equation. Cases 1/7 share one sigma; 2/8 allow regional sigma; other additions compare to 2/8. |
| 13-18 | Which predictors need to change the vertical shape? | Original equation, MLD-only varying effect, distance-only varying effect, region-specific MLD varying effects, region-specific MLD and distance varying effects, and a shared MLD-by-distance varying effect. |
| 19-22 | Does more regional depth flexibility matter after changing the noise model? | Regional depth k=7 versus 10/15; numeric-by MLD/distance k stays 5. |
| 23-28 | Do coefficient or group/smooth priors change predictions and decomposition? | Change only coefficient scale, only group/smooth SD prior family, or only mean-model priors to brms defaults. |
| 29-46 | What sampler settings buy accuracy or speed? | For each control: delta .90/.95/.97; tree 12/14; warmup 500/1,000/2,000; retained draws 1,000/2,000; 1/2/4 threads; repeat seed. Each contrast changes one setting from its control. |
| 47-50 | Does a regional noise model improve the learning curve? | The same equations at 10%, 20% and 50% training, common validation and fixed predictor scales. |
| 51-52 | Does the model generalize to unseen nested date groups? | Two additional 10%-training fits with a >=10% date-blocked holdout. Their cohorts differ from primary tests. |
| 53-58 | Are the competing proper priors plausible on the response scale? | Prior-only simulations for the two controls and three proper prior specifications. No accuracy ranking. |
| 59-60 | Does residual dependence matter after modelling the mean? | Gaussian AR(1) within cast, constant sigma, compared with cases 62/63. |
| 61 | Does fitting down to 250 m change the 10-150 m curves? | Original raw-predictor equation, original user priors and historical sampler settings; matched to the existing 150 m baseline case 1 on the same 108 training casts. |
| 62-63 | Provide a fair computation reference for the AR models. | Same mean equation and scalar sigma as cases 59/60, same one CPU per chain and no compiled within-chain threading, with the AR term omitted. |

For a further basis-budget comparison, new case 7 uses `(k_depth,k_MLD,k_distance)=(7,5,5)` and existing extensive150 case 49 uses `(7,7,7)` with the same standardized mean equation, custom priors, constant sigma and sampler settings. The common rescoring verifies training identity before making this comparison. Cases 19-22 then change depth k separately with a fixed effect basis. Existing full tensor cases remain relevant; unfinished cases are not declared inferior because they take longer.

The cast-curve model uses four fixed natural-spline basis columns, knots 45/80/115 m and boundaries 10/150 m, plus a cast intercept. Their independent random coefficients permit flexible between-cast deviations and predictable new-cast simulation without estimating a separate unconstrained fixed curve for every cast. It is compared with the corresponding no-cast model; compare its cast-intercept counterpart as well to isolate the additional curve slopes.

The AR comparison uses `ar(time=depth, gr=cast_key, p=1, cov=TRUE)` and no cast random effect. It describes adjacent observed rows after sorting by depth. It does not estimate a continuous correlation length in metres; irregular spacing is recorded. Validation predictive simulation uses the covariance formulation and does not condition on held-out chlorophyll responses to forecast a later row.

## Priors

On standardized covariates, the custom mean priors are Intercept student_t(3,0.73,5.448), coefficients normal(0,2), group/smooth SDs student_t(3,0,1.5), and constant residual sigma student_t(3,0,2.5). These match the stated historical values for the raw equation in case 61. Standardizing a predictor changes the meaning of a fixed coefficient prior, so raw-versus-standardized comparisons can reflect prior scaling as well as numerical conditioning.

For regional sigma, each log-sigma coefficient has normal(0,1). Depth-varying log-sigma smooth SDs have half-normal(0,0.5). Cast coefficient SDs have half-normal(0,1). The AR coefficient has normal(0,0.3) within its model-defined bounds.

`b_tight` changes only mean coefficients to normal(0,1). `scales_halfnormal` changes only mean group/smooth SDs to half-normal(0,1.5), keeping mean coefficients and sigma priors fixed. `mu_default` uses default brms mean-model priors while retaining the same explicit log-sigma priors. It is not labelled an all-default-prior model. Fully default flat coefficient priors do not define the proper prior-predictive simulations requested here.

## Compare every completed fit

Read convergence diagnostics first, but keep flagged fits in the tables and pairwise comparisons. Then inspect cast-balanced error and bias, demeaned profile shape, peak-depth and peak-height errors, integrated positive enhancement, regional performance, predictive coverage, interval width and the proper 90% interval score. The interval score penalizes both missed observations and unnecessarily wide intervals. Raw observed maxima are noisy; they do not alone determine a preferred equation.

Use all-pairs matched-cast differences to see how improvements are distributed. Inspect seen versus unseen dates/ENSO groups and depth/MLD/distance strata. DCM occurrence here means the observed or expected maximum is deeper than 15 m within 10-150 m; maxima at the 150 m boundary are flagged. It does not describe maxima deeper than the analysis window.

Use R-hat, bulk/tail ESS, MCSE/SD, divergences, E-BFMI, tree-hit fraction, ESS per sampling second and ESS per requested core-hour to compare sampler settings. A larger dataset does not by itself justify fewer posterior draws. Four chains are held fixed so shorter warmup and longer retained sampling have interpretable comparisons. Tree efficiency and convergence are reported separately.

The primary holdout has already guided these targeted choices. It is a development validation set, not an untouched final test. After selecting a small set of defensible models, freeze the specification and use a new whole-cast/date holdout from the remaining data or an outer grouped cross-validation design for the final performance claim. The new date-blocked panel is an additional generalization stress test, not a substitute for that final evaluation.

References: [Stan diagnostics](https://mc-stan.org/learn-stan/diagnostics-warnings.html), [brms distributional models](https://paulbuerkner.com/brms/articles/brms_distreg.html), [brms AR terms](https://paulbuerkner.com/brms/reference/ar.html), [brms threading](https://paulbuerkner.com/brms/articles/brms_threading.html).
