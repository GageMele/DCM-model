# Experiment plan

## Scientific target and interpretation

Predict chlorophyll enhancement (`chl_norm`) over observed depths 10–150 m.
Region, distance from shore and mixed-layer depth may change profile shape, level,
or both. Every candidate retains exactly:

```r
(1 | Ocean:Region:date) + (1 | Ocean:Region:enso)
```

These are intercepts for the nested date and ENSO grouping combinations. They
do not add separate basin-level intercepts or temporal correlation. We keep the
requested grouping structure fixed so the other comparisons are interpretable.

The existing results justify testing simpler structures before increasing every
sampler limit. In the latest completed equation campaign, the additive four-chain
case had maximum R-hat 1.0066 and minimum bulk ESS about 707, with zero sampling
divergences/tree hits. The four-chain MLD-varying-plus-Region model had maximum
R-hat 1.0064, minimum bulk ESS about 751, zero divergences, and 3 tree hits in 4,000
sampling draws. These are useful candidates; fit quality still needs assessment.

## Data and fair comparisons

The frozen 10%, 20% and 50% subsets contain 108/216/540 complete casts and
15,228/30,455/76,139 observations, respectively. All eight Regions are represented.
The subsets are nested and copied unchanged from the baseline campaign. There is
no added smoothing, binning, row subsampling, or change to chlorophyll normalization.
All eligible observed depths within each selected cast are retained.

Validation uses **108 additional whole casts (15,226 observations)**, selected
reproducibly by Region from outside the 50% training subset. The same validation
casts are used for every new posterior case, including the larger-data fits.
The validation observations are never used for fitting, predictor centering,
scaling or prior construction. No new empirical prior is estimated from full-data
chlorophyll outcomes. Historical custom prior constants remain the user-supplied
constants; their earlier provenance cannot be made independent retroactively.

The 10% training casts are also the fixed in-sample assessment set for every model.
Thus differences in those reconstruction scores are not caused by evaluating
different casts. Training error is reported separately from validation error.

Validation is for model selection. After comparing many candidates on it, its best
score is not an unbiased final-test estimate. Some validation casts may share
date/ENSO groups with training. We report seen/new group strata, but this is not
equivalent to leaving entire cruises or dates out. A later date/cruise-blocked
evaluation is needed before claiming transfer to independent sampling occasions.

## Equation and basis comparisons: cases 1–49

Every row below also has an intercept, Region fixed effects and the two random
intercepts. `L` means the original `log10_distance`; `M` means original
`MLD_temp_MIX`; `zL` and `zM` are explicit training-only standardized versions.
`R(depth)` means `s(depth, by=Region, bs='cr', k=k_depth)`.

| Cases | Equation terms beyond the common terms | Question |
|---|---|---|
| 1–4 | `L + R(depth)` | Original equation reference |
| 5–8 | `L + M + R(depth)` | Do additive level shifts suffice? |
| 9–12 | `L + s(depth, by=M) + R(depth)` | Does MLD alter vertical shape? |
| 13–16 | `M + s(depth, by=L) + R(depth)` | Does distance alter vertical shape? |
| 17–20 | `s(depth, by=M) + s(depth, by=L) + R(depth)` | Both shape effects with Region curves |
| 21–24 | `s(depth) + s(depth, by=M) + s(depth, by=L)` | Shared shape relationships; Region offsets |
| 25–28 | `L + t2(depth, M, by=Region)` | Nonlinear MLD–depth surfaces within Region |
| 29–32 | `M + t2(depth, L, by=Region)` | Nonlinear distance–depth surfaces within Region |
| 33–36 | `L + t2(depth, M)` | Cheaper shared nonlinear MLD–depth surface |
| 37–40 | `M + t2(depth, L)` | Cheaper shared nonlinear distance–depth surface |
| 41–44 | `t2(depth, L, M, by=Region)` | Small three-dimensional interactions |
| 45–48 | Both varying effects + Region curves | Unequal k for the three smooths |
| 49 | `s(depth, by=zM) + s(depth, by=zL) + R(depth)` | Centered/scaled varying effects |

Cases 1–40 each use depth k = 5, 7, 10, 15. All one-dimensional splines use cubic
regression bases. For varying-coefficient models, each numeric-by spline initially
has the same depth k as the Region spline. For two-dimensional tensors, the other
marginal k is 3. No claim is made that k=15 is necessary simply because it was used
at 250 m; k controls capacity, while smoothing priors also control realized wiggle.

The three-dimensional k triples, in order **depth, distance, MLD**, are
(5,3,3), (7,3,3), (5,4,3), (5,3,4). The mixed one-dimensional k triples, in order
**Region-depth, MLD-by-depth, distance-by-depth**, are (10,5,5), (15,5,5),
(10,10,5), (10,5,10). The order is explicit in `cases.tsv` and generated formulas.

A numeric-by spline is `M * g(depth)`, not an arbitrary nonlinear surface in MLD.
Its effect is linear in MLD at any fixed depth, with the coefficient changing
smoothly with depth. It already includes a constant-in-depth coefficient for MLD;
the formula therefore avoids adding a duplicate standalone M term. The same
applies to the distance-by smooth. The tensor cases test richer nonlinear
relationships. brms supports `s()` and `t2()`; this package does not assume `ti()`
is a supported brms formula term.

Example candidate, case 49:

```r
bf(chl_norm ~ 1 + Region +
     s(depth, by = Region, bs = "cr", k = 7) +
     s(depth, by = mld_z, bs = "cr", k = 7) +
     s(depth, by = distance_z, bs = "cr", k = 7) +
     (1 | Ocean:Region:date) + (1 | Ocean:Region:enso),
   family = gaussian())
```

Actual-data construction gives this candidate 50 penalized spline columns,
compared with 104 for the original Region-only k=15 smooth. These are spline
columns, not total parameter counts; fixed/null-space and group parameters also
contribute. Small 3D cases still have 296–440 penalized spline columns. The
`basis_counts.csv` file gives every distinct tested construction.

Scaling uses one observation per 10%-training cast to calculate MLD/distance
means and SDs. Those same constants are reused at 20%/50%. `chl_norm` and depth
are unchanged. Giving standardized coefficients the same numeric priors changes
the implied prior in raw units; case 49 is a parameterization-and-prior comparison,
not a pure computational invariance test.

## Sampling and cost comparisons: cases 50–73 and 80–83

Two anchors cover different complexity: additive MLD/distance + Region-depth
at k=7 (case 6), and both standardized varying effects + Region-depth at k=7
(case 49). The reference settings are four chains, 1,000 warmup + 1,000 retained
draws per chain (`iter=2000`), adapt_delta=.95, max_treedepth=12, two threads per
chain and seed 1234.

| Factor | Tested values including the reference | What stays fixed |
|---|---|---|
| adapt_delta | .90, .95, .97, .99 | Tree cap 12, equation and draw budget |
| max_treedepth | 10, 12, 14 | adapt_delta .95 and draw budget |
| Warmup per chain | 250, 500, 1,000, 2,000 | 1,000 retained draws per chain |
| Retained draws per chain | 500, 1,000, 2,000 | Warmup 1,000 |
| Chains, equal total retained draws | 2 × 2,000 vs 4 × 1,000 | 4,000 total retained draws; warmup per chain 1,000 |
| Chains, equal per-chain length | 2 × 1,000 vs 4 × 1,000 | Same warmup and per-chain iterations |
| Threads per chain | 1, 2, 4 | Four chains and the same statistical settings |

These are matched one-factor comparisons except the explicitly labelled
equal-total-draw chain comparison, which changes chain count and per-chain length
together. Warmup cost across chains necessarily differs in that comparison.
This first wave does not estimate every interaction among all sampler settings.
If a particular equation needs both delta=.99 and depth=14, that is a targeted
follow-up rather than treating a higher cap as a universal improvement.

Larger data sets do not provide extra MCMC draws. They can make likelihood
evaluation slower and the posterior more concentrated. We choose run length from
R-hat, bulk/tail ESS, Monte Carlo error and reproducibility of profile predictions.
Four chains with 1,000 retained draws each is a starting comparison, not an assumed
sufficient final run. Warmup divergences/tree hits are reported separately from
sampling diagnostics. Increasing adapt_delta may increase leapfrog work and tree
hits; increasing tree depth can make an iteration much more expensive.

## Priors, data sizes and repeatability

Cases 74–79 add three prior choices to each anchor's historical custom prior:

| Component | Historical custom | Tighter | Half-normal SD sensitivity |
|---|---|---|---|
| Intercept | student_t(3, .73, 5.448) | same | same |
| b | normal(0, 2) | normal(0, 1) | normal(0, 2) |
| sd and sds | positive student_t(3, 0, 1.5) | positive normal(0, .75) | positive normal(0, 1.5) |
| sigma | positive student_t(3, 0, 2.5) | same | positive normal(0, 2.5) |

The third alternative is `prior=NULL`, allowing the installed brms defaults.
Resolved priors are recorded, rather than describing defaults as if they were
constant across versions/data. Cases 97–102 check the three proper custom prior
sets on both anchors. Default flat coefficient priors do not define a proper
prior-predictive distribution, so no default-prior-only run is submitted.
Prior scale values are sensitivity experiments, not established biological truth.

Cases 84–91 repeat additive, MLD-varying, both-scaled-varying, and Region-specific
MLD tensor models at 20% and 50%. These use the same validation casts as the 10%
fits and keep k and sampling settings fixed. Your previously queued historical
baseline runs retain the original 5,000/1,000 settings; they are not resubmitted
by this package. The new original-equation cases supply references under the new
matched budget and prediction protocol.

Cases 92–93 repeat both-scaled-varying with seeds 20261016 and 20261017; case 96
repeats additive with seed 20261016. Seeds are kept fixed elsewhere, but a shared
seed does not make trajectories or hardware timing noise identical.

Cases 94–95 compare a Student-t likelihood with the Gaussian likelihood on the
two anchors. This tests sensitivity to heavy residual tails/outliers while
keeping the mean equation fixed. The Student-t degrees-of-freedom prior is the
resolved brms default. There is no transformation of the outcome.

## How to select a model

1. Inspect completed sampling diagnostics. The convergence screen requires
   maximum varying-parameter R-hat <=1.01, minimum bulk/tail ESS >=400, no sampling
   divergences and per-chain E-BFMI >=0.3. Parameter mean MCSE/SD is reported.
   This is a screening threshold; important predictions may need more precision.
2. Inspect tree-depth hit rates separately, alongside leapfrog work and ESS per
   second/core-hour. A few hits do not carry the same interpretation as divergent
   transitions or poor mixing. R-hat/ESS passing does not erase a saturation issue.
3. Compare cast-balanced validation RMSE/MAE, bias and **profile-shape RMSE**.
   Shape error removes each cast's observed/predicted mean before comparing
   depths, distinguishing profile reconstruction from merely predicting its level.
4. Inspect errors by Region, depth, MLD and distance, plus seen/new date/ENSO
   groups. Check predictive interval coverage and width together. Broad intervals
   should not win simply because they cover most observations.
5. Review observed/predicted profiles, raw peak-depth/peak-height errors, and
   residual lag-1 correlation. Observed maxima can be noisy; they are secondary
   descriptive checks. Cast-balanced row errors are not joint profile likelihoods.
6. Compare models that give similar predictive performance by cost and stability.
   Check seed repeats and whether higher k materially changes profiles. Repeat
   the smallest sufficient-looking sampler setting on the selected equation.
7. Confirm finalists with date/cruise-blocked validation and adequate chains/draws
   before a final full-data fit. Estimate final posterior uncertainty only after
   checking residual dependence and the uncertainty calibration.

Conditional predictions include the fitted date/ENSO effects where available and
sample effects for unseen groups. Population predictions set group effects to
zero. Their residual predictive intervals are conditional on zero effects and do
not represent the full predictive spread for a new date. Both are labelled.
Predictions use 400 posterior draws for this screen; individual latent-parameter
diagnostics use all retained draws. The numerical prediction files allow direct
repeat-seed profile comparisons. Differences close to Monte Carlo noise should
be re-evaluated with more predictive draws before model selection.

One concern beyond spline k is correlation among depths within a cast. Date and
ENSO random intercepts do not model depth-dependent residual correlation. This
campaign records residual correlation and holds out whole casts, but all likelihoods
still assume conditionally independent observations. If correlation remains strong,
the next model revision should test a justified cast-level structure or residual
correlation. An AR(1) specification should respect actual depth spacing; it is not
added blindly to irregular observations.

## Sources

- Stan, [diagnostics and convergence guidance](https://mc-stan.org/learn-stan/diagnostics-warnings.html).
- brms, [smooth specification](https://paulbuerkner.com/brms/reference/s.html).
- mgcv, [numeric-by smooths](https://stat.ethz.ch/R-manual/R-patched/library/mgcv/html/gam.models.html)
  and [tensor construction](https://stat.ethz.ch/R-manual/R-patched/library/mgcv/html/t2.html).
- Pedersen et al. (2019), *Hierarchical generalized additive models in ecology*,
  [PeerJ 7:e6876](https://doi.org/10.7717/peerj.6876).
- Roberts et al. (2017), *Cross-validation strategies for data with temporal,
  spatial, hierarchical, or phylogenetic structure*,
  [Ecography](https://doi.org/10.1111/ecog.02881).
