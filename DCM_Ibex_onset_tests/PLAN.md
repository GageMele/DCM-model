# MLD-onset tests — 19 October 2026

Question: does mixed-layer depth locate the beginning of subsurface chlorophyll
enhancement, rather than merely improving whole-profile error or peak depth?

This package registers 34 cases: **24 new posterior fits, 4 prior-predictive fits,
and 6 evaluations of existing fits**. No old model is refitted or cancelled.

## Equations

All retain Region, nested Ocean:Region:date and Ocean:Region:enso intercepts,
Gaussian observations and sigma ~ 0 + Region. Custom standardized priors are
copied from the existing follow-up. Distance is standardized log10_distance.

1. Absolute-depth varying MLD: existing s(depth,by=mld_z,k=5) alongside Region
   depth curves and a linear distance effect. This is a matched simpler MLD control.
2. MLD shift: replace depth in Region curves by depth - MLD, k=7. This tests
   translation with MLD; the enhancement can start above or below MLD.
3. Half shift: Region curves use depth - 0.5*MLD, k=7. This tests a weaker
   relationship instead of assuming a one-metre translation per metre of MLD.
4. MLD ratio: Region curves use depth / MLD, k=7. This tests stretching with MLD.
5. Hybrid: retain absolute-depth Region curves k=7 and add a shared smooth of
   depth - MLD, k=5. This allows regional absolute-depth structure and a separate
   MLD-relative contribution to coexist; check posterior dependence/diagnostics.
6. Shift plus distance: shifted Region curves k=7 plus a shared coefficient
   s(depth - MLD,by=distance_z,k=5), allowing distance to change shifted shape.

Every equation is run on the exact existing primary 108-cast training set, primary
216-cast training set, and 108-cast date-blocked training set: 18 posterior fits.
Six additional fits isolate repeat seed, robust sampler settings, and k=10 for
the shift and hybrid models. Four prior-only fits check those two new forms.
Controls reuse follow-up cases 1/5,13/15,22/24 (no-MLD and product interaction),
with raw training identity checks before scoring. Their original sampler settings
are preserved. New fits use four threads per chain; old controls use two, so
equation comparisons are statistical comparisons, not controlled speed benchmarks.

## Sampler settings

Default: four chains, 1000 warmup + 2000 retained draws per chain (iter=3000),
adapt_delta=.95, max_treedepth=14, four threads per chain (16 CPUs). Robust checks
use warmup=2000 and adapt_delta=.99, keeping 2000 retained draws. Prior-only fits
use 500 warmup + 500 retained and two threads (8 CPUs). Memory is 32 GB for 10%
fits and 64 GB for 20%. Wall limit 7 days is a scheduler limit, not a runtime estimate.
All 28 new fitting jobs together request 416 CPUs and 1088 GB if simultaneously
allocated; Slurm decides actual concurrency. Evaluation/preparation jobs are extra.

## Data and honest comparisons

No new split selection, response smoothing, binning, cast deletion or MLD imputation.
At least 10% (108 casts) is used to train every posterior model. Validation is
another >=10%, plus the existing disjoint 108-cast reserve. Primary validation may
share dates; date-blocked validation excludes training Ocean/Region/date groups.
Copied raw rows must be identical. New deterministic predictor columns change RDS
hashes, so a separate source-identity audit is retained. Reserve is now a development
cohort because we have already inspected it; it is not an untouched final test.
Nonpositive MLD stops preparation rather than dropping rows or silently changing
the ratio formula. Every actual formula, requested prior and design matrix is
constructed on Ibex before fits can start; rank deficiency fails preparation.

## New onset diagnostics

The response is (Chl-surface_chl)/surface_chl. Fixed thresholds .25,.5,1 therefore
mean 25%,50%,100% enhancement relative to the original 5-15 m surface baseline.
The primary onset threshold is .5; .25 and 1 are prespecified sensitivities.
Onset is the first depth at/after 15 m above threshold for five consecutive
one-metre grid points. For scoring only, both observed and expected profiles use
a five-metre running median. Fit inputs remain unsmoothed. Grid interpolation
stays inside observed 10-150 m support and never bridges gaps >3 m.

Threshold already reached at the upper boundary is left-censored. No threshold
crossing is labelled not reached in the observed window, never replaced by 150 m.
Gaps before first detection are flagged. Export detection agreement, false
positives/negatives, censor counts, conditional onset MAE, onset minus MLD,
and posterior onset probabilities/intervals from up to 400 expected-curve draws.
Paired onset-MAE contrasts use the same jointly detectable casts. Always inspect
false positives/negatives and the sample count alongside conditional MAE; a model
must not appear better just because it fails to predict difficult onsets.

Also export steepest-rise depth (positive centred four-metre slope >=.02 chl_norm
per metre, profile amplitude >=.5), error within +/-20 m of MLD, full-profile RMSE,
shape error, peak depth, interval calibration, regional metrics and MCMC diagnostics.
The steepest-rise location is a complementary operational metric, not a biological
definition of DCM onset. Do not force onset or peak depth to equal MLD. No physical
truncation or removal of large Red Sea errors is applied.

## Execution

Upload queue_onset_dcm.sh and run it from /ibex/user/meleg/DCM_MODEL. Repeating the
same installer resumes recorded submissions. ONSET_NEW=1 explicitly creates a new
campaign; leave it unset for recovery. Accepted old or new fit jobs are not replaced.
An ambiguous submission stops rather than risking duplicate jobs. Unique job names
permit recovery after interruption. All helper checks run before queueing; preparation
is an afterok dependency of new fits and reused evaluations. Source evaluations also
wait for active source fits/checks. The final collector includes all recorded earlier
campaigns and retains diagnostic failures and incomplete states.

Refresh: bash DCM_Ibex_onset_tests/collect_results.sh
Monitor all: bash DCM_Ibex_onset_tests/monitor_all_dcm_v7.sh --diagnostics
Export all: bash DCM_Ibex_onset_tests/export_results.sh

Scientific motivation: Mignot et al. 2014, Global Biogeochemical Cycles,
doi:10.1002/2013GB004781, describes mixed-layer and DCM structure while highlighting
light-driven DCM movements. It motivates testing the relationship, not fixing it.
Smooth syntax: https://paulbuerkner.com/brms/reference/s.html
