# Next logical steps

**Updated 15 November 2026.** The full menu of untried ideas is in [TRIED_AND_UNTRIED.md](TRIED_AND_UNTRIED.md). The ranked list below is the recommended order. The 25 October list follows it for history; its items 1–2 (finish the existing campaign and rebuild the all-model comparison) remain open for shape150.

## Recommended order (15 November 2026)

1. **Let round3 case 11 finish, or cancel it** (`scancel 55404452`). It is the Gaussian pooled regional MLD fit on seen dates. Its unseen-date twin (case 12) already shows that the model does not earn its cost, so the result cannot change a decision.
2. **Add standard errors before more ranking.** Compute paired per-date elpd differences and bootstrap them over dates (or cruises), from each run's `checks/validation_conditional_date_elpd.csv`. No refits are needed. Treat |Δelpd| < 4 per date, or an interval spanning zero, as a tie and prefer the simpler model.
3. **Fix unseen-date calibration at depth,** starting from the lognormal half-shift on 5 m bins:
   - (a) a date effect that varies with depth, for example a date-level slope on depth relative to MLD;
   - (b) sigma as a function of depth relative to MLD, small inside the mixed layer;
   - (c) AR(1) within cast.

   Score them on the date-blocked split first.
4. **Make the peak profile converge on unseen dates.** Options: tighter priors on `zm`/`lwid`, initial values from next150 case 13, 5 m bins, or fitting only DCM-present casts (max chl_norm > 1). It is the most interpretable model and the best on seen dates.
5. **Check the Equatorial deep values near chl_norm = −1** and the six boundary-maximum casts in the raw data before modelling around them.
6. **Collect the shape150 results after 27 October and diagnose the 48 fit failures.** This is reporting only. The Gaussian 40–100% results still answer the sample-size question.
7. **Then scale up:** fit the chosen lognormal setup at 40–100% (on 5 m bins), and run date-grouped K-fold on the 2–3 finalists.
8. **Obtain the extra data** for the tests that need it: T/S profiles for a pycnocline coordinate, cast time for quenching, raw Chl for the ratio assumption, salinity for the Red Sea GAIW. Also define a fresh, untouched grouped test set before any final claim.

## Previous plan (25 October 2026)

### 1. Finish the already authorized campaign

Read the latest pointer, queue, accounting, dispatcher ledger, and timestamp. Collect completed fit/check outputs. Repair reporting without refitting. Do not assume the October 25 incomplete cases are still incomplete today.

Prioritize the outstanding regional tensors, 40/50/75/100% fits, unseen-date tests, higher regional k, and longer-chain controls. These directly answer the open questions and are already in the design.

### 2. Rebuild the all-model comparison

Use the full registry included here and obtain the cast-level scoring tables from the named source archive. Include every registered model and its state. Keep prior-only, reused, failed, and invalid fits visible. Separate posterior ranking, prior plausibility, computational failure, and reconstruction.

Use exact shared cohorts and prediction modes. Report regional and individual-cast changes, influence, common detected-onset counts, transition width, slope ratios, peak errors, intervals, residual correlation, and cost.

### 3. Confirm the strongest current candidates

Compare no-MLD, half-shift, full-shift, regional MLD k7, and Central/South RS k15/others k7. Include regional tensors as they complete. Confirm the promising delta .95 and tighter-b prior variants at larger sizes and across cast selections.

Require acceptable diagnostics before interpreting small performance differences. Resolve divergences through geometry/prior/parameterization analysis rather than automatically maximizing every sampler control.

### 4. Separate the scientific targets

A model that improves average chlorophyll offsets may not improve onset or peak location. A translated curve may align onset while missing steepness. A mean curve may be appropriately smoother than an individual noisy profile, but systematic underprediction of slopes needs examination.

Use MLD slices within observed support and compare learned regional deformations. Test the user's hypothesis that deeper MLD can sharpen the rise without imposing that relation universally unless evidence supports it.

### 5. Establish a final evaluation design

Primary and reserve sets are development data now. Define a fresh grouped evaluation before final model selection claims: ideally unseen dates/cruises or an external dataset, consistent with the intended deployment target. Confirm raw data lineage and cast/date independence.

If fitting 100% for the final descriptive model, report reconstruction and obtain independent performance evidence elsewhere. Do not rename reconstruction as validation.

### 6. Resolve upstream methods and release readiness

Document the MLD derivation/reference depth/threshold, chlorophyll calibration/units, distance computation, cruise/station provenance, original normalization differences, and ENSO labeling rationale. Establish data/code licensing and attribution with the owner.

Record an explicit environment specification for future fits without rewriting historical versions. Preserve original models, data, source hashes, selections, and logs. Update the decision record after each new evidence snapshot.
