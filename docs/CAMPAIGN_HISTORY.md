# Campaign history and registration counts

Dates identify actual frozen campaign directories, not merely when a plan was discussed. Status counts come from the October 25 archive. They may have changed on Ibex.

| Campaign | Registered | Complete | Incomplete | Main purpose |
|---|---:|---:|---:|---|
| baseline150_20261015T110518Z_1736107 | 12 | 12 | 0 | Original equation, custom/default priors, 10/20/50%, long runs |
| equation150_20261015T104335Z_1647064 | 18 | 18 | 0 | Additive, shared MLD varying, MLD-only; six sampler variants each |
| extensive150_20261015T160120Z_2161146 | 102 | 102 | 0 | Broad equation, likelihood, priors, k, sampler and compute screening |
| fast150_20261014T135159Z_1121880 | 24 | 2 | 22 | Early short sampler/prior/subset/thread screen |
| followup150_20261018T150024Z_3898749 | 34 | 34 | 0 | Controlled follow-up, uncertainty/dependence/data-domain questions |
| k150_20261014T182008Z_1910449 | 12 | 12 | 0 | Lower k and larger tree limits |
| onset150_20261020T121808Z_3759601 | 34 | 34 | 0 | Shift/ratio/hybrid onset models and reused controls |
| pilot | 1 | 0 | 1 | Initial pilot registration |
| shape150_20261022T102851Z_3426481 | 282 | 79 | 203 | Learned MLD shape, larger samples, regional k and repeats |
| targeted150_20261016T104550Z_296126 | 63 | 63 | 0 | Controlled equations, variance, original-fit bridge and common rescoring |
| threading_20261014T094438Z_3994764 | 4 | 0 | 4 | Threading registrations |
| **Total** | **586** | **356** | **230** | All states retained |

Complete means completed checking/recorded output under the collector's definitions, not passing every diagnostic. The 356 complete registrations comprise 322 posterior registrations (including six reused controls) and 34 prior-only runs. They are not 356 independent newly fitted posterior models.

The [complete case catalogue](CASE_CATALOGUE.md) and [581-row matrix](../provenance/all_cases.tsv) preserve all recovered designs. The [586-entry registry](../results/20261025/all_586_registrations.csv) retains all states.

## Before the Ibex campaigns

The original final notebook fitted regional depth splines through 250 m. Its saved model contains 255,369 observations, no MLD, four chains, 1000 warmup/4000 retained per chain, delta .97, cap 14. Original notebook alternatives include smoothing/binning and varying normalization windows; do not infer execution history from their order.

Earlier discussions also considered global/regional smooth decompositions, linear Region-depth terms, robust likelihoods, and INLA/RW2 alternatives. Those proposals lack equivalent verified matrices and are not counted as completed experiments. October 13 discussions considered a regional depth-distance tensor plus depth-varying MLD, separate regional depth-MLD surfaces, and three-way interactions.

A subsequent complex MLD/distance/depth tensor motivated numerical and computational screening. The main Ibex work standardized unsmoothed preparation and focused on 10–150 m.

## Fast and k campaigns

Fast cases 1–12 varied three prior families, delta .90/.95, and tree cap 10/12 on 10% data. Cases 13–18 explored 5% and 20% with those prior families. Case 19 was short (iter 500/warmup 250), 20 four chains, 21 two threads, 22 four threads, and 23/24 prior-only checks.

The pasted intermediate tree diagnostic excluded cases 23/24 deliberately. It reported warmup only at that time. Case 24's fast prior check was complete but had R-hat around 1.047, and was not a posterior accuracy benchmark.

The k campaign tested lower spline dimensions and higher tree limits. Some completed runs were severely unconverged; retain them but do not use their error values to select k.

## Equation and baseline campaigns

Equation A added MLD as a linear offset, B used a numeric MLD-by-depth smooth plus regional depth smooths, and C omitted regional depth smooths. Each had delta/tree/warmup/longer-chain variants. The posted eighteen-row output recorded all as complete; short-run R-hat and tree saturation remained concerns.

Baseline tests restored the original equation and long settings with custom versus default priors at 10%, 20%, and 50%, plus controlled setting changes. They are a common scientific reference, not proof that a similarly named original Windows fit used those subsets.

## Extensive and targeted campaigns

Extensive has 96 posterior and six prior cases. It broadened mean equations, distributional variance, likelihood, dependence, k, priors, run length, and resources. Review all results, including unsuccessful tensors and invalid short runs.

Targeted has 57 posterior and six prior cases. It froze common cohorts, required at least 10% training, compared standardized equations, variance/priors/compute changes, and audited old fitted objects for evaluation without refitting. The original-model review included direct saved-fit diagnostic extraction and descriptive profile overlays.

On October 16, interrupted old-fit evaluation scheduling was repaired: all 63 original fit/check registrations were preserved, 67 previous evaluations plus 50 additional evaluations made 117. Collector 53867752 waited for 164 then-active recorded jobs. No model fits were submitted or cancelled by that repair.

## Follow-up and later rescoring

Follow-up has 30 posterior and four prior cases. Its prepared 250 m source became the frozen source for later onset/shape campaigns.

On October 18, rescoring added predictions for newly completed k cases 5/9 and extensive cases 41/42, jobs 53988446–53988449. Collector 53988451 waited for eight active recorded jobs. Those were prediction jobs reusing fits, not four new models.

## Onset campaign and repair

The plan was delivered October 19; the actual campaign directory dates to October 20. It contains 24 new posterior fits, four prior checks, and six reused-control evaluations. Shift, half-shift, ratio, hybrid, and shifted-distance structures were compared at primary 10/20% and date-blocked 10%, with seed/sampler/basis and prior controls.

The first exported table marked most checks failed with “Insufficient vertical support,” although new fits often had R-hat about 1.002–1.005 and large ESS. A singleton validation/reserve support case caused the scoring failure. The supplied support repair reran evaluations only. Preserve both the original and corrected methods.

## Shape campaign

The 282-case matrix includes 270 posterior fits and 12 prior checks, eleven full-data fits, 10/20/40/50/75/100% samples, three extra 20% cast selections, two date-blocked selections, regional k allocation, MLD tensor basis tests, fixed warps, prior sensitivity, and sampler controls.

Preparation 54250664 completed in 2m48s. Initial dispatcher 54250665 ran the scheduler. The early October 22 console showed 45 fitting, three checking, and 234 not submitted, with 732 requested CPUs. Later archived status at 2026-10-25T11:14:41.409930Z showed 79 complete, 48 fitting, 155 not submitted, and 768 requested CPUs.

The early all-campaign count 583 was caused by nonnumeric threading IDs being coerced to NA and collapsed in the registry. Correct string-preserving keys restore 586. This is a reporting repair, not three new fits.

## What is next in the historical sequence

Collect the existing shape campaign and resolve any check/scheduler failures before adding another broad grid. The current evidence cannot tell us which queued fits completed after October 25. Update this table from a new timestamped export, preserving the previous snapshot.
