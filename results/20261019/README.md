# DCM comparison data and figures — 19 October 2026

The audited tables cover all 270 registered runs, including 207 completed
posterior fits and 18 completed prior-predictive fits. Current scores and
diagnostic summaries come from the 19 October result export. Earlier frozen
registration and parameter sources are explicitly labelled where used.

The companion workbook and HTML report summarize these tables. The CSVs retain
individual models, regions, casts, prediction targets, and diagnostic failures.
Models are not discarded because they have poor errors or sampler warnings.

## Recreate standard scientific plots in R

1. Extract this package and open the `DCM_review_20261019` directory in RStudio.
2. Install `ggplot2`, `dplyr`, and `tidyr` if they are not already available.
3. Run `source("reproduce_figures.R")` from that directory.

The annotated script creates five PDF figures from the audited CSVs. It does
not load fitted models, submit jobs, or rerun Stan. The PDFs are not included
in this package: the reproduction script was provided without local R execution.
Five standalone SVG figures are included and embedded in the HTML report.

## Comparison basis

- 204 posterior fits use the same 108-cast primary validation data hash.
- Three additional fits use a separate 117-cast date-blocked validation set.
- Sixteen follow-up fits have reserve scores on the same 108 reserve casts;
  15 share one saved data hash and the date-blocked model uses a separate hash.
- Conditional predictions include learned date/ENSO effects. Population
  predictions set those effects to zero.
- RMSE is the square root of the mean of individual-cast MSEs. Every cast has
  equal weight; regions have unequal numbers of casts.
- All 860 aggregate RMSEs and all 20,814 saved paired comparisons were checked.
- The focused 5,000-replicate bootstrap resamples Ocean/Region/date groups
  within regions. Its intervals are exploratory and do not adjust for model search.

## Important files

- `all_registered_runs.csv`: all cases, statuses, settings, and available metrics.
- `completed_posterior.csv`: all 207 completed posterior fits.
- `focused_date_bootstrap_contrasts.csv`: direct scientific comparisons.
- `cast_influence.csv`: contributions of individual casts to squared error.
- `training_composition_date_contributions.csv`: dates driving the effect of
  changing training casts.
- `paired_influence_sensitivity.csv`: the same reference-defined casts omitted
  from both members of a pair for influence assessment only.
- `predictor_support.csv`: regional predictor ranges and unseen-date coverage.
- `selected_observed_predicted_profiles.csv`: saved representative examples,
  not observation-level predictions for the full dataset.
- `selected_population_reference_profiles.csv`: expected profiles at saved
  regional predictor settings, with date/ENSO effects set to zero.
- `audit.json`, `rmse_audit.csv`, and `deliverable_validation.json`: checks.

No observations were clipped, transformed, smoothed, binned, or removed in
this review. Influence diagnostics do not replace the full-data scores.

Methods references: Stan diagnostic guidance
(https://mc-stan.org/learn-stan/diagnostics-warnings.html), brms expected
prediction documentation (https://paulbuerkner.com/brms/reference/posterior_epred.brmsfit.html),
and Roberts et al. (2017), Ecography, https://doi.org/10.1111/ecog.02881.
