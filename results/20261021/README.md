# Complete DCM comparison — 21 October 2026

Open `DCM_complete_review_20261021.html` for the illustrated scientific review,
or `DCM_all_models_comparison_20261021.xlsx` for the complete comparison tables.

The registry retains all 304 cases, including bad diagnostics, prior checks
and entries without completed results. Six onset controls reuse older fits;
the 255 completed posterior evaluations represent 249 distinct model fits.
Twenty-two prior-predictive fits are separate from prediction accuracy.
The original full-data 250 m model is a separate historical reference.

Every current metric and status comes from the uploaded 21 October compact
archive. The review independently recomputed all 1,136 profile RMSE summaries
and all 360 onset summaries from individual-cast tables. The optional separate
audit of every saved model pair was stopped by automatic approval review
after an unexpected telemetry request and was not retried. The saved all-pair
file is retained as source output, without a new verification claim.

Important files:

- `all_registered_runs.csv`: settings, diagnostics, every available score.
- `onset_all_threshold_scores.csv`: onset detection, censoring and errors.
- `onset_all_cast_metrics.csv`: every scored profile, threshold and mode.
- `focused_onset_contrasts.csv`: identical-cast onset comparisons and
  exploratory date-cluster bootstrap intervals.
- `focused_profile_contrasts.csv`: matched profile comparisons.
- `focused_regional_contrasts.csv`: regional differences for each comparison.
- `all_*_cast_metrics.csv`: full training/validation/reserve cast-level errors.
- `paired_influence_sensitivity.csv`: same reference-defined omitted casts in
  both models, for influence checks only; full-data scores remain primary.
- `review_findings.csv`, `review_methods.csv`: conclusions and interpretation.
- `rmse_audit.csv`, `onset_score_audit.csv`, `audit_summary.json`: verification.
- `historical_original_model.csv`: original-model diagnostics, not a held-out
  accuracy comparison with the newer fits.

Run the annotated `reproduce_figures.R` in RStudio to generate three PDF plots
from these CSVs, using ggplot2, dplyr and tidyr. It does not load fitted models
or rerun Stan. The supplied SVGs are already embedded in the HTML report.

The Python analysis scripts are included as an audit trail. Re-running the
full numerical analysis additionally requires the source archive unpacked
with its directory structure and the earlier configuration metadata noted
in the script; the standalone R plotting script requires only this folder.
No model fitting, prediction generation, job submission or observation
filtering was performed in this review.

The source archive SHA256 is
`aa21a660aaa005ec74fd259e5459c8c0fae464210eba20195e4f699780204643`.
