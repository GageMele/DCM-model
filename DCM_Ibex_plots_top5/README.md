# Posterior plots for the top 5 setups (15 Nov 2026)

[09_plot_top5.sh](09_plot_top5.sh) is read-only. It submits:
- one Slurm array task per fitted case: job name `dcm-plot-top5`, 2 CPUs, 64 GB, 12 h;
- one dependent comparison job, `dcm-plot-compare`.

It loads each `fit.rds` and the frozen data and writes everything to `plots/top5_<UTC stamp>/`. No experiment directory is changed. The latest output folder is recorded in `logs/latest_plots_top5.txt`. The 15 Nov run was array 55624665 plus comparison job 55624666, writing to `plots/top5_20261115T155349Z`. All jobs finished within about 35 minutes.

| Tag | Experiment / case | Model |
|---|---|---|
| LN_half_seen / _unseen | next150 1 / 2 | Lognormal half-shift |
| LN_sigdepth_seen / _unseen | next150 3 / 4 | Lognormal half-shift + sigma by depth |
| LN_distvary_seen / _unseen | next150 11 / 12 | Lognormal half-shift + sigma by depth + distance by depth |
| LN_peak_seen | next150 13 | Lognormal peak profile (case 14 did not converge, so it is not plotted) |
| LN_hist1m_seen / _unseen | bins1m 7 / 8 | Lognormal 1 m histogram (bin means) |

## What `plot_case.R` makes

[plot_case.R](plot_case.R) produces these outputs for each case, separately for training casts (in-sample) and for held-out validation casts:

- `*_profiles.png`: a typical cast and a strong-DCM cast per region. Each shows observed points, the posterior mean, and 50% and 90% posterior predictive bands.
- `*_ppc_density_region.png` and `*_ppc_density_depth.png`: the observed distribution of Chl/surface against 40 posterior predictive replicates, by region and by depth band, on a log axis.
- `*_quantile_profiles.png/.csv`: the across-cast 10/50/90% quantiles per region and 5 m band, observed against posterior predictive. The bands show the 90% uncertainty of the predicted quantiles.
- `*_pit.png`: a PIT histogram, P(replicate < observed), by region. A flat histogram means the model is calibrated, a U shape means the intervals are too narrow, and a hump means they are too wide.
- `*_obs_vs_pred.png`: observed against posterior mean, on log axes.
- `parameters.png` and `parameter_summary.csv`: posteriors of population coefficients, smooth/group SDs and sigma terms, with R-hat and bulk ESS.

Predictions use 300 posterior draws. Unseen date/ENSO levels are drawn from their fitted distribution, once per date group, matching the campaign scoring. Lognormal predictions are shifted back to chl_norm.

## What `combine_plots.R` makes (`compare/`)

[combine_plots.R](combine_plots.R) puts the models side by side:

- `heldout_summary.csv`: RMSE, MAE, 50% and 90% coverage, median 90% width and PIT tail fractions.
- `pit_by_model.png` and `coverage90_by_depth.png`.
- `ppc_density_by_model.png`.
- `profiles_<split>.png`: one strong-DCM held-out cast per region, with every model's mean.
- `quantile_profiles_<split>.png`.

The figures and a reading of them are in [results/20261115](../results/20261115/README.md).

## Running it again

```bash
cd /ibex/user/meleg/DCM_MODEL && bash 09_plot_top5.sh 2>&1 | tee logs/plot_top5.txt
squeue -u $USER -o "%.18i %.20j %.9T %.10M" | grep dcm-plot        # empty = finished
cd plots && tar czf top5_png.tar.gz --exclude='*.rds' --exclude='*_rows.csv' top5_<stamp>
```

To plot other cases, edit the `printf` model list in the script. Each row has four fields: experiment, case, tag and display name.
