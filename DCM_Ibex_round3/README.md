# round3: Gaussian versions and depth-bin ("histogram") models (Nov 2026)

[06_create_round3.sh](06_create_round3.sh) installed this package on 9 November 2026.

Gage asked for two things:
- **(A)** Gaussian versions of the models he had already tried. The shape150 campaign was already all Gaussian, so the non-Gaussian models in question are the lognormal next150 cases. Round 3 refits each of them with a Gaussian likelihood on chl_norm and keeps every other setting the same.
- **(B)** New Bayesian models that treat depth as separate bins ("a histogram rather than a density plot"). These keep the same covariates and random effects.

| Item | Value |
|---|---|
| Experiment | `/ibex/user/meleg/DCM_MODEL/experiments/round3_20261109T161506Z_3488274` |
| Pointer | `logs/latest_round3_experiment.txt` |
| Token / job prefix | `4cb3655563` / `dcm-s150-4cb3655563-<stage>` |
| Controller | `controller_001`, job 55404217; fits are jobs 55404439–55404457 |
| Status, 15 Nov 2026 | 15 of 16 cases checked. Case 11 (Gaussian pooled regional MLD, primary) was still fitting at 5.8 days, against a 7-day fit limit |
| Results | [results/20261115/round3_bins1m_2026-11-12.md](../results/20261115/round3_bins1m_2026-11-12.md) |

## Settings

Every case uses p20 training data, both splits (primary, and date-blocked with seed 20261022), 4 × (1,000 + 2,000) draws, `adapt_delta = 0.95` and treedepth 14. Priors are as in next150. The 5 m bin models use 32 GB; the DCM-depth models use 4 CPUs and 8 GB.

## Cases

| Cases | Model | Data fitted | Notes |
|---|---|---|---|
| 1–2 | `bin_means`: Gaussian mean enhancement per 5 m bin | 5 m bin medians (28 bins, 6,048 rows) | `(1 + mld_z ‖ depth_bin) + (1 + mld_z ‖ Region:depth_bin)` + distance_z + mld_z + date/ENSO; `sigma ~ 0 + Region + depth_bin` (fixed); scored on raw rows |
| 3–4 | `bin_class`: cumulative-logit enhancement class per bin | 5 m bin medians | Classes <+10%, +10–50%, +50–100%, +100–200%, ≥+200%; same bin terms; gives P(bin reaches each class) |
| 5–6 | `dcm_depth`: cumulative logit over which 10 m bin holds each cast's DCM | One row per cast (216) | `mld_z + (0 + mld_z \| Region)` + distance_z + date/ENSO; `thres(13)` |
| 7–8 | Gaussian version of next150 3–4 | Raw | Half-shift + sigma by depth |
| 9–10 | Gaussian version of next150 7–8 | Raw | Pooled GI shape + sigma by depth |
| 11–12 | Gaussian version of next150 9–10 | Raw | Pooled regional MLD + sigma by depth |
| 13–14 | Gaussian version of next150 11–12 | Raw | Distance by depth + sigma by depth |
| 15–16 | Gaussian version of next150 13–14 | Raw | Peak profile: `chl_norm = background + peak − 1`, modelled on its own scale |

Population predictions in the bin models keep the bin terms and drop only the date/ENSO effects (`population_re()`). The bin terms are group-level effects, but they are the model's depth structure. The `bin_means` sigma uses fixed bin effects because `re_formula` would drop group-level sigma terms.

The ordinal models are scored with:
- held-out log predictive density, per row and jointly per date;
- the ranked probability score (RPS);
- Brier scores for each "at least class k";
- modal accuracy and median class MAE;
- calibration by Region × bin.

The DCM-depth model also reports DCM depth MAE in metres and 80% coverage. Its plots are `bin_probability_median_mld.png`, `bin_probability_mld_slices.png` and `dcm_depth_distribution.png` in each run's `checks/` folder.

Class counts on real data:
- Enhancement classes (case 3): 2,545 / 1,262 / 661 / 693 / 887.
- DCM 10 m classes (case 5): 9 8 12 8 19 18 7 20 31 21 18 17 12 16.

## Results (15 Nov 2026)

| Model | Seen elpd | Unseen elpd | Unseen RMSE | Unseen cov90 | Fit h |
|---|---|---|---|---|---|
| Gaussian + sigma by depth (7/8) | −13,550 | −16,965 | 1.975 | 0.880 | 17.5 / 30.5 |
| Gaussian GI (9/10) | −13,484 | −17,007 | 1.975 | 0.880 | 93.3 / 79.8 |
| Gaussian pooled regional MLD (11/12) | still fitting | −17,005 | 1.755 | 0.905 | >139 / 92.3 |
| Gaussian distance by depth (13/14) | −13,511 | −16,965 | 1.976 | 0.872 | 23.3 / 31.2 |
| Gaussian peak profile (15/16) | −15,533 | −17,234 | 2.016 | 0.849 | 14.0 / 8.9; R-hat 1.32 / 1.10 (not converged) |
| 5 m bin means, Gaussian (1/2) | −15,036 | −18,992 | 2.154 | 0.853 | 0.5 / 0.5 |
| 5 m enhancement class (3/4) | lpd −0.803, RPS 0.072 | lpd −0.839, RPS 0.077 | – | – | 0.4 |
| 10 m DCM depth (5/6) | MAE 13.9 m | MAE 10.9 m | – | – | ~0 |

Reference: Gaussian half-shift is −17,134 / −19,107, and lognormal half-shift is −11,554 / −14,573.

Making sigma vary with depth recovers most of the Gaussian deficit, but every Gaussian model still trails its lognormal counterpart. GI pooling and distance by depth add cost and no gain.
