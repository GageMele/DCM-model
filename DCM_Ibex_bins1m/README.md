# bins1m: depth-bin ("histogram") models at 1 m resolution (Nov 2026)

[07_create_bins1m.sh](07_create_bins1m.sh) installed this package on 9 November 2026. Gage asked for the histogram approach at 1 m rather than the 5 m used in round3. Here every measured depth from 10 to 150 m is its own bin: the factor `depth_m` has 141 levels, `m010`–`m150`. The raw 1 m rows are fitted directly, with no bin medians. The code is a generalisation of round3, and it builds identical formulas for the round3 cases.

| Item | Value |
|---|---|
| Experiment | `/ibex/user/meleg/DCM_MODEL/experiments/bins1m_20261109T175123Z_3650341` |
| Pointer | `logs/latest_bins1m_experiment.txt` |
| Token / job prefix | `dbfdd2cd23` / `dcm-s150-dbfdd2cd23-<stage>` |
| Controller | `controller_001`, job 55406551, COMPLETED after 1-04:14 |
| Status | All 8 cases checked by 12 Nov 2026 |
| Results | [results/20261115/round3_bins1m_2026-11-12.md](../results/20261115/round3_bins1m_2026-11-12.md) |

## Settings

Every case uses p20 training data, both splits, 4 × (1,000 + 2,000) draws, `adapt_delta = 0.95` and treedepth 14. The bin models use 16 CPUs and 64 GB; the DCM-depth models use 4 CPUs and 16 GB. The data are 30,455 / 30,456 raw rows from 216 casts.

| Cases | Model | Data |
|---|---|---|
| 1–2 | `bin_means`, Gaussian: `(1 + mld_z ‖ depth_m) + (1 + mld_z ‖ Region:depth_m)`, `sigma ~ 0 + Region + depth_m` | Raw rows |
| 3–4 | `bin_class`, ordinal, with the same 5 enhancement classes | Raw rows (class counts 12,561 / 6,497 / 3,397 / 3,490 / 4,510) |
| 5–6 | `dcm_depth`, ordinal over 141 1 m classes; DCM = depth of the cast's maximum raw chl_norm | 216 casts |
| 7–8 | `bin_means`, lognormal: response `ratio = chl_norm + 1`, Intercept `student_t(3, 0, 1)` | Raw rows (10 and 4 rows floored) |

Cases 7–8 were added because the lognormal likelihood won in next150.

## Caveats

- **Adjacent 1 m rows are highly correlated (about 0.97),** so 1 m bins give more overconfident intervals than 5 m bins.
- **The 1 m DCM model spreads 216 casts over 141 classes,** so its histogram is spiky, and many thresholds rest on the prior.
- **Boundary maxima:** 4 casts have their maximum at 150 m and 2 at 10 m. These maxima sit at the edge of the depth range: the DCM is deeper than 150 m, or the cast has no subsurface maximum. The same casts affect the 10 m DCM model in round3. Flag or exclude them when interpreting DCM-depth results.

## Results (12 Nov 2026)

| Model | Seen elpd | Unseen elpd | Seen / unseen RMSE | Unseen cov90 | Fit h | Diagnostics |
|---|---|---|---|---|---|---|
| Gaussian bin means (1/2) | −15,812 | −19,559 | 1.627 / 2.116 | 0.946 | 11.5 / 24.2 | R-hat 1.03; 22 / 62 divergences |
| Lognormal bin means (7/8) | **−10,466** | −18,181 | **1.401 / 1.857** | 0.809 | 13.6 / 19.2 | R-hat 1.010 / 1.024; 4 / 0 divergences |
| Enhancement class (3/4) | lpd −0.848, RPS 0.076 | lpd −0.859, RPS 0.079 | – | – | 4.7 / 12.5 | Clean |
| DCM depth (5/6) | MAE 15.6 m | MAE 10.8 m | – | – | ~0 | Clean |

The lognormal 1 m histogram has the lowest RMSE of any model and the second-best seen-date elpd, behind the lognormal peak profile at −9,001. On unseen dates it falls to −18,181, against −14,573 for the lognormal half-shift smooth. This means the per-bin effects do not transfer to new dates.

Compared with 5 m bins, 1 m bins add 10–25 h per fit and give no gain in accuracy or DCM depth. The one exception is coverage for the Gaussian bin means.
