# next150 results (all 23 cases checked, read 2026-11-09)

Scores are on whole held-out casts. "elpd (date)" is the joint log predictive density summed over held-out date groups (higher is better; comparable only within a split). lpd/row is per raw row on the chl_norm scale. Onset is the MAE of the +50% onset depth (m). Hours are fit time; for the re-scored anchors (15-21) it's their original shape150 fit time. LN = lognormal, G = Gaussian.

## Seen dates (primary p20), sorted by elpd

| case | model | family | h | val RMSE | res RMSE | val cov90 | lpd/row | elpd (date) | onset m |
|---|---|---|---|---|---|---|---|---|---|
| 13 | peak profile | LN | 8.6 | 1.419 | 2.641 | 0.849 | -0.801 | **-9001** | 11.3 |
| 11 | half + distance-by-depth | LN | 12.6 | 1.358 | 2.738 | 0.859 | -0.797 | -10823 | 10.6 |
| 7 | half, pooled shape (GI) | LN | 58.3 | 1.356 | 2.746 | 0.857 | -0.802 | -10913 | 9.9 |
| 3 | half + sigma by depth | LN | 18.5 | 1.357 | 2.745 | 0.856 | -0.804 | -10981 | 9.9 |
| 9 | half + pooled regional MLD | LN | 123.6 | 1.359 | 2.707 | 0.846 | -0.840 | -11110 | 9.0 |
| 1 | half | LN | 16.5 | 1.377 | 2.679 | 0.863 | -0.843 | -11554 | 10.1 |
| 22 | half, 5 m bins | G | **0.8** | 1.418 | 2.767 | 0.919 | -1.212 | -16235 | 12.3 |
| 18 | regional MLD (shape150 129) | G | 33.9 | 1.401 | **2.506** | 0.913 | -1.230 | -16789 | 14.0 |
| 17 | regional MLD (shape150 34) | G | 33.2 | 1.409 | 2.559 | 0.913 | -1.234 | -16901 | 13.6 |
| 16 | half (shape150 25) | G | 27.2 | 1.489 | 2.648 | 0.922 | -1.236 | -17134 | 13.1 |
| 5 | half, date-only groups | G | 5.3 | 1.501 | 2.650 | 0.922 | -1.256 | -17147 | 12.9 |
| 15 | no MLD (shape150 19) | G | 19.7 | 1.540 | 2.739 | 0.927 | -1.236 | -17298 | 14.6 |

## Unseen dates (date-blocked p20, seed 20261022), sorted by elpd

| case | model | family | h | val RMSE | res RMSE | val cov90 | lpd/row | elpd (date) | onset m |
|---|---|---|---|---|---|---|---|---|---|
| 2 | half | LN | 12.4 | 1.845 | 2.557 | 0.822 | -1.147 | **-14573** | 11.7 |
| 12 | half + distance-by-depth (5 div) | LN | 13.1 | 1.839 | 2.649 | 0.805 | -1.171 | -14852 | 11.4 |
| 4 | half + sigma by depth | LN | 14.5 | 1.842 | 2.651 | 0.808 | -1.174 | -15097 | 11.0 |
| 8 | half, pooled shape (GI) | LN | 59.7 | 1.838 | 2.653 | 0.808 | -1.180 | -15107 | 11.0 |
| 23 | half, 5 m bins | G | **0.9** | 1.943 | 2.756 | 0.925 | -1.363 | -18619 | **10.7** |
| 10 | half + pooled regional MLD (2 div) | LN | 77.7 | 1.950 | 2.558 | 0.843 | -1.220 | -18964 | 13.3 |
| 20 | half (shape150 267) | G | 29.1 | 1.963 | 2.637 | 0.949 | -1.429 | -19107 | 11.2 |
| 6 | half, date-only groups | G | 6.2 | 2.002 | 2.652 | 0.935 | -1.515 | -19148 | 18.7 |
| 21 | regional MLD (shape150 270, 7 div) | G | 52.6 | 2.121 | 2.721 | 0.944 | -1.562 | -19817 | 12.0 |
| 19 | no MLD (shape150 265) | G | 32.6 | 2.002 | 2.722 | 0.950 | -1.449 | -19904 | 12.8 |
| 14 | peak profile | LN | 53.9 | not converged: R-hat 1.68, 948 divergences | | | | | |

## What it says
1. **Lognormal beats Gaussian clearly.** Same half-shift equation: elpd -11554 vs -17134 on seen dates and -14573 vs -19107 on unseen dates. RMSE and onset error are also lower. The catch: lognormal 90% intervals are too narrow (80-86% coverage, against 91-95% for Gaussian).
2. **On unseen dates the plain lognormal half-shift is best.** Depth-varying sigma, the pooled GI shape and distance-by-depth don't add anything there. Pooled regional MLD is clearly worse (-18964) and takes 3-5 days. That matches shape150: regional MLD helps seen dates but doesn't carry to new dates.
3. **The peak profile is the best on seen dates by a wide margin (-9001, 8.6 h), but its unseen-date fit did not converge.** It needs fixing (tighter priors or inits from the primary fit) before it can be trusted.
4. **5 m binning is a large speed-up at no cost.** 0.8-0.9 h vs 27-29 h for the same Gaussian model, with similar or better validation RMSE, elpd and onset; reserve RMSE is ~0.1 worse. That makes full-data and lognormal binned fits cheap.
5. **Dropping the ENSO term** leaves elpd unchanged but worsens unseen-date onset error (18.7 vs 11.2 m).
6. Reserve RMSE (~2.5-2.8) barely separates models; the regional-MLD Gaussian fits are lowest there.

Caveat: elpd differences have no standard errors yet. The per-date elpd CSVs in each run's checks/ folder allow paired differences with SEs; that's the next scoring step.
