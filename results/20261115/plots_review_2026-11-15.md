# Top-5 posterior plots: what they show (15 Nov 2026)

Source: plots/top5_20261115T155349Z (job 55624665/66). All panels show held-out casts. Seen dates cover 108 casts and 15,226 rows; unseen dates cover 117 casts and 16,497 rows. The predictions use 300 posterior draws, and dates the model never saw get freshly drawn date effects.

## Summary table (held-out rows)
| Model | Split | RMSE | MAE | cov50 | cov90 | 90% width | PIT<0.05 | PIT>0.95 |
|---|---|---|---|---|---|---|---|---|
| LN half-shift | seen | 1.378 | 0.655 | 0.523 | **0.859** | 1.33 | 0.063 | 0.072 |
| LN + sigma by depth | seen | 1.358 | 0.645 | 0.514 | 0.852 | 1.10 | 0.070 | 0.072 |
| LN + distance by depth | seen | **1.357** | 0.644 | 0.514 | 0.852 | 1.10 | 0.071 | 0.071 |
| LN 1 m histogram | seen | 1.402 | 0.672 | 0.516 | **0.873** | 1.23 | 0.056 | 0.063 |
| LN peak profile | seen | 1.412 | **0.634** | 0.513 | 0.845 | 0.95 | 0.067 | 0.084 |
| LN half-shift | unseen | 1.848 | 0.780 | 0.538 | **0.818** | 1.82 | 0.145 | 0.033 |
| LN + sigma by depth | unseen | 1.839 | **0.766** | 0.531 | 0.805 | 1.55 | 0.152 | 0.039 |
| LN + distance by depth | unseen | **1.837** | 0.768 | 0.521 | 0.801 | 1.56 | 0.155 | 0.039 |
| LN 1 m histogram | unseen | 1.855 | 0.782 | 0.546 | 0.807 | 1.73 | 0.156 | 0.034 |

The nominal targets are cov50 = 0.50 and cov90 = 0.90. For a calibrated model, the two PIT tail columns should each be about 0.05.

## What the plots say
1. **The four smooth lognormal models give almost the same predictions.** Their mean profiles nearly overlap in every region. The 1 m histogram follows the same curves but is noisy from metre to metre. The differences in elpd come mainly from interval width, not from the predicted shape.
2. **The shape and depth of the DCM are predicted well, but the size of the peak is underestimated.** In the example casts, the depth of the peak is about right in all 8 regions. Strong peaks are predicted at about half their observed size: for example, Central RS observed 22 vs 9–10 predicted, Mariana 5 vs 2.5 and Hawaiian 4 vs 2. Peak size depends on the date, so on a new date the model reverts to the regional average.
3. **On seen dates, calibration is close to correct.** The PIT histogram is flat apart from a mild U shape (6–8% in each tail instead of 5%). That pattern means the intervals are slightly too narrow.
4. **On unseen dates, coverage breaks down below about 60 m.** Above 50 m, 90% coverage is 0.90–0.95. Below 60 m it falls to 0.63–0.78. About 15% of rows fall below the model's 5th percentile, compared with only about 3.5% above the 95th. So on new dates the model often predicts too much enhancement at depth.
5. **The plain half-shift model has the best coverage on unseen dates** because its intervals are wider at depth. That is why it has the best unseen elpd even though its mean predictions are no better.
6. **The observed values have a low tail the models cannot produce.** On unseen dates the observed density has a bump at 3–30% of surface Chl, and no model generates any values there. The clearest example is Equatorial on unseen dates: below about 80 m the observed median drops to chl_norm ≈ −0.9, while the models predict −0.5 to −0.75. Hawaiian and Mariana also come out lower at depth on unseen dates than predicted. The Equatorial values near −1 are worth checking in the raw data before modelling around them.
7. **There is a sharp spike at Chl = surface (ratio 1), which every model smooths over.** Those rows are in the mixed layer, where Chl is close to its surface value by construction. Making sigma depend on depth in metres does not capture this. The peak-profile model gets closest on seen dates, but it also has the narrowest intervals and the most U-shaped PIT.
8. **The peak profile on seen dates** fits the overall distribution best, but it underestimates the regional medians at the DCM in North RS and South RS.

## Suggested next steps
- **Let the date effect vary with depth.** For example, use a date-level random slope or smooth on depth relative to MLD, rather than a date intercept only. This targets points 2 and 4: new dates would get wider intervals at depth, and the PIT spike near 0 would shrink.
- **Model sigma as a function of depth relative to the MLD,** with very small values inside the mixed layer. This targets point 7.
- **Fix the unseen-date peak-profile fit** using tighter priors or initial values taken from the seen-date fit.
- **Check the Equatorial casts** that drop to near-zero Chl below 80 m.
