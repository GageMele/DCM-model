# Overall results so far (15 Nov 2026)

These results combine shape150, next150, round3 and bins1m. Elpd is the joint per-date log predictive density on held-out casts, and higher is better. Seen = held-out casts from dates that were in training (primary_p20). Unseen = held-out casts from dates that were never trained on (date-blocked p20). Unseen is the real generalisation test. Elpd differences have no standard errors yet, so treat gaps of a few hundred as uncertain.

## Ranking (continuous chl_norm models)
| Rank | Setup | Seen elpd | Unseen elpd | Unseen cov90 | Fit h | Verdict |
|---|---|---|---|---|---|---|
| 1 | **Lognormal half-shift** (next150 1/2) | -11554 | **-14573** | 0.82 | 12–17 | Best on unseen dates. Current recommended setup |
| 2 | LN half + sigma by depth + distance by depth (next150 11/12) | -10823 | -14852 | 0.81 | 13 | Slightly better on seen dates, slightly worse on unseen |
| 3 | LN half + sigma by depth (next150 3/4) | -10981 | -15097 | 0.81 | 15–19 | Same pattern |
| 4 | LN pooled GI shape + sigma by depth (next150 7/8) | -10913 | -15107 | 0.81 | 58–60 | No better than 3 at 3× the cost |
| 5 | **LN peak profile** (next150 13/14) | **-9001** | not converged | – | 9 / 54 | Best on seen dates by far; the unseen fit needs fixing |
| 6 | LN 1 m histogram bin means (bins1m 7/8) | -10466 | -18181 | 0.81 | 14–19 | Good on seen dates, poor on unseen |
| 7 | Gaussian half + sigma by depth (round3 7/8) | -13550 | -16965 | 0.88 | 18–31 | Best Gaussian |
| – | Gaussian half (shape150 baseline) | -17134 | -19107 | 0.95 | 27–29 | Reference |

## What made a substantial difference
- **Lognormal likelihood (largest gain).** Elpd improves by about 5,600 on seen dates and 4,500 on unseen dates, and RMSE and onset error also improve. The cost is that the 90% intervals are too narrow, with only 81–86% coverage.
- **Sigma varying by depth.** This is large for Gaussian models (about +3,600 seen and +2,100 unseen). For lognormal models it helps seen dates (+570) but hurts unseen dates (−520).
- **Peak-profile equation.** It gives +2,550 over the LN half-shift on seen dates and has the best fitted-profile shape (lowest residual RMSE among LN fits). It has not been shown to work on unseen dates yet.
- **5 m binning of the training data.** It is about 35× faster with the same accuracy. This is a speed gain, not an accuracy gain.

## What did not help
Pooled GI shape, pooled regional MLD, distance by depth, dropping ENSO, 1 m bins (compared with 5 m) and the Gaussian peak profile. Each of these adds cost and gives at best a small gain on seen dates, with no gain on unseen dates. Regional and bin-level detail tends to fit the training dates better without carrying over to new dates.

Note: next150 cases 7–12 and round3 cases 7–14 all include sigma by depth (`sigma ~ 0 + Region + s(depth, by = Ocean, k = 5)`), so they are compared against cases 3/4, not 1/2.

## Histogram (ordinal) models
The class and DCM-depth models are scored on different targets, so they can't be ranked against the table above. The DCM-depth error is about 11–16 m whether the bins are 5 m or 1 m.
