# Current best model setups

**As of 15 November 2026.** This page names the setups to start from and gives their exact formulas. The evidence is in [results/20261115](../results/20261115/README.md) and [results/20261109](../results/20261109/README.md); the full history is in [RESULTS_AND_INTERPRETATION.md](RESULTS_AND_INTERPRETATION.md).

No final production model has been chosen. The rankings below come from one 20% training fraction, one primary split and one date-blocked split, with no standard errors on elpd yet.

## Which setup for which purpose

| Purpose | Use | Why |
|---|---|---|
| Predicting casts on dates the model has not seen (generalisation) | **Lognormal half-shift**, next150 case 1/2 | Best unseen-date elpd (−14,573), the widest honest intervals at depth, and simple |
| Describing profiles on sampled dates | **Lognormal peak profile**, next150 case 13 | Best seen-date elpd (−9,001), 8.6 h, interpretable parameters. The unseen-date fit has not converged yet |
| Lowest point error on sampled dates | Lognormal 1 m histogram, bins1m case 7 | RMSE 1.401. It generalises poorly (unseen elpd −18,181) |
| Fast screening of new structure | Any of the above fitted on **5 m bin medians** | About 35× faster with the same accuracy (next150 22/23) |
| Best Gaussian, if a Gaussian is required | Gaussian half-shift + sigma by depth, round3 7/8 | −13,550 / −16,965 |

The regional MLD-varying Gaussian (shape150 cases 34/129) was the October RMSE frontrunner on seen dates. It is not a recommended setup. On unseen dates its validation RMSE was worse than no MLD in both date-blocked seeds: 2.121 vs 2.002 (cases 270 vs 265, with 7 divergences) and 2.135 vs 2.127 (279 vs 274). Its reserve RMSE was equal or better. The half-shift beat no-MLD in all 8 date-blocked comparisons.

## 1. Lognormal half-shift (recommended default)

```r
ratio <- pmax(chl_norm + 1, 0.01)          # Chl / surface Chl; 10 of 30,455 rows floored
u     <- depth - 0.5 * MLD_temp_MIX        # "depth_mld_half"

bf(ratio ~ 1 + Region + distance_z +
     s(depth_mld_half, by = Region, bs = "cr", k = 7) +
     (1 | Ocean:Region:date) + (1 | Ocean:Region:enso),
   sigma ~ 0 + Region,
   family = lognormal())
```

Priors:
- Intercept `student_t(3, 0, 1)`;
- `b ~ normal(0, 2)`;
- `sd` and `sds ~ student_t(3, 0, 1.5)`;
- sigma coefficients `normal(0, 1)`.

Sampler: 4 chains × (1,000 warmup + 2,000 retained), `adapt_delta = 0.95`, `max_treedepth = 14`, 4 threads per chain (16 CPUs), 96 GB. Fit time is 12–17 h on 20% of casts (216 casts, about 30k raw rows). On 5 m bins it is expected to take well under 1 h; this has not been run yet.

Scores (seen / unseen dates):
- elpd −11,554 / −14,573;
- RMSE 1.377 / 1.845;
- cov90 0.86 / 0.82;
- onset MAE 10.1 / 11.7 m.

Known weaknesses:
- The intervals are too narrow overall.
- On unseen dates, coverage below 60 m falls to about 0.65–0.75, and the model over-predicts deep enhancement.
- Strong DCM peaks are predicted at about half their size on new dates.
- The model cannot reproduce the sharp ratio = 1 spike from mixed-layer rows.

Code: [DCM_Ibex_next150/models.R](../DCM_Ibex_next150/models.R) (`campaign_spec`, equation `half`, family `lognormal`).

## 2. Lognormal peak profile (most interpretable; needs an unseen-date fix)

```r
# z = depth_hm = (depth - 0.5 * MLD) / 100
bf(ratio ~ log(exp(lbg + slope * depth_hm) +
               exp(lamp) * exp(-0.5 * ((depth_hm - zm) / exp(lwid))^2)),
   lbg   ~ 0 + Region,
   slope ~ 0 + Region,
   lamp  ~ 0 + Region + mld_z + distance_z + (1 | g | Ocean:Region:date),
   zm    ~ 0 + Region + mld_z + (1 | g | Ocean:Region:date),
   lwid  ~ 0 + Region + mld_z,
   sigma ~ 0 + Region,
   nl = TRUE, family = lognormal())
```

Priors:
- `lbg ~ normal(0, 0.5)`;
- `slope ~ normal(-1, 1)`;
- `lamp ~ normal(0, 1)`;
- `zm ~ normal(0.6, 0.3)`;
- `lwid ~ normal(-1.4, 0.5)`;
- `sd(lamp) ~ normal(0, 0.5)`;
- `sd(zm) ~ normal(0, 0.3)`.

The fit uses `init = 0`.

Scores: seen elpd −9,001, RMSE 1.419, cov90 0.85, 8.6 h. The date-blocked fit (case 14) took 53.9 h and did **not** converge (R-hat 1.68, 948 divergences).

Posterior (seen dates):
- MLD adds no displacement beyond the half-shift (`zm ~ mld_z` 0.005 [−0.002, 0.012]).
- MLD lowers peak amplitude (−7.5% per SD) and slightly narrows the peak (−3% per SD).
- Date-to-date amplitude variation is large (sd 0.68 on the log scale).

See [results/20261115](../results/20261115/README.md#peak-profile-posterior-how-mld-acts-on-the-dcm).

## 3. Lognormal 1 m histogram (bin means)

```r
bf(ratio ~ 1 + Region + distance_z + mld_z +
     (1 + mld_z || depth_m) + (1 + mld_z || Region:depth_m) +
     (1 | Ocean:Region:date) + (1 | Ocean:Region:enso),
   sigma ~ 0 + Region + depth_m,
   family = lognormal())
```

Here `depth_m` is a 141-level factor, one level per metre from 10 to 150 m. Scores (seen / unseen): elpd −10,466 / −18,181, RMSE 1.401 / 1.857, 14–19 h. Code: [DCM_Ibex_bins1m/models.R](../DCM_Ibex_bins1m/models.R).

## Choices that are settled

- **Likelihood:** lognormal on the ratio beats Gaussian on chl_norm everywhere tested.
- **Depth coordinate:** the half-shift `depth − 0.5·MLD` is the best-supported MLD representation on unseen dates. The peak model independently finds no extra MLD displacement beyond 0.5.
- **Regional curves:** keep `s(u, by = Region)` (Pedersen model "I"). GI pooling neither helps nor hurts.
- **Basis size:** k = 7 is sufficient at 20%. More k was not monotonically better in shape150.
- **Draws:** 2,000 retained per chain is enough. Shape150 case 244, with 4,000 draws, was identical to case 34.
- **Training data:** 5 m bin medians are acceptable for screening and much faster. Score on raw rows.
- **Group effects:** keep date and ENSO. Dropping ENSO did not change elpd but worsened unseen-date onset error.

## Choices that are still open

- Whether sigma should vary with depth, and how. It helps seen dates and hurts unseen dates. A sigma tied to depth relative to MLD (small inside the mixed layer) has not been tried.
- How to give unseen dates honest uncertainty at depth: a depth-varying date effect, or AR(1).
- Whether the peak profile can be made to converge on the date-blocked split.
- Whether a pycnocline or thermocline coordinate beats MLD. This needs T/S profiles; it is test T11 in the [test design](research/TEST_DESIGN_20261027.md).
- Performance at larger training fractions and on a fresh, final grouped test set.

See [TRIED_AND_UNTRIED.md](TRIED_AND_UNTRIED.md) for the full list.
