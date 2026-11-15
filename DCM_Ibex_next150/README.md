# next150: likelihood, variance, pooling and peak-profile tests (Oct–Nov 2026)

This package was installed on Ibex by [04_create_next150.sh](04_create_next150.sh) on 28 October 2026. Cases 22–23 were added the same day by [05_add_binned_pair.sh](05_add_binned_pair.sh). The design follows the literature review in [docs/research/TEST_DESIGN_20261027.md](../docs/research/TEST_DESIGN_20261027.md), whose test labels (T1–T14) appear below.

| Item | Value |
|---|---|
| Experiment | `/ibex/user/meleg/DCM_MODEL/experiments/next150_20261028T114047Z_653459` |
| Pointer | `logs/latest_next150_experiment.txt` |
| Token / job prefix | `f6c43cea00` / `dcm-s150-f6c43cea00-<stage>` |
| Controller | `controller_001`, job 54722393 |
| Data | Frozen shape150 data, symlinked. Primary p20: 216 training casts, 30,455 rows. Date-blocked p20 (seed 20261022): 216 casts, 30,456 rows |
| Status | All 23 cases checked by 9 Nov 2026 |
| Results | [results/20261109](../results/20261109/README.md), with an overall ranking in [results/20261115](../results/20261115/README.md) |

## Files

The standalone `models.R`, `scoring.R`, `02_fit_case.R`, `preflight.R` and `cases.tsv` are the live experiment files after `05_add_binned_pair.sh` ran. That script added `bin_depth()` and cases 22–23. The installer `04_create_next150.sh` embeds the earlier version (cases 1–21) and checks its md5s. Scheduler, worker, data and onset/profile helpers are copied from the shape150 experiment (see [DCM_Ibex_shape_tests](../DCM_Ibex_shape_tests/)).

## Shared settings

All new fits train on 20% of casts (p20) with depth basis k = 7, 4 chains × (1,000 warmup + 2,000 retained), `adapt_delta = 0.95`, `max_treedepth = 14`, 16 CPUs (4 chains × 4 threads) and 96 GB. Cases 22–23 use 48 GB. Each test is run twice: on the primary split (held-out casts from seen dates) and on the date-blocked split, seed 20261022 (held-out casts from unseen dates).

Priors: Intercept `student_t(3, 0, 1)` for lognormal or `student_t(3, 0.73, 5.448)` for Gaussian; `b ~ normal(0, 2)`; `sd` and `sds ~ student_t(3, 0, 1.5)`; sigma coefficients `normal(0, 1)`; sigma smooth `sds ~ normal(0, 0.5)`.

The lognormal cases model `ratio = chl_norm + 1 = Chl / surface`. Ratios ≤ 0 are floored at 0.01: 10 of 30,455 primary rows and 4 of 30,456 date-blocked rows. Scores are always computed on the unmodified chl_norm. Because `ratio` is a pure shift of chl_norm, the log density needs no Jacobian, and elpd is directly comparable between lognormal and Gaussian fits.

## Cases

The anchor is the Gaussian half-shift model, shape150 case 25 (primary) and case 267 (date-blocked): `s(depth - 0.5*MLD, by = Region, k = 7)` + distance_z + date and ENSO intercepts, with `sigma ~ 0 + Region`.

| Cases | Test (report label) | Change from the anchor |
|---|---|---|
| 1–2 | Lognormal likelihood (T3) | `ratio ~ lognormal()`; no predictive draws below chl_norm = −1 |
| 3–4 | + depth-varying sigma (T3) | Lognormal + `sigma ~ 0 + Region + s(depth, by = Ocean, k = 5)` |
| 5–6 | Date-only groups (T2) | Gaussian; drops `(1 \| Ocean:Region:enso)` |
| 7–8 | Pooled shape, Pedersen "GI" (T5) | `s(u, by = Ocean) + s(u, by = Region, m = 1)`; lognormal + sigma by depth |
| 9–10 | Pooled regional MLD (T6) | + `s(u, by = mld_z)` + per-region `s(u, by = region_mld_i, k = 5, m = 1)`; lognormal + sigma by depth |
| 11–12 | Island mass effect (T7) | Distance effect varies with depth: `s(u, by = distance_z, k = 5)`; lognormal + sigma by depth |
| 13–14 | Gaussian-peak profile, nonlinear (T9) | See below; lognormal; date groups only |
| 15–21 | Rescore shape150 anchors 19, 25, 34, 129 (primary) and 265, 267, 270 (date-blocked) | No refit; the saved fits are linked and rescored to add elpd |
| 22–23 | 5 m binned half-shift (T1) | Same as anchor 25/267, fitted on each cast's 5 m bin medians (6,048 rows); scored on raw 1 m rows |

Here `u = depth - 0.5 * MLD`.

The peak profile (cases 13–14) is, with `z = (depth - 0.5*MLD)/100`:

```
median ratio(z) = exp(lbg + slope*z) + exp(lamp) * exp(-0.5*((z - zm)/exp(lwid))^2)
lbg   ~ 0 + Region            slope ~ 0 + Region
lamp  ~ 0 + Region + mld_z + distance_z + (1|g|Ocean:Region:date)
zm    ~ 0 + Region + mld_z + (1|g|Ocean:Region:date)
lwid  ~ 0 + Region + mld_z
```

MLD enters peak depth, amplitude and width separately. This makes "MLD moves the DCM" distinguishable from "MLD weakens the DCM". The fits start at `init = 0`.

## Headline results

See [results/20261109](../results/20261109/README.md) for the full tables.

- **Lognormal clearly beats Gaussian.** For the same half-shift equation, joint per-date elpd is −11,554 vs −17,134 on seen dates and −14,573 vs −19,107 on unseen dates. The lognormal 90% intervals are too narrow: coverage is 0.80–0.86.
- **On unseen dates the plain lognormal half-shift (case 2) is best.** Sigma by depth, GI pooling and distance by depth each improve seen dates slightly but make unseen dates slightly worse.
- **The peak profile is the best model on seen dates** (−9,001, 8.6 h). Its date-blocked fit (case 14) did not converge: R-hat 1.68 and 948 divergences.
- **5 m binning was about 35× faster** (0.8–0.9 h vs 27–29 h) with equal accuracy.
- **Pooled regional MLD was expensive** (78–124 h) and clearly worse on unseen dates (−18,964).
- **Dropping ENSO left elpd unchanged** but worsened unseen-date onset error (18.7 vs 11.2 m).
