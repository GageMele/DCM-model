# What has been tried, and what has not

**As of 15 November 2026.** This page is one ledger of every model idea, setting and data choice that was tested, with its verdict, followed by every idea that is still untested. Case IDs are campaign-qualified. Numbers before 27 October come from the [25 October review](../results/20261025/README.md) and [RESULTS_AND_INTERPRETATION.md](RESULTS_AND_INTERPRETATION.md). Later numbers come from [results/20261109](../results/20261109/README.md) and [results/20261115](../results/20261115/README.md).

Verdicts: **Adopt** = part of the current best setup. **Keep testing** = promising but unresolved. **Rejected** = tested, no gain or worse. **Invalid** = runs that did not converge or had bugs, so no conclusion.

## Tried

### Response, likelihood and data treatment

| Idea | Where | Result | Verdict |
|---|---|---|---|
| Gaussian on chl_norm | All campaigns to Oct 2026 | 6–7% of predictive draws below −1 (impossible negative Chl) | Superseded |
| **Lognormal on ratio = chl_norm + 1** | next150 1–4, 7–14; bins1m 7–8 | Elpd +5,600 (seen) / +4,500 (unseen) vs Gaussian; no draws below −1; 90% intervals too narrow (0.80–0.86) | **Adopt** |
| Student-t, estimated nu | extensive / targeted | nu ≈ 1, unstable | Rejected |
| Student-t, nu = 4 | extensive / targeted | Worse error; about 82% coverage | Rejected (Gaussian era) |
| Raw 1 m rows | All campaigns | Baseline | Standard for scoring |
| **5 m bin medians for fitting** | next150 22–23; round3 1–4 | About 35× faster with equal RMSE and onset error; lpd not like-for-like | **Adopt for screening** |
| 1 m bins as factor levels (histogram) | bins1m | Same accuracy as 5 m, 10–25 h more per fit | Rejected (no gain) |
| Depth domain 10–250 m vs 10–150 m | followup | 250 m alters the fitted shape above 150 m | Keep 150 m |

### Variance

| Idea | Where | Result | Verdict |
|---|---|---|---|
| Single sigma | Original model | – | Superseded |
| Regional sigma, `sigma ~ 0 + Region` | targeted onward | Interval score 7.34 → 5.49 | **Adopt** |
| Sigma by depth and Region (Gaussian, k = 5) | targeted / followup | Better interval score; worse RMSE and coverage; 2× time | Mixed |
| Sigma by depth by Ocean, Gaussian | round3 7/8 | +3,600 (seen) / +2,100 (unseen) over Gaussian half | Adopt *if Gaussian* |
| Sigma by depth by Ocean, lognormal | next150 3/4 | +570 (seen) / −520 (unseen) | Keep testing (helps seen, hurts unseen) |
| Sigma by 5 m or 1 m bin (fixed) | round3 1–2; bins1m 1–2, 7–8 | Part of the bin models | – |

### Depth / MLD structure

| Idea | Where | Result | Verdict |
|---|---|---|---|
| Regional depth smooths, no MLD (original_z) | All | Baseline, val/res RMSE 1.540/2.739 | Control |
| MLD-only, no regional smooths | equation C | RMSE 1.925 vs 1.463 | Rejected |
| Additive linear MLD | equation A, targeted | Small gain | Rejected |
| Shared depth-varying MLD, `s(depth, by = MLD)` | equation B, shape150 | 1.537/2.738, about the same as no MLD | Rejected |
| **Half-shift, `s(depth − 0.5·MLD, by = Region)`** | onset, shape150 25/267, next150 | Beat no-MLD in all 8 date-blocked comparisons; best unseen elpd in lognormal form | **Adopt** |
| Full-shift, `depth − MLD` | onset, shape150 28 | 1.505/2.624; 1 divergence; slightly worse than half on unseen dates | Rejected in favour of half |
| Ratio coordinate depth/MLD; hybrid | onset | No consistent gain | Rejected |
| Shift + distance | onset | 1.390, but 6 divergences | Invalid |
| Regional depth-varying MLD (mld_vary_region) | shape150 34, 129, 241, 249 | Best seen-date RMSE (1.401–1.409); worse or equal on unseen dates; gain concentrated in a few casts | Rejected for generalisation |
| Regional tensors t2(depth, MLD, by = Region) and shift/half-shift variants | shape150 142–196, 203–220, 271–282 | 72–81 h at 20%; mostly worse (e.g. 151: 1.604/2.608); 154 and 219 not converged | Rejected |
| Large regional 3D tensor (depth × distance × MLD) | extensive | RMSE 4.69–8.88; no convergence | Invalid |
| Fixed warps (±0.25, +0.5) on half-shift | shape150 221–232 | No warp improved both cohorts | Rejected |
| Pooled shape, Pedersen GI (Ocean curve + regional m = 1 deviations) | next150 7/8; round3 9/10 | Equal to the sigma-by-depth reference at 3–5× the cost | Rejected (no gain) |
| Pooled regional MLD (global by-MLD + shrunk regional deviations) | next150 9/10; round3 11/12 | LN: −18,964 on unseen dates at 78–124 h; Gaussian: −17,005 | Rejected |
| **Gaussian-peak profile, nonlinear (MLD on peak depth, amplitude, width)** | next150 13/14; round3 15/16 | LN: best seen elpd (−9,001); date-blocked not converged (R-hat 1.68). Gaussian: not converged | **Keep testing** (needs a convergence fix) |
| Depth as independent 5 m bins (bin_means, bin_class) | round3 1–4 | Gaussian bin means −15,036/−18,992; class RPS 0.072/0.077 | Not better than smooths |
| Depth as independent 1 m bins | bins1m 1–4, 7–8 | LN bin means: best RMSE (1.401) but unseen −18,181 | Rejected for generalisation |
| DCM depth as an ordinal outcome (10 m and 1 m classes) | round3 5–6; bins1m 5–6 | DCM depth MAE 10.8–15.6 m; 6 casts have boundary maxima | Descriptive only |

### Distance, groups and dependence

| Idea | Where | Result | Verdict |
|---|---|---|---|
| Linear standardised log-distance | All | Standard | Adopt |
| Distance varying with depth (island-mass effect) | next150 11/12; round3 13/14 | LN: +160 / +245 over the sigma-by-depth reference; −280 vs plain LN half on unseen | Keep testing (near-neutral) |
| Date + ENSO nested intercepts | All | Standard | **Adopt** |
| Date only (drop ENSO) | next150 5/6 | Same elpd; unseen onset error 18.7 vs 11.2 m | Rejected |
| Cast random intercepts / cast spline curves | extensive / followup | No held-out gain; curves worsened prediction | Rejected |
| AR(1) along depth within cast (Gaussian, scalar sigma) | followup | RMSE 1.464 → 1.443; interval score not improved | Not tested on the current best |

### Sampler, size and basis

| Idea | Result | Verdict |
|---|---|---|
| 2 chains, short runs | R-hat 1.01–1.03 | Rejected |
| 4 × (1,000 + 2,000), treedepth 14 | Clean for most models | **Adopt** |
| More draws (4,000/chain) or longer warmup | Identical results (case 244 vs 34) | Not needed |
| adapt_delta 0.99 | More tree hits, no gain | Rejected; 0.95–0.97 is fine |
| k = 7 / 11 / 15; k = 15 only in Central and South RS | More k not monotonically better; k15 helped South RS | k = 7 default |
| Training fraction 10% vs 20% | Regional MLD at 10% was worse than no-MLD on validation | 20% is used for all comparisons |
| Four different 20% cast selections | Same direction of results; some divergences | Done |
| Tighter b priors normal(0, 1) | Clean; similar error (case 249) | Optional |
| Custom vs default priors | baseline campaign | Custom used |

## Not yet tried

Each item is marked with the reason it matters. Items 1–8 can run on the existing data.

| # | Idea | Why | Notes |
|---|---|---|---|
| 1 | **Date effect that varies with depth**, e.g. `(1 + u_z \| Ocean:Region:date)` or a date-level smooth on depth relative to MLD | On unseen dates, coverage below 60 m is 0.63–0.78 and 15% of rows fall below the 5th percentile; peak size varies about ×2 between dates | Highest priority from the 15 Nov plots |
| 2 | **Sigma as a function of depth relative to MLD**, small inside the mixed layer | Every model smooths over the ratio = 1 spike from mixed-layer rows | Could use `s(depth - MLD)` in sigma, or a mixed-layer indicator |
| 3 | **Date-blocked peak profile that converges** | Best seen-date model, untested on unseen dates | Tighter priors; initial values from case 13; fit only DCM-present casts (max chl_norm > 1); 5 m bins |
| 4 | Lognormal half-shift on 5 m bins | Should match case 1/2 at under 1 h | Makes K-fold and larger fractions affordable |
| 5 | **Paired elpd differences with date-bootstrap SEs** | All current rankings lack uncertainty | Use `checks/*_date_elpd.csv`; no refits needed |
| 6 | AR(1) within cast on the lognormal / binned model (T8) | Residual lag-1 correlation of about 0.97 makes intervals too narrow | Check `stancode()` compatibility with distributional sigma and threading |
| 7 | Lognormal half + distance by depth without sigma by depth | Distance helped on top of sigma by depth; the plain LN half is the best unseen model | Cheap on 5 m bins |
| 8 | Student-t or Gamma(log) on the ratio (T4, T3 alternative) | Central Red Sea tail; lognormal intervals too narrow | Only if the PPC still shows tail problems |
| 9 | Pedersen "GS" factor-smooth pooling, `s(u, Region, bs = "fs", m = 1)` | GI pooling did not help; GS shares one smoothing parameter | Low priority |
| 10 | Smooth background + parametric peak (T10) | If the peak model converges but misfits the background | Optional |
| 11 | mgcv `bam()` pre-screen (Phase 0): shift-coefficient grid 0–1, I vs GI vs GS, global t2 | Cheap structure screen before brms | Not done |
| 12 | Lognormal models at 40–100% training and on a final full-data fit | All current comparisons use 20% | Full data is reconstruction, not validation |
| 13 | Collect the remaining shape150 results (40–100% Gaussian fits, regional k at 50/75%, tensors, warps) | 195 fits complete as of 15 Nov; results after 27 Oct not yet collected; 48 fit failures not investigated | Reporting only; no new fits |
| 14 | Equatorial deep values near chl_norm = −1 | Largest unseen-date misfit | Check the raw data first |
| 15 | Cruise-blocked and ENSO-phase-blocked stress tests | ENSO labels are nearly a function of cruise | Needs a cruise ID |
| 16 | A fresh, untouched grouped test set (or external data) | Primary and reserve cohorts have been inspected repeatedly | Needed before any final claim |

### Needs data not in the current bundle

| # | Idea | Data needed |
|---|---|---|
| T11 | Pycnocline/thermocline depth coordinate (z(N²max), z20/z25, isopycnals); density-criterion MLD (0.03 and 0.005 kg m⁻³) | Temperature and salinity profiles |
| T12 | Daytime fluorescence quenching in the 5–15 m denominator (night flag) | Cast UTC time and position |
| T13 | Is deep Chl one-for-one with surface Chl? Model log Chl with a log-surface slope | Raw Chl and raw surface Chl |
| T14 | Red Sea GAIW salinity tracer (50–75 m) | Salinity profiles |
| – | Island traits for the island-mass effect (atoll vs island, slope, reef area, bottom depth) | Island metadata |

The reasoning and literature for T1–T14 are in [research/TEST_DESIGN_20261027.md](research/TEST_DESIGN_20261027.md), with background notes in [research/](research/).
