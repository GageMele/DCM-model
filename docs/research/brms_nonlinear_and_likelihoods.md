# brms nonlinear profile models and likelihoods for a relative-enhancement (chl_norm) response

Scope: concrete brms 2.23 / cmdstanr 0.9 / CmdStan 2.37 patterns for (a) parametric Gaussian-peak-on-background DCM profile models, (b) likelihoods for chl_norm = (Chl - surface_chl)/surface_chl, (c) distributional (sigma) models, (d) along-depth autocorrelation, and (e) compute tricks. Each section separates what sources say (Cited Findings) from my reasoning and code (Inferences). The code was written against the documented brms syntax but was **not executed in this session**. Run `stancode()` / `standata()` on it before any long fit.

Data context assumed in the code: data frame `d` with columns `chl` (raw fluorescence-derived Chl), `surface_chl`, `chl_norm`, `depth` (m, 10-150), `Region` (factor, 8 levels), `cast` (factor, 1080 levels), `MLD`, `distance_z`, `date`, `enso`, `Ocean`.

---

## Q1. How to write a Gaussian-peak-on-background (Lewis 1983 / Platt) profile as a brms nonlinear model

### Takeaway
brms fits this model directly with `bf(y ~ <peak formula>, <nlpar> ~ <linear predictor>, nl = TRUE)`. Each nonlinear parameter (zm, log-amplitude, log-width, background) gets its own full linear predictor, which can include Region, MLD, and correlated cast-level random effects through the `|ID|` syntax. Priors on every nonlinear parameter are mandatory, and weak priors cause multimodality across chains. Parameterize positive quantities on the log scale, scale depth to O(1), and use peak amplitude rather than integrated biomass h.

### Cited Findings
- The standard parametric form is P(z) = P_max·exp[-(z - z_m)²/(2σ²)] with P_max = h/(σ√(2π)). Here h is the vertically integrated Chl, z_m the depth of the maximum, and σ the width, so the layer thickness is 2σ and the layer spans z_m ± σ. The paper attributes this form to Lewis et al. (1983) and notes background terms (constant or depth-dependent, citing Richardson et al. 2003, Uitz et al. 2006, Mignot et al. 2011) for profiles with nonzero surface Chl. It also reports that Gaussian fits reproduced DCM depth and thickness reasonably well at SEATS/HOT/BATS but estimated intensity poorly. — [Gong et al. 2015, Biogeosciences 12:905](https://bg.copernicus.org/articles/12/905/2015/bg-12-905-2015.pdf)
- brms nonlinear syntax: `bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE)`, with priors given per nonlinear parameter via `nlpar`, e.g. `prior(normal(1, 2), nlpar = "b1")`. Bounded priors use `lb`/`ub`, e.g. `prior(beta(1, 1), nlpar = "guess", lb = 0, ub = 1)`. — [brms_nonlinear vignette](https://paulbuerkner.com/brms/articles/brms_nonlinear.html)
- The vignette's loss-curve example puts a random intercept on one nonlinear parameter (`ult ~ 1 + (1|AY), omega ~ 1, theta ~ 1`) and uses `control = list(adapt_delta = 0.9)`. — [brms_nonlinear vignette](https://paulbuerkner.com/brms/articles/brms_nonlinear.html)
- The vignette states that priors on population-level effects are often mandatory to identify nonlinear models. Chains converging to different regions signals an identification problem, and the suggested fix is narrower priors. If the nonlinear formula already includes the link, set the family link to identity. — [brms_nonlinear vignette](https://paulbuerkner.com/brms/articles/brms_nonlinear.html)
- Bürkner (2018) says to "think of non-linear parameters as placeholders for linear predictor terms rather than as parameters themselves." Each has its own linear predictor that is "compatible with all terms supported by brms", which includes splines and GPs. The paper shows correlated varying effects across nonlinear parameters with `b1 ~ z + (1|ID|g); b2 ~ (1|ID|g); b3 ~ (1|ID|g)`. It warns that "non- or weakly informative priors in non-linear models often leads to problems", that "Stan works best when parameters are roughly on the same order of magnitude", and suggests rescaling. — [Bürkner 2018, R Journal / arXiv 1705.11123](https://arxiv.org/pdf/1705.11123)
- The `brmsformula` reference says: priors on nonlinear parameters are mandatory. Terms sharing an ID are correlated even across different nonlinear-parameter formulas, e.g. `a1 ~ 1 + (1 |2| g), a2 ~ x + (x |2| g)` or equivalently `(1 | gr(g, id = 2))`. `center` defaults to FALSE for nonlinear parameters. `loop = TRUE` (default) evaluates the nonlinear formula inside a loop over observations, and `FALSE` evaluates it outside. `nlf()`/`lf()` combine nonlinear and distributional formulas, e.g. `bf(y ~ eta, nl = TRUE) + lf(eta ~ 1 + x) + nlf(sigma ~ tau * sqrt(eta)) + lf(tau ~ 1)`. — [brmsformula reference](https://paulbuerkner.com/brms/reference/brmsformula.html)
- brms 2.20.0 added matrix data as nonlinear covariates. 2.22.0 made internal functions accessible when evaluating nonlinear formulas. 2.23 fixed `log_sum_exp()` in `nlf()` formulas. — [brms NEWS](https://paul-buerkner.r-universe.dev/brms/NEWS)
- A Stan-forum user who tried to combine a smooth `s()` on the main formula with separate nonlinear parameters in a single `bf()` hit errors. The only reply said brms spline handling made that particular construction unworkable. No working code or developer statement was given. — [Stan Discourse thread 19519](https://discourse.mc-stan.org/t/is-it-possible-to-have-a-linear-smoothing-spline-and-non-linear-components-in-a-single-model-in-brms/19519)

### Inferences
**Reparameterization choices (my reasoning):**
- **Use amplitude, not area.** Use log-amplitude `logA` (= log P_max) instead of h. h = A·σ·√(2π) is strongly negatively correlated with σ in the posterior. Peak height and width are closer to orthogonal. Derive h as a generated quantity afterwards.
- **Scale depth to O(1).** Use `z = depth/100`, giving 0.1-1.5, as recommended above.
- **Keep z_m inside the observed window.** Either use a tight prior, or map it through `zm = 0.1 + 1.4*inv_logit(zmr)`. `lb`/`ub` on an nlpar prior bounds every coefficient of that nlpar. That only makes sense for a cell-means `0 + Region` predictor with no slopes, so it breaks once MLD slopes are added.
- **Watch for weakly identified profiles.** If a profile has no peak (mixed or monotone), `logA → -∞` and `zm`/`logw` become unidentified. Cast-level random effects with partial pooling plus informative priors on `logw` handle this better than per-cast fits. Expect divergences if many casts lack a DCM; a mixture (peak vs no-peak) is a heavier alternative.
- **Response scale.** Fit on the log-ratio scale, `lratio = log(chl/surface_chl) = log1p(chl_norm)` (see Q2). The background is then `b0 + slope*z`, with b0 ≈ 0 near the reference depth. The peak becomes multiplicative enhancement, exp(peak).

**Code (written to documented syntax; not run):**

```r
library(brms)      # 2.23
library(cmdstanr)  # 0.9

d <- d |>
  dplyr::filter(chl > 0, surface_chl > 0) |>       # see Q2 for alternatives to dropping
  dplyr::mutate(
    lratio = log(chl / surface_chl),               # == log1p(chl_norm)
    z      = depth / 100,                          # 0.1 - 1.5
    mld_z  = as.numeric(scale(MLD)),
    cast   = factor(cast)
  )

f_peak <- bf(
  lratio ~ b0 + slope * z +
           exp(logA) * exp(-0.5 * ((z - zm) / exp(logw))^2),
  b0    ~ 0 + Region,
  slope ~ 0 + Region,
  logA  ~ 0 + Region + mld_z + (1 | c | cast),
  zm    ~ 0 + Region + mld_z + (1 | c | cast),
  logw  ~ 0 + Region + (1 | c | cast),
  nl = TRUE
)
# Optional: add (1|Ocean:Region:date) to zm/logA as in the current model.
# Cast is nested in date, so a cast-level RE probably makes the date RE redundant
# for these parameters.

pri_peak <- c(
  prior(normal(0, 0.5),   nlpar = "b0"),
  prior(normal(0, 1),     nlpar = "slope"),
  prior(normal(0, 1),     nlpar = "logA"),                  # exp(0)=1 log-unit peak ≈ 2.7x enhancement
  prior(normal(0, 0.5),   nlpar = "logA", coef = "mld_z"),
  prior(normal(0.7, 0.3), nlpar = "zm"),                    # ~70 m (z = depth/100)
  prior(normal(0, 0.2),   nlpar = "zm",   coef = "mld_z"),
  prior(normal(-1.6, 0.5), nlpar = "logw"),                 # exp(-1.6)*100 ≈ 20 m sigma-width
  prior(exponential(2),   class = "sd", nlpar = "logA"),
  prior(exponential(5),   class = "sd", nlpar = "zm"),
  prior(exponential(2),   class = "sd", nlpar = "logw"),
  prior(lkj(2),           class = "cor"),
  prior(exponential(2),   class = "sigma")
)

# Inspect generated code / parameter names before fitting
stancode(f_peak, data = d, prior = pri_peak)
sd0 <- standata(f_peak, data = d, prior = pri_peak)   # gives K_zm, K_logA, ... used below

nR <- nlevels(d$Region)
init_fun <- function() list(
  b_b0    = rep(0, sd0$K_b0),
  b_slope = rep(0, sd0$K_slope),
  b_logA  = c(rep(0, nR), 0)[seq_len(sd0$K_logA)],
  b_zm    = c(rep(0.7, nR), 0)[seq_len(sd0$K_zm)],
  b_logw  = rep(-1.6, sd0$K_logw)
)  # unspecified parameters get CmdStan's default random inits

fit_peak <- brm(
  f_peak, data = d, family = gaussian(), prior = pri_peak,
  backend = "cmdstanr", chains = 4, cores = 4, threads = threading(4),
  init = init_fun, normalize = FALSE,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  iter = 1500, warmup = 750, seed = 1
)
```

**Bounded-z_m variant:** replace `zm` with `(0.1 + 1.4 * inv_logit(zmr))` in the formula, give `zmr ~ 0 + Region + mld_z + (1|c|cast)`, and put the prior on `nlpar = "zmr"` with `normal(0, 1.5)`.

**Hybrid (parametric peak + smooth background):** put the smooth into the background nonlinear parameter. Per Bürkner (2018), nlpar formulas accept all brms terms. Keep k small so the smooth cannot absorb the peak.

```r
f_hybrid <- bf(
  lratio ~ bg + exp(logA) * exp(-0.5 * ((z - zm) / exp(logw))^2),
  bg   ~ 0 + Region + s(z, by = Region, k = 4),
  logA ~ 0 + Region + mld_z + (1 | c | cast),
  zm   ~ 0 + Region + mld_z + (1 | c | cast),
  logw ~ 0 + Region + (1 | c | cast),
  nl = TRUE
)
# priors: as above, plus prior(normal(0,0.5), nlpar="bg") and
# prior(exponential(2), class="sds", nlpar="bg")
```

- The failed forum attempt (thread 19519) put `s()` on the main left-side formula. The documented route is `s()` inside an nlpar formula, as shown here. I could not find a published brms example of this exact use, so confirm with `stancode()` first.
- Identifiability: a flexible smooth background and a Gaussian peak trade off against each other. Diagnose by checking the posterior correlation between `sds_bg` and `logA`/`logw` and by comparing across chains.
- Sanity check: fit the pure-parametric and hybrid models on a 50-100-cast subsample, with `pathfinder` or short chains, before the full fit.

### Gaps
- I found no published paper that fits a Lewis/Platt DCM profile in brms or Stan. The oceanographic literature I found fits it with nonlinear least squares per profile, and Gong et al. (2015) do not report their fitting procedure.
- I did not confirm whether `loop = FALSE` auto-vectorizes `*`/`/` in the nonlinear formula. Test it with `stancode(update(f_peak, loop = FALSE), ...)`. It is a likely speed win if it compiles.
- No runtime numbers exist for brms nonlinear models of this size (~150k rows, ~4k random effects).

---

## Q2. Likelihoods for a ratio-to-surface response

### Takeaway
The current Gaussian on chl_norm puts mass below -1 because nothing enforces Chl > 0. Three better options:
- **Log-ratio model:** Gaussian or Student-t on `log(chl/surface_chl) = log1p(chl_norm)`. This is the same as `lognormal()` on chl with `offset(log(surface_chl))`.
- **Gamma(log link):** on chl with the same offset. Here μ is the mean ratio directly.
- **Skew-normal:** on chl_norm. This keeps the original scale and a mean-parameterized μ, but still does not enforce > -1.

The log-scale models give enhancement as exp(·) - 1, and you must choose whether to report the median or the mean ratio.

### Cited Findings
- brms lognormal density: f(y) = 1/(√(2π)σy)·exp(-½((log y - μ)/σ)²). σ is "the residual standard deviation on the log-scale", so μ is on the log scale. Default link is identity for μ and log for σ. — [brms_families vignette](https://paulbuerkner.com/brms/articles/brms_families.html); [brmsfamily reference](https://paulbuerkner.com/brms/reference/brmsfamily.html)
- brms Gamma density: f(y) = (α/μ)^α/Γ(α)·y^(α-1)·exp(-αy/μ), so μ is the mean and α the shape. brms's default Gamma link is log, unlike stats::Gamma which defaults to inverse. — [brms_families vignette](https://paulbuerkner.com/brms/articles/brms_families.html); [brmsfamily reference](https://paulbuerkner.com/brms/reference/brmsfamily.html)
- brms parameterizes skew_normal "in terms of the mean": ω = σ/√(1 - (2/π)·α²/(1+α²)) and ξ = μ - ω·α/√(1+α²)·√(2/π). μ is the mean, σ the SD, α the skewness, and α = 0 reduces to the Gaussian. Default links are identity (μ), log (σ), and identity (α). — [brms_families vignette](https://paulbuerkner.com/brms/articles/brms_families.html); [brmsfamily reference](https://paulbuerkner.com/brms/reference/brmsfamily.html)
- brms `student()` "can be used for robust linear regression that is less influenced by outliers". `exgaussian` is parameterized so that μ = ξ + β is the mean, giving a right-skewed Gaussian plus exponential. — [brmsfamily reference](https://paulbuerkner.com/brms/reference/brmsfamily.html); [brms_families vignette](https://paulbuerkner.com/brms/articles/brms_families.html)
- Offsets go in brms formulas as `offset(log(x))`. The reference shows `sigma ~ 1 + offset(log(xdisp))`, and for Poisson `rate(denom)` is equivalent to `offset(log(denom))`. — [brmsformula reference](https://paulbuerkner.com/brms/reference/brmsformula.html)

### Inferences
**Equivalence of the log-ratio model and lognormal-with-offset (my derivation):**
- `lognormal()` on `chl` with `offset(log(surface_chl))` in the μ formula is the same model as `gaussian()` on `log(chl/surface_chl)`. Both say log(chl) - log(surface) ~ N(η, σ).
- Their log-likelihoods differ only by the Jacobian -Σlog(chl), which is constant in the parameters. Posteriors are identical, but raw LOO/ELPD values are not comparable with the current chl_norm-scale Gaussian without a Jacobian correction.
- For a model on y' = log1p(chl_norm), add -log(1 + chl_norm_i) = -y'_i to each pointwise log-lik to put it on the chl_norm scale. This is a standard change-of-variables step. Do it whenever any of these models are compared with `loo_compare()` against the current fit.

**Back-transforming enhancement (my derivation, standard lognormal algebra):**
- For a log-scale Gaussian with mean η and SD σ:
  - Median ratio = exp(η), so median enhancement = exp(η) - 1.
  - Mean ratio = exp(η + σ²/2), so mean enhancement = exp(η + σ²/2) - 1.
- With depth- or Region-specific σ (Q3), the mean back-transform must use the matching σ. `posterior_epred()` for `lognormal` returns the mean, exp(μ + σ²/2). For a Gaussian on the transformed response it returns η, so transform the draws by hand.
- The Gaussian-peak parameters then describe the log-enhancement profile: exp(logA) is the peak's log-ratio excess, and the DCM depth z_m is unchanged by monotone transforms of the peak term.
- **Gamma(log) with offset** models the *mean* ratio directly: E[chl/surface] = exp(η). The coefficients are therefore on a mean-ratio scale, which is closer to the current chl_norm interpretation. The variance is ∝ μ², so the coefficient of variation is constant (unless shape is modelled).

**Handling zero/negative fluorescence:** the log and Gamma options require chl > 0 and surface_chl > 0. Options, from least to most work:
1. Floor chl at a small positive value, such as the sensor's dark/blank uncertainty, and run a sensitivity analysis.
2. Treat values at or below a detection limit as left-censored, using brms's documented censoring addition term: `chl | cens(cens_flag) ~ ...` with `lognormal()`, where `cens_flag = "left"` for floored values.
3. Drop non-positive values. This is fine only if they are rare and not concentrated at depth in oligotrophic regions, where they would bias the deep tail.

**Small-denominator risk:** surface_chl in the denominator is itself a noisy single-profile estimate. Very small surface_chl, as in oligotrophic Red Sea surface waters, inflates the ratio. This is a plausible cause of the Central Red Sea right tail. On the log scale it becomes an additive cast-level shift, which a cast random intercept on `b0` absorbs. Defining surface_chl as a mean over, say, 10-20 m instead of a single value would also help.

**Code sketches (not run):**

```r
# (i) Log-ratio, Gaussian or Student-t, smooth mean as currently used
f_lr <- bf(lratio ~ 1 + Region + distance_z + s(depth, by = Region, k = 7) +
             (1 | Ocean:Region:date) + (1 | Ocean:Region:enso),
           sigma ~ 0 + Region)
fit_lr <- brm(f_lr, data = d, family = student(), backend = "cmdstanr", ...)

# (ii) Same likelihood written as lognormal on chl with offset (identical posterior)
f_ln <- bf(chl ~ 1 + Region + distance_z + s(depth, by = Region, k = 7) +
             offset(log(surface_chl)) +
             (1 | Ocean:Region:date) + (1 | Ocean:Region:enso),
           sigma ~ 0 + Region)
fit_ln <- brm(f_ln, data = d, family = lognormal(), backend = "cmdstanr", ...)

# (iii) Gamma with log link: mu = mean Chl; exp(eta) = mean ratio to surface
f_ga <- bf(chl ~ 1 + Region + distance_z + s(depth, by = Region, k = 7) +
             offset(log(surface_chl)) +
             (1 | Ocean:Region:date) + (1 | Ocean:Region:enso),
           shape ~ 0 + Region)
fit_ga <- brm(f_ga, data = d, family = Gamma(link = "log"), backend = "cmdstanr", ...)

# (iv) Skew-normal on the original chl_norm scale (mu = mean; alpha by Region)
f_sn <- bf(chl_norm ~ 1 + Region + distance_z + s(depth, by = Region, k = 7) +
             (1 | Ocean:Region:date) + (1 | Ocean:Region:enso),
           sigma ~ 0 + Region, alpha ~ 0 + Region)
fit_sn <- brm(f_sn, data = d, family = skew_normal(), backend = "cmdstanr", ...)

# Left-censoring variant for floored values
# d$cens_flag <- ifelse(d$chl_raw <= dl, "left", "none"); d$chl <- pmax(d$chl_raw, dl)
# bf(chl | cens(cens_flag) ~ ... + offset(log(surface_chl)), sigma ~ ...), family = lognormal()
```

- **Ranking for testing:** (i)/(ii) first, since they are the cheapest and fix the < -1 problem by construction. Then Student-t vs Gaussian on the log scale for the Red Sea tail. Then Gamma(log) if a mean-ratio interpretation matters.
- **Skew-normal** does not fix the bound violation. It is worth testing only if staying on the chl_norm scale is a hard requirement.

### Gaps
- I found no oceanographic source that specifically compares Gaussian vs lognormal vs Gamma likelihoods for DCM enhancement ratios.
- I found no source giving a standard detection-limit floor for CTD fluorometers. It is sensor- and calibration-specific.

---

## Q3. Distributional regression: sigma as a function of depth/Region

### Takeaway
brms models sigma (log link) with any linear predictor, including smooths and random effects. `sigma ~ 0 + Region + s(depth, by = Region, k = 5)` is valid syntax. On the log-ratio scale, depth-dependent sigma is likely needed because variability concentrates around the DCM and near the base of the profile. I found no profile-specific literature quantifying heteroscedasticity.

### Cited Findings
- Distributional regression examples include `bf(y ~ group, sigma ~ group)`, where sigma coefficients are contrasts of log residual SDs and are back-transformed with `hypothesis(fit, "exp(sigma_Intercept) = 0")`. Smooths and random effects in sigma also work: `bf(y ~ s(x1) + s(x2) + (1|fac), sigma ~ s(x0) + (1|fac))` with `control = list(adapt_delta = 0.95)`. — [brms_distreg vignette](https://paulbuerkner.com/brms/articles/brms_distreg.html)
- Positive distributional parameters default to a log link. The identity link is available but "much more likely to lead to problems". Other distributional parameters (`nu`, `shape`, `alpha`, `zi`, `hu`) can be predicted the same way. — [brmsformula reference](https://paulbuerkner.com/brms/reference/brmsformula.html)
- Nonlinear and distributional parts combine with `+ lf(...)` / `+ nlf(...)`, e.g. `nlf(sigma ~ tau * sqrt(eta))` to tie sigma to the mean. — [brmsformula reference](https://paulbuerkner.com/brms/reference/brmsformula.html)
- In a palaeo time-series GAM, Simpson (2018) handled heteroscedasticity with observation weights (gamm/gam) after dropping an unidentifiable CAR(1) term. — [Simpson 2018, Front. Ecol. Evol. (bioRxiv full text)](https://www.biorxiv.org/content/10.1101/322248v2.full)

### Inferences
Code (not run):

```r
# Smooth mean + depth/Region-varying sigma on the log-ratio scale
f_ds <- bf(lratio ~ 1 + Region + distance_z + s(depth, by = Region, k = 7) +
             (1 | Ocean:Region:date) + (1 | Ocean:Region:enso),
           sigma ~ 0 + Region + s(depth, by = Region, k = 5))

# Nonlinear peak model + distributional sigma
f_peak_ds <- f_peak + lf(sigma ~ 0 + Region + s(z, by = Region, k = 5))
pri_peak_ds <- c(pri_peak[pri_peak$class != "sigma", ],
                 prior(normal(-1, 1), class = "b",   dpar = "sigma"),
                 prior(exponential(2), class = "sds", dpar = "sigma"))
```

- For Student-t, add `nu ~ 0 + Region` to let the Red Sea have heavier tails. For skew-normal, add `alpha ~ 0 + Region`.
- **Cost:** a smooth in sigma adds roughly 8 × (k-1) spline coefficients plus 8 sds parameters. That is a small cost compared with the 150k-row likelihood.
- Heteroscedasticity and residual autocorrelation interact: an AR(1) term on raw residuals with non-constant sigma can be misleading. Fit sigma first, then check standardized residuals.

### Gaps
- I did not find oceanographic or vertical-profile literature on heteroscedastic modelling of Chl profiles specifically. The evidence base is general GAM/GAMLSS practice.

---

## Q4. Within-profile autocorrelation (residual lag-1 ≈ 0.97)

### Takeaway
brms `ar(time, gr, p = 1)` with `cov = FALSE` is the cheap regression formulation for Gaussian-type families. It uses observation order within each group, so it is defensible on regular 1-m bins. Reported overheads range from about +50% runtime (panel data) to impractical runtimes when `cov = TRUE` builds full covariance matrices.

Ignoring lag-1 ≈ 0.97 makes intervals anti-conservative, because effective sample size is far below 152k. Wiggly smooths and strong autocorrelation also compete to explain the same pattern.

The most practical options, in my judgment:
1. Model profile-level structure explicitly: cast-level random effects on the nonlinear peak parameters, or cast-level factor smooths.
2. Fit on depth-binned/thinned data, e.g. 5-m bins.
3. Add `ar(p = 1)` on top only if residual autocorrelation remains.

### Cited Findings
- `ar(time = NA, gr = NA, p = 1, cov = FALSE)`: `time` sets the observation order (default: data row order), and `gr` restricts correlation to within groups. `cov = TRUE` uses residual covariance matrices and supports only stationary order-1 ARMA; it adds latent residuals if the family lacks natural residuals. `cov = FALSE`, the default, uses a regression formulation that is "much faster", allows higher orders, and is available only for Gaussian and some of its generalizations. — [brms ar() reference](https://paulbuerkner.com/brms/reference/ar.html)
- Paul Bürkner, on a panel model with `cor_ar(~ year | firm, p = 1)` that ran about 50% slower: he was "surprised it does not add more than 50%". He said the default ARMA formulation doesn't require stationarity, and that `gr` should be the finest grouping (e.g. `country:firm`). — [Stan Discourse 10118](https://discourse.mc-stan.org/t/ar-model-in-brms-how-to-speed-it-up-and-for-unbalanced-panel-structure/10118)
- A ZINB GAM with `cor_ar(..., cov = TRUE)` did not finish after about a week on a 24-core workstation. Respondents noted that the covariance approach builds a full covariance matrix and suggested approximating AR(1) with an additional smooth, or coding the model directly in Stan. — [Stan Discourse 11965](https://discourse.mc-stan.org/t/computational-time-of-a-fairly-complex-gam-with-arma-structure-in-brms/11965)
- For unevenly spaced series, a forum answer (N. Clark) says a continuous-time AR(1) is not available in brms. It is available in mvgam via `trend_model = CAR()`, where correlation decays as φ^Δt. — [Stan Discourse 10001](https://discourse.mc-stan.org/t/autocorrelation-for-unevenly-spaced-time-series/10001/8)
- Simpson (2018) makes four relevant points:
  - "A wiggly trend and an autocorrelation function with large ϕ are two ways to describe the same pattern." The two can be separated only when they act at different scales.
  - Ignoring residual dependence gives narrower intervals and coverage "lower than the anticipated 1 – α".
  - Uneven spacing needs CAR(1), which he fits with mgcv `gamm()`. In one series GAM + CAR(1) was non-identifiable (φ CI 0-1), so he dropped it and used k = 40 with REML.
  - GP smooths with Matérn/power-exponential kernels are an alternative.
  — [Simpson 2018 (bioRxiv full text)](https://www.biorxiv.org/content/10.1101/322248v2.full)
- Pedersen et al. (2019):
  - `bs = "fs"` creates a copy of the basis per group level with **one** shared smoothing parameter (model "GS"/"S").
  - Group-level smoothers increase cost roughly with the number of groups squared.
  - GS/S are much faster than GI/I with many groups.
  - `bam()` was upward of an order of magnitude faster than `gam()` beyond about 16 groups.
  - brms can fit factor-smooth interactions and `t2()` but not `te()`.
  - The paper does not treat autocorrelation as a modelling issue.
  — [Pedersen et al. 2019, PeerJ 7:e6876](https://peerj.com/articles/6876)
- `gp(..., k = NA, c = 5/4, cov = "exp_quad" | "matern32" | "matern52" | "exponential", gr = TRUE, by = NA)`: setting `k` gives a Hilbert-space approximate GP, with c = 5/4 as a suggested default. `gr = TRUE` shares latent values among identical predictor values, which "can greatly speed sampling". `by = factor` fits one GP per level. Kernel choice was added in brms 2.22.0. — [brms gp() reference (2.22)](https://www.rdocumentation.org/packages/brms/versions/2.22.0/topics/gp); [brms NEWS](https://paul-buerkner.r-universe.dev/brms/NEWS)

### Inferences
**Effective sample size (my arithmetic):**
- With AR(1)-like residuals and ρ ≈ 0.97 per 1-m step, the effective number of independent residuals per cast is roughly n(1-ρ)/(1+ρ) ≈ 141 × 0.03/1.97 ≈ 2.
- So the 152k rows carry information closer to a few thousand independent residuals. Intervals from the current fit are likely much too narrow, and smoothing parameters likely too wiggly.

**Thinning or binning:**
- Averaging into 5-m bins would cut N to about 30k and lag-1 to about 0.97^5 ≈ 0.86, if AR(1)-like. 10-m bins would cut N to about 15k with lag-1 ≈ 0.74.
- This should cut runtime roughly in proportion to N, since likelihood cost dominates.
- The DCM (peak σ-width ~10-30 m) is still resolved at 5-m spacing.
- This is the highest value-per-effort test.

**Cast-level structure is the physically motivated fix:**
- Lag-1 ≈ 0.97 mostly reflects each cast deviating smoothly from the regional mean curve. Date-level intercepts only shift whole profiles.
- Two options capture the shape deviations:
  1. Nonlinear peak model: cast-level correlated random effects on `zm`, `logA`, `logw` (Q1). That is about 3 × 1080 parameters and is interpretable.
  2. Smooth model: `s(depth, cast, bs = "fs", k = 5, m = 1)`, about 5 × 1080 coefficients with one shared penalty, per Pedersen et al. This is likely expensive in brms; no runtimes were found.
- After either, recompute residual lag-1. If still above about 0.5, add `ar(time = depth_idx, gr = cast, p = 1)` to the Gaussian/Student log-ratio model.

**`ar()` code:** it needs an integer index on a regular grid. Gaps in depth are not distance-aware: they are treated as adjacent.

```r
d <- d |> dplyr::arrange(cast, depth) |>
  dplyr::mutate(depth_idx = as.integer(round(depth)))   # 1-m bins; or depth/5 for 5-m bins
f_ar <- bf(lratio ~ 1 + Region + distance_z + s(depth, by = Region, k = 7) +
             (1 | Ocean:Region:date) + (1 | Ocean:Region:enso) +
             ar(time = depth_idx, gr = cast, p = 1),
           sigma ~ 0 + Region)
```

- **ar() and nonlinear models:** I did not verify whether `ar()` is supported in nonlinear (`nl = TRUE`) formulas.
- **Per-cast GPs:** `gp(depth, by = cast, k = ...)` would create 1080 GPs with separate hyperparameters, which is not sensible. A GP replacing the regional smooth (`gp(z, by = Region, k = 15, cov = "matern32")`) is possible, but it changes the mean model rather than the residual correlation.
- **Coverage:** none of these fixes guarantee nominal coverage. Validate by refitting on held-out casts with grouped `kfold(fit, group = "cast")`, or by simulation-based checks on a subsample.

### Gaps
- I did not find a brms-specific benchmark of `ar()` on ~10^5 rows.
- I did not find documentation on whether `ar()` combines with `threads = threading()` or with nonlinear formulas.
- I did not find runtimes for brms `bs = "fs"` smooths with ~1000 levels.
- I found no oceanographic application of brms to vertical profiles with along-depth autocorrelation.

---

## Q5. Computational tricks for large brms GAMs

### Takeaway
The biggest gains come from reducing N (depth binning/thinning) and fixing the likelihood. After that: cmdstanr with within-chain threading (4 chains × 4 threads on 16 cores), `normalize = FALSE`, informative priors and good inits, and testing on subsamples with `algorithm = "pathfinder"` before full NUTS. Threading speedups for cheap Gaussian likelihoods are modest, around 1.8x on 2-4 cores in the vignette's Poisson example and less for normal likelihoods.

### Cited Findings
- **Threading setup and speedups:**
  - Within-chain threading needs Stan ≥ 2.26: `update(fit, chains = 2, cores = 2, backend = "cmdstanr", threads = threading(2))`.
  - Use no more threads in total than physical cores; hyper-threading did not help.
  - Example speedups (Poisson, 10^4 obs, 1000 groups): about 1.8x with 2 cores and about 1.8-1.9x with 4.
- **Threading caveats:**
  - Only the likelihood is parallelized.
  - Expensive likelihoods (Poisson) parallelize better than cheap ones (normal).
  - Multilevel models with many parameters carry more overhead.
  - Rows should be randomly shuffled so chunks cost about the same.
- **Threading tuning:** grainsize can be tuned by benchmarking `threading(cores, grainsize = ceiling(N/chunks), static = TRUE)`.
- — [brms_threading vignette](https://paulbuerkner.com/brms/articles/brms_threading.html)
- brms 2.21.0 added experimental `pathfinder` and `laplace` algorithms (cmdstanr backend). 2.22.0 added `threading(force = TRUE)` to force threading activation. — [brms NEWS](https://paul-buerkner.r-universe.dev/brms/NEWS)
- `decomp = "QR"` (default NULL) "helps in fitting models with highly correlated predictors". Centering the population-level design matrix improves sampling, and it is on by default for linear and distributional parameters but off for nonlinear parameters. — [brmsformula reference](https://paulbuerkner.com/brms/reference/brmsformula.html)
- Approximate GPs via `gp(..., k = , c = 5/4)` and `gr = TRUE` auto-grouping of identical predictor values (e.g. repeated depth levels) "can greatly speed sampling". — [brms gp() reference](https://www.rdocumentation.org/packages/brms/versions/2.22.0/topics/gp)
- In mgcv, `bam()` is up to an order of magnitude faster than `gam()` with many groups. Shared-penalty factor smooths (GS/S) are much cheaper than per-level penalties (GI/I). — [Pedersen et al. 2019](https://peerj.com/articles/6876)

### Inferences
**Recommended run template (not run):**

```r
set.seed(1)
d_fit <- d[sample(nrow(d)), ]          # random row order helps reduce_sum load balance
# BUT: for ar() models keep rows sorted by cast/depth.
# Do not shuffle when ar() is used (my inference; brms sorts by time within gr,
# but verify via standata()).

fit <- brm(
  f_lr, data = d_fit, family = student(),
  backend = "cmdstanr", chains = 4, cores = 4,
  threads = threading(4, grainsize = ceiling(nrow(d_fit) / 64)),
  normalize = FALSE,                  # drop normalizing constants (faster; loo still works on log_lik)
  init = 0,                           # or init_fun as in Q1 for nl models
  iter = 1500, warmup = 750,
  control = list(adapt_delta = 0.9, max_treedepth = 11)
)

# Fast structural screening before a 20-80 h run:
fit_pf <- brm(f_lr, data = d_sub, family = student(), backend = "cmdstanr",
              algorithm = "pathfinder")   # experimental since brms 2.21
```

**Further speedups (my reasoning):**
- **Bin depths** (Q4): about 5x fewer rows. This is likely the largest speedup available.
- **Drop redundant random effects.** If cast-level effects are added, the `(1|Ocean:Region:date)` intercept is probably near-collinear with a cast intercept (same-day casts), which slows sampling.
- **Use mgcv for structure screening.** Run `bam(..., discrete = TRUE)` with `rho`/`AR.start` (AR1) for fast comparison of mean structures, k, and the AR1 effect, then fit only the chosen structures in brms. This comes from my own knowledge of mgcv's documented `bam` arguments, not verified in this session.
- **Use pathfinder draws for inits.** Pass a pathfinder fit to cmdstanr's `$sample(init = ...)` by running the brms-generated model directly. Get `stancode()`/`standata()`, compile with `cmdstan_model()`, run `$pathfinder()` and then `$sample(init = pf)`. This is a cmdstanr ≥ 0.7 feature as I recall, not verified here.

### Gaps
- I found no published runtimes for brms GAMs at ~150k rows with 8 by-smooths and crossed random effects.
- I did not verify:
  - Whether `normalize = FALSE` affects `loo()`. I believe it does not, because log_lik is computed separately in generated quantities.
  - Whether brms's OpenCL/GPU support helps this model class.
  - Whether `threads` is compatible with `ar()` terms.
