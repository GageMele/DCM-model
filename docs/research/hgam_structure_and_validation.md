# Hierarchical GAM structure, smoothness priors, and validation/comparison for hierarchical CTD profile data (brms 2.23 / mgcv 1.9)

Context assumed: response `chl_norm` (chlorophyll relative to surface), ~152k obs, 1080 casts, 8 Regions in 2 Oceans, 156 Region-date groups (`RegDate`), ENSO nearly confounded with cruise. Current models: `s(depth, by=Region)`, depth shifted by 0.5*MLD, numeric-by MLD smooths per region, `t2(depth, mld_z, by=Region)`, and random intercepts `(1|Ocean:Region:date) + (1|Ocean:Region:enso)`. Observed result: regional MLD models win on held-out casts from seen dates, the depth-shift model wins on unseen dates.

## Q1. Pedersen et al. 2019 HGAM models (G, GS, GI, S, I): definitions, syntax, when global + group smooths beat separate by-factor smooths, partial pooling for regions with few casts

### Takeaway
The current `s(depth, by=Region)` with no global smooth is Pedersen's **model I**: each region has its own smoothing penalty and nothing is shared about curve shape. **GS** (global smooth + `bs="fs"` group smooths with a shared penalty) or **GI** (global smooth + `by=` group smooths with `m=1`) partially pool shape toward a common profile. That pooling is the main protection for thin regions such as the Red Sea (58-84 casts), and only models with a global smooth can predict a group that was not in training. Pedersen et al. say the choice between global and no-global is scientific, warn against picking it by AIC alone, and recommend checking predictions on held-out data.

### Cited Findings
- **Definitions and mgcv syntax** (CO2 example; `x` = covariate, `g` = group factor) — [Pedersen et al. 2019, PeerJ 7:e6876](https://peerj.com/articles/6876):
  - **G**, one shared smooth plus a group random intercept: `s(x, k=5, bs="tp") + s(g, k=12, bs="re")`
  - **GS**, a global smooth plus group smooths with **one shared** smoothing parameter: `s(x, k=5, m=2) + s(x, g, k=5, bs="fs", m=2)`
  - **GI**, a global smooth plus group smooths with **separate** smoothing parameters, a lower-order penalty and an explicit random intercept: `s(x, k=5, m=2, bs="tp") + s(x, by=g, k=5, m=1, bs="tp") + s(g, bs="re", k=12)`
  - **S**, shared-penalty group smooths with no global: `s(x, g, k=5, bs="fs", m=2)`
  - **I**, separate-penalty group smooths with no global: `s(x, by=g, k=5, bs="tp", m=2) + s(g, bs="re", k=12)`
- "Similar wiggliness" (shared penalty) and "similar shape" (a global smooth) are separate assumptions. The two can be combined in either way — [Pedersen et al. 2019](https://peerj.com/articles/6876).
- "The choice to include a global term should be made based on scientific considerations." Models S and I "cannot be used to predict … for [groups] that are not part of the training set". GS "can be used to simulate functional variation for unobserved group levels" — [Pedersen et al. 2019](https://peerj.com/articles/6876).
- A shared penalty (GS/S) can over-smooth wiggly groups or under-smooth smooth ones when true wiggliness differs between groups. Separate penalties (GI/I) can over-smooth in some simulations. The authors suggest fitting both and comparing. With too few data per group, GS or S may still be preferable — [Pedersen et al. 2019](https://peerj.com/articles/6876).
- "There is no way for the model to share information on function shape between grouping levels without the global smoother." Model S estimates are expected to vary more than GS estimates when groups have little data — [Pedersen et al. 2019](https://peerj.com/articles/6876).
- **Concurvity:** "the global term is either entirely or almost entirely concurved with the group-level smoothers". This can inflate the uncertainty of the global smooth, or flatten it. Two mitigations are given: use `m=1` (a first-derivative penalty) for the group-level smooths and `m=2` for the global smooth, or make the group k clearly different from the global k — [Pedersen et al. 2019](https://peerj.com/articles/6876).
- **Cost:** GI and I add one smoothing parameter per group, so they are usually slower than GS and S when there are many groups, though the bird example was an exception. `bam()` "can be upward of an order of magnitude faster" when there are more than 16 groups — [Pedersen et al. 2019](https://peerj.com/articles/6876).
- **Empirical results:** On the CO2 data, AIC favoured I and S (ΔAIC: G 101, GS 22, GI 4, S 1, I 0). On the bird data, AIC slightly preferred S even though GS generated the data, which the authors blame on AIC's handling of global terms. On the Daphnia data, "Model GI is able to predict as well or better than model G or GS for all lakes". In the zooplankton out-of-sample test, S beat I — [Pedersen et al. 2019](https://peerj.com/articles/6876).
- "We strongly recommend not selecting models based purely on AIC." GLRTs (generalized likelihood ratio tests) are not recommended. For predictive accuracy, use held-out data or k-fold CV, and consider model averaging — [Pedersen et al. 2019](https://peerj.com/articles/6876).
- brms can fit all of these models through `fs` and `t2` terms, but "the brms package does not currently support te() tensor products" — [Pedersen et al. 2019](https://peerj.com/articles/6876).
- In brms, an `fs` smooth is a fully penalized spline per level. It bundles a random intercept, a random linear slope and a random wiggly part, so the fit reports **three sds variance terms** for the term. All levels share one smoothness parameter set — [Stan Discourse, ucfagls / jd_c](https://discourse.mc-stan.org/t/how-does-one-interpret-the-output-of-using-a-s-term-in-brms/22426).

### Inferences
- In Pedersen's terms, the current regional model is **I**. Its region curves are partially pooled only through the region intercept, and curve shape is not pooled at all. With 58-84 casts per Red Sea region and roughly 100-200 obs per cast, there is plenty of within-region data for the *mean* curve. Shape pooling matters more for the *date-to-date* variation and for MLD interactions, where the effective sample is the number of dates per region, not the number of observations.
- **Test GI next to I.** GI keeps separate per-region penalties and adds a global depth smooth, with `m=1` on the regional deviations. If the regional deviations are small, the penalty shrinks them toward the global curve. This is principled partial pooling of shape.
- **Test GS as well.** With 8 regions, one shared smoothing parameter for the regional deviations is estimable, and it is cheaper than GI.
- The unseen-date result is consistent with this theory. More region-specific flexibility (the regional MLD model) helps interpolation within known dates but generalizes less well. GI/GS-style shrinkage of MLD effects toward a global MLD effect is the obvious middle ground to test (see Q2).
- Runnable mgcv prescreen with consistent bases. `depth_s` is the shifted depth, and `RegDate` and `Region` are factors.
```r
library(mgcv)
ctrl <- list(nthreads = 8)
m_I  <- bam(chl_norm ~ Region + s(depth_s, by = Region, k = 12) + s(RegDate, bs = "re"),
            data = d, method = "fREML", discrete = TRUE, control = ctrl)
m_GI <- bam(chl_norm ~ Region + s(depth_s, k = 12) +
              s(depth_s, by = Region, k = 12, m = 1) + s(RegDate, bs = "re"),
            data = d, method = "fREML", discrete = TRUE, control = ctrl)
m_GS <- bam(chl_norm ~ s(depth_s, k = 12) +
              s(depth_s, Region, bs = "fs", k = 12, m = 1) + s(RegDate, bs = "re"),
            data = d, method = "fREML", discrete = TRUE, control = ctrl)
# bs="fs" already contains region-specific intercepts, so no Region fixed effect in GS
k.check(m_GI); concurvity(m_GI, full = TRUE); summary(m_GI)
```
- The brms equivalents (`s()` in brms is built with mgcv's `smoothCon`, so `m=` and `bs=` carry over):
```r
library(brms)
f_GI <- bf(chl_norm ~ Region + s(depth_s, k = 12) + s(depth_s, by = Region, k = 12, m = 1) +
             (1 | RegDate))
f_GS <- bf(chl_norm ~ s(depth_s, k = 12) + s(depth_s, Region, bs = "fs", k = 12, m = 1) +
             (1 | RegDate))
```

### Gaps
- Pedersen et al. give no rule for "how few casts per group" makes GS better than GI. Their guidance is qualitative: compare the two, and prefer the shared penalty when groups are data-poor.
- I found no published HGAM comparison specifically for oceanographic profiles that are nested by cast and date.

## Q2. Ocean-level nesting (global + Ocean + Region deviation) and how to add MLD without regional tensors

### Takeaway
Nesting follows directly from the HGAM framework: global depth smooth + Ocean-level deviation + Region-level deviation (`m=1`). Ocean has only 2 levels, so an Ocean "random" smooth (`bs="fs"`) cannot estimate its own variance well. Treat Ocean as a fixed `by=` effect and keep partial pooling for Region within Ocean. For MLD, a single global `t2(depth, mld)` interaction plus regional depth-only deviations is far cheaper and better constrained than `t2(..., by=Region)`.

### Cited Findings
- For multi-dimensional group smooths, Pedersen et al. use a global `te(...)` plus a group term `t2(..., full=TRUE)` with `bs="re"` on the group margin. In brms this has to be written with `t2` throughout, because brms lacks `te()` — [Pedersen et al. 2019](https://peerj.com/articles/6876).
- mgcv drops some basis functions from tensor-product group terms, so the range of group deviations can differ between groups. This is a bigger problem "when few basis functions or few groups are used" — [Pedersen et al. 2019](https://peerj.com/articles/6876).
- Gelman's guidance: a uniform prior on a group-level SD is a reasonable default "unless the number of groups J is low (below 5, say)". With J = 1 or 2 the uniform prior gives an improper posterior, and with J = 3 it gives a very long right tail. Use a weakly informative half-t or half-Cauchy prior with a sensible scale instead — [Gelman 2006, Bayesian Analysis](https://econwpa.ub.uni-muenchen.de/econ-wp/em/papers/0404/0404001.pdf).

### Inferences
- With J = 2 oceans, any random-effect variance for Ocean (`(1|Ocean)` or `s(depth, Ocean, bs="fs")`) is identified almost entirely by its prior. Model Ocean as a **fixed** `by=` smooth: `Ocean + s(depth_s, by = Ocean, k = 12)`. Then pool regions within Ocean with a GS-type `fs` term (8 levels) or a GI-type `by=Region, m=1` term.
- Suggested nested formula, the HGAM analogue of global, then Ocean, then Region:
```r
f_nest <- bf(chl_norm ~ Ocean + s(depth_s, by = Ocean, k = 12) +                  # ocean-level curves
               s(depth_s, Region, bs = "fs", k = 10, m = 1) +                    # pooled region deviations
               (1 | RegDate))
```
  Two versions are equivalent here: an Ocean `by=` smooth with no separate global smooth, or a global smooth plus an `Ocean` deviation `s(depth_s, by = Ocean, m = 1)`. The second is more concurved, so prefer the first.
- **MLD with partial pooling instead of 8 regional tensors.** Use one global interaction surface plus regional depth deviations. Optionally add a pooled regional MLD slope-by-depth term:
```r
f_mld_pool <- bf(chl_norm ~ Ocean + s(depth_s, by = Ocean, k = 12) +
                   t2(depth_s, mld_z, k = c(10, 5)) +               # one shared depth x MLD surface
                   s(depth_s, Region, bs = "fs", k = 10, m = 1) +
                   (1 + mld_z | Region) + (1 | RegDate))
# Alternative "varying-coefficient" form: global linear-in-MLD effect at each depth
#   s(depth_s, by = mld_z, k = 10) + s(depth_s, by = interaction(Region, ...)) ... only if needed
```
  A global `s(depth_s, by = mld_z)` is the pooled counterpart of the current "numeric-by MLD per region" model. Adding regional deviations from it with `m = 1` gives a GI-style test of whether MLD effects really differ by region.
- Consider making the depth-shift itself hierarchical rather than fixed at 0.5. That requires a nonlinear brms formula (`nl = TRUE`) with a shift coefficient, global or per Ocean, and is harder to sample. Treat it as an optional later test.

### Gaps
- I found no source that directly compares fixed-by-Ocean against random Ocean smooths for J = 2. The recommendation follows from Gelman's few-groups result and is not an HGAM-specific result.

## Q3. Random-effect structure when ENSO is confounded with cruise/date; few-level groups; sd priors

### Takeaway
`(1|Ocean:Region:enso)` sits above `(1|Ocean:Region:date)`, and ENSO is nearly a function of cruise. The two variances are therefore weakly identified, and the posterior split between them is prior-driven. For prediction to unseen dates, the date-level random intercept is integrated out as new-level noise. What matters is that the total between-date variance is estimated well, not how it is divided between date and ENSO. Use one date-level random effect, and add ENSO only as a fixed contrast (if at all), with an informative half-normal or exponential prior on the SDs.

### Cited Findings
- Few groups: with J < ~5, use a weakly informative half-t/half-Cauchy on the SD, not a flat or inverse-gamma(ε, ε) prior. Gelman shows inverse-gamma(ε, ε) is not noninformative and its results depend on ε. In regression settings "more than 3 groups" may be needed to estimate σ_α from a noninformative prior — [Gelman 2006](https://econwpa.ub.uni-muenchen.de/econ-wp/em/papers/0404/0404001.pdf).
- brms default: group-level SDs have "a half student-t prior with 3 degrees of freedom" whose scale "depends on the standard deviation of the response after applying the link function", with a minimum of 2.5. This is meant to be weakly informative while giving "at least some regularization". Syntax: `set_prior("<prior>", class = "sd", group = "<group>", coef = "<coef>")` — [brms set_prior docs](https://rdrr.io/cran/brms/man/set_prior.html).
- Penalised-complexity (PC) priors shrink toward a simpler "base model", such as zero random-effect variance. For a Gaussian random-effect SD they reduce to an exponential prior on σ, scaled by a statement P(σ > U) = α — [Simpson, Rue, Riebler, Martins & Sørbye 2017, Statistical Science 32(1)](https://projecteuclid.org:443/euclid.ss/1491465621). *I could not retrieve the full text in this session. The exponential-on-SD form is the standard PC-prior result from that paper, but the reader should check the exact wording.*
- Random CV and mixed models still underestimate predictive error when the target is new groups. Mixed models can give unbiased parameters, but that does not fix the optimism of random splits — [Roberts et al. 2017, Ecography 40:913-929](https://www.biometrie.uni-freiburg.de/mitarbeiter/dormann/roberts-et-al-2017-ecography.pdf).

### Inferences
- **Options to compare**, all with `(1 | RegDate)` as the base:
  1. `(1 | RegDate)` alone (ENSO absorbed into date variance). This is the simplest and is honest about confounding.
  2. `enso_f:Ocean` as a **fixed** contrast + `(1 | RegDate)`. The ENSO effect is then identified only from between-cruise contrasts, so read it as "cruise-cluster difference" and not causal ENSO.
  3. The current `(1|Region:enso) + (1|Region:date)`. Keep it only if the posterior for `sd(Region:enso)` moves clearly away from its prior. Compare `sample_prior = "yes"` draws with posterior draws.
- For unseen-date prediction, `RegDate` is a new level. brms draws new-level effects with `sample_new_levels = "gaussian"` by default (see Q5). The predictive interval then scales with `sd(RegDate)`, so this SD drives the elpd on date-blocked folds.
- Suggested priors. chl_norm is O(1), so check the scale against `sd(d$chl_norm)`:
```r
pri <- c(
  prior(normal(0, 1), class = b),                     # fixed effects incl. linear null-space parts of smooths
  prior(exponential(2), class = sd),                  # PC-style: P(sd > 1.5) ~ 0.05
  prior(exponential(2), class = sds),                 # smooth wiggliness (see Q4)
  prior(exponential(1), class = sigma)
)
# few-level group (e.g. if Ocean or ENSO kept as random): tighter, explicit
# prior(normal(0, 0.5), class = sd, group = "Ocean:Region:enso")
```
- Check identifiability directly: `pairs(fit, variable = c("sd_Ocean:Region:date__Intercept", "sd_Ocean:Region:enso__Intercept"))`. A strong negative ridge means the two variances trade off and one should be dropped.

### Gaps
- I found no source that addresses exactly the case of a grouping factor nested and confounded with the date grouping in brms. The recommendation is derived from general identifiability arguments.

## Q4. Smoothness priors in brms (`sds`), PC priors, whether k matters, mgcv k checks, and prescreening with mgcv/bam

### Takeaway
brms fits each smooth as a mixed model. The null-space part (the linear term, named e.g. `sdepth_s_1`) is a fixed effect in class `b` with a flat default prior. The wiggly part is a random effect whose SD (`sds`) plays the role of the inverse smoothing parameter, so a prior on `sds` is a smoothness prior. Once a term is penalized, k only needs to be "large enough", and mgcv's `k.check`/`gam.check` gives a cheap way to confirm that. Prescreening structures and k in mgcv `bam()` is well justified: brms uses mgcv to build its bases, and REML is an empirical-Bayes version of the same model. I found no formal statement that this is "standard practice", though.

### Cited Findings
- "Smoothing splines are implemented in brms using the 'random effects' formulation". The SD class "is called `sds`", is set with `set_prior("<prior>", class = "sds", coef = "<term label>")`, and its "default prior is the same as for standard deviations of group-level effects" — [brms set_prior](https://rdrr.io/cran/brms/man/set_prior.html).
- "The functions in the penalty null space get treated as fixed effects terms in the model". "The variance of the random effect is inversely proportional to the smoothing parameter". "Smaller variances = larger smoothing parameters = less wiggly smooth functions". Priors go on the SDs, not on λ — [Stan Discourse, ucfagls (G. Simpson)](https://discourse.mc-stan.org/t/better-priors-non-flat-for-gams-brms/23012).
- The default `get_prior()` for `s(wt)` lists a flat prior on class `b` (including `swt_1`) and `student_t(3, 0, 2.5)` on `sds`. The suggested workflow is a Gaussian or t prior on the linear term, a separate prior on `sds`, then `sample_prior = "only"` with `conditional_smooths()` to check the implied wiggliness — [Stan Discourse](https://discourse.mc-stan.org/t/better-priors-non-flat-for-gams-brms/23012).
- "sds(...) is the sd of the smooth weights (spline coefficients)". It "determines the amount of 'wiggliness'". Its numeric value "is not very practically interpretable". A value near 0 "is like 'complete-pooling' of the basis functions" — [Stan Discourse, jd_c](https://discourse.mc-stan.org/t/how-does-one-interpret-the-output-of-using-a-s-term-in-brms/22426).
- k: "exact choice of k is not generally critical". It should be large enough to represent the truth and small enough to compute. When the EDF is well below k-1, checking is less important. Check by smoothing residuals with a larger k, by refitting with a larger k, or with `gam.check`'s simulation-based test — [mgcv choose.k](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/choose.k.html).
- Pedersen et al.: when the EDF is well below k, increasing k has little effect. Use `gam.check` / the k-index, and use REML/ML rather than GCV ("We strongly recommend using either REML or marginal likelihood (ML)"). The penalty is a prior precision matrix, REML is empirical Bayes, and mgcv intervals are Bayesian credible intervals — [Pedersen et al. 2019](https://peerj.com/articles/6876).
- PC priors (exponential on SD toward a zero-variance base model) are the principled default for "shrink toward simpler" components — [Simpson et al. 2017](https://projecteuclid.org:443/euclid.ss/1491465621). *Full text not re-checked here; see the caveat in Q3.*

### Inferences
- **Effect of sds priors in this setting.** With ~152k observations, the likelihood dominates the global smooths and the prior on their `sds` barely matters. The prior matters for terms identified by few units: the regional deviation smooths, per-region MLD terms, Ocean terms, and anything with a small effective group count. The default `student_t(3,0,2.5)` with a minimum scale of 2.5 is wide relative to an O(1) chl_norm. An `exponential(2)` to `exponential(5)` prior, or `normal(0, 0.5)` half-normal, on the deviation-smooth `sds` gives PC-style shrinkage toward "no regional deviation". It also tends to help convergence of the `t2(..., by=Region)` terms.
- Target priors at specific terms:
```r
get_prior(f_GI, data = d)   # list exact coef names
pri_s <- c(
  prior(exponential(1), class = sds, coef = "s(depth_s, k = 12)"),                 # global: weak
  prior(exponential(4), class = sds),                                               # deviations: shrink
  prior(normal(0, 1), class = b)
)
fit_prior <- brm(f_GI, data = d, prior = pri_s, sample_prior = "only",
                 chains = 2, iter = 1000, backend = "cmdstanr")
plot(conditional_smooths(fit_prior), ask = FALSE)  # prior predictive check of wiggliness
```
  The coef label must match `get_prior()` output exactly. The string above is illustrative.
- **k once penalized.** Raising k from 7 to 15 should change little *if* the EDF stays well below k-1. If the EDF approaches k-1 at k = 7, the earlier comparisons over k may have been driven by basis truncation and not by smoothness. Check this in mgcv:
```r
k.check(m_GI)         # k', edf, k-index, p-value per smooth
m_GI_k20 <- update(m_GI, . ~ . - s(depth_s, by = Region, k = 12, m = 1) +
                     s(depth_s, by = Region, k = 20, m = 1))
AIC(m_GI, m_GI_k20)   # or compare EDF / fitted curves
```
- **Prescreening in mgcv first.** This is defensible for three reasons: (1) brms smooths are built from mgcv bases, so the structures are identical; (2) REML smoothing parameters are an empirical-Bayes estimate of the same model (Pedersen); (3) `bam(discrete = TRUE)` fits a 150k-row HGAM in seconds to minutes, against hours in Stan. Use mgcv to (a) choose k, (b) check concurvity between global and group smooths, (c) run cheap date-blocked CV over many candidate structures, and then (d) refit only the 2-4 finalists in brms for full uncertainty, elpd and stacking.
- Date-blocked CV in mgcv. A new `RegDate` level is zeroed out with `exclude=`, which matches the marginal prediction):
```r
set.seed(1)
dates  <- unique(d$RegDate)
fold_of_date <- sample(rep(1:10, length.out = length(dates)))
d$fold <- fold_of_date[match(d$RegDate, dates)]
cv_bam <- function(form, d) {
  out <- vector("list", 10)
  for (k in 1:10) {
    tr <- droplevels(d[d$fold != k, ]); te <- d[d$fold == k, ]
    m  <- bam(form, data = tr, method = "fREML", discrete = TRUE)
    te$RegDate_orig <- te$RegDate
    te$RegDate <- factor(levels(tr$RegDate)[1], levels = levels(tr$RegDate))  # placeholder level, zeroed by exclude
    te$pred <- as.numeric(predict(m, te, exclude = "s(RegDate)"))
    out[[k]] <- te
  }
  r <- do.call(rbind, out)
  per_date <- tapply((r$chl_norm - r$pred)^2, r$RegDate_orig, mean)
  list(rmse = sqrt(mean((r$chl_norm - r$pred)^2)), per_date_mse = per_date)
}
```
  Use the same `d$fold` vector for every candidate so the comparisons are paired.

### Gaps
- I found no authoritative statement that "prescreen in mgcv, then fit in brms" is standard practice. It is a pragmatic workflow supported by the shared basis construction and the REML/empirical-Bayes equivalence.
- I found no published guidance on calibrated `sds` prior scales for spline deviations in brms. The thread notes that sds values are not directly interpretable, so prior predictive simulation is the recommended route.

## Q5. Model comparison for hierarchical, autocorrelated profile data: PSIS-LOO vs grouped K-fold / LOGO / LFO / block CV; elpd vs RMSE; marginal vs conditional prediction

### Takeaway
Pointwise LOO is not *invalid* under dependence. It answers the wrong question: it estimates prediction of a single missing depth bin inside a cast that is otherwise observed, which is near-interpolation along an autocorrelated profile. The CV split must match the prediction task. Held-out casts from seen dates correspond to **conditional** prediction (known date effect). Held-out dates correspond to **marginal / leave-one-cluster-out** prediction (date effect integrated out). Theory says these can rank models differently, and random or within-structure splits favour overly complex models. The observed result is therefore what theory predicts. Pick the target (most likely unseen cruises/dates), use grouped K-fold by `RegDate` in brms with `joint="group"`, and report both elpd and RMSE.

### Cited Findings
- "LOO is valid if the focus is in the conditional data (observation) model". "Leave-one-group-out (LOGO) cross-validation usually doesn't work with PSIS approach". "LOO and cross-validation in general do not require independence or conditional independence". "It depends which level of the model is interesting" — [Vehtari, Cross-validation FAQ](https://users.aalto.fi/~ave/CV-FAQ.html).
- "It is often natural to match the actual prediction task and how the data is partitioned in cross-validation" — [CV FAQ §3.1](https://users.aalto.fi/~ave/CV-FAQ.html).
- Random splits leave dependent data on both sides of the split, so error estimates become "too optimistic, and model selection will favour too complex models". Block by group or individual when predicting to new groups. Blocks should be at least as large as the residual autocorrelation range, and the spatial simulation showed they needed to be "substantially larger". Blocked CV and buffered LOO did well for both selection and error estimation. In the phylogenetic simulation, random CV and unbuffered LOO "more often chose overly complex models". Blocking can force extrapolation and be conservative, but the authors judge that cost small next to the overconfidence of random CV. The usual 5-10 fold advice is "perhaps more appropriate for random splitting" — [Roberts et al. 2017, Ecography](https://www.biometrie.uni-freiburg.de/mitarbeiter/dormann/roberts-et-al-2017-ecography.pdf).
- Conditional vs marginal: whether the model should predict clusters already in the data or new clusters. "Marginal WAIC corresponds to leave-one-cluster-out (LOcO) cross-validation; conditional WAIC corresponds to leave-one-unit-out" — [Merkle, Furr & Rabe-Hesketh 2019, Psychometrika 84:802-829 / arXiv:1802.04452](https://arxiv.org/abs/1802.04452v2).
- `kfold.brmsfit`, from the [brms reference](https://paulbuerkner.com/brms/reference/kfold.brmsfit.html):
  - `folds = "grouped"` with `group=` uses `loo::kfold_split_grouped`. `folds = NULL` with `group=` gives leave-one-group-out, and K is ignored. A numeric vector assigns folds by hand.
  - `joint = FALSE/"obs"` gives a pointwise elpd that can be compared with `loo`. `joint = "fold"` and `joint = "group"` use joint log-likelihoods per fold or per group.
  - `Ksub` evaluates a subset of folds. `save_fits = TRUE` keeps the refits.
  - Held-out new levels are drawn with `sample_new_levels = "gaussian"` by default.
  - `future_args` and `future::plan()` parallelize the refits.
- `joint = "group"` sums the log-likelihoods within a group before `log_mean_exp`. Vehtari: `joint=FALSE` fits predicting new individuals, while `joint="group"` fits predicting a new group "accounting for within-group dependence". It is at least as sensitive as `joint=FALSE`, and a difference between the two signals strong dependence in the joint predictive. Bürkner confirmed the implementation — [Stan Discourse](https://discourse.mc-stan.org/t/difference-between-kfold-joint-group-and-kfold-joint-false-for-brms-model/39081).
- LFO-CV: LOO "uses future observations to predict past ones", which makes it optimistic for forecasting. LFO respects time order and has a PSIS approximation with diagnostics — [Bürkner, Gabry & Vehtari 2020, JSCS / arXiv:1902.06281](https://arxiv.org/abs/1902.06281).
- Scores: "We recommend log predictive density (log score) for model comparison in general … as it measures the goodness of the whole predictive distribution". Also use "application specific utility and loss functions". "It is possible that one model is better than others, but still not useful for practice" — [CV FAQ §11](https://users.aalto.fi/~ave/CV-FAQ.html).

### Inferences
- **The two cohorts answer different questions.** The seen-date cohort is a *conditional* test: the date intercept is learned from sibling casts. The unseen-date cohort is a *marginal / LOcO* test. Under Merkle et al. and Roberts et al., more region-specific flexibility is expected to win the conditional test and to lose or tie the marginal one. If the scientific use is describing the DCM response to MLD for new cruises or climatology, the **unseen-date (marginal) result should be primary**.
- **Choose the block unit by residual dependence.** Profiles are autocorrelated in depth within a cast, casts within a cruise share conditions, and ENSO is tied to cruise. Following Roberts et al., the blocks should be at least as large as the dependence. That points to **whole `RegDate` groups**, or even whole cruises across regions if one cruise spans several regions. Cast-level blocks are fine for a conditional estimand only.
- **LFO** applies only if the use case is forecasting later cruises. Because ENSO is confounded with cruise, LFO would mostly test extrapolation to a new ENSO state. Grouped K-fold over dates is the better primary design. LFO, or a leave-one-ENSO-phase/period-out split, can be an extrapolation stress test.
- **brms grouped K-fold** (K = 10 over 156 dates). Use the same folds for all models:
```r
library(loo)
set.seed(2024)
folds_date <- loo::kfold_split_grouped(K = 10, x = d$RegDate)   # identical for every model
kf_shift <- kfold(fit_shift, folds = folds_date, group = "RegDate", joint = "group",
                  save_fits = TRUE, chains = 2, cores = 2,
                  future_args = list(seed = TRUE))
kf_mld   <- kfold(fit_mld,   folds = folds_date, group = "RegDate", joint = "group",
                  save_fits = TRUE, chains = 2, cores = 2)
loo_compare(kf_shift, kf_mld)
# observation-level elpd too (comparable with loo), for diagnostics:
kf_shift_obs <- kfold(fit_shift, folds = folds_date, joint = "obs", chains = 2)
# cheaper: evaluate a subset of folds first
kf_quick <- kfold(fit_mld, folds = folds_date, Ksub = 3, chains = 1)
```
  With 152k rows, each refit costs as much as the full fit. Use `Ksub`, `future::plan(multisession)`, `chains = 2`, and consider thinning the profiles (for example every 2 m, or 1 m bins) for CV only. Binning also reduces within-cast autocorrelation.
- **RMSE from the saved fits.** This gives the application loss alongside elpd:
```r
rmse_kfold <- function(kf, d) {
  out <- lapply(seq_along(kf$fits[, "fit"]), function(i) {
    idx <- kf$fits[, "omitted"][[i]]
    pr  <- posterior_epred(kf$fits[, "fit"][[i]], newdata = d[idx, ],
                           allow_new_levels = TRUE, sample_new_levels = "gaussian")
    data.frame(i = idx, pred = colMeans(pr))
  })
  r <- do.call(rbind, out); r <- r[order(r$i), ]
  r$err <- d$chl_norm[r$i] - r$pred; r$RegDate <- d$RegDate[r$i]; r
}
```
  For the marginal estimand, also compute RMSE with `re_formula = NA`. This is the "population" curve. `sample_new_levels = "gaussian"` matters for the predictive distribution (elpd) but averages out of the posterior-mean prediction.

### Gaps
- No source gives an autocorrelation-range rule for vertical profiles specifically. Roberts et al. deal with spatial, temporal, phylogenetic and individual blocking.
- I did not verify whether the brms `kfold` output with `joint = "group"` returns one pointwise row per group. Check `dim(kf$pointwise)` before using it in stacking.

## Q6. Stacking / model averaging: is stacking the half-shift and regional-MLD models a reasonable final model?

### Takeaway
Yes, if the weights are estimated with the CV design that matches the target (date-blocked) and the weights are not degenerate. Stacking is built for the M-open case where no candidate is true. It optimizes the combined predictive density, and it handles near-duplicate models gracefully. Report the weights as predictive-combination weights, not model probabilities. A single pooled model (GI/GS with a global MLD surface) is the structural alternative, and it should be one of the stacking candidates or a competitor to the stack.

### Cited Findings
- "Bayesian model averaging is flawed in the M-open setting". Stacking "generalize[s] to the combination of predictive distributions" and works with any proper scoring rule. The authors "recommend stacking of predictive distributions", with Pseudo-BMA+ (Bayesian bootstrap) as a cheaper approximation — [Yao, Vehtari, Simpson & Gelman 2018, Bayesian Analysis](https://projecteuclid.org/euclid.ba/1516093227); [loo stacking vignette](https://mc-stan.org/loo/articles/loo2-weights.html).
- `loo_model_weights()` takes loo objects or pointwise log-likelihood matrices. `stacking_weights()` and `pseudobma_weights()` take a **pointwise elpd matrix** with one column per model. Pseudo-BMA weights get diluted by near-duplicate models, while stacking shares weight between similar models and leaves distinct models' weights intact — [loo vignette](https://mc-stan.org/loo/articles/loo2-weights.html).
- "If one of the models in the model selection is not clearly the best, it may be better to average over many models". "Bayesian stacking has better performance than LOO weights and LOO-BB weights" — [CV FAQ §6](https://users.aalto.fi/~ave/CV-FAQ.html).
- Pedersen et al. also note that averaging over multiple HGAMs may help predictive accuracy — [Pedersen et al. 2019](https://peerj.com/articles/6876).

### Inferences
- Because `stacking_weights()` takes any pointwise elpd matrix, weights can (and should) come from **date-blocked K-fold** elpd, not from PSIS-LOO:
```r
lpd <- cbind(shift = kf_shift$pointwise[, "elpd_kfold"],
             mld   = kf_mld$pointwise[,   "elpd_kfold"],
             gi    = kf_gi$pointwise[,    "elpd_kfold"])
w <- loo::stacking_weights(lpd); w
# predictions from the stack (posterior draws mixed by weight):
pp <- pp_average(fit_shift, fit_mld, fit_gi, weights = as.numeric(w),
                 newdata = nd, allow_new_levels = TRUE, summary = TRUE)
```
  Use group-level (joint) elpd rows if available, so the 152k correlated rows do not over-weight within-cast detail.
- Interpretation: weights near (0.5, 0.5) mean the models carry complementary predictive information. A plausible pattern is regional MLD helping on dates whose MLD sits inside the regional training range and the shift model helping elsewhere. Weights near (1, 0) mean pick the winner. Weights from cast-blocked CV will lean toward the regional-MLD model and should not be used for the unseen-date target.
- A stack is harder to interpret scientifically, because there is no single MLD-effect estimate. If interpretation matters, the cleaner theoretical answer is a hierarchical model that *contains* both ideas: a depth shift plus a pooled (GI/GS) regional MLD deviation. Let shrinkage decide how much regional MLD structure is supported.

### Gaps
- I found no worked example of stacking with grouped K-fold elpd in brms. The mechanics (a matrix input to `stacking_weights`) are documented, but group-level stacking is my extrapolation.

## Q7. Simplicity, projection and paired comparisons: when differences are within noise, prefer the simpler model? How to test fairly?

### Takeaway
Yes. Vehtari's guidance is to choose a smaller model with similar predictive performance, and to read |Δelpd| < 4 as small. With dependent data the normal-approximation SE is unreliable, so build the uncertainty from **paired differences aggregated to the independent unit** (date, or cruise), with a bootstrap over dates. Projection predictive selection (projpred) is the formal "simplest model that matches the reference" method, but its support for brms GAMs with many smooths and random effects is limited at this scale.

### Cited Findings
- "If elpd difference … is less than 4, the difference is small". If it is larger, compare it with the SE. "SE tends to be underestimated especially if the number of observations is small or the models are badly misspecified". "The whole normal approximation tends to fail if the models are very similar" — [CV FAQ §12, §16](https://users.aalto.fi/~ave/CV-FAQ.html).
- "Choose a smaller model with similar predictive performance as the biggest model". More elaborate approaches, "such as projection predictive variable selection, are recommended". "If cross-validation prefers the simpler model, this is not necessarily evidence that the effect is exactly zero". "In general there is no need to do any model selection" (averaging or encompassing models are alternatives) — [CV FAQ §1.3, §5](https://users.aalto.fi/~ave/CV-FAQ.html).
- Pedersen et al.: do not select purely on AIC. Global terms are penalized unfairly because group smooths can reproduce them — [Pedersen et al. 2019](https://peerj.com/articles/6876).
- Random-split and LOO-type CV more often pick overly complex models when data are structured — [Roberts et al. 2017](https://www.biometrie.uni-freiburg.de/mitarbeiter/dormann/roberts-et-al-2017-ecography.pdf).

### Inferences
- **Paired, group-level comparison.** Per-observation SEs treat 152k correlated rows as independent and overstate certainty. Aggregate to dates (156) or casts (1080) first:
```r
# elpd pointwise (obs-level) from identical folds
dd <- data.frame(RegDate = d$RegDate,
                 diff = kf_mld_obs$pointwise[, "elpd_kfold"] - kf_shift_obs$pointwise[, "elpd_kfold"])
g  <- tapply(dd$diff, dd$RegDate, sum)           # per-date elpd difference
c(total = sum(g), se = sd(g) * sqrt(length(g)), frac_dates_mld_better = mean(g > 0))
# cluster bootstrap over dates (also for RMSE differences)
set.seed(1)
boot <- replicate(4000, sum(sample(g, replace = TRUE)))
quantile(boot, c(.025, .5, .975))
# same for RMSE: per-date MSE difference
r1 <- rmse_kfold(kf_shift, d); r2 <- rmse_kfold(kf_mld, d)
gm <- tapply(r2$err^2 - r1$err^2, r1$RegDate, mean)
quantile(replicate(4000, mean(sample(gm, replace = TRUE))), c(.025, .975))
```
  If one cruise spans several regions, bootstrap over **cruise** rather than `RegDate`, so that cruise-level (ENSO) dependence is respected.
- **Break the comparison down by Region and Ocean** with `tapply(g, region_of_date, sum)`. Theory (Q1) predicts that regional MLD terms help most where a region has many dates spanning a wide MLD range, and hurt where a region has few dates (Red Sea).
- **Decision rule grounded in the sources:**
  1. Use date-blocked elpd (`joint = "group"`) as the primary criterion and RMSE as the application check.
  2. If |Δelpd| < 4, or the date-bootstrap interval spans 0, prefer the simpler or more pooled model (depth shift, or GI/GS) rather than the regional-MLD model.
  3. If the models are complementary (weights well away from 0/1), either stack or build the encompassing hierarchical model.
  4. Do not use PSIS-LOO or cast-blocked results to choose between structures meant for new dates.
- **Recommended next test set**, in theory order:
  - (a) I vs GI vs GS for depth, with and without the 0.5·MLD shift.
  - (b) Nested Ocean-`by` + Region-`fs` model.
  - (c) A pooled MLD model: global `t2(depth, mld)` or `s(depth, by = mld)` plus `m=1` regional deviations, as the theoretically preferred replacement for the per-region `t2` terms.
  - (d) Random-effect variants: `(1|RegDate)` only vs adding fixed ENSO vs the current ENSO random effect.
  - Prescreen all of these in `bam` with date-blocked folds. Refit the top 2-4 in brms with PC-style `exponential` priors on the deviation `sds` and `sd`. Compare by grouped K-fold elpd and RMSE, paired by date. Stack only if no model dominates.

### Gaps
- I did not research whether projpred currently supports brms models with `s()`/`t2()` smooths plus several random-effect terms at 150k rows. Its feasibility for this workflow is unverified.
- I found no source giving a threshold for an "important" RMSE difference in chlorophyll-profile modelling. That should be set from the application, for example DCM depth or magnitude error tolerances.
