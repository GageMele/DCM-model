# Priors, sampling, and computation

## Original user-defined configuration

~~~r
priors <- c(
  prior(student_t(3, 0.73, 5.448), class = "Intercept"),
  prior(normal(0, 2), class = "b"),
  prior(student_t(3, 0, 1.5), class = "sd"),
  prior(student_t(3, 0, 1.5), class = "sds"),
  prior(student_t(3, 0, 2.5), class = "sigma")
)
fit <- brm(
  andrew_form, data = CTD, prior = priors,
  cores = 4, chains = 4, iter = 5000, warmup = 1000,
  threads = threading(2), refresh = 50,
  backend = "cmdstanr", save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.97, max_treedepth = 14)
)
~~~

5.448 is exactly 2 × 2.724. Scale parameters have positive support. The brms Intercept prior refers to its centered parameterization; inspect the design before interpreting it as an ecological baseline.

The original Windows prefix was C:/Users/gagem/OneDrive/Desktop/brms_fits/my_model_50_FINAL. The filename does not establish a 50% sample; the saved fit is authoritative.

## Prior families

| Name/idea | Change |
|---|---|
| Current/custom raw | Original mean, group/smooth, and residual priors above |
| Default brms | Omit priors for the relevant original specification; defaults depend on version/data |
| Empirical | Response-derived location/scale from full prepared data with cast weighting |
| Standardized custom | Numerical coefficient scales applied to standardized predictors |
| Tighter extensive | b N(0,1), sd/sds positive N(0,0.75), recorded residual specification |
| Half-normal scales | Positive normal instead of Student-t scale priors; values depend on campaign |
| b_tight in shape | Mean b N(0,1); retain other main priors |
| scales_halfnormal in shape | Mean sd/sds positive N(0,1.5); retain b N(0,2) |
| mu_default | Mean defaults while retaining proper variance priors; not all-default |

The early empirical calculation uses cast-weighted response median q50 and s = max(0.5, 1.4826 × weighted MAD, IQR/1.349). Intercept is Student-t(3,q50,2s), b N(0,s), and positive sd/sds/sigma Student-t(3,0,s). It uses full prepared response information. It is not an independently elicited prior or a procedure that never used held-out outcomes.

Regional sigma uses N(0,1) on class b, dpar sigma: a log-SD prior replacing the direct scalar sigma prior. Depth-varying sigma uses smooth-scale priors 0.5 or 0.25 in recorded tests. Cast and AR priors also differ by role. Resolve every run's prior table rather than applying a mean prior to every distributional component.

The shape baseline uses Student-t(3,0.73,5.448) Intercept, N(0,2) mean b, Student-t(3,0,1.5) mean sd/sds, and N(0,1) regional log-sigma coefficients.

## Iteration accounting

iter includes warmup. Four chains with iter 3000 and warmup 1000 produce 8,000 retained draws before prediction subsampling, not necessarily 8,000 effective samples.

More data does not remove the need for sufficient MCMC. More observations can tighten the posterior and alter geometry. Effective sample size, mixing, and Monte Carlo precision remain necessary. Warmup adapts step size and metric; retained draws estimate the posterior.

| Configuration | Chains | Warmup/chain | Retained/chain | Delta | Tree cap |
|---|---:|---:|---:|---:|---:|
| Original/baseline main | 4 | 1000 | 4000 | .97 | 14 |
| Baseline longer warmup | 4 | 2500 | 4000 | case table | case table |
| Fast main | 2 | 500 | 500 | .90/.95 | 10/12 |
| Fast short | 2 | 250 | 250 | case table | case table |
| k tests | 2 | 500 | 500 | .95 | 12/14 |
| Equation main | 2 | 500 | 500 | .95/.99 | 10/12/14 |
| Equation long | 4 | 1000 | 1000 | .95 | 12 |
| Targeted/extensive usual | 4 | 1000 | 1000 | .95 | 12 |
| Follow-up usual | 4 | 1000 | 2000 | .95 | 14 |
| Shape main | 4 | 1000 | 2000 | .97 | 14 |
| Shape longer warmup | 4 | 2000 | 2000 | .97 | 14 |
| Shape longer sampling | 4 | 1000 | 4000 | .97 | 14 |
| Shape delta controls | 4 | 1000 | 2000 | .95/.99 | 14 |
| Shape prior checks | 4 | 500 | 500 | case table | case table |

These are summaries; individual matrices retain exceptions. Warmup 2500 plus retained 4000 requires iter 6500.

## Computation diagnostics

Use post-warmup R-hat, bulk/tail ESS, divergences, treedepth hits, energy diagnostics, and relevant parameter/predictive precision. Keep warmup separate. A completed file or Slurm COMPLETED state is not proof of diagnostic validity.

The early fast output showed high warmup tree saturation. It did not establish final posterior saturation. Later equation B at cap 10 had 997/1000 post-warmup hits; C had 419/1000. Raising the cap resolved many hits, but short chains still mixed inadequately.

Higher adapt_delta generally requests smaller steps and can reduce divergences, while increasing cost and tree saturation. Higher max_treedepth permits longer trajectories, without repairing parameterization. Depth 14 permits far more leapfrog work than depth 10; use observed n_leapfrog, not a theoretical maximum, to explain runtime.

Across 67 completed shape posterior fits, maximum R-hat was at most 1.00531, minimum bulk ESS about 1332, minimum tail ESS about 1808, and minimum E-BFMI about .622. Eleven fits still diverged. Fifty-six passed the review convergence screen, 27 the stricter zero-tree-hit screen. Retain raw diagnostics and both flags.

## Provisional settings

Four chains, 1000 warmup plus 2000 retained, delta .97 and cap 14 are a useful controlled baseline, not a proved minimum.

Regional MLD at delta .95 (case 241) was clean at about 28.09 hours versus 33.23 hours for the .97 anchor, with almost identical errors. Tight mean b N(0,1), case 249, was clean at 24.11 hours and retained the improvement. Confirm on larger cohorts before selecting an optimum.

For half-shift, delta .99 (case 234) increased hits to 136 from two in the anchor without useful error improvement. Increasing controls is not uniformly beneficial.

## Resource and timing interpretation

Shape fits request four chains × four threads = 16 CPUs through 50% data; four × eight = 32 CPUs at 75/100%. Memory is 64–256 GiB by case. Checks request 4 CPUs/64 GiB, controller 1 CPU/4 GiB, preparation 1 CPU/96 GiB.

Historical threading comparisons include unthreaded controls. The environment caps BLAS/OpenMP helper threads at one to avoid oversubscription while Stan uses requested within-chain threads.

Report compilation, warmup, sampling, fit walltime, allocation, CPU model, and CPU-hours where available. Selected current runs used Xeon Gold 6248 at 2.5 GHz, but different nodes/load limit speed claims. Original workstation totals near 10–11 days per chain are not direct benchmarks against smaller 150 m Ibex subsets.
