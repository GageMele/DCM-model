# Validation performed on 15 October 2026

- Parsed every R file successfully.
- Constructed and checked the full 102-case matrix: 96 posterior + 6 prior-only;
  no duplicate configurations; CPUs match chains × threads; all iteration
  budgets equal warmup plus retained draws.
- Verified that warmup comparisons retain the same 1,000 sampling draws per
  chain and the two-chain equal-budget comparisons retain 4,000 draws total.
- Built Stan data and generated Stan code for all **75 distinct equation/data/
  prior/family/thread designs** covering the 102 cases, using the actual prepared
  data and local brms 2.23.0. Resolved all priors, including the defaults and
  Student-t family. `basis_counts.csv` records the constructions.
- Checked that copied 10%/20%/50% input bundles have unchanged MD5 hashes and
  contain the original rows, grouping labels and model values. Their casts remain
  nested. All fitted depths are 10–150 m.
- Selected 108 validation casts (15,226 observations), verified no overlap with
  the 540-cast training set, and checked representation of all eight Regions.
- Checked cast-balanced RMSE/MAE and profile-shape error against known values on
  an unequal-profile-length fixture.
- Executed collector regressions for zero completed runs, only completed
  prior-predictive runs, mixed partial and complete results, different available
  summary columns, matched contrasts, and paired cast comparisons. Prior-only
  and failed-check cases are excluded from predictive rankings.
- Reviewed the submission manifest, resource fields, job dependencies, frozen
  data inputs and monitor discovery patterns.

The new posterior fits have **not** been run locally or submitted to Ibex from
this session. There is no authenticated Ibex connection here. The local Windows
environment cannot start the Rtools Bash runtime, so the Linux mock submission
test could not execute locally. The submitter runs that test automatically on
Ibex before any real `sbatch` submission. It verifies 102 fit jobs, per-case
resources, dependencies and 206 total jobs. Shell syntax is also checked there.

Stan C++ compilation, sampling, and end-to-end posterior predictions will be
exercised by the actual Ibex jobs. Generating Stan code/data is not the same as
fitting a model. Preparation repeats all construction checks with Ibex's installed
package versions and gates fitting through `afterok` dependencies. A failed
preparation prevents fitting; failed checks preserve diagnostic summaries.

`preflight_session.txt` records the local validation environment. No local
validation claims successful MCMC convergence or a preferred scientific model.
