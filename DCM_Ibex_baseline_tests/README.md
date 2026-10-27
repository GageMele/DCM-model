# Historical DCM model baseline on Ibex

This campaign uses the exact requested Gaussian formula:

```r
# Region-specific depth profiles; a shared linear distance effect.
andrew_form <- brms::bf(
  chl_norm ~ 1 + Region + log10_distance +
    s(depth, by = Region, bs = "cr", k = 15) +
    (1 | Ocean:Region:date) + (1 | Ocean:Region:enso),
  family = gaussian()
)
```

The formula has no MLD effect. Region changes the depth profile; distance has
a shared linear effect on the Gaussian mean. The depth range remains 10–150 m
to match the current campaigns; this is not a reproduction of the old 250 m
data range. Chlorophyll normalization, observation values and cast membership
are preserved. No smoothing or binning is applied.

## Twelve cases

Every case uses four chains, `cores = 4`, `threading(2)`, Gaussian likelihood,
`k = 15`, `refresh = 50`, `save_pars(all = TRUE)`, CmdStanR, and a fixed seed
1234. Four chains × two threads request eight physical CPUs per fit.

| Cases | Whole-cast subset | Priors | Iterations | Warmup | adapt_delta | Maximum tree depth |
|---|---|---|---:|---:|---:|---:|
| 1 / 2 | 10% | Custom / default | 5000 | 1000 | 0.97 | 14 |
| 3 / 4 | 20% | Custom / default | 5000 | 1000 | 0.97 | 14 |
| 5 / 6 | 50% | Custom / default | 5000 | 1000 | 0.97 | 14 |
| 7 / 8 | 10% | Custom / default | 5000 | 1000 | 0.97 | 12 |
| 9 / 10 | 10% | Custom / default | 5000 | 1000 | 0.99 | 14 |
| 11 / 12 | 10% | Custom / default | 6500 | 2500 | 0.97 | 14 |

The actual `warmup = 1000`, `adapt_delta = 0.97`, `max_treedepth = 14` values
in the supplied code define the baseline. Alternative values mentioned only
in its comments do not replace them. Longer warmup adds 1500 total iterations
to retain the same 4000 post-warmup draws per chain (16000 across four chains).

Custom priors are exactly:

```r
# Historical constants; 2 * 2.724 is exactly 5.448.
priors <- c(
  brms::set_prior("student_t(3, 0.73, 5.448)", class = "Intercept"),
  brms::set_prior("normal(0, 2)", class = "b"),
  brms::set_prior("student_t(3, 0, 1.5)", class = "sd"),
  brms::set_prior("student_t(3, 0, 1.5)", class = "sds"),
  brms::set_prior("student_t(3, 0, 2.5)", class = "sigma")
)
```

Default cases pass `prior = NULL` and save the defaults resolved by the
installed brms version for that subset. They do not retain a custom prior on
any class. Defaults may differ between data subsets.

## Inputs and nesting

The source is the frozen first fast150 campaign recorded in
`logs/latest_fast150_experiment.txt`. Its full input and 10%/20% RDS files are
copied and compared byte-for-byte. Preparation also verifies that the copied
10% and 20% data reproduce the parent data's original deterministic selection.
It extends that same regional cast order to 50%, keeping every eligible
observation from selected casts. Thus 10% is nested in 20%, which is nested
in 50%. Fractions refer to casts, not rows.

The actual prepared input validates to 108 / 216 / 540 casts and 15228 /
30455 / 76139 observations at 10% / 20% / 50%, respectively. All eight regions
remain represented. If the source differs, preparation checks the real data
before starting any fit; see `subset_report.csv` for its authoritative counts.

## Upload and submit

Download `DCM_Ibex_baseline_tests_150m.zip`. In local Windows PowerShell:

```powershell
# Upload the downloaded package to the existing Ibex project.
scp "$env:USERPROFILE\Downloads\DCM_Ibex_baseline_tests_150m.zip" meleg@ilogin.ibex.kaust.edu.sa:/ibex/user/meleg/DCM_MODEL/
```

In the Ibex terminal:

```bash
# Extract the new package and queue one baseline campaign.
cd /ibex/user/meleg/DCM_MODEL
unzip -oq DCM_Ibex_baseline_tests_150m.zip
bash DCM_Ibex_baseline_tests/submit_tests.sh

# Read scheduler states, printed progress and phase-separated CSV diagnostics.
bash DCM_Ibex_baseline_tests/monitor_all_dcm_v3.sh --diagnostics
```

The submitter first runs an isolated mock-Slurm self-test, checks shell/R
syntax, validates all 12 case settings, and checks Slurm's largest resource
request with `sbatch --test-only`. It then submits a preparation job, 12 fits,
12 checks, and one collector: 26 actual jobs. Preparation validates all six
subset/prior combinations with real brms design matrices and Stan code.
Fits require successful preparation. Checks run after their fit ends; a
completed check job can still mean the fit was incomplete, so read the report.

Each fit has a seven-day limit by default. The 10%/20% cases request 32 GB;
50% cases request 64 GB. Checks request 32 GB and two hours. If all fits start
together, requests total 96 CPUs and 448 GiB; Slurm controls actual concurrency.
These are full-length baseline runs, not short timing pilots. A time limit is
not an estimate of completion time.

To change only the scheduler time allowance, set `BASE_TIME` before submission.
For a different frozen source campaign, set `SOURCE_EXPERIMENT` to its absolute
path. Existing runs and campaigns are not modified. Running the submitter again
creates a second campaign, so do not rerun it to monitor an existing submission.

## Monitoring and reports

The updated monitor includes the original pilot, threading benchmarks, fast150,
k150, equation150 and baseline150 campaigns. The basic form (without
`--diagnostics`) is cheaper for frequent checks; full CSV scans are best run
occasionally. Sampling diagnostics are separated from warmup.

```bash
# Show status and progress every five minutes; Ctrl-C stops only the monitor.
watch -n 300 'bash DCM_Ibex_baseline_tests/monitor_all_dcm_v3.sh'

# Locate this submission's frozen files and later regenerate its report.
experiment=$(cat logs/latest_baseline150_experiment.txt)
source "$experiment/ibex_environment.sh"
Rscript --vanilla "$experiment/04_collect_results.R" "$experiment"
```

Each fit saves `run_config.rds`, requested/resolved priors, model code,
`fit.rds`, and chain CSV files including warmup. Check jobs save coefficient
R-hat, bulk/tail ESS, divergences, tree-hit counts and fractions, leapfrog
counts, E-BFMI, time/CPU efficiency, trace plots and in-sample predictive checks.
The collector produces `case_status.csv`, `all_results.csv`, `comparison.txt`,
and `matched_baseline_pairs.csv`. Zero completed fits are handled explicitly.

Tree saturation is an efficiency concern: interpret it with R-hat, ESS,
divergences and predictive behavior. The conservative diagnostic screen is not
a finding of scientific validity. Training predictive checks do not measure
independent predictive performance. Compare equations on the same data and
account for their different draw/threading budgets; raw walltime alone is not
a fair comparison.

## Local validation limits

The package was checked with real prepared data and brms 2.23.0. It generated
Stan data/code for both prior choices at all three fractions and reproduced
the exact requested formula and prior constants. No MCMC was run locally and
no Slurm job has been submitted by this local preparation. Local Bash execution
is restricted; shell syntax, the mock submission self-test, and cluster resource
acceptance are checked on Ibex by the submission command above.

Reference: https://paulbuerkner.com/brms/reference/brm.html
