# Fast exploratory chlorophyll model tests on Ibex

This package starts a separate experiment under `/ibex/user/meleg/DCM_MODEL/experiments/fast150_*`.
It leaves the existing pilot and threading benchmarks running with their original 250 m data.
It does not edit the existing preparation, fit, plotting, or environment scripts.

The new tests use the existing model formula and original priors extracted from a frozen copy
of `02_fit_model.R`. The Gaussian family, region-specific three-dimensional tensor, random
effects, response definition, and spline dimensions `k = c(15, 5, 5)` are retained.
The new depth range is **10–150 m**, inclusive. The original 5–15 m surface baselines are
retained. No chlorophyll smoothing, depth binning, thinning, clipping, or averaging is added.

## Experiment design

The baseline is case 3: 10% of full eligible casts, two parallel unthreaded chains,
1,000 total iterations per chain, including 500 warmup, `adapt_delta = 0.95`,
`max_treedepth = 10`, and the original priors. This retains 1,000 posterior draws in total.

| Cases | Purpose | Settings |
|---|---|---|
| 1–12 | Sampler and prior comparisons on the same 10% subset | Three prior choices × delta 0.90/0.95 × tree depth 10/12; two chains; 1,000/500 iterations/warmup |
| 13–15 | Smaller dataset | 5% casts; three prior choices; baseline sampler settings |
| 16–18 | Larger exploratory dataset | 20% casts; three prior choices; baseline sampler settings |
| 19 | Shorter fit | 10% casts; current priors; two chains; 500/250 iterations/warmup |
| 20 | Matched chain-count comparison | 10% casts; current priors; four parallel chains; 1,000/500 |
| 21–22 | Matched threading comparison | 10% casts; current priors; two chains × two/four threads |
| 23–24 | Prior predictive checks | 10% casts; current/empirical proper priors; two chains; 500/250; likelihood disabled |

There are **22 posterior fits and two prior-predictive fits**. Each is independently eligible
to run after preparation; Slurm decides actual concurrency. Each gets its own output directory
and a separate checking/plotting job. All files are snapshotted before submission.

The fit jobs request **58 CPUs and 768 GiB RAM in total** if all 24 run simultaneously.
Each fit requests 32 GiB and 12 hours by default. This is separate from storage allocation and
from resources requested by your existing jobs. These are time limits, not runtime predictions.
The preparation job requests one CPU, 16 GiB, and two hours; checking jobs use one CPU and 8 GiB.

Two chains reduce total work and resource requests compared with four. Because the existing
four chains already run simultaneously, halving the number of chains does not automatically
halve the time needed by an individual chain. Smaller datasets and lower tree-depth limits
are likely to have a larger effect. Lower tree limits can reduce exploration and effective
sample size, so speed is evaluated alongside divergences, tree hits, R-hat, ESS, and E-BFMI.

## Whole-cast subsets

Subsets are drawn from **the full prepared bundle**, not the existing 25% pilot.
First, observations deeper than 150 m are removed. Casts are then selected within regions,
retaining every remaining observation from each selected cast. All regions are retained.
The 5%, 10%, and 20% selections are nested and fixed across prior/sampler comparisons.

Selection preserves at least five distinct distance and MLD values per region, as required
by the existing basis. Counts are allocated proportionally subject to this support constraint.
If a requested sample is smaller than the minimum support requirements, its achieved fraction
is increased and explicitly reported. Subsets are never made smaller by reducing the basis.
Rows are shuffled once for balanced threading; measurement values and cast identities are unchanged.

Actual-data validation on the local full bundle produced:

| Dataset | Casts | Observations | Depth range |
|---|---:|---:|---|
| Original full prepared data | 1,080 | 258,984 | 10–250 m |
| Full after depth restriction | 1,080 | 151,962 | 10–150 m |
| 5% casts | 54 | See `subset_report.csv` | 10–150 m |
| 10% casts | 108 | See `subset_report.csv` | 10–150 m |
| 20% casts | 216 | See `subset_report.csv` | 10–150 m |

Cast percentages do not guarantee identical observation percentages when profile lengths
vary. The reports written by the Ibex preparation job are authoritative for that experiment.

## Priors

The original script currently supplies:

| Parameter class | Current prior |
|---|---|
| Intercept | `student_t(3, 0.73, 5.448)` |
| Population coefficients `b` | `normal(0, 2)` |
| Group SD `sd` | Half `student_t(3, 0, 1.5)` |
| Spline SD `sds` | Half `student_t(3, 0, 1.5)` |
| Residual SD `sigma` | Half `student_t(3, 0, 2.5)` |

The `current` cases read these definitions from your frozen original script instead of
hard-coding another copy. The `default` cases supply `prior = NULL`, letting the installed
brms version determine its defaults for that dataset. These are not all flat: coefficient
defaults are improper flat priors, while other classes have distributional defaults.
The exact resolved priors are saved as `resolved_priors.csv` and `.txt` after every successful fit.
Default-prior-only fits are excluded because improper coefficient priors do not define a
proper prior-predictive distribution.

The `empirical` cases use the **full prepared response after restriction to 150 m**.
Each cast has equal total weight when calculating response location and scale, preventing
longer profiles from dominating. Let:

- `m` = cast-balanced response median;
- `s` = maximum of 0.5, cast-balanced MAD × 1.4826, and cast-balanced IQR / 1.349.

The custom priors are:

| Parameter class | Empirical exploratory prior |
|---|---|
| Intercept | `student_t(3, m, 2*s)` |
| `b` | `normal(0, s)` |
| `sd`, `sds`, `sigma` | Half `student_t(3, 0, s)` |

Local actual-data statistics gave `m ≈ 0.2373924` and `s ≈ 0.8767800`;
the observation-weighted response mean was 1.102393 and SD was 2.806814.
These numbers are recomputed on Ibex and saved in `full_response_statistics.csv`.
The robust scale reduces sensitivity to the extreme tail, without modifying observations.

These are **empirical Bayes sensitivity choices**, not independent prior evidence.
They deliberately reuse response information from the full dataset. Their spline-scale
choices are heuristic because spline coefficients depend on the basis parameterization.
Inspect the prior-predictive checks and compare against the current/default cases before
considering a production prior. Do not present these tests as independent validation.

## Upload and submit

Download `DCM_Ibex_fast_tests_150m.zip`. In Windows PowerShell, change `$fastZip` to its
actual downloaded location, then upload it:

```powershell
# Use the actual location where you saved the downloaded archive.
$fastZip = "C:\Users\gagem\Downloads\DCM_Ibex_fast_tests_150m.zip"
if (-not (Test-Path -LiteralPath $fastZip -PathType Leaf)) {
    throw "File not found: $fastZip"
}
scp "$fastZip" "meleg@ilogin.ibex.kaust.edu.sa:/ibex/user/meleg/DCM_MODEL/"
ssh meleg@ilogin.ibex.kaust.edu.sa
```

On Ibex:

```bash
cd /ibex/user/meleg/DCM_MODEL

# Extract only this isolated package; existing experiment snapshots stay intact.
unzip -o DCM_Ibex_fast_tests_150m.zip

# Validate shell syntax, freeze inputs, and queue the independent experiments.
bash -n DCM_Ibex_fast_tests/submit_tests.sh
bash DCM_Ibex_fast_tests/submit_tests.sh
```

The helper checks R and Slurm script syntax before submission. It freezes the full input,
current model source, environment, case manifest, and experiment scripts. Preparation,
24 fitting jobs, 24 checking jobs, and one collector are queued with dependencies.
All model/data preparation occurs on compute nodes. Your old jobs are not canceled.

To choose different requested resources for a **new** campaign:

```bash
# Example: request 24-hour limits and 48 GiB per fit for a new campaign.
TEST_TIME=24:00:00 TEST_MEMORY=48G \
  bash DCM_Ibex_fast_tests/submit_tests.sh
```

Every invocation creates and submits a new campaign; use it once for the desired settings.
The default CPU constraint is `intel&cascadelake`; actual CPU models are recorded.

## Monitor preparation and fits

```bash
cd /ibex/user/meleg/DCM_MODEL
experiment=$(cat logs/latest_fast150_experiment.txt)
prepare_job=$(cat "$experiment/preparation_job_id.txt")

# Confirm preparation succeeded before interpreting pending dependent fit jobs.
sacct -j "$prepare_job" --format=JobID,State,ExitCode,Elapsed,MaxRSS
tail -n 40 "logs/dcm-f150-prepare-${prepare_job}.out"
tail -n 40 "logs/dcm-f150-prepare-${prepare_job}.err"

# The three reports document the actual observation/depth/subset choices.
cat "$experiment/depth_trim_report.csv"
cat "$experiment/subset_report.csv"
cat "$experiment/full_response_statistics.csv"

# Check all fit jobs with readable names and elapsed times.
fit_ids=$(awk 'NR>1 {print $3}' "$experiment/manifest.tsv" | paste -sd, -)
squeue -j "$fit_ids" -o '%.12i %.22j %.10M %.12l %.2t %.25R'

# Show the baseline (case 3) log once its fit job starts.
baseline_job=$(awk -F '\t' '$1==3 {print $3}' "$experiment/manifest.tsv")
if [[ -f "logs/dcm-f150-c3-${baseline_job}.out" ]]; then
    tail -n 40 "logs/dcm-f150-c3-${baseline_job}.out"
fi
```

Preparation failures leave fits pending on dependency; inspect the preparation error first.
For partial submission, `job_events.tsv` records accepted jobs even if a later submission fails.
Do not launch another campaign merely because the first jobs are pending.

## Results

The collector writes:

- `case_status.csv`: success, failure, incomplete cases, and job IDs;
- `all_results.csv`: runtime, retained draws, ESS throughput, diagnostics, and training PPC;
- `diagnostic_passes.csv`: posterior fits that pass the conservative diagnostic screen;
- `matched_threading_results.csv`: comparable cases 3, 21, and 22 with CPU-model checks;
- `comparison.txt`: readable compact table;
- `READ_RESULTS_FIRST.txt`: interpretation limits.

Each case is under `runs/case_XX/`. Look in its `checks/` directory for trace plots,
predictive density, population profiles at real cast distance/MLD pairs, parameter diagnostics,
and sampler diagnostics **including warmup**. Population profiles stop at 150 m.

```bash
# These files appear after the dependent collector completes.
cat "$experiment/comparison.txt"
cat "$experiment/case_status.csv"
cat "$experiment/matched_threading_results.csv"

# Check observed runtime/memory and any failures or timeouts.
sacct -j "$fit_ids" --format=JobID,State,ExitCode,Elapsed,AllocCPUS,MaxRSS
```

Throughput is median bulk ESS of population/scale parameters per sampling second or per
requested CPU-hour. Fit runtime includes compilation, basis construction, sampling, and
serialization; per-chain Stan warmup/sampling times are recorded separately. Checking jobs
are separate and are not included in fit runtime. Queue waiting is excluded.

The diagnostic screen is R-hat <= 1.01, minimum bulk/tail ESS >= 400, zero post-warmup
divergences and tree-depth hits, and minimum E-BFMI >= 0.3. Prior-only runs are not ranked
as posterior fits. Two-chain short runs can fail these screens; that is useful information,
not a reason to lower the criteria. The package preserves runtime diagnostics if plotting fails.

Predictive checks use up to 60 existing observations per region and up to 100 posterior draws.
Training RMSE and interval coverage are **in-sample diagnostics**, not test-set accuracy.
Different subset sizes and different priors define different posterior problems, so do not
choose a production setting by sorting all rows on runtime or training RMSE. Compare matched
cases, inspect diagnostic and predictive plots, then confirm the candidate with a longer fit.

## Existing-job time extensions

Ibex's documented maximum walltime is 14 days for eligible accounts, but Slurm normally
requires an administrator to increase the time limit of an already running job.
The current 48-hour limits remain in effect unless an update is actually accepted.

```bash
# Read-only inspection of your original pilot and four existing benchmarks.
bash DCM_Ibex_fast_tests/inspect_or_extend_current_jobs.sh

# Optional: attempt the authorized extension. Permission rejection is expected
# for an ordinary user; the helper leaves jobs running and reports the result.
bash DCM_Ibex_fast_tests/inspect_or_extend_current_jobs.sh --attempt
```

If rejected, request an extension from Ibex Support for jobs
`53749517, 53759873, 53759874, 53759876, 53759877`, to a **14-day total walltime**.
No support message is sent by this package. Extending time alone does not improve sampling.
Changing a script's `#SBATCH --time` does not change the limit of an already submitted job.

## Validation performed here

R syntax and case consistency were checked. The subset workflow was executed against the
actual full prepared data; whole-cast preservation, unchanged response/predictor measurements,
the 150 m limit, nested fractions, regional basis support, robust weighting, and diagnostic CSV
parsing were verified. The scientific formula/prior assignments were verified in the current
source. brms/CmdStan fits and Slurm scheduling have not been executed locally; those require
your installed Ibex environment. Submission includes Linux shell/R syntax checks.

## Official references

- [brms prior definitions](https://paulbuerkner.com/brms/reference/set_prior.html)
- [brms threading guidance](https://paulbuerkner.com/brms/articles/brms_threading.html)
- [Stan warmup and NUTS configuration](https://mc-stan.org/docs/reference-manual/mcmc.html)
- [Slurm running-job time-limit restrictions](https://slurm.schedmd.com/scontrol.html)
- [Ibex scheduling policy](https://docs.hpc.kaust.edu.sa/policy/ibex.html)
