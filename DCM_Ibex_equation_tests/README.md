# DCM: three exact equations on the same frozen 10% subset

This campaign tests the three equations requested on 2026-10-15. It creates new
jobs and directories. Existing fits and queued jobs are unchanged.

## Exact model definitions

All formulas are explicit in `equation.R`, use `family = gaussian()`, and retain
`1 + Region + log10_distance + (1|Ocean:Region:date) + (1|Ocean:Region:enso)`.
All depth terms use cubic regression splines with **k=15**.

- **A / additive:** `MLD_temp_MIX + s(depth, by=Region, bs="cr", k=15)`.
- **B / mld_varying:** `s(depth, by=MLD_temp_MIX, bs="cr", k=15) +
  s(depth, by=Region, bs="cr", k=15)`.
- **C / mld_only:** `s(depth, by=MLD_temp_MIX, bs="cr", k=15)`.

In A, MLD has a shared constant linear coefficient, while each region has its
own depth smooth. In B, the shared MLD coefficient varies smoothly with depth,
in addition to each region's depth smooth. In C, that MLD-dependent depth
component is shared across regions; regional intercepts provide vertical
offsets at a fixed MLD and distance. A numeric by variable multiplies the depth
smooth: its contribution is `MLD_temp_MIX * g(depth)`. It is linear in MLD at
each depth. No tensor or additional main MLD term is silently inserted.

MLD and log10_distance retain their original numeric units. The source of the
numerical priors is the frozen original fit script from the previous campaign.
The new models keep the current priors, not default or empirical alternatives.
Identical numerical prior constants can imply different function priors across
different formulas; these are scientific model comparisons, not pure CPU tests.

## Frozen input

The default source is the campaign in `logs/latest_fast150_experiment.txt`,
regardless of whether its fits timed out. Its preparation must have completed.
The script copies `data/p10.rds` byte-for-byte: the known input has **108 casts,
15,228 observations, eight regions, and depths 10-150 m**. No casts are
reselected, no normalization is recomputed, and no smoothing or binning is used.
`SOURCE_EXPERIMENT=/absolute/path` overrides that source when needed.

## Eighteen matched cases

The following six settings are repeated for each equation. Thus cases 1/7/13,
2/8/14, and so forth differ only in formula. All use four threads per chain.

| A case | B case | C case | Chains | Total iterations/chain | Warmup/chain | Delta | Max tree depth |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 7 | 13 | 2 | 1000 | 500 | .95 | 10 |
| 2 | 8 | 14 | 2 | 1000 | 500 | .95 | 12 |
| 3 | 9 | 15 | 2 | 1000 | 500 | .95 | 14 |
| 4 | 10 | 16 | 2 | 1000 | 500 | .99 | 12 |
| 5 | 11 | 17 | 2 | 1500 | 1000 | .95 | 12 |
| 6 | 12 | 18 | 4 | 2000 | 1000 | .95 | 12 |

The 1500/1000 cases isolate longer warmup while retaining 500 draws per chain.
The four-chain cases produce 4000 retained draws overall and have more warmup
work. They are longer confirmation runs, not isolated chain-count speed tests.

Each fit requests **32 GiB** and a **seven-day walltime limit**. Two-chain fits
request eight CPUs; four-chain fits request 16 CPUs. If all 18 run together they
request **168 CPUs and 576 GiB RAM**. Slurm determines their actual start times.
The walltime is an allocation limit, not a completion forecast. Override with
`ETEST_TIME`, `ETEST_MEMORY`, and `ETEST_CONSTRAINT` when appropriate. Hardware
defaults to `intel&cascadelake`, as in the preceding campaigns.

## Upload and submit once

Windows PowerShell:

```powershell
# Upload the complete new package, then connect to Ibex.
$equationZip = "$env:USERPROFILE\Downloads\DCM_Ibex_equation_tests_150m.zip"
if (-not (Test-Path -LiteralPath $equationZip)) { throw "Archive not found: $equationZip" }
scp "$equationZip" "meleg@ilogin.ibex.kaust.edu.sa:/ibex/user/meleg/DCM_MODEL/"
ssh meleg@ilogin.ibex.kaust.edu.sa
```

Ibex:

```bash
cd /ibex/user/meleg/DCM_MODEL
unzip -o DCM_Ibex_equation_tests_150m.zip

# Check shell syntax, freeze inputs, and submit exactly one new campaign.
bash -n DCM_Ibex_equation_tests/submit_tests.sh
bash DCM_Ibex_equation_tests/submit_tests.sh
```

The helper validates Slurm resource acceptance before any submission. It then
submits one model-construction/preparation job, 18 fits, 18 checks, and one
collector: **38 jobs total**. Fits wait for successful preparation. That job
constructs the complete brms data and Stan code for all three formulas on the
actual frozen data before releasing any fit.

Each submission-helper invocation creates a new campaign. Do not repeat it
because jobs are pending. If submission fails partway, inspect the recorded
`manifest.tsv` and `job_events.tsv` before submitting another campaign.

## Monitor all model campaigns

```bash
cd /ibex/user/meleg/DCM_MODEL

# All recorded DCM campaigns, including the new equation comparisons.
bash DCM_Ibex_equation_tests/monitor_all_dcm_v2.sh

# Add saved warmup/post-warmup tree-depth and divergence diagnostics.
bash DCM_Ibex_equation_tests/monitor_all_dcm_v2.sh --diagnostics

# Locate this campaign and inspect preparation output when available.
experiment=$(cat logs/latest_equation150_experiment.txt)
prepare_job=$(cat "$experiment/preparation_job_id.txt")
for extension in out err; do
    path="logs/dcm-e150-prepare-${prepare_job}.${extension}"
    [[ ! -f "$path" ]] || tail -n 40 "$path"
done
[[ ! -f "$experiment/equation_validation.csv" ]] || cat "$experiment/equation_validation.csv"
```

The new monitor is separate from the old monitor. Its scans include the
original pilot, threading benchmarks, fast150, k150 and equation150. Saved CSV
counts may lag printed iteration checkpoints because of buffering. Warmup
diagnostics are reported separately from posterior diagnostics.

## Repair the failed first-campaign report

The confirmed failure occurs when only prior-predictive results exist: a scalar
`FALSE` or `NA_real_` cannot be assigned as a new column of a zero-row table.
The repair uses vectors of length `nrow(hardware)` instead.

```bash
# Create a separate repaired reporting script and regenerate the old reports.
# This does not restart models or change the original frozen reporting script.
bash DCM_Ibex_equation_tests/repair_fast150_report.sh
```

This repair targets the campaign referenced by `logs/latest_fast150_experiment.txt`.
An explicit campaign path can be passed as its first argument. It preserves the
original `04_collect_results.R` and writes `04_collect_results_repaired.R`.
Only reporting outputs are regenerated. The timed-out fits remain incomplete;
the repaired report does not recover a completed fit from partial CSVs.

The new equation collector also handles zero completed results and retains
Slurm states when accounting is available. A completed reporting job means its
reports were written; inspect `case_status.csv` to see whether models finished.

## Results

Each completed case has the fit, exact formula, numerical priors, configuration,
hardware information, warmup-inclusive sampler CSVs, R-hat/ESS diagnostics,
E-BFMI, predictive density, traces and example population profiles. White image
backgrounds keep exported plots legible in light and dark viewers.

Campaign reports include `all_results.csv`, `comparison.txt`, `case_status.csv`,
`diagnostic_passes.csv`, and `matched_equation_pairs.csv`. Reported predictive
checks use training observations; they are not held-out predictive validation.
Compare profile behavior and the scientific assumptions alongside runtime.
Lower tree-hit percentages or faster sampling alone do not establish a better
ecological model. Numeric priors and Gaussian support remain sensitivity issues.

## Documentation

- Numeric/factor by smooths: https://stat.ethz.ch/R-manual/R-patched/RHOME/library/mgcv/help/s.html
- brms smooths: https://paulbuerkner.com/brms/reference/s.html
- Stan diagnostics: https://mc-stan.org/learn-stan/diagnostics-warnings.html
