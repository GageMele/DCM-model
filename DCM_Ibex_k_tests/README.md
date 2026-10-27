# DCM: smaller tensor bases and higher tree-depth limits

This isolated campaign submits 12 posterior fits. It does not cancel, modify,
extend or restart any existing job. Run its submission helper once.

## Frozen scientific input

The source defaults to the campaign in `logs/latest_fast150_experiment.txt`.
Override it with `SOURCE_EXPERIMENT=/absolute/path/to/previous/campaign`.
The helper copies that campaign's `data/p10.rds` byte-for-byte, together with its
original model script and environment. No new casts are selected. There is no
smoothing, depth binning, or change in surface normalization. Data remain at
10-150 m. The known input has 108 casts and 15,228 observations across eight regions.

The model retains its three-dimensional region-specific depth-distance-MLD
interaction, Gaussian likelihood, group effects and existing numerical priors.
The k vector is ordered **depth, log10 distance, MLD**. Lower k still allows
smooth interactions but reduces their representational capacity. Changing the
basis also changes its induced function prior, even with identical numerical priors.

## Cases

| Case | k | Max tree depth | Total iterations/chain | Warmup/chain | Delta |
|---|---|---:|---:|---:|---:|
| 1 | 15,5,5 | 12 | 1000 | 500 | .95 |
| 2 | 10,5,5 | 12 | 1000 | 500 | .95 |
| 3 | 10,5,5 | 14 | 1000 | 500 | .95 |
| 4 | 7,5,5 | 12 | 1000 | 500 | .95 |
| 5 | 7,5,5 | 14 | 1000 | 500 | .95 |
| 6 | 5,5,5 | 12 | 1000 | 500 | .95 |
| 7 | 5,5,5 | 14 | 1000 | 500 | .95 |
| 8 | 7,4,4 | 12 | 1000 | 500 | .95 |
| 9 | 7,4,4 | 14 | 1000 | 500 | .95 |
| 10 | 5,3,3 | 12 | 1000 | 500 | .95 |
| 11 | 7,4,4 | 12 | 1500 | 1000 | .95 |
| 12 | 7,4,4 | 12 | 1000 | 500 | .99 |

Every case uses two concurrent chains and four Stan threads per chain: 8 CPUs.
Each chain retains 500 posterior draws; case 11 isolates longer adaptation
without increasing its retained-draw budget.

Cases 1/2/4/6 isolate depth k, keeping distance and MLD k at five. Cases 2/3,
4/5, 6/7 and 8/9 isolate the tree-depth limit. Cases 8/11 isolate warmup length;
cases 8/12 isolate adapt_delta. Smaller depth range motivates investigating k,
but does not prove a smaller k is adequate for DCM profiles.

## Upload and submit

From Windows PowerShell, after downloading the ZIP:

```powershell
# Upload the separate package, then connect to Ibex.
$kZip = "$env:USERPROFILE\Downloads\DCM_Ibex_k_tests_150m.zip"
if (-not (Test-Path -LiteralPath $kZip)) { throw "Archive not found: $kZip" }
scp "$kZip" "meleg@ilogin.ibex.kaust.edu.sa:/ibex/user/meleg/DCM_MODEL/"
ssh meleg@ilogin.ibex.kaust.edu.sa
```

In the Ibex shell:

```bash
cd /ibex/user/meleg/DCM_MODEL
unzip -o DCM_Ibex_k_tests_150m.zip

# Validate shell syntax and submit one new campaign.
bash -n DCM_Ibex_k_tests/submit_tests.sh
bash DCM_Ibex_k_tests/submit_tests.sh
```

The helper first checks Slurm resource acceptance with `sbatch --test-only`.
It then submits one preparation job, 12 independent fits, 12 checks/plot jobs,
and one results job: 26 jobs in total. Fits wait for successful preparation,
which constructs every proposed brms basis on the actual frozen data.

Each fit requests 32 GiB RAM and a **seven-day** limit. This is an allocation
limit, not a completion forecast. If all 12 fits run together they request
96 CPUs and 384 GiB. Slurm controls actual start times. Increase of an already
running job's limit is not attempted. Default hardware constraint is
`intel&cascadelake`, consistent with the earlier campaign's request.

Resource overrides, when needed, are `KTEST_TIME`, `KTEST_MEMORY` and
`KTEST_CONSTRAINT`. For example, `KTEST_TIME=2-00:00:00 bash
DCM_Ibex_k_tests/submit_tests.sh` creates a new campaign with a two-day limit.
If the initial resource check rejects the request, it submits no jobs.
Seven-day requests require the account to be eligible under Ibex policies.

Never repeat the submission helper simply because jobs are pending: each call
creates a new experiment. If a submission fails partway through, inspect
`job_events.tsv` and `manifest.tsv` before starting another campaign.

## Monitoring

```bash
cd /ibex/user/meleg/DCM_MODEL

# Print accounting plus current per-chain warmup/sampling diagnostics.
bash DCM_Ibex_k_tests/monitor.sh

experiment=$(cat logs/latest_k150_experiment.txt)
cat "$experiment/manifest.tsv"

# Inspect preparation and its basis report after it runs.
prepare_job=$(cat "$experiment/preparation_job_id.txt")
for extension in out err; do
    path="logs/dcm-k150-prepare-${prepare_job}.${extension}"
    [[ ! -f "$path" ]] || tail -n 30 "$path"
done
[[ ! -f "$experiment/basis_validation.csv" ]] || cat "$experiment/basis_validation.csv"

# Show each chain's latest printed iteration.
while IFS=$'\t' read -r case_id label fit_job check_job cpus run_dir; do
    [[ "$case_id" =~ ^[0-9]+$ ]] || continue
    path="logs/dcm-k150-c${case_id}-${fit_job}.out"
    [[ -f "$path" ]] || continue
    echo "Case $case_id: $label"
    awk '/^Chain [1-4] Iteration:/ {latest[$2]=$0}
         END {for (i=1;i<=4;i++) if (i in latest) print latest[i]}' "$path"
done < "$experiment/manifest.tsv"
```

Live saved draws can lag the printed iteration because of file buffering.
The live report separates warmup from sampling, counts `treedepth__ >= cap`,
and reports divergence counts and average leapfrog steps. Its `recent%` is
the tree-hit percentage among the last 50 saved draws per phase (or fewer,
when fewer are available). Warmup diagnostics are provisional, not a verdict
on the eventual posterior. Compare the late-warmup and post-warmup rows.

The monitor also works on a specified older campaign:
`bash DCM_Ibex_k_tests/monitor.sh /absolute/path/to/experiment` if its manifest
has the same six columns and its directory contains `live_diagnostics.awk`.
The default use above requires no changes to earlier campaign directories.

## Results and interpretation

After checks finish, the new campaign has:

- `comparison.txt`, `all_results.csv`, `case_status.csv`, and `diagnostic_passes.csv`.
- `matched_treedepth_pairs.csv` for completed pairs at tree depth 12 versus 14.
- Per-case `fit.rds`, `model.stan`, input copy, exact formula, k, priors and hardware.
- Per-case sampler CSVs, including warmup.
- `checks/sampler_summary_by_phase.csv`, parameter R-hat and bulk/tail ESS.
- Trace plots, predictive density, and expected population profiles.

The automated screen uses R-hat <=1.01, bulk/tail ESS >=400, no post-warmup
divergences or tree hits, and E-BFMI >=.3. This conservative exploratory screen
does not prove scientific adequacy. Assess profile shape and sensitivity to k
alongside sampling diagnostics. Training predictive checks are not independent
predictive validation. Lower k changes the model; compare accuracy and the
scientific quantities of interest rather than ranking only by runtime or
parameter ESS, whose parameter set changes with k.

Higher max tree depth allows substantially more leapfrog work; a run at 14
can be slower than one at 12 even when it explores more effectively. A larger
adapt_delta may reduce post-warmup divergences while increasing runtime and
tree-depth pressure. Longer warmup is tested separately, not assumed to fix
all sampling problems.

Incomplete/timed-out fits remain listed. Check their actual Slurm states with
`sacct`; preserved CSVs are not a completed fit or a general-purpose restart.

## Source documentation

- mgcv tensor bases: https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/t2.html
- Basis choice: https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/choose.k.html
- Stan diagnostics: https://mc-stan.org/learn-stan/diagnostics-warnings.html
- Ibex limits: https://docs.hpc.kaust.edu.sa/policy/ibex.html
