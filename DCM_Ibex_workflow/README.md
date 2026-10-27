# Chlorophyll enhancement model: three scripts and an Ibex run plan

Prepared for Gage Mele, 14 October 2026, from the supplied 3,298-line notebook and the existing repository's Ibex launchers. This package contains exactly three scientific R scripts. The shell files only manage the software environment and Slurm execution.

| File | Responsibility | Main output |
| --- | --- | --- |
| `01_prepare_data.R` | Read CTD and MLD inputs; join, normalize original observations, select casts, audit exclusions | `model_data.rds` and QC CSVs |
| `02_fit_model.R` | Define the current model and priors; construct the spline basis; run CmdStan | `fit.rds`, persistent chain CSVs, configuration, exact data copy |
| `03_check_and_plot.R` | Read the saved fit; assess sampling and predictions; plot profiles and probabilities | Diagnostic CSVs, profile/probability tables, PNG/PDF figures |
| `ibex_environment.sh` | Load the same R module and user packages for every stage | Environment settings |
| `setup_environment.sh` | Install missing packages and a pinned CmdStan version on a compute node | Software and version records |
| `run_stage.slurm` | Launch preparation, fitting, or checks with separately chosen resources | Slurm logs |

## The scientific model retained

```r
# The current Gaussian model, including distance x MLD x depth interactions.
andrew_form <- brms::bf(
  chl_norm ~ 1 + Region +
    t2(depth, log10_distance, MLD_temp_MIX,
       by = Region, bs = c("cr", "cr", "cr"), k = c(15, 5, 5)) +
    (1 | Ocean:Region:date) +
    (1 | Ocean:Region:enso),
  family = gaussian()
)
```

The response is `(Chl - surface_chl) / surface_chl`, where `surface_chl` is each cast's median raw chlorophyll from 5–15 m. Thus `chl_norm = 0` means the surface baseline, `chl_norm = 0.5` means 1.5 times that baseline, and `chl_norm = 1` means twice the baseline. Adding 1 converts the response to a chlorophyll ratio. It is not a percentage unless multiplied by 100.

The final active priors in the attachment are retained: Intercept `student_t(3, 0.73, 5.448)`; fixed coefficients `normal(0, 2)`; group SDs and smooth SDs `student_t(3, 0, 1.5)`; residual SD `student_t(3, 0, 2.5)`. Historical alternative formulas and priors are removed. Review the included prior predictive run before treating these historical constants as scientifically justified for the new tensor model.

The tensor has up to `15 * 5 * 5 = 375` marginal basis products per region before constraints and brms/mgcv reparameterization. Its computational cost cannot be inferred from the previous one-dimensional smooth. Raising all three k values to 15 would create 3,375 products per region. The pilot is for measuring feasibility and sampling behavior, not choosing k solely by R2.

## Explicit cleanup choices

1. **Sampling occurs once.** Literal sequential execution of the notebook's 10% then 25% sampling blocks would retain approximately 2.5% of casts; its later 100% sampling block cannot restore excluded casts. The new default is all eligible casts. A pilot fraction is sampled once within each region, with seed 1234. This is an intentional change; it does not promise the exact same sample as a previous interactive session.
2. **Use unsmoothed, unbinned observations.** Every eligible row retains its input `Chl` and depth. All running-median and depth-binning code has been removed, with no option to enable it. The 5–15 m median is used only to calculate each cast's surface normalization baseline; the depth profile itself is not aggregated. The model's fitted tensor smooth remains part of the regression formula.
3. **Keep all join fields and standardize metadata.** The four-column join remains `LocationCode`, `Region`, `cast_id`, `Ocean`. Known aliases are standardized: `Hawaiian` to the existing CTD/model label `Hawiian`, and `Red Sea` to `RedSea`. Missing, empty and `Not specified` location codes are treated as unknown locations; an unknown location matches another unknown location only within the same Region/cast_id/Ocean tuple. Required cast, region, ocean and date identifiers are still checked, as is uniqueness across dates. Exact duplicate right-table rows are collapsed; conflicting MLD values per key stop execution. Unmatched CTD and unused MLD casts are listed. Rows without usable MLD are excluded explicitly because the formula requires MLD. `identifier_qc.csv` records missing-location counts.
4. **Cast identity uses multiple fields.** Region, location, ocean, cast ID and date keep identically named casts in different sources separate. If the four requested MLD join keys identify several dates, the workflow stops and asks you to add date to both tables/keys.
5. **Dates remain readable grouping labels.** The formula uses a date random intercept, not a continuous date trend. Replacing scaled date labels with ISO date factors retains date group membership.
6. **The retained regions follow active code.** KAUST is excluded; Equatorial is retained. The notebook's active ENSO year lookup is retained (2009 El Nino; 2010/2011 La Nina; other years neutral). This is not a newly validated ENSO classification.
7. **The fitted data and provenance travel with the fit.** Output folders are unique; existing prepared data and run folders are not overwritten. Each fit keeps its input-data copy, configuration, priors, source scripts and software versions. CSV chains remain outside `/tmp`.
8. **All prediction grids contain MLD.** Reference predictions use an actual cast-level distance/MLD pair nearest the regional medians, rather than fixing distance to 1 km in every region. The chosen cast/pair is recorded. Additional plots use observed pairs near low/high distance and low/high MLD targets; differences between those curves do not isolate one predictor's causal effect.
9. **Probabilities have correct labels.** `mean(pred > 0.5)` is a probability of exceeding 1.5 times the surface baseline, not the probability of belonging to the interval 1.5–2 times. The plotting script exports exceedance probabilities for both expected profiles and simulated observations, with legends matching the thresholds. Red Sea figures use >7x through >12x; Pacific figures use >1.5x through >6x. Both sets remain in the tables for all regions.

## 1. Inspect Ibex storage and create the project

Run this in Windows PowerShell:

```powershell
# Connect to the CPU login service; log in as usual.
ssh meleg@ilogin.ibex.kaust.edu.sa
```

Then run this in the Ibex Bash shell:

```bash
# Your 1.5 TB allocation is storage capacity, separate from job RAM.
df -h /ibex/user/meleg
mkdir -p /ibex/user/meleg/DCM_MODEL/{data/raw,data/processed,logs,runs,software}
cd /ibex/user/meleg/DCM_MODEL

# Discover the LIVE R module names. The supplied default is the verified R 4.5.0 module.
module avail R
```

This directory shares your personal allocation with CROCO and other projects. Check its actual free space before large model runs. KAUST documents `/ibex/user/$USER` as the fast 1.5 TB personal HPC filesystem. A Slurm request such as `--mem=64G` reserves job memory; it does not change your disk quota.

## 2. Upload the scripts and two data files

Download and extract the ZIP on your computer. In a NEW Windows PowerShell window:

```powershell
# Replace this path with the folder containing the extracted .R/.sh/.slurm files.
$dcmCode = "C:\Users\gagem\Downloads\DCM_Ibex_workflow"
scp "$dcmCode\*.R" "$dcmCode\*.sh" "$dcmCode\*.slurm" "$dcmCode\README.md" "meleg@ilogin.ibex.kaust.edu.sa:/ibex/user/meleg/DCM_MODEL/"

# This CTD location is from your existing workflow; adjust if your current CSV differs.
scp "K:\Gage Mele\chla\CTD_NOAA_WHOI_KAUST.csv" "meleg@ilogin.ibex.kaust.edu.sa:/ibex/user/meleg/DCM_MODEL/data/raw/"

# Change this example to the ACTUAL location of your exported MLD table.
$dcmMld = "C:\REPLACE_WITH_ACTUAL_FOLDER\mld_full.csv"
scp "$dcmMld" "meleg@ilogin.ibex.kaust.edu.sa:/ibex/user/meleg/DCM_MODEL/data/raw/"
```

The CTD needs `LocationCode`, `Chl`, `distance`, and original columns `DeepCTDID`, `Region.x`, `UTCDate.x`, `Depth_m`, `Temp_DegC` (standardized aliases are also accepted). The MLD table needs `LocationCode`, `Region`, `cast_id`, `Ocean`, `MLD_temp_MIX`. Ocean labels use `RedSea`/`Pacific`/`Other`; `Red Sea` is accepted as an alias. CTD dates must be ISO `yyyy-mm-dd` or begin with that ISO date. Confirm that `distance` uses the intended units: the original workflow treats it as kilometres, but the script does not guess or convert units.

Back on Ibex:

```bash
cd /ibex/user/meleg/DCM_MODEL

# Normalize Windows line endings in uploaded text files; raw CSVs are unchanged.
sed -i 's/\r$//' *.R *.sh *.slurm README.md
bash -n ibex_environment.sh
bash -n setup_environment.sh
bash -n run_stage.slurm

# If module avail shows another name, edit the R_MODULE default in this file.
nano ibex_environment.sh
```

## 3. One-time software setup on a compute node

```bash
cd /ibex/user/meleg/DCM_MODEL
# Allocate CPUs and RAM for package/CmdStan compilation.
srun --partition=batch --nodes=1 --ntasks=1 --cpus-per-task=4 \
  --mem=32G --time=04:00:00 --pty bash -l

# These commands run inside the compute allocation.
cd /ibex/user/meleg/DCM_MODEL
bash setup_environment.sh
exit
```

The setup uses CmdStan 2.37.0 as an explicit reproducibility choice; it is not a claim that this is the latest release. If you choose another supported version, set `CMDSTAN_VERSION` before setup and record it. Installed R package versions are recorded. This is a version inventory, not a complete lockfile; retain this environment for matched reruns. Use the same R/compiler module on setup and all jobs. No Windows binaries or Windows R libraries are copied to Ibex. If package downloading is blocked on compute nodes or a system library is missing, retain the installation error and use KSL's supported installation route rather than compiling on the login node.

## 4. Prepare full data and a pilot subset

If you used the earlier package to prepare data, rerun this stage from the original CTD CSV with a new `MODEL_DATA` filename (for example `full_unsmoothed.rds`). Then point the fit job at that new file. Previously smoothed/binned observations cannot be recovered from the old bundle; the updated fit script rejects bundles that do not explicitly record unsmoothed, unbinned preparation. Original observations may increase the number of model rows substantially, so reassess runtime and RAM using the new pilot.

```bash
cd /ibex/user/meleg/DCM_MODEL

# Submit full preparation using original observations without smoothing or binning.
sbatch --job-name=dcm-prep-full --cpus-per-task=1 --mem=32G --time=02:00:00 \
  --export=ALL,CAST_FRACTION=1,MODEL_DATA=/ibex/user/meleg/DCM_MODEL/data/processed/full.rds \
  run_stage.slurm prepare

# Prepare an independent 25% cast pilot, sampled once within each region.
sbatch --job-name=dcm-prep-pilot --cpus-per-task=1 --mem=32G --time=02:00:00 \
  --export=ALL,CAST_FRACTION=0.25,MODEL_DATA=/ibex/user/meleg/DCM_MODEL/data/processed/pilot.rds \
  run_stage.slurm prepare
```

Wait for preparation to complete. Inspect `logs/dcm-prep-*.out`, the matching `.err` logs, and these QC tables before fitting:

```bash
cat data/processed/full_qc/region_summary.csv
cat data/processed/full_qc/row_exclusions_by_region.csv
cat data/processed/full_qc/unmatched_ctd_casts.csv
```

The tensor's conservative support check requires at least 15 distinct depths, five distinct distances, and five distinct MLDs in every retained region. A small pilot can fail this check. Increase the pilot fraction and use a new output filename rather than silently dropping regions or changing k. For a region with inadequate support even in full data, review the model specification before proceeding. Full-dataset preparation may need more than 32 GB if the raw CSV is very large; inspect measured `MaxRSS` and adjust.

## 5. Run a prior predictive check and a posterior pilot

```bash
cd /ibex/user/meleg/DCM_MODEL

# Optional but useful: assess the retained priors on the actual pilot design.
sbatch --job-name=dcm-prior --time=2-00:00:00 \
  --export=ALL,FIT_MODE=prior,MODEL_DATA=/ibex/user/meleg/DCM_MODEL/data/processed/pilot.rds \
  run_stage.slurm fit

# Posterior pilot: 4 parallel chains; 2 threads each; 2,000 iterations,
# including 1,000 warmup iterations per chain; delta 0.97; tree depth 14.
sbatch --job-name=dcm-pilot --time=2-00:00:00 \
  --export=ALL,FIT_MODE=pilot,MODEL_DATA=/ibex/user/meleg/DCM_MODEL/data/processed/pilot.rds \
  run_stage.slurm fit
```

The default fit request is 8 CPUs and 64 GB on one node. These are initial requests, not measured resource requirements. More cores only help if the algorithm uses them; the code explicitly checks parallel chains times threads against the allocation. It is a shared-memory CPU workload and does not request a GPU or MPI ranks.

Monitor the job IDs printed by `sbatch`:

```bash
# Replace 123456 with your actual job ID.
squeue -u meleg
tail -f logs/dcm-pilot-123456.out
# Ctrl+C stops following the log; the submitted job continues.
sacct -j 123456 --format=JobID,JobName,State,Elapsed,AllocCPUS,ReqMem,MaxRSS
```

The pilot's run folder is `runs/pilot_JOBID`; the prior folder is `runs/prior_JOBID`. A job finishing successfully only confirms sampling completed. Run diagnostics next. If a job reaches its time or RAM limit, inspect logs and accounting before changing the request. An interrupted run may retain partial CSV files, but these scripts do not resume an unfinished NUTS chain; use a new run directory for a new run.

## 6. Check and plot the completed pilot

```bash
cd /ibex/user/meleg/DCM_MODEL
# Replace the example ID with the fit job ID, not this new checking job's ID.
sbatch --job-name=dcm-check --cpus-per-task=2 --mem=64G --time=08:00:00 \
  --export=ALL,RUN_DIR=/ibex/user/meleg/DCM_MODEL/runs/pilot_123456,PLOT_DRAWS=1000 \
  run_stage.slurm check
```

For a prior-only run, submit the same check stage with `RUN_DIR=.../runs/prior_JOBID`. It produces prior predictive density checks and avoids treating that run as a fitted posterior.

Review `checks_and_plots/diagnostic_status.txt`, `parameter_diagnostics.csv`, `flagged_parameters.csv`, and `sampler_summary.csv`. The screen uses R-hat <=1.01, bulk/tail ESS >=400, zero divergences and tree-depth hits, and E-BFMI >=0.3. Inspect trace plots and predictive/residual plots as well. A high in-sample R2 is not sufficient evidence of convergence or out-of-sample performance.

`posterior_epred()` describes uncertainty about an expected response. `posterior_predict()` also adds residual variability. Both types are exported separately. `re_formula=NA` sets group effects to zero; these probabilities therefore describe the specified reference groups, not full variability among new dates or ENSO groups. The Gaussian observation model can generate `chl_norm < -1`, implying a physically impossible negative chlorophyll ratio; the probability table records that fraction so you can evaluate the adequacy of the retained family.

Peak depth, half-maximum width and positive-enhancement-weighted spread are computed from reference expected profiles. A peak is flagged as subsurface only when positive enhancement peaks below 15 m; peaks at the grid boundary are flagged. Width uses the contiguous half-maximum component around the peak and is NA when truncated by a depth boundary. These are enhancement-profile summaries, not independent evidence of a chlorophyll/biomass DCM. Ocean curves average region profiles with equal region weights within each posterior draw; they are descriptive standardized summaries, not Ocean coefficients.

## 7. Submit production after the pilot is acceptable

```bash
cd /ibex/user/meleg/DCM_MODEL
# Start from 4 chains x 2 threads, 5,000 iterations with 2,500 warmup.
# Keep the current delta 0.97/tree depth 14 rather than loosening them on full data.
# Set memory using the pilot and full-data basis size; 128G is an initial example.
sbatch --job-name=dcm-full --cpus-per-task=8 --mem=128G --time=14-00:00:00 \
  --export=ALL,FIT_MODE=production,MODEL_DATA=/ibex/user/meleg/DCM_MODEL/data/processed/full.rds \
  run_stage.slurm fit
```

KAUST documents up to 14 days for ordinary associated Ibex users; unassociated accounts have tighter limits. Check your actual account limits. The historical notebook's roughly 10.9-day execution comment is not a measured runtime for this new three-way MLD tensor. This production model may need a different strategy if the pilot indicates it will exceed the available limit. Do not infer that requesting 14 days guarantees completion.

Once production finishes, submit checks again with `RUN_DIR=.../runs/production_JOBID`. Fit and checking jobs are deliberately separate so adjusting figures never reruns sampling. Submitted jobs continue after you close SSH.

## 8. Retrieve results

Run in Windows PowerShell:

```powershell
# Change the ID and destination. Download the WHOLE run folder, retaining Stan CSVs.
scp -r "meleg@ilogin.ibex.kaust.edu.sa:/ibex/user/meleg/DCM_MODEL/runs/production_123456" "C:\Users\gagem\Downloads\"
```

For local plotting in RStudio, install the same R packages, then run:

```r
# Use the path to your downloaded completed run.
Sys.setenv(RUN_DIR = "C:/Users/gagem/Downloads/production_123456", PLOT_DRAWS = "1000")
source("C:/Users/gagem/Downloads/DCM_Ibex_workflow/03_check_and_plot.R")
check_and_plot()  # Reads fit.rds and the matching model_data.rds; does not refit.
```

## Validation performed and remaining checks

All three R scripts were parsed successfully. Synthetic checks passed for original chlorophyll/depth retention, irregular observations, normalization, one-time sampling, the mgcv tensor basis and profile metrics. Additional tests passed for unknown-location matching, region/ocean aliases, conflicting MLD rejection, date ambiguity rejection and required cast-ID validation. Preparation was then executed on your actual local CSV copies: all 1,081 input cast identifiers matched MLD metadata before depth filtering; full preparation retained 258,984 observations from 1,080 casts, and the 25% pilot retained 64,283 observations from 267 casts across eight regions. All pilot regions passed the marginal predictor-support check, and retained chlorophyll/depth values exactly matched the original input. Your Ibex log confirms installation of R 4.5.0, brms 2.23.0 and CmdStan 2.37.0, with successful shell checks and R parsing. The Bayesian fit and production figures remain untested on the actual data. No Ibex job was submitted by the assistant.

For scientific validation beyond these implementation checks, use cast-blocked holdouts (rather than random depth-row splits), check residual structure within casts, examine prior and posterior predictions, and assess basis adequacy at the chosen k. The current formula retains date/ENSO grouping but has no cast-specific intercept or residual autocorrelation term. Whether this adequately captures within-cast dependence needs assessment on your actual data; this cleanup does not add those terms.

## Primary documentation

- KAUST filesystems: https://docs.hpc.kaust.edu.sa/systems/ibex/filesystem.html
- KAUST scheduling/storage policy: https://docs.hpc.kaust.edu.sa/policy/ibex.html
- KAUST software catalogue: https://docs.hpc.kaust.edu.sa/apps_catalogue/ibex.html
- Existing repository launchers: https://github.com/GageMele/pelagic-subsidies-chlorophyll/tree/main/scripts/ibex
- brms fitting: https://paulbuerkner.com/brms/reference/brm.html
- brms threading: https://paulbuerkner.com/brms/reference/threading.html
- brms expected predictions: https://paulbuerkner.com/brms/reference/posterior_epred.brmsfit.html
- brms predictive observations: https://paulbuerkner.com/brms/reference/posterior_predict.brmsfit.html
- CmdStan installation: https://mc-stan.org/cmdstanr/reference/install_cmdstan.html
- Stan diagnostic guidance: https://mc-stan.org/learn-stan/diagnostics-warnings.html
