# Ibex jobs: verified capture and historical inventory

## Latest verified capture: 27 October 2026

The uploaded archive records the **entire visible queue for user meleg** at **2026-10-27 10:47:40 UTC / 13:47:40 Saudi time (UTC+3)**. The query ran from 10:47:38.260894 to 10:47:40.353886 UTC. This is a dated observation, not a continuously refreshed status page.

| Measure | Verified count |
|---|---:|
| All account queue entries | 335 |
| RUNNING jobs, entire account | 20 |
| RUNNING DCM fits, exact manifest matches | 19 |
| RUNNING jobs without a DCM manifest match | 1 |
| PENDING jobs, entire account | 315 |
| PENDING jobs matched to DCM | 0 |
| Visible DCM controllers / checking jobs | 0 / 0 |
| CPUs recorded for running DCM fits | 304 |
| CPUs recorded for all running account jobs | 400 |

All 19 DCM fits belong to `shape150_20261022T102851Z_3426481`. Each requests 16 CPUs and a seven-day walltime. Eighteen record 96G minimum memory; case 111 records 64G. These are scheduler resource fields, not measurements of CPU use or peak memory.

The [original complete CSV](../operations/jobs/20261027T104740Z/raw/all_queued_jobs.csv) lists **every one of the 335 entries**. The [original running CSV](../operations/jobs/20261027T104740Z/raw/running_jobs.csv) contains **all 20 running jobs**, with exact names, starts, elapsed time, limits, resources, nodes and reasons. The [readable original capture](../operations/jobs/20261027T104740Z/raw/README.md) lists all account entries.

## All 19 running DCM fits

Every row below has an exact job-ID association in the [captured manifest matches](../operations/jobs/20261027T104740Z/raw/dcm_manifest_matches.csv). Job names are `dcm-s150-77caaba437-fit_NNN`, where NNN is the three-digit case ID. All are RUNNING on partition `batch`, with one requested node and 16 CPUs each. Elapsed times are `days-hours:minutes:seconds`.

The basis column is `k_depth / k_mld` from the frozen plan. For the numeric MLD-varying equation, the generic `k_mld` column is present in the matrix but is not an independent tensor dimension.

| Case | Job ID | Equation | Training casts | k depth / MLD | Elapsed | Memory | Node |
|---|---|---|---|---|---|---|---|
| 111 | 54564846 | mld_vary_region | 20% | 13/4 | 1-12:52:45 | 64G | cn512-24-l |
| 142 | 54566189 | tensor_region | 10% | 7/4 | 1-12:13:11 | 96G | cn512-03-l |
| 143 | 54567669 | tensor_region | 10% | 11/4 | 1-11:28:38 | 96G | cn512-04-l |
| 144 | 54577971 | tensor_region | 10% | 15/4 | 1-04:19:36 | 96G | cn512-21-r |
| 189 | 54250973 | tensor_region | 20% | 7/5 | 4-22:59:41 | 96G | cn509-28-r |
| 190 | 54250974 | tensor_region | 20% | 7/7 | 4-22:45:36 | 96G | cn512-22-l |
| 193 | 54250975 | tensor_halfshift_region | 20% | 7/5 | 4-22:38:46 | 96G | cn512-14-r |
| 194 | 54250976 | tensor_halfshift_region | 20% | 7/7 | 4-22:38:46 | 96G | cn509-23-r |
| 211 | 54251087 | tensor_region | 20% | 7/4 | 4-19:45:46 | 96G | cn512-22-r |
| 212 | 54251089 | tensor_halfshift_region | 20% | 7/4 | 4-19:45:46 | 96G | cn512-20-r |
| 246 | 54407221 | tensor_halfshift_region | 20% | 7/4 | 3-16:41:52 | 96G | cn512-27-r |
| 247 | 54412325 | tensor_halfshift_region | 20% | 7/4 | 3-16:09:18 | 96G | cn512-02-r |
| 248 | 54414270 | tensor_halfshift_region | 20% | 7/4 | 3-15:50:16 | 96G | cn509-05-l |
| 271 | 54443754 | tensor_region | 20% | 7/4 | 3-08:30:51 | 96G | cn512-24-r |
| 272 | 54444266 | tensor_shift_region | 20% | 7/4 | 3-08:11:53 | 96G | cn512-24-l |
| 273 | 54445934 | tensor_halfshift_region | 20% | 7/4 | 3-07:48:55 | 96G | cn509-14-r |
| 280 | 54474462 | tensor_region | 20% | 7/4 | 2-23:54:10 | 96G | cn512-20-r |
| 281 | 54486850 | tensor_shift_region | 20% | 7/4 | 2-21:53:40 | 96G | cn512-04-r |
| 282 | 54489283 | tensor_halfshift_region | 20% | 7/4 | 2-19:48:06 | 96G | cn509-14-l |

The [enriched DCM CSV](../operations/jobs/20261027T104740Z/dcm_running_with_settings.csv) preserves exact job names, full labels, manifest/run paths, all Slurm fields, and every frozen case setting. Its `plan_*` fields come from the [archived case matrix](../evidence/shape_20261025/campaign/cases.tsv); they are documented settings, not newly inspected run objects.

### Settings and scientific comparisons in progress

All 19 planned fits use four chains, four parallel cores, four threads per chain, sampler seed 1234, `max_treedepth=14`, `custom_z` priors, Gaussian observations, regional residual SD, and the 150 m maximum depth. `iter` includes warmup; retained draws are per chain.

| Case | Comparison | Training selection | iter | Warmup | Retained / chain | adapt_delta |
|---|---|---|---|---|---|---|
| 111 | intermediate_k | primary / frozen | 3000 | 1000 | 2000 | 0.97 |
| 142 | mld_surface | primary / frozen | 3000 | 1000 | 2000 | 0.97 |
| 143 | mld_surface | primary / frozen | 3000 | 1000 | 2000 | 0.97 |
| 144 | mld_surface | primary / frozen | 3000 | 1000 | 2000 | 0.97 |
| 189 | mld_basis | primary / frozen | 3000 | 1000 | 2000 | 0.97 |
| 190 | mld_basis | primary / frozen | 3000 | 1000 | 2000 | 0.97 |
| 193 | mld_basis | primary / frozen | 3000 | 1000 | 2000 | 0.97 |
| 194 | mld_basis | primary / frozen | 3000 | 1000 | 2000 | 0.97 |
| 211 | cast_selection | primary / selection seed 20261024 | 3000 | 1000 | 2000 | 0.97 |
| 212 | cast_selection | primary / selection seed 20261024 | 3000 | 1000 | 2000 | 0.97 |
| 246 | delta99 | primary / frozen | 3000 | 1000 | 2000 | 0.99 |
| 247 | warmup2000 | primary / frozen | 4000 | 2000 | 2000 | 0.97 |
| 248 | samples4000 | primary / frozen | 5000 | 1000 | 4000 | 0.97 |
| 271 | date_blocked | date_blocked / selection seed 20261022 | 3000 | 1000 | 2000 | 0.97 |
| 272 | date_blocked | date_blocked / selection seed 20261022 | 3000 | 1000 | 2000 | 0.97 |
| 273 | date_blocked | date_blocked / selection seed 20261022 | 3000 | 1000 | 2000 | 0.97 |
| 280 | date_blocked | date_blocked / selection seed 20261023 | 3000 | 1000 | 2000 | 0.97 |
| 281 | date_blocked | date_blocked / selection seed 20261023 | 3000 | 1000 | 2000 | 0.97 |
| 282 | date_blocked | date_blocked / selection seed 20261023 | 3000 | 1000 | 2000 | 0.97 |

For frozen training cohorts, `selection_seed=0` is a plan sentinel; the saved membership defines the data. It does not mean those casts were newly drawn with seed zero.

Eight running fits use `tensor_region`, eight use `tensor_halfshift_region`, two use `tensor_shift_region`, and one uses `mld_vary_region`. The tensor models allow the profile shape to vary with MLD; using a shifted coordinate inside a tensor still allows shape change. Three fits train on 10% of casts and sixteen on 20%. No 40-100% DCM fit is visible in this capture; their completion/submission states require the refreshed campaign records. See [equations](EQUATIONS.md) for the exact implementations.

The longest elapsed fit is case 189 at **4 days, 22:59:41**. Elapsed job time includes all work within the allocation; this snapshot does not report chain iterations, diagnostics, fit accuracy, or completed sampling time.

## The other running job

| Job ID | Exact job name | State | CPUs | Memory | Elapsed | Time limit | Node |
|---|---|---|---|---|---|---|---|
| 53515610 | C_D1_202202_run | RUNNING | 96 | 64G | 7:26:55 | 1-00:00:00 | cn504-05 |

This job has no captured DCM manifest association and does not use the DCM job-name prefix. Its project ownership is not established by this archive; it is retained to satisfy the complete account inventory.

## All pending jobs

All 315 pending entries lack a captured DCM manifest match and have non-DCM job names. The [pending CSV](../operations/jobs/20261027T104740Z/pending_jobs.csv) retains every ID and field.

| Recorded pending reason | Count |
|---|---|
| (Dependency) | 310 |
| (DependencyNeverSatisfied) | 4 |
| (Priority) | 1 |

Name groups provide a descriptive cross-check, not a project-ownership claim:

| Name prefix or exact name | RUNNING | PENDING |
|---|---|---|
| C_D1_ | 1 | 136 |
| VR_ | 0 | 5 |
| V_ | 0 | 173 |
| dcm- | 19 | 0 |
| testB_spinup | 0 | 1 |

The four entries reporting `DependencyNeverSatisfied` are:

| Job ID | Exact captured name |
|---|---|
| 53593830 | V_compare_D0_2021-02_53593829_5359318 |
| 53593828 | V_compare_D0_2021-01_53593827_5359318 |
| 53593270 | V_compare_D0_2021-02_53593269_5359318 |
| 53593265 | V_compare_D0_2021-01_53593264_5359318 |

The one priority-pending job is `54656619`, `testB_spinup`. No dependency graph or accounting history was supplied, so this snapshot does not establish the cause of the four unsatisfied dependencies.

## Operational interpretation and next checks

**No DCM controller or check job appears in the captured queue.** The saved campaign dispatcher was designed to submit checks and further fits as slots became available. Inspect its latest controller ledger, accounting state and logs before assuming the remaining cases will launch automatically. The snapshot does not establish a controller failure, its cause, or whether another controller appeared after capture.

Use the [runbook](IBEX_RUNBOOK.md) to read the current queue, controller accounting, saved dispatcher timestamp, fit/check markers, and latest manifest. A queue entry disappearing does not prove scientific completion. In particular, the October 25 count of 79 complete cases must not be replaced with `282 - 19`.

The repository's scientific results remain the **October 25** export: 67 completed new posterior fits and 12 prior checks. This queue upload contains no new RMSE, onset scores, R-hat, ESS, divergences, or treedepth results. Refresh and export existing results to update those comparisons.

## Evidence and audit

All six original archive members are retained byte for byte under the [dated snapshot](../operations/jobs/20261027T104740Z/README.md). The [archive metadata](../operations/jobs/20261027T104740Z/raw/metadata.json) matches the earlier [reported terminal metadata](../operations/jobs/20261027T104740Z/reported_capture_metadata.json). The [validation report](../operations/jobs/20261027T104740Z/validation.json) verifies unique IDs, raw/CSV agreement, running subset, manifest associations, and case labels/resources. [Provenance](../operations/jobs/20261027T104740Z/provenance.json) records archive and member hashes.

The capture expands array tasks; IDs remain strings. Matching searches the capture script's known manifest locations, so an unmatched job is not automatically proof of another project's ownership. In this snapshot all `dcm-` names exactly match the 19 manifest-associated fits. `start_time` is preserved as Slurm emitted it without a timezone offset; only the query timestamps explicitly record UTC. For running jobs, the captured `%R` reason field contains the node list.

No job was submitted, cancelled, or changed to assemble this inventory.

## Last verified shape dispatcher snapshot

At 2026-10-25T11:14:41.409930Z the dispatcher recorded 79 complete, 48 fitting, and 155 not submitted, with 768 requested CPUs. The [archived dispatcher status](../evidence/shape_20261025/campaign/dispatch_status.json) and [all 282 job registrations](../operations/jobs/20261025/shape_all_282_job_registrations.csv) are included.

The table below lists all 48 accepted shape fitting jobs whose checks were incomplete in that archive, matching the dispatcher count. **This is not a fresh squeue RUNNING list:** accepted fitting jobs may have been pending or running, and may have finished since. Other DCM campaigns and unrelated jobs on the account are not inferable from this one dispatcher.

| Case | Accepted fit job | CPUs | Label |
|---:|---|---:|---|
| 14 | 54508550 | 16 | s014_size_basis_mld_vary_shared_p10_k11_uniform |
| 15 | 54508552 | 16 | s015_size_basis_mld_vary_shared_p10_k15_uniform |
| 16 | 54509510 | 16 | s016_size_basis_mld_vary_region_p10_k7_uniform |
| 17 | 54510248 | 16 | s017_size_basis_mld_vary_region_p10_k11_uniform |
| 18 | 54510297 | 16 | s018_size_basis_mld_vary_region_p10_k15_uniform |
| 20 | 54512543 | 16 | s020_size_basis_original_z_p20_k11_uniform |
| 21 | 54512583 | 16 | s021_size_basis_original_z_p20_k15_uniform |
| 23 | 54512589 | 16 | s023_size_basis_interaction_z_p20_k11_uniform |
| 24 | 54513421 | 16 | s024_size_basis_interaction_z_p20_k15_uniform |
| 26 | 54514245 | 16 | s026_size_basis_mld_halfshift_p20_k11_uniform |
| 27 | 54514502 | 16 | s027_size_basis_mld_halfshift_p20_k15_uniform |
| 29 | 54515138 | 16 | s029_size_basis_mld_shift_p20_k11_uniform |
| 30 | 54515183 | 16 | s030_size_basis_mld_shift_p20_k15_uniform |
| 32 | 54515327 | 16 | s032_size_basis_mld_vary_shared_p20_k11_uniform |
| 33 | 54516249 | 16 | s033_size_basis_mld_vary_shared_p20_k15_uniform |
| 35 | 54517221 | 16 | s035_size_basis_mld_vary_region_p20_k11_uniform |
| 36 | 54523817 | 16 | s036_size_basis_mld_vary_region_p20_k15_uniform |
| 100 | 54528646 | 16 | s100_intermediate_k_original_z_p20_k9_uniform |
| 151 | 54250970 | 16 | s151_mld_surface_tensor_region_p20_k7_uniform |
| 154 | 54250971 | 16 | s154_mld_surface_tensor_shift_region_p20_k7_uniform |
| 157 | 54250972 | 16 | s157_mld_surface_tensor_halfshift_region_p20_k7_uniform |
| 189 | 54250973 | 16 | s189_mld_basis_tensor_region_p20_k7_uniform |
| 190 | 54250974 | 16 | s190_mld_basis_tensor_region_p20_k7_uniform |
| 193 | 54250975 | 16 | s193_mld_basis_tensor_halfshift_region_p20_k7_uniform |
| 194 | 54250976 | 16 | s194_mld_basis_tensor_halfshift_region_p20_k7_uniform |
| 203 | 54251071 | 16 | s203_cast_selection_tensor_region_p20_k7_uniform |
| 204 | 54251074 | 16 | s204_cast_selection_tensor_halfshift_region_p20_k7_uniform |
| 211 | 54251087 | 16 | s211_cast_selection_tensor_region_p20_k7_uniform |
| 212 | 54251089 | 16 | s212_cast_selection_tensor_halfshift_region_p20_k7_uniform |
| 219 | 54251604 | 16 | s219_cast_selection_tensor_region_p20_k7_uniform |
| 220 | 54251605 | 16 | s220_cast_selection_tensor_halfshift_region_p20_k7_uniform |
| 242 | 54382620 | 16 | s242_delta99_mld_vary_region_p20_k7_uniform |
| 244 | 54395444 | 16 | s244_samples4000_mld_vary_region_p20_k7_uniform |
| 245 | 54399246 | 16 | s245_delta95_tensor_halfshift_region_p20_k7_uniform |
| 246 | 54407221 | 16 | s246_delta99_tensor_halfshift_region_p20_k7_uniform |
| 247 | 54412325 | 16 | s247_warmup2000_tensor_halfshift_region_p20_k7_uniform |
| 248 | 54414270 | 16 | s248_samples4000_tensor_halfshift_region_p20_k7_uniform |
| 251 | 54425774 | 16 | s251_prior_sensitivity_tensor_halfshift_region_p20_k7_uniform |
| 252 | 54425805 | 16 | s252_prior_sensitivity_tensor_halfshift_region_p20_k7_uniform |
| 270 | 54442291 | 16 | s270_date_blocked_mld_vary_region_p20_k7_uniform |
| 271 | 54443754 | 16 | s271_date_blocked_tensor_region_p20_k7_uniform |
| 272 | 54444266 | 16 | s272_date_blocked_tensor_shift_region_p20_k7_uniform |
| 273 | 54445934 | 16 | s273_date_blocked_tensor_halfshift_region_p20_k7_uniform |
| 275 | 54447135 | 16 | s275_date_blocked_interaction_z_p20_k7_uniform |
| 279 | 54469872 | 16 | s279_date_blocked_mld_vary_region_p20_k7_uniform |
| 280 | 54474462 | 16 | s280_date_blocked_tensor_region_p20_k7_uniform |
| 281 | 54486850 | 16 | s281_date_blocked_tensor_shift_region_p20_k7_uniform |
| 282 | 54489283 | 16 | s282_date_blocked_tensor_halfshift_region_p20_k7_uniform |

## Capture every currently running and queued job

The [read-only capture script](../scripts/capture_ibex_jobs.py) records all visible queue states for the authenticated user, expands array tasks, saves a separate RUNNING-only CSV, and maps known DCM manifest registrations. It writes timestamped reports and a small archive; it does not submit, cancel or update jobs.

On an Ibex login node, from a clone of this repository:

~~~bash
python3 scripts/capture_ibex_jobs.py   --user meleg   --project-root /ibex/user/meleg/DCM_MODEL
~~~

If this repository is not yet on Ibex, download only the capture script and run it:

~~~bash
cd /ibex/user/meleg/DCM_MODEL
curl -fL https://raw.githubusercontent.com/GageMele/DCM-model/main/scripts/capture_ibex_jobs.py   -o logs/capture_ibex_jobs.py
python3 logs/capture_ibex_jobs.py --user meleg --project-root "$PWD"
~~~

The script prints the archive path. Preserve each new snapshot under a new timestamp and update this inventory from its contents. It includes job IDs/names, states, partitions, elapsed/limit/start times, CPUs, memory requests, nodes and reasons. It preserves all account jobs, with DCM matches in a separate file, and records query start/end UTC.

For a quick terminal list:

~~~bash
squeue --user meleg --array --states=RUNNING   --format="%.18i %.60j %.12T %.12M %.6C %N"
~~~

A successful empty query means no matching RUNNING jobs at capture time. A failed connection or failed squeue command does not establish an empty queue.

## Update policy

Keep this dated archive. Add each new live capture under a new timestamp and link it here. Do not relabel an old result as current, infer RUNNING from a fit registration alone, or remove completed/failed jobs from the scientific history.
