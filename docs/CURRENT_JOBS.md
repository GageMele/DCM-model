# Ibex jobs: verified snapshot and current-status boundary

**Live queue status on 27 October 2026 could not be obtained.** The fresh SSH attempt from this workspace failed with “connect to host ... port 22: Permission denied.” This is a connection failure, not a Slurm result. The repository cannot claim which jobs are running now.

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

The script prints the archive path. Return that snapshot to populate a verified current inventory. It includes job IDs/names, states, partitions, elapsed/limit/start times, CPUs, memory requests, nodes and reasons. It preserves all account jobs, with DCM matches in a separate file, and records query start/end UTC.

For a quick terminal list:

~~~bash
squeue --user meleg --array --states=RUNNING   --format="%.18i %.60j %.12T %.12M %.6C %N"
~~~

A successful empty query means no matching RUNNING jobs at capture time. A failed connection or failed squeue command does not establish an empty queue.

## Update policy

Keep this dated archive. Add each new live capture under a new timestamp and link it here. Do not relabel an old result as current, infer RUNNING from a fit registration alone, or remove completed/failed jobs from the scientific history.
