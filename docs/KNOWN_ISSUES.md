# Known issues, repairs, and interpretation traps

| Issue | Evidence | Resolution/boundary |
|---|---|---|
| Empty fast report | same_cpu_model assignment of one value to zero rows; job 53779690 | Guard zero-row tables/vector lengths; do not refit |
| Wrong monitor path | baseline_tests/monitor_all_dcm_v4.sh missing | v4 belonged to extensive; current shape uses v8 |
| Interrupted targeted dependencies | 67 of 117 old-fit evaluations registered | Repair added 50 evaluations and collector, preserving 63 fit/check pairs |
| Onset checks all failed | “Insufficient vertical support” in reused and new fits | Singleton-depth support handling; rerun checks from saved fits |
| All-campaign count 583 | integer coercion of C4_T0/T2/T4/T8 | Preserve string IDs, restore 586 distinct entries |
| Original rstan DLL unavailable locally | Saved final fit could not load its old DLL | Read hydrated sample slots and temporary in-memory prediction adapter; no fit/data mutation |
| SSH disconnect | client_loop connection reset | Reconnect and query Slurm; not evidence jobs died |
| ZIP download difficulty | User could not retrieve an earlier deliverable | Repository preserves readable source; named large assets remain separately transferable |
| shape150 controller died 26 Oct | `squeue --me` failed during an Ibex Slurm outage; the scheduler treats any failed Slurm query as fatal | New controller_002 (job 54668377) on 27 Oct. Making Slurm query failures retryable would prevent a repeat |
| Account at `QOSMaxCpuPerUserLimit` | Many DCM fits plus other projects' jobs (V_, C_, S5*, D0B, D1B, V5_) | New fits wait for slots. The non-DCM jobs belong to other projects; leave them alone |
| Population predictions dropped a group-level sigma | `(1 \| depth_bin)` in sigma is removed by `re_formula` | round3/bins1m use fixed `sigma ~ 0 + Region + depth_bin` |
| ggplot `geom_errorbarh` and `if/else` precedence in plot code | Local tests of round3 plotting | Fixed before launch |
| Lognormal ratio ≤ 0 | 10 / 4 rows with Chl ≤ 0 | Floored at ratio 0.01 for the likelihood only |
| 48 shape150 fit failures (by 12 Nov) | `FIT_FAILED.txt` count from 08_status_all.sh | Not yet diagnosed |
| next150 peak profile, date-blocked (case 14) | R-hat 1.68, 948 divergences, 53.9 h | Not converged; do not interpret. Gaussian peak (round3 15/16) also unconverged (R-hat 1.32 / 1.10) |
| Boundary DCM maxima | 4 casts peak at 150 m, 2 at 10 m | Flag or exclude when interpreting DCM-depth models |

The included targeted dependency patch is archived separately from the original targeted source. The onset repair contains original and patched checksum manifests: a mismatch against original_SHA256SUMS for a deliberately patched file is expected. The top-level repair SHA256SUMS identifies the delivered repair package.

The original shape package and embedded launcher predate the collector string-ID fix. Do not silently call them the patched live cluster version. Keep original bytes, record any patch, and inspect scheduler source-identity checks before resuming.

## Statistical issues still open

- Some fits have divergences despite otherwise good R-hat and ESS.
- Tree limits and diagnostic-screen thresholds differ across campaigns.
- Residual neighboring-depth correlation remains high.
- Gaussian support permits chl_norm below −1. The lognormal models fix this.
- Lognormal intervals are too narrow: 90% coverage is 0.80–0.86 overall and 0.63–0.78 below 60 m on unseen dates.
- Strong DCM peaks are under-predicted by about half on new dates, because between-date amplitude variance is large.
- Elpd differences do not yet have standard errors. Totals can be compared only within one split. Binned-fit lpd is not like-for-like with raw-row fits.
- High-k and tensor results are incomplete or, for older 3D trials, invalid.
- MLD shape improvements may concentrate in influential casts.
- Conditional gains need confirmation for unseen dates.
- Repeatedly inspected reserve data is no longer an untouched final test.
- Full-data fits cannot use former holdouts as independent validation.
- Point estimates from differing cohorts, depth domains, draw counts, or prediction modes are not interchangeable.

## Source and environment limitations

Original notebooks are interactive histories with alternative definitions, paths, and installation calls. A saved model is stronger evidence of its fitted settings than a comment. The recovered original preparation package and notebook are included here. Raw/prepared data and fitted objects remain external. Current code depends on prepared data and the existing environment.

No complete historical lockfile/container was supplied. Original source calibration, MLD algorithm, distance construction, ENSO index provenance, and redistribution permissions require upstream documentation.
