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

The included targeted dependency patch is archived separately from the original targeted source. The onset repair contains original and patched checksum manifests: a mismatch against original_SHA256SUMS for a deliberately patched file is expected. The top-level repair SHA256SUMS identifies the delivered repair package.

The original shape package and embedded launcher predate the collector string-ID fix. Do not silently call them the patched live cluster version. Keep original bytes, record any patch, and inspect scheduler source-identity checks before resuming.

## Statistical issues still open

- Some fits have divergences despite otherwise good R-hat and ESS.
- Tree limits and diagnostic-screen thresholds differ across campaigns.
- Residual neighboring-depth correlation remains high.
- Gaussian support permits chl_norm below −1.
- High-k and tensor results are incomplete or, for older 3D trials, invalid.
- MLD shape improvements may concentrate in influential casts.
- Conditional gains need confirmation for unseen dates.
- Repeatedly inspected reserve data is no longer an untouched final test.
- Full-data fits cannot use former holdouts as independent validation.
- Point estimates from differing cohorts, depth domains, draw counts, or prediction modes are not interchangeable.

## Source and environment limitations

Original notebooks are interactive histories with alternative definitions, paths, and installation calls. A saved model is stronger evidence of its fitted settings than a comment. The recovered original preparation package and notebook are included here. Raw/prepared data and fitted objects remain external. Current code depends on prepared data and the existing environment.

No complete historical lockfile/container was supplied. Original source calibration, MLD algorithm, distance construction, ENSO index provenance, and redistribution permissions require upstream documentation.
