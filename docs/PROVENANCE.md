# Provenance and reproducibility

## Evidence layers

This handover was assembled on 27 October 2026 from the user's project workspace, original notebooks, delivered campaign packages, earlier analyses, the October 25 shape comparison archive, and its audited review. The scientific status cutoff is the archive timestamp, not the handover date.

The [source inventory](../provenance/source_inventory.csv) records repository path, original relative source path, byte count, SHA256 and role for every copied file. The [archive checksums](../provenance/archive_checksums.json) identify the source archive and review ZIP. The [repository checksums](../provenance/repository_checksums.json) verify the delivered tree, excluding the hash list itself and its generated validation report.

Paths beginning shape/ or review/ in the inventory refer to archive extractions. Other source paths refer to the local project mirror. These labels preserve provenance and are not runtime path requirements.

## Delivered versus executed source

The DCM_Ibex_* directories preserve original deliveries and their validation notes/checksums. The evidence/shape_20261025/campaign directory preserves the actual cluster export, including the corrected collector. These are intentionally separate. A later repair can make an executed file differ from its original package checksum.

Root queue/repair scripts are self-contained deployment artifacts, often with embedded payloads. Editing an extracted package does not rewrite its launcher. Archived builders preserve how those packages were constructed; many contain original absolute paths and must be adapted before reuse.

The source archive excludes data/, stan_csv/, fit objects, shell/Python source and certain large files. Its R/configuration files are direct executed evidence; the retained scheduler, worker and environment source comes from delivered packages. Do not describe the compact export as a complete filesystem backup.

## Numerical validation already performed

- The October 21 review independently recomputed 1,136 profile RMSE summaries and 360 onset aggregates. Its optional all-pairs audit was not completed; saved all-pair tables are not described as independently verified.
- The October 25 review independently recomputed 2,412 aggregate values, with a maximum discrepancy of 5.51e-14, and 344 planned paired contrasts. It verified identical holdout membership across four 20% selections.
- The original saved fit was inspected without refitting. Hydrated posterior draws and temporary prediction adapters supported historical plots.
- This repository handover validates file integrity, local documentation links, and case/registry counts. It does not rerun the statistical models or claim a new full posterior audit.

See [validation.json](../provenance/validation.json) for the checks run on this repository.

## Status and deduplication

The 586 registrations are not 586 completed fits. Reused controls and prediction-only evaluations need source-fit identity. Numeric case IDs repeat across campaigns. The stable key combines campaign name and case identifier. Nonnumeric threading IDs must remain strings: coercing them to integers caused the documented 583-to-586 correction.

A saved marker, Slurm state and diagnostic flag answer different questions. Preserve all three where available. Snapshot collection can happen while files are changing; retain the timestamp and archive checksum.

## Reproducibility limits

Exact upstream MLD derivation, raw CTD calibration and final data permissions are not established by the campaign packages. Earlier interactive notebook execution cannot be inferred solely from the order of chunks. The original saved fit lacks some cast/date information. Historical package versions are recorded but not completely locked.

Some current limitations are scientific rather than software bugs: regional heterogeneity, residual dependence, physical predictive support, incomplete generalization evidence, and censoring. Do not change source outputs to conceal them.

## External methods references

The handover checked these official references on 27 October 2026:

- brms smooth terms: https://paulbuerkner.com/brms/reference/s.html
- mgcv tensor smooths: https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/t2.html
- Stan diagnostics: https://mc-stan.org/learn-stan/diagnostics-warnings.html
- Slurm submission/dependencies: https://slurm.schedmd.com/sbatch.html

Other primary links in archived packages are retained as historical references. Current online versions do not establish the version used in an earlier fit. Obtain formal software citations from the installed packages and journal requirements before publication.

## Updating the record

Append dated snapshots and decisions. Do not replace a historical result CSV with newer rows under the same date. Record corrections alongside old versions. Changes to the response, MLD construction, depth range, excluded casts, scaling, priors or prediction target must remain visible as changes to the scientific comparison.

No copyright or licence grant has been invented. Source ownership, data permissions and publication credits need the owner's review before a formal release.
