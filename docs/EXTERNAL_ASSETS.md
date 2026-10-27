# External assets needed to continue

This repository contains source, configurations, compact summaries and selected evidence. Fresh fitting and full numerical reproduction require external research assets.

## Data and fitted objects

| Asset | Location or identity | Purpose |
|---|---|---|
| Raw CTD CSV | `DCM_MODEL/data/raw/CTD_NOAA_WHOI_KAUST.csv` | Rebuild preparation from original observations |
| Raw MLD CSV | `DCM_MODEL/data/raw/mld_full.csv` | Cast-level MLD join; obtain upstream derivation |
| Full prepared 250 m bundle | Frozen `full_input_250m.rds` or audited processed equivalent | Preserve preparation and eligible cast pool |
| Shape source bundle | `followup150_20261018T150024Z_3898749/full_input_250m.rds` | MD5 `cb5348ec5d9474d3f10570508e06d70d` in the latest source audit |
| Frozen split RDS files | Each `experiment/data/` | Exact observations and row ordering |
| Factor levels and scaling | `region_levels.rds`, `scaling.rds` | Saved factor order and training transforms |
| Fitted object and input | Each run's `fit.rds`, `model_data.rds`, `run_config.rds` | Posterior, model specification and exact training data |
| CmdStan parameter CSVs | Each run's `stan_csv/*.csv` | Persistent posterior and saved warmup draws |
| Original saved fit | `my_model_50_FINAL.rds`, `fit_brms_model_50_FINAL.RData` | Historical pre-MLD model; not two independent fits |
| Full prediction tables | Run `checks/` directories | Detailed numerical reproduction |

Obtain the actual files from the project owner or authorized Ibex storage and verify their identities. Do not silently substitute a newly prepared dataset into an existing frozen campaign.

## Archives

The latest reviewed source is `dcm_shape_comparison_20261025T111507Z_160796.tar.gz`, size 77,804,086 bytes. Its SHA256 and the review ZIP checksum are recorded in [archive_checksums.json](../provenance/archive_checksums.json).

The archive contains current shape checks and combined legacy summaries. It excludes full fitted objects, parameter chains, prepared input data and complete older-campaign prediction evidence.

The extracted contents of `DCM_review_20261025.zip` (918,260 bytes) are in [results/20261025](../results/20261025/). Its analysis script accepts an extracted source experiment and an output directory; it needs Python, numpy, pandas and matplotlib. It does not fit or submit models. Full reproduction needs the source archive, not just the compact review tables.

Earlier October 16/18/19/21 snapshots support historical analyses. Selected reports, findings, contrasts and source code are included. Full spreadsheet/HTML deliveries and large cast tables remain external where absent from the inventory.

## Software and current state

Retain the existing R user library, CmdStan path/version, environment script and per-run version inventories. Some packages are development builds. Installing today's CRAN packages does not reproduce their historical versions.

Live manifests and scheduler records are also needed to continue safely. Accepted job IDs and source-fit mappings prevent duplicate submissions. The shape export excluded shell and Python source, so this repository retains delivered scheduler/worker/environment code alongside directly exported R/configuration evidence.

## Permissions and missing methods

The public repository does not establish permission to redistribute all underlying observations. Licence selection and data permissions remain with the owner. Ibex requires separate authorized access.

Acquire upstream chlorophyll calibration, MLD algorithm and thresholds, distance derivation, and cruise/station provenance before writing final methods. The historical ENSO lookup is preserved for reproducibility, not presented as a newly validated monthly index.

## Before a new production run

Confirm input identity, response definition, original observation retention, cast membership, depth range, factor levels, saved scaling, formula, priors, software versions, sampler settings and seed, prediction target, compute allocation, and output paths. Explicitly distinguish reconstruction from held-out evaluation.
