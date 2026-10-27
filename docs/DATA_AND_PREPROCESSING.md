# Data, response and preprocessing

## Authoritative current preparation

Read [01_prepare_data.R](../DCM_Ibex_workflow/01_prepare_data.R) and the [workflow README](../DCM_Ibex_workflow/README.md). Later campaigns freeze its prepared bundle and explicitly verify that observations are unchanged. The current raw-data files are expected at:

- `/ibex/user/meleg/DCM_MODEL/data/raw/CTD_NOAA_WHOI_KAUST.csv`
- `/ibex/user/meleg/DCM_MODEL/data/raw/mld_full.csv`

The original CTD upload example used `K:\Gage Mele\chla\CTD_NOAA_WHOI_KAUST.csv` on Windows. That is a historical workstation path, not a portable input dependency.

## Field dictionary

| Field | Meaning and treatment |
|---|---|
| `Chl` | Original chlorophyll observation; preserve input values. Absolute calibration/units need upstream confirmation. |
| `depth` / `Depth_m` | Depth in metres; original observation depths retained. |
| `surface_chl` | Within-cast median finite raw Chl from 5–15 m, computed before model-depth filtering. |
| `chl_norm` | `(Chl - surface_chl) / surface_chl`, fractional enhancement. |
| `distance` | Positive distance from shore; notebook treats units as km. No automatic unit conversion. |
| `log10_distance` | `log10(distance)`; nonfinite values excluded explicitly. |
| `MLD_temp_MIX` | Cast-level mixed-layer depth in metres from the external MLD table. Derivation is not performed by these campaign scripts. |
| `Region` | Factor retaining the verified source labels. `Hawiian` is the historical canonical spelling. |
| `Ocean` | Basin grouping label, including `RedSea` and `Pacific`. |
| `date` / `date_raw` | ISO date group label in the current workflow; not a continuous time covariate. |
| `enso` | Historical year lookup: 2009 el_nino; 2010/2011 la_nina; otherwise neutral. Not a newly validated ENSO index. |
| `cast_id` | Source cast identifier; alone may not be globally unique. |
| `cast_key` | Compound identity distinguishing location, region, ocean, cast and date. Preserve exactly. |
| `LocationCode` | Join metadata; unknown location values standardized. |
| `mld_z` | MLD centered/scaled using saved training-cast moments. |
| `distance_z` | Log 10-distance centered/scaled using saved training-cast moments. |
| `mld_distance_z` | Product `mld_z * distance_z`; not standardization of a raw product. |
| `depth_mld` | `depth - MLD_temp_MIX`. |
| `depth_mld_half` | `depth - 0.5 * MLD_temp_MIX`. |

The eight-region inventory is recorded in the [frozen cast metadata](../evidence/shape_20261025/campaign/full_cast_metadata.csv). Do not silently correct spelling in an existing fit: factor levels and model matrices depend on it.

## Response interpretation

| chl_norm | Chlorophyll relative to surface | Enhancement |
|---:|---:|---:|
| -1 | 0 times surface | -100% |
| -0.5 | 0.5 times surface | -50% |
| 0 | 1 times surface | 0% |
| 0.25 | 1.25 times surface | +25% |
| 0.5 | 1.5 times surface | +50% |
| 1 | 2 times surface | +100% |

Add 1 for a ratio; multiply enhancement by 100 for percent. An exceedance probability `P(chl_norm > 0.5)` is not the probability of falling between 1.5 and 2 times surface concentration. Negative normalized observations are permitted. Gaussian draws below -1 imply negative chlorophyll on this scale and are recorded as a model-support diagnostic.

Because the baseline uses a surface median, individual observations at 10 m are not forced to exactly zero. Neither a Region intercept nor centering a spline imposes a zero-value boundary condition.

## Join and quality rules

1. Read CTD aliases: cast from DeepCTDID/cast_id, Region from Region.x/Region, date from UTCDate.x/date, depth from Depth_m/depth, and temperature from Temp_DegC/temperature.
2. Exclude KAUST according to active notebook logic; retain Equatorial despite an older contrary comment.
3. Standardize metadata aliases: Hawaiian → Hawiian, Red Sea → Red Sea, and missing/blank/Not specified location codes to a common unknown-location representation.
4. Join MLD using LocationCode, Region, cast_id and Ocean. Unknown locations only match within the other identifying fields. Exact duplicate right-hand rows collapse; conflicting MLD values for a key stop preparation.
5. Verify required identity fields and ambiguity across dates. If one four-column join key represents multiple dates, stop rather than silently merge casts.
6. Compute the 5–15 m baseline from finite raw observations before restricting depth.
7. Exclude invalid/nonpositive surface baselines, depths outside 10–250 m, nonfinite response, invalid log-distance and missing/invalid MLD. Base preparation accepts zero MLD; the later ratio campaign requires strictly positive MLD and stops on nonpositive values.
8. Sort by cast/depth; preserve original chlorophyll/depth values; write unmatched-cast, unused-MLD and exclusion QC.
9. Sample entire casts once within Region. Do not successively sample a subset and mislabel the retained fraction.

The prepared 250 m source contains 258,984 observations from 1,080 casts after quality/depth filtering. Current shape fitting restricts that eligible source to 10–150 m, giving 151,962 observations. These counts differ from the original saved fit's 255,369 rows; that older fit is not the identical prepared dataset.

## No fitting-data smoothing or binning

The user explicitly required unsmoothed, unbinned inputs in October 2026. Original values/depths are retained. The surface median is a normalization baseline, not profile smoothing. Five-metre median operations used later for onset/sharpness **scoring** do not replace model-fitting observations.

Historical notebooks contain rolling medians, five-metre bins, alternative windows and sequential sample blocks. Their presence does not establish which blocks preceded a saved interactive fit. Do not execute those notebooks from beginning to end as the current pipeline.

## Exact frozen shape datasets

| Dataset | Casts | Observations | Role |
|---|---:|---:|---|
| Primary 10% | 108 | 15,228 | Training |
| Primary 20% | 216 | 30,455 | Training |
| Primary 40% | 432 | 60,911 | Training |
| Primary 50% | 540 | 76,139 | Training |
| Primary 75% | 810 | 114,136 | Training |
| Primary 100% | 1,080 | 151,962 | Full training; no independent holdout |
| Each of three alternate 20% selections | 216 | 30,456 | Seeds 20261023/24/25 |
| Each of two date-blocked 20% selections | 216 | 30,456 | Seeds 20261022/23 |
| Primary validation | 108 | 15,226 | Fixed development holdout |
| Date-blocked validation | 117 | 16,497 | Fixed unseen nested-date cohort |
| Reserve | 108 | See frozen membership | Separate cast cohort, now used for development |

The [split counts](../evidence/shape_20261025/campaign/split_counts.csv) and [membership](../evidence/shape_20261025/campaign/split_membership.csv) are authoritative. The reserve has one cast with a single retained depth; ordinary profile errors keep it, while onset cannot be located.

The primary 10/20/50% sets are frozen from earlier campaigns. The new 40% set lies between the existing 20/50% sets. The 75% set contains 50%, and all primary subsets below 100% exclude validation and reserve. Together those holdouts occupy 20%; the remaining 80% permits a 75% training sample. Unused casts are not silently reassigned.

Date-blocked training excludes every Ocean–Region–date group represented in date validation. Reserve is a separate cast set, not guaranteed to consist solely of unseen dates. Selection uses predictor/Region metadata and fixed seeds, not chlorophyll responses.

Every new posterior fit trains on at least 10% of full eligible casts. Historical 5% entries remain in the record but are not a design to repeat. Do not split a 10% subset into 5% training and 5% validation.

## Scaling and identity

Primary scaling comes from the earlier primary 10% training casts and is kept fixed across larger and repeat subsets. Date-blocked scaling is separate and frozen. Do not recalculate centers/scales on validation or a prediction grid. Changing scaling changes the meaning of coefficient priors as well as numerical conditioning.

Use data MD5, exact raw-row comparisons, selected cast membership, factor levels, and source script hashes together. Adding deterministic features changes serialized RDS hashes without necessarily changing raw observations, so a hash mismatch alone does not prove different training data.

The latest source identity points to `followup150_20261018T150024Z_3898749/full_input_250m.rds`, MD5 `cb5348ec5d9474d3f10570508e06d70d`. The shape preparation independently verifies raw identity after its depth restriction.

## Remaining upstream documentation

Obtain the MLD algorithm/reference/threshold, absolute chlorophyll calibration and units, precise distance derivation, cruise/station provenance, full exclusion tallies for each raw-data version, and permissions for data redistribution. Preserve the historical ENSO labels for reproducibility; assess an updated monthly index only as a separately versioned scientific change.

See the [region and split inventory](REGION_AND_SPLIT_INVENTORY.md) for exact regional cast counts and predictor ranges.
