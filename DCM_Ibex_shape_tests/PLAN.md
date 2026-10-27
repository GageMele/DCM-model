# MLD shape, training size and basis-dimension campaign

This campaign contains **282 fits: 270 posterior and 12 prior predictive**.
There are **11 posterior fits on 100% of the eligible prepared data**.
The fitting depth interval remains **10–150 m**. “100%” means all casts and
observations in that prepared interval, retaining the established response,
quality filtering and surface-normalization workflow.

## Scientific questions

1. Can MLD change the steepness, width and peak of a regional profile?
2. Is a shift plus shape change preferable to absolute-depth shape change?
3. Do findings survive different 20% cast selections and withheld dates?
4. Which basis dimensions are supported, particularly in Central/South RS?
5. Do results transfer as training grows to 40%, 50%, 75% and 100%?

Every equation includes Region, nested Ocean:Region:date and Ocean:Region:enso
random intercepts and regional residual sigma. `distance_z` and `mld_z` use
the previously frozen, training-derived scaling. No new scaling is estimated
from a validation response or fitted separately on prediction data.

## Equations

For the explanations below, z is depth in metres, m is MLD in metres, M is
standardized MLD, D is standardized log10 distance, and R indexes Region.
The table shows additional mean terms after the Region and random effects.

| Equation | Mean-profile structure | What MLD can change |
|---|---|---|
| original_z | D + f_R(z) | Benchmark without MLD |
| interaction_z | f_R(z) + M g(z) + D h(z) + MD j(z) | Depth-dependent shared MLD/distance effects |
| mld_halfshift | D + f_R(z − 0.5m) | Profile position |
| mld_shift | D + f_R(z − m) | Profile position |
| mld_vary_shared | D + f_R(z) + M g(z) | Shared, depth-dependent shape deformation |
| mld_vary_region | D + f_R(z) + M g_R(z) | Region-specific, depth-dependent shape deformation |
| tensor_region | D + F_R(z, M) | Nonlinear changes in shape, height, width and position |
| tensor_shift_region | D + F_R(z − m, M) | Shape changes around an MLD-relative coordinate |
| tensor_halfshift_region | D + F_R(z − 0.5m, M) | Shape changes around the half-shift coordinate |
| warp_half_pos025 | D + f_R((z − 0.5m) exp(0.25M)) | Explicit narrowing/steepening as MLD increases |
| warp_half_pos050 | D + f_R((z − 0.5m) exp(0.50M)) | Stronger prescribed narrowing/steepening |
| warp_half_neg025 | D + f_R((z − 0.5m) exp(−0.25M)) | Opposite width response as a control |

The width-warp coefficients are fixed hypotheses, not estimated physical
constants. The regional tensors and varying-coefficient models learn shape
changes from data. Tensor models use `t2(..., by = Region, k = c(k_depth,
k_mld))` without adding a redundant copy of the same main-depth smooth.

## Data and experimental design

- Fractions are 10%, 20%, 40%, 50%, 75% and 100% of all eligible whole casts.
- Existing primary 10%, 20% and 50% cast selections are retained exactly.
  The new 40% selection is nested between the existing 20/50% sets; the 75%
  selection contains the existing 50% set.
- The same primary 10% validation and 10% reserve cohorts are excluded from
  every primary training subset below 100%. The 75% subset is feasible within
  the remaining 80% pool. Unused casts remain outside those subsets.
- Three additional 20% selections use seeds 20261023/20261024/20261025 while
  holding the MCMC seed at 1234. Thus cast-selection and sampler randomness
  are different experimental factors.
- Eighteen date-blocked fits use two new 20% selections. Every nested
  Ocean/Region/date group represented in the fixed date validation set is
  excluded from those training sets. Reserve remains an independent cast
  cohort; it is not claimed to contain exclusively unseen date groups.
- Full-data fits include former validation/reserve casts. Their outputs are
  named `training_full`, `reconstruction_primary` and `reconstruction_reserve`.
  They never receive held-out validation scores or held-out rankings.
- Previously inspected validation/reserve sets remain development cohorts.
  The campaign does not create an untouched final test by relabelling them.

The full frozen case matrix, including every exception, is `cases.tsv`.

## Basis dimensions and controls

Regional depth bases are tested at **7, 9, 11, 13 and 15**. The broad
size-by-equation comparison uses 7/11/15; extra 9/13 tests occur at 20/50%.
Eighteen additional fits use **k=15 in Central and South RS and k=7 elsewhere**.
Tensor MLD margins use k=4, with additional 5/7 checks at 20%; this is a
separate basis dimension from depth and must be supported by distinct MLDs
within each region. Requested dimensions are never silently reduced.

The original MLD–distance interaction retains its established shared smooth
bases at 5 as a control, while its regional depth basis varies. New shared
MLD-varying depth smooths use the tested depth k. Region-specific MLD-varying
smooths follow the regional depth-basis setting.

The main posterior configuration is four chains, 1,000 warmup and 2,000
retained draws per chain, adapt_delta=0.97 and max_treedepth=14. Sixteen
targeted fits change one of adapt_delta (0.95/0.99), warmup (2,000) or retained
draws (4,000) at a time. Four prior-sensitivity fits and 12 prior-predictive
checks retain their explicit labels.

Fits request four threads/chain (16 CPUs) through 50% and eight threads/chain
(32 CPUs) at 75/100%. Memory is 64–256 GiB according to size/complexity.
Walltime is capped at the partition's reported MaxTime. Every fit saves its
CPU model, allocation, sampler CSVs, warmup diagnostics, priors and hydrated fit.

## Shape and prediction evaluation

All raw eligible observations stay in profile RMSE, MAE, shape RMSE, peak
depth/height, positive enhancement area, predictive coverage and interval
score. Conditional and zero-group-effect population predictions are separate.
Regional results and individual-cast errors remain available for every model.

Additional questions are assessed with:

- Sustained +25%, +50% and +100% enhancement onsets, with missing support,
  non-detection and left/gap censoring explicit.
- Maximum positive 4 m centred slope and the depth of the strongest rise.
- Transition width between 25% and 75% of each profile's positive enhancement
  peak, with inadequate enhancement, boundary peaks and censoring flagged.
- Model-based MLD quantile slices at fixed regional distance, with expected
  profile uncertainty and pairwise training-hull support flags saved.
- Planned pairwise differences using the same casts. Onset/width differences
  use casts comparable in both models and report the number retained.

Onset and sharpness diagnostics use the established 1 m grid/5 m running
median on supported observations only. These scoring operations never modify
the data used to fit the model. A singleton cast remains in ordinary errors.
The new width metric is an operational profile-shape diagnostic, not a new
universal DCM definition.

## Submission and recovery

`queue_shape_dcm.sh` installs a versioned package, freezes a new campaign once,
and submits preparation and a lightweight dispatcher. Formula/data checks run
on a compute node before fitting begins. A failed individual preflight is
retained and skipped; it does not cause all valid models to be abandoned.

By default the dispatcher allows **48 fitting/checking jobs and 1,024 requested
CPUs** at once. Slurm still determines actual execution. Explicit submission
quota rejections are retried as capacity becomes available; statistical
settings and resource requests are not silently changed after other errors.
Other DCM/CROCO jobs are not cancelled or modified.

The dispatcher prioritizes 20% anchors and quick prior checks, submits later
cases as slots open, checks completed fits, and refreshes every model's summary.
It arranges a successor before its walltime ends; if a submission quota prevents
this it attempts to requeue only the lightweight dispatcher. If the cluster
disallows that action it records a clear restart instruction. Accepted model
fits are never automatically repeated. Re-running the launcher resumes the
same frozen campaign. An uncertain sbatch acknowledgement is recovered from
its unique job name or stopped for inspection rather than duplicated.

Change operational concurrency before launching, or re-run with updated limits:

```bash
# Optional: use more concurrent slots if the account has sufficient capacity.
SHAPE_MAX_ACTIVE=64 SHAPE_MAX_CPUS=1536 bash queue_shape_dcm.sh
```

Preparation/runtime code uses the R/CmdStan environment from the existing
prepared follow-up campaign. No packages are installed automatically. This
delivery prepares submission; it cannot establish new Ibex job IDs until the
launcher is executed through the user's authenticated Ibex session.

## Methods references

- brms smooth support: https://paulbuerkner.com/brms/reference/s.html
- mgcv tensor specification: https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/t2.html
- Slurm dependencies and submitted-job behaviour: https://slurm.schedmd.com/sbatch.html

Array throttling alone does not bypass MaxSubmitJobs; this campaign instead
limits actual submissions. No new fitting is inferred from a dependency ID
or a case registration before Slurm returns an accepted job ID.
