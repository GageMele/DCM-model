# Scientific background and objectives

## Observational target

A deep chlorophyll maximum (DCM) is a subsurface maximum in chlorophyll. This project's measured target is the vertical pattern of **chlorophyll enhancement relative to a cast-specific surface baseline**. A fluorescence/chlorophyll maximum is not automatically a maximum in phytoplankton biomass, productivity, or a demonstrated causal response to mixing.

The immediate modelling question is practical and ecological: how well can a model recover the water-column profile from depth, Region, distance from shore and mixed-layer depth, while controlling for date and ENSO state?

The supplied CTD workflow combines NOAA/WHOI/KAUST-labelled source material and Red Sea/Pacific regional data. Current fitting excludes KAUST and retains eight other regions. The original input filename is not evidence that every row represents an independent instrument, cruise or calibrated product. Upstream instrument calibration, MLD derivation, and sampling provenance need their own source documentation; they are not reconstructed by this repository.

## Why these predictors matter to the question

Depth describes the vertical position of the response. Region permits very different background conditions and vertical patterns. Distance from shore is log10 transformed because the input distances cover a broad range; it is a predictor of association, not an experimental treatment.

MLD was introduced because the user expects it to help locate where chlorophyll begins to rise relative to surface water. The later, more specific hypothesis is that MLD can change **shape**: a deeper mixed layer may be associated with a sharper increase toward the DCM, not just a translated copy of a curve.

Date and ENSO are nuisance/context grouping variables in the current models. Their nested random intercepts account for shared mean offsets among observations in the same Ocean–Region–date or Ocean–Region–ENSO group. They do not fit a continuous seasonal trend, a mechanistic ENSO response, or group-specific profile shape. Date and ENSO are separate grouping terms, not one nested chain of date within ENSO.

## Questions the experiments separate

1. Does MLD only change average enhancement, or does its effect depend on depth?
2. Is the MLD effect common across regions or region-specific?
3. Can an MLD-relative coordinate align onset better?
4. Does a nonlinear depth–MLD surface improve on a linear-in-MLD varying coefficient?
5. Does distance need its own depth-varying effect or an interaction with MLD?
6. Is poor calibration mainly a shared residual-variance problem, a missing mean structure, or within-cast dependence?
7. How much spline complexity is supported in each region?
8. Do conclusions persist with different cast selections, larger samples and unseen dates?
9. Which sampler settings achieve reliable posterior computation at acceptable walltime/core cost?

These questions require different outcomes. Whole-profile RMSE, demeaned shape error, onset localization, peak depth, transition width, steepness, and predictive intervals are not interchangeable.

## Scientific interpretation of the latest evidence

Region-specific MLD coefficients improve completed 20% conditional profile errors, particularly in Central/South Red Sea. Half-shift can improve onset localization while making peak-depth predictions worse. Fixed universal sharpening warps have not consistently helped. Together, these results support keeping the shape-change hypothesis open and comparing regional representations.

The evidence does **not** yet establish a universal relation of the form “deeper MLD always sharpens the DCM.” Covariate distributions differ among regions; extreme fixed-distance MLD slices can be outside observed support. A fitted association also does not isolate MLD from all unmeasured environmental differences.

The large residual correlations and strong underprediction of steep rises show that a model can have excellent MCMC diagnostics yet still miss ecologically important structure.

## Historical trajectory

The recovered original final fit used a Gaussian regional depth spline to 250 m, raw log-distance, and nested date/ENSO intercepts. It sampled well, but fitting was expensive. The October 2026 work shortened the main analysis domain to 150 m, standardized and audited preparation, tested compute settings, and progressively added MLD.

The initial broad three-dimensional regional tensor was expensive and often poorly diagnosed. Simpler additive and varying-coefficient models made testing feasible. Regional residual SD improved calibration. Shifted coordinates then addressed onset. The newest campaign returns to shape changes with regional MLD coefficients, two-dimensional regional tensors, targeted basis sizes, repeated samples, and larger data.

See [CAMPAIGN_HISTORY.md](CAMPAIGN_HISTORY.md) for dated evidence and [EQUATIONS.md](EQUATIONS.md) for precise mathematical meanings.

## Limits to the background record

This is a project handover, not a completed literature review. It preserves the stated scientific motivation and empirical findings. Original acquisition/collaboration details, MLD algorithm/reference depth/threshold, chlorophyll calibration and units, and a definitive causal interpretation remain external or unresolved. Request source methods before writing publication-level methods for those upstream stages.
