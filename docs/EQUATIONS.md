# Equation catalogue and interpretation

The actual constructors are [targeted campaign.R](../DCM_Ibex_targeted_tests/campaign.R), [follow-up campaign.R](../DCM_Ibex_followup_tests/campaign.R), [onset campaign.R](../DCM_Ibex_onset_tests/campaign.R), and [shape models.R](../DCM_Ibex_shape_tests/models.R). Read them with the case matrices: a label does not encode every setting.

## Notation and shared structure

Let y be chl_norm, z physical depth in metres, m raw MLD in metres, L raw log10 distance, M standardized MLD, D standardized log10 distance, and R Region. Standardization uses frozen training-derived moments. f/g/h/j denote penalized smooth functions; F denotes a tensor surface.

Unless a named experiment changes it, the mean also contains:

~~~r
1 + Region +
  (1 | Ocean:Region:date) +
  (1 | Ocean:Region:enso)
~~~

The colon creates a combined grouping factor. This does not automatically add separate Ocean or Region random intercepts. Date and ENSO are separate nested group effects; date is not nested inside ENSO. Their random intercepts shift the mean without changing profile shape.

The original Gaussian model used a single residual SD. Most later leading models use:

~~~r
bf(mean_formula, sigma ~ 0 + Region, family = gaussian())
~~~

The sigma submodel uses a log link: its coefficients are log residual SDs, not chlorophyll slopes. Region in the mean and Region in sigma serve different purposes.

## Original model and the three early MLD equations

The user's original mean formula was:

~~~r
andrew_form <- bf(
  chl_norm ~ 1 + Region + log10_distance +
    s(depth, by = Region, bs = "cr", k = 15) +
    (1 | Ocean:Region:date) +
    (1 | Ocean:Region:enso),
  family = gaussian()
)
~~~

The commented-out global depth smooth was not active. The recovered saved final fit has this structure, no MLD, and a 10–250 m fitting domain.

| Equation campaign | Additional terms beyond intercept/Region/groups | Interpretation |
|---|---|---|
| A, additive | L + m + f_R(z) | MLD changes a depth-independent mean offset |
| B, mld_varying | L + m g(z) + f_R(z) | Shared depth-dependent MLD deformation plus regional baseline curves |
| C, mld_only | L + m g(z) | No regional depth smooth; Region retains intercept differences |

The B mean terms requested by the user were:

~~~r
log10_distance +
  s(depth, by = MLD_temp_MIX, bs = "cr", k = 15) +
  s(depth, by = Region, bs = "cr", k = 15)
~~~

C removes the final regional smooth. A adds MLD_temp_MIX as an ordinary linear predictor instead.

### Numeric versus factor by variables

A numeric by-MLD smooth contributes m g(z). It can change profile shape, not merely its intercept. At fixed depth it is linear in MLD. Its depth derivative is m g'(z), so slopes can change with MLD.

A factor by variable such as Region creates separate curves. A numeric-by smooth is not a general bivariate surface and does not mean “move depth by MLD.” Numeric-by terms can contain a constant-in-depth component, so redundant linear terms require care. Frozen constructors are authoritative.

## Extensive and targeted families

The extensive campaign screened:

| Family | Additional mean structure |
|---|---|
| Original | L + f_R(z) |
| Additive MLD | L + m + f_R(z) |
| Shared MLD varying | L + m g(z) + f_R(z) |
| Distance varying | m + L h(z) + f_R(z) |
| Both varying | m g(z) + L h(z) + f_R(z) |
| Both standardized | M g(z) + D h(z) + f_R(z) |
| Global baseline | f(z) + m g(z) + L h(z) |
| Regional/global MLD tensor | L + F_R(z,m), or L + F(z,m) |
| Regional/global distance tensor | m + F_R(z,L), or m + F(z,L) |
| Regional three-dimensional tensor | F_R(z,L,m) |

A 15 × 5 × 5 tensor requests 375 product-basis components per Region before constraints and penalty decomposition. This helps explain the cost of the early regional 3D trials. Their poor computation does not establish that every smaller 2D tensor is unsuitable.

Targeted variants include:

| Label | Additional mean structure |
|---|---|
| original_raw | L + f_R(z) |
| original_z | D + f_R(z) |
| additive_z | D + M + f_R(z) |
| mld_vary_z | D + M g(z) + f_R(z) |
| distance_vary_z | M + D h(z) + f_R(z) |
| both_z | M g(z) + D h(z) + f_R(z) |
| region_mld_z | D + M g_R(z) + f_R(z) |
| region_both_z | M g_R(z) + D h_R(z) + f_R(z) |
| interaction_z | f_R(z) + M g(z) + D h(z) + MD j(z) |

Standardization changes prior interpretation as well as computation. It need not preserve the posterior under unchanged numerical priors. MD is the product of the separately standardized predictors.

## Half-shift and full-shift

Full shift defines u = z − m; half shift defines u = z − 0.5m. They replace the coordinate inside the regional curve:

~~~r
CTD$depth_mld <- CTD$depth - CTD$MLD_temp_MIX
CTD$depth_mld_half <- CTD$depth - 0.5 * CTD$MLD_temp_MIX

full_shift <- bf(
  chl_norm ~ 1 + Region + distance_z +
    s(depth_mld, by = Region, bs = "cr", k = 7) +
    (1 | Ocean:Region:date) + (1 | Ocean:Region:enso),
  sigma ~ 0 + Region, family = gaussian()
)
half_shift <- bf(
  chl_norm ~ 1 + Region + distance_z +
    s(depth_mld_half, by = Region, bs = "cr", k = 7) +
    (1 | Ocean:Region:date) + (1 | Ocean:Region:enso),
  sigma ~ 0 + Region, family = gaussian()
)
~~~

Treating MLD as a vertical coordinate means locating observations relative to the mixed layer before evaluating the smooth. If MLD increases from 30 to 50 m, full shift translates a fixed feature 20 m deeper; half shift translates it 10 m deeper.

Neither forces onset or maximum to occur at MLD. A feature at u = 15 occurs at z = m + 15 for full shift. The 0.5 coefficient is a fixed hypothesis, not an estimated mixing constant.

Holding other predictors fixed, translation preserves intrinsic width and steepness: the derivative is f'_R(z − a m). It changes which part of the curve appears at a physical depth. The finite observation window can truncate a translated curve, but this is not learned shape deformation.

Other onset models used D + f_R(z/m), requiring positive MLD; a hybrid D + f_R(z) + g(z−m); and shifted distance effects f_R(z−m) + D h(z−m). The ratio stretches depth proportionally to MLD. The hybrid permits more than translation, but overlapping smooth explanations can be hard to identify.

## Current shape campaign: twelve families

| Label | Additional mean structure | Flexibility |
|---|---|---|
| original_z | D + f_R(z) | No MLD |
| interaction_z | f_R(z) + M g(z) + D h(z) + MD j(z) | Shared MLD/distance shape effects |
| mld_halfshift | D + f_R(z−0.5m) | Translation |
| mld_shift | D + f_R(z−m) | Translation |
| mld_vary_shared | D + f_R(z) + M g(z) | Shared linear-in-MLD deformation |
| mld_vary_region | D + f_R(z) + M g_R(z) | Regional deformation |
| tensor_region | D + F_R(z,M) | Nonlinear regional MLD/depth shape |
| tensor_shift_region | D + F_R(z−m,M) | Shift plus nonlinear deformation |
| tensor_halfshift_region | D + F_R(z−0.5m,M) | Half-shift plus nonlinear deformation |
| warp_half_pos025 | D + f_R((z−0.5m) exp(0.25M)) | Prescribed narrowing with deeper MLD |
| warp_half_pos050 | D + f_R((z−0.5m) exp(0.50M)) | Stronger prescribed narrowing |
| warp_half_neg025 | D + f_R((z−0.5m) exp(−0.25M)) | Opposite width response |

For regional varying effects, the constructor makes numeric columns region_mld_i = I(Region == region_i) × mld_z and gives each its own depth smooth alongside the regional baseline. Each Region learns its own MLD deformation; at fixed depth it remains linear in MLD.

A representative two-dimensional tensor is:

~~~r
bf(
  chl_norm ~ 1 + Region + distance_z +
    t2(depth, mld_z, by = Region,
       bs = c("cr", "cr"), k = c(7, 4)) +
    (1 | Ocean:Region:date) + (1 | Ocean:Region:enso),
  sigma ~ 0 + Region, family = gaussian()
)
~~~

The tensor replaces the regional baseline depth term; no redundant copy is added. MLD dimensions 4/5/7 differ from depth dimensions 7–15. Distinct predictor values within each Region must support them. Preflight does not silently reduce k.

For u = (z−0.5m) exp(lambda M), the depth derivative is exp(lambda M) f'_R(u). Positive lambda imposes narrowing/steepening. These fixed coefficients are hypotheses, not inferred physical laws.

## Basis dimension and regional allocation

k limits basis dimension; penalties determine effective wiggliness. Raising k allows complexity without forcing it. An unconverged fit cannot identify an adequate k from RMSE.

The broad grid uses 7/11/15, with selected 9/13 tests. The regional allocation uses 15 for Central RS and South RS, 7 elsewhere. For unequal bases, ordered off/on factors create each intended baseline smooth without duplicating Region intercepts. Numeric regional MLD smooths follow that allocation.

The interaction_z control keeps its shared MLD, distance, and product bases at 5 while regional k changes. New mld_vary_shared smooths use the case's depth k. Not every smooth in a case has the same k.

## Variance, likelihood, and dependence experiments

- Shared Gaussian sigma: one observation SD.
- Regional sigma: sigma ~ 0 + Region.
- Regional depth-varying sigma: sigma ~ 0 + Region + s(depth, by=Region, k=5), with proper log-scale and smooth-SD priors.
- Student likelihood with estimated degrees of freedom and fixed nu = 4.
- Cast random intercepts for cast offsets.
- Cast-specific natural spline coefficients: (1 + cast_b1 + cast_b2 + cast_b3 + cast_b4 || cast_key), knots 45/80/115 m, boundaries 10/150 m. Double bars omit correlations among these coefficients.
- AR(1): ar(time=depth, gr=cast_key, p=1, cov=TRUE), in the recorded scalar-sigma configuration with matched unthreaded controls. This is dependence along ordered rows, not a continuous-distance correlation per metre.

Cast effects change new-cast prediction. Do not judge them only by in-sample fit.

## Primary software references

[brms smooths](https://paulbuerkner.com/brms/reference/s.html), [mgcv t2](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/t2.html), [Stan diagnostics](https://mc-stan.org/learn-stan/diagnostics-warnings.html). Historical versions and frozen constructors take precedence over assuming current documentation implies identical old behavior.
