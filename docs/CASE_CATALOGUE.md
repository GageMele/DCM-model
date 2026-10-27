# Complete case catalogue

Every recovered case matrix is included. A case number is unique only within its campaign. The machine-readable [all_cases.tsv](../provenance/all_cases.tsv) is the union of all original columns; blanks mean that field was not present in that matrix, not that the setting was zero.

## DCM_Ibex_baseline_tests

12 cases. [Exact case matrix](../DCM_Ibex_baseline_tests/cases.tsv).

- **equation**: `original_region_depth`.
- **mode**: `posterior`.
- **subset**: `p10`, `p20`, `p50`.
- **prior**: `custom`, `default`.
- **chains**: `4`.
- **warmup**: `1000`, `2500`.
- **adapt_delta**: `0.97`, `0.99`.
- **treedepth**: `12`, `14`.
- **threads**: `2`.
- **k_depth**: `15`.

## DCM_Ibex_equation_tests

18 cases. [Exact case matrix](../DCM_Ibex_equation_tests/cases.tsv).

- **equation**: `additive`, `mld_only`, `mld_varying`.
- **mode**: `posterior`.
- **subset**: `p10`.
- **prior**: `current`.
- **chains**: `2`, `4`.
- **warmup**: `1000`, `500`.
- **adapt_delta**: `0.95`, `0.99`.
- **treedepth**: `10`, `12`, `14`.
- **threads**: `4`.
- **k_depth**: `15`.

## DCM_Ibex_extensive_tests

102 cases. [Exact case matrix](../DCM_Ibex_extensive_tests/cases.tsv).

- **block**: `chains_half_total`, `chains_same_total`, `data_size`, `delta`, `draws`, `equation_k`, `likelihood`, `mixed_k`, `prior_predictive`, `priors`, `repeat_seed`, `scaling`, `small_3d`, `threads`, `tree`, `warmup`.
- **equation**: `additive`, `both_scaled`, `both_vary`, `distance_vary`, `global_both`, `global_tensor_distance`, `global_tensor_mld`, `mld_vary`, `original`, `region_tensor_3d`, `region_tensor_distance`, `region_tensor_mld`.
- **mode**: `posterior`, `prior`.
- **subset**: `p10`, `p20`, `p50`.
- **prior**: `custom`, `default`, `half_normal`, `tighter`.
- **family**: `gaussian`, `student`.
- **chains**: `2`, `4`.
- **warmup**: `1000`, `2000`, `250`, `500`.
- **samples**: `1000`, `2000`, `500`.
- **adapt_delta**: `0.9`, `0.95`, `0.97`, `0.99`.
- **treedepth**: `10`, `12`, `14`.
- **threads**: `1`, `2`, `4`.
- **k_depth**: `10`, `15`, `5`, `7`.
- **k_mld**: `10`, `15`, `3`, `4`, `5`, `7`.
- **k_distance**: `10`, `15`, `3`, `4`, `5`, `7`.

## DCM_Ibex_fast_tests

24 cases. [Exact case matrix](../DCM_Ibex_fast_tests/cases.tsv).

- **mode**: `posterior`, `prior`.
- **subset**: `p05`, `p10`, `p20`.
- **prior**: `current`, `default`, `empirical`.
- **chains**: `2`, `4`.
- **warmup**: `250`, `500`.
- **adapt_delta**: `0.90`, `0.95`.
- **treedepth**: `10`, `12`.
- **threads**: `0`, `2`, `4`.

## DCM_Ibex_followup_tests

34 cases. [Exact case matrix](../DCM_Ibex_followup_tests/cases.tsv).

- **block**: `date_blocked`, `delta`, `depth_basis`, `effect_basis_high`, `effect_basis_low`, `prior_sigma`, `sigma_control`, `sigma_depth`, `sigma_prior`, `stability`, `threads`, `training_composition`, `training_size`.
- **equation**: `both_z`, `interaction_z`, `original_z`.
- **mode**: `posterior`, `prior`.
- **subset**: `p10`, `p10_alt`, `p20`.
- **prior**: `custom_z`.
- **family**: `gaussian`.
- **variance**: `region`, `region_depth`.
- **cast**: `none`.
- **autocor**: `none`.
- **chains**: `4`.
- **warmup**: `1000`, `2000`, `500`.
- **samples**: `2000`, `500`.
- **adapt_delta**: `0.9`, `0.95`, `0.99`.
- **treedepth**: `14`.
- **threads**: `2`, `4`.
- **k_depth**: `10`, `7`.
- **k_mld**: `3`, `5`, `7`.
- **k_distance**: `3`, `5`, `7`.

## DCM_Ibex_k_tests

12 cases. [Exact case matrix](../DCM_Ibex_k_tests/cases.tsv).

- **mode**: `posterior`.
- **subset**: `p10`.
- **prior**: `current`.
- **chains**: `2`.
- **warmup**: `1000`, `500`.
- **adapt_delta**: `0.95`, `0.99`.
- **treedepth**: `12`, `14`.
- **threads**: `4`.
- **k_depth**: `10`, `15`, `5`, `7`.
- **k_mld**: `3`, `4`, `5`.
- **k_distance**: `3`, `4`, `5`.

## DCM_Ibex_onset_tests

34 cases. [Exact case matrix](../DCM_Ibex_onset_tests/cases.tsv).

- **block**: `depth_basis`, `equation`, `existing_control`, `prior_predictive`, `repeat_seed`, `robust_sampler`.
- **equation**: `interaction_z`, `mld_halfshift`, `mld_hybrid`, `mld_ratio`, `mld_shift`, `mld_shift_distance`, `mld_vary_z`, `original_z`.
- **mode**: `posterior`, `prior`.
- **subset**: `p10`, `p20`.
- **prior**: `custom_z`.
- **family**: `gaussian`.
- **variance**: `region`.
- **cast**: `none`.
- **autocor**: `none`.
- **chains**: `4`.
- **warmup**: `1000`, `2000`, `500`.
- **samples**: `2000`, `500`.
- **adapt_delta**: `0.95`, `0.99`.
- **treedepth**: `14`.
- **threads**: `2`, `4`.
- **k_depth**: `10`, `7`.
- **k_mld**: `5`.
- **k_distance**: `5`.

## DCM_Ibex_shape_tests

282 cases. [Exact case matrix](../DCM_Ibex_shape_tests/cases.tsv).

- **block**: `cast_selection`, `date_blocked`, `delta95`, `delta99`, `full_data`, `full_data_high_k`, `full_data_surface`, `intermediate_k`, `mld_basis`, `mld_surface`, `mld_width`, `prior_predictive`, `prior_sensitivity`, `regional_k`, `samples4000`, `size_basis`, `warmup2000`.
- **equation**: `interaction_z`, `mld_halfshift`, `mld_shift`, `mld_vary_region`, `mld_vary_shared`, `original_z`, `tensor_halfshift_region`, `tensor_region`, `tensor_shift_region`, `warp_half_neg025`, `warp_half_pos025`, `warp_half_pos050`.
- **mode**: `posterior`, `prior`.
- **subset**: `p10`, `p100`, `p20`, `p40`, `p50`, `p75`.
- **prior**: `b_tight`, `custom_z`, `scales_halfnormal`.
- **family**: `gaussian`.
- **variance**: `region`.
- **cast**: `none`.
- **autocor**: `none`.
- **chains**: `4`.
- **warmup**: `1000`, `2000`, `500`.
- **samples**: `2000`, `4000`, `500`.
- **adapt_delta**: `0.95`, `0.97`, `0.99`.
- **treedepth**: `14`.
- **threads**: `4`, `8`.
- **k_depth**: `11`, `13`, `15`, `7`, `9`.
- **k_mld**: `4`, `5`, `7`.
- **k_distance**: `5`.
- **k_profile**: `rs15_other7`, `uniform`.
- **selection_seed**: `0`, `20261022`, `20261023`, `20261024`, `20261025`.

## DCM_Ibex_targeted_tests

63 cases. [Exact case matrix](../DCM_Ibex_targeted_tests/cases.tsv).

- **block**: `ar_reference`, `data_size`, `date_blocked`, `delta`, `depth_k`, `draws`, `equation`, `historical_depth_bridge`, `likelihood`, `prior_predictive`, `priors`, `repeat_seed`, `residual_ar`, `structure`, `threads`, `tree`, `warmup`.
- **equation**: `additive_z`, `both_z`, `distance_vary_z`, `interaction_z`, `mld_vary_z`, `original_raw`, `original_z`, `region_both_z`, `region_mld_z`.
- **mode**: `posterior`, `prior`.
- **subset**: `p10`, `p20`, `p50`.
- **prior**: `b_tight`, `custom_z`, `mu_default`, `scales_halfnormal`.
- **family**: `gaussian`, `student4`.
- **variance**: `constant`, `region`, `region_depth`.
- **cast**: `curve`, `intercept`, `none`.
- **autocor**: `ar1`, `none`.
- **chains**: `4`.
- **warmup**: `1000`, `2000`, `500`.
- **samples**: `1000`, `2000`, `4000`.
- **adapt_delta**: `0.9`, `0.95`, `0.97`.
- **treedepth**: `12`, `14`.
- **threads**: `1`, `2`, `4`.
- **k_depth**: `10`, `15`, `7`.
- **k_mld**: `5`.
- **k_distance**: `5`.

## Field meanings

| Field | Meaning |
|---|---|
| package / campaign | Delivery family versus actual timestamped experiment identity |
| case_id / label | Frozen within-campaign identifier and descriptive label |
| equation | Constructor switch label; qualify by campaign |
| mode / execution | Posterior or prior; new fit versus reused control where present |
| subset / fraction | Whole-cast training fraction, not row fraction |
| split / dataset_id | Cohort scheme and saved dataset identity |
| selection_seed / seed | Cast selection versus MCMC randomness |
| k_depth / k_mld / k_distance | Marginal basis budgets or effect-specific k, as defined by constructor |
| k_profile | Uniform versus targeted regional k |
| variance / cast / autocor / family | Residual structure, cast terms, dependence and likelihood |
| prior | Requested prior family; inspect resolved priors |
| chains / cores / threads / cpus | Sampling chains, parallel-chain limit, within-chain threads and allocation |
| iter / warmup / samples | Total, adaptation and retained iterations per chain |
| adapt_delta / treedepth | NUTS adaptation target and maximum tree-depth cap |
| memory_gb / refresh | Requested memory and progress-print interval |
| reference_case / *_reference_case | Direct planned comparison within the campaign |
| evaluation_role / heldout_eligible | Reconstruction versus genuinely excluded evaluation data |

The registry's fitted result fields also include observed diagnostics, errors and timings. A case matrix states intent; run_config, resolved priors and saved outputs establish execution.
