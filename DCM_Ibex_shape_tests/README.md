# Queue the MLD shape campaign on Ibex

Run the single delivered launcher in `/ibex/user/meleg/DCM_MODEL`:

```bash
# Install the versioned package and submit/resume its prepared campaign.
bash queue_shape_dcm.sh

# Monitor ALL older DCM campaigns and this new one.
bash DCM_Ibex_shape_tests/monitor_all_dcm_v8.sh --diagnostics

# Refresh current summary tables; no models are refitted.
bash DCM_Ibex_shape_tests/collect_results.sh

# Export compact comparison tables and supporting evidence.
bash DCM_Ibex_shape_tests/export_results.sh
```

There are 282 fits (270 posterior, 12 prior predictive), including 11 on 100%
of the eligible prepared 10–150 m data. Read `PLAN.md` and `cases.tsv` for all
equations, training fractions, basis dimensions and controls.

The launcher defaults to resuming the same campaign. Set `SHAPE_NEW=1` only
when intentionally requesting an additional independent campaign. Source data
are found through `logs/latest_followup150_experiment.txt`; `SHAPE_SOURCE`
can explicitly select another compatible prepared follow-up directory.

The source datasets, existing jobs and older campaign directories are retained.
Full-data reconstruction is labelled separately from held-out validation.
Actual fitting begins only after the compute-node preparation checks pass.
