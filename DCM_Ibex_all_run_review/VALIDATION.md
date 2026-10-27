# Validation, 16 October 2026

The collector passed local tests using R 4.5.2 with base R only:

```text
PASS: all-run retention, poor diagnostics, partial checks, prior separation, mixed schemas,
legacy evaluation labels, shared validation identity, cast/Region evidence, raw-draw exclusions,
unresolved benchmarks, empty campaigns, overwrite protection and all 102 cases without truncation.
```

The fixtures include a completed fit with 100 divergences and high prediction error, partial diagnostic output after failed checks, prior-only fits, original runs outside manifests, marker-only completed runs, and a benchmark without a resolvable run path. A separate 102-case fixture verifies that every case survives, including the highest-error case. Empty inputs produce a report without the historical zero-row assignment error. Existing output directories are rejected.

Reproduce on Ibex:

```bash
Rscript --vanilla DCM_Ibex_all_run_review/test_collection.R DCM_Ibex_all_run_review
bash -n DCM_Ibex_all_run_review/run_review.sh
```

The Linux launcher has not been executed on Ibex from this session. No real run results were modified or models submitted. The package exports available scores; it does not generate missing old-run held-out predictions or perform the final scientific comparison without the result files.
