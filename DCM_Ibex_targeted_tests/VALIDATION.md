# Validation record

Completed offline checks on 16 October 2026:

- Every packaged R file parses.
- The generated case table has 63 unique cases, 57 posterior and 6 prior, with four chains and no training fraction below 10%.
- A scheduler stub verifies CPU/memory requests, preparation-to-fit and fit-to-check dependencies, unique recorded IDs, old-fit evaluation scheduling and collection. No actual Slurm jobs were submitted by these local tests.
- Empty-result collection retains all 63 statuses without a zero-row assignment error. Synthetic completed models with R-hat 2 and 100 divergences remain in reports and pairwise comparisons. Different validation checksums are not paired.
- Independent offline checks verify the actual 108+108 primary and 108+117 date-blocked cast counts, no cast overlap, no shared nested date groups in the blocked panel, and unchanged training responses/depths.
- Profile/peak/AUC metric calculations, incremental old-fit discovery, duplicate prevention and prediction-only scheduling pass offline checks.
- The historical posterior diagnostics and read-only predictions completed. Three figures were produced using base R from saved CSVs. Matched regional distance references, depth-ordered profile lines and exclusive probability-bin sums were checked.
- Local actual-data model construction wrote complete Stan code, formula and resolved-prior files for cases 1-44. Separate earlier prototypes constructed the fixed-nu Student-t and AR covariance formulations.

The local full model-construction validation was stopped by automatic approval review after it detected an HTTPS telemetry connection with an unspecified payload that could expose model/data information. It was not retried or bypassed. All subsequent checks used base R or other local offline file operations. Consequently **the full 63-case construction and actual fit/prediction execution are not claimed to have passed locally**.

The Ibex preparation job must construct all 63 cases successfully before any fit can start. It verifies raw-data identity, whole-cast training/validation independence, proper prior-only specifications and predictor support through brms construction, then writes `PREPARATION_COMPLETE.txt`. A construction error writes `PREPARATION_FAILED.txt`; dependent fits are prevented from running. Standata/code construction does not establish that sampling will converge or that every downstream prediction will complete. Fit/check failures remain visible in the monitor and result tables.

The shell scripts could not be executed in the local Windows environment. Before any real submission, the installer and launcher run `bash -n` over their shell scripts on Ibex and execute the offline scheduler/reporting tests there. The requested maximum resource allocation is checked using `sbatch --test-only` before the first real submission.

No new Ibex job IDs exist until the user runs the installer in their authenticated Ibex terminal. This package does not contain an SSH credential or an authenticated cluster connection.
