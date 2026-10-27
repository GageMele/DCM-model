# Validation status

Passed locally without starting R/Stan: independent 34-row design generation;
24 new posterior/4 prior fits and six reused source fits; exact source-control
parameter matching; unique labels/IDs; CPU/memory/iteration arithmetic; reference
ordering; balanced R strings/delimiters; no fitting call in source reuse code;
Unix newlines; installer payload paths and SHA256 equality with source files.

R parser/runtime, model construction, Slurm and Bash are not claimed as locally
verified. Local R was previously blocked by automatic approval review because of
an attempted external telemetry connection; it is not retried here. This package
runs its base-R tests and Bash syntax checks on Ibex before submitting anything.
The tests cover onset shifts, absence, censoring, persistence, gaps, posterior
intervals, exact new-job counts, six reused fits, dependency/resources, repeats,
orphan recovery, uncertain-submission handling and failed-preparation blocking.

An Ibex preparation job then verifies the actual raw training/validation/reserve
rows, >=10% training, source-control settings, positive MLD, all requested brms
formulas/priors/Stan code and fixed-design rank before any model starts.
Failure stops fitting through afterok dependencies. New spline combinations
and posterior sampling therefore remain subject to real cluster preflight and
diagnostics; successful static checks do not guarantee model convergence.
