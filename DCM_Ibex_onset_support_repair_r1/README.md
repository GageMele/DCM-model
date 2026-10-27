# Onset support repair R1

The old scorer aborted a whole posterior check when a profile had fewer than
seven one-metre grid points. Unsupported onset metrics are now explicitly flagged
insufficient_support, including singletons and insufficient support after 15 m.
Every observation remains in ordinary prediction errors; missing onset is neither
imputed nor treated as absence. Summaries count unsupported and comparable casts.
Standard scores are saved before onset-specific calculations.

Only onset_metrics.R and predictive.R in the selected experiment are patched,
with exact source hashes and backups. Data, cases, formulas, priors, fits, sampler
settings and original job registrations are preserved. Recovery workers archive
failed checking outputs and reuse the existing diagnostic/prediction scripts.
They wait for active original jobs, including unfinished case 18, and skip completed
checks and prior-predictive cases. No fitting jobs are submitted or cancelled.

Local checks: design/source identities, R delimiters, permitted patch targets,
safe archive paths, payload equality and SHA256. Actual R parsing and tests run
on Ibex before queueing. Tests cover short/empty/singleton, shifted/absent/censored/
gapped onsets, posterior uncertainty, denominators and submission recovery. All
real validation/reserve casts are audited, with counts and cast IDs exported.
Local R is not retried after the earlier automatic-review block.

Run bash repair_onset_checks.sh after upload. Repeating resumes accepted checking
registrations. A final collector refreshes all registered campaigns. Later use
bash DCM_Ibex_onset_tests/collect_results.sh to refresh completed results manually.
