# A new researcher's handover

## What this project is trying to accomplish

Predict the enhancement of chlorophyll with depth relative to each CTD cast's near-surface baseline. Reconstruct the observed profile, including when enhancement begins, how rapidly it increases, and the depth/height of the maximum. Determine whether mixed-layer depth (MLD), distance from shore, and Region explain that vertical structure without making the Bayesian fit impractical.

The project is in model development. It has a substantial comparison history and a large partially completed campaign. It has not selected a final production equation or established final out-of-sample performance.

## Reading order

Read the root README, then [data definitions](DATA_AND_PREPROCESSING.md), [equations](EQUATIONS.md), [results](RESULTS_AND_INTERPRETATION.md), and [evaluation](EVALUATION.md). Read the [runbook](IBEX_RUNBOOK.md) before interacting with Ibex. Use [campaign history](CAMPAIGN_HISTORY.md) to interpret old labels.

The authoritative configuration hierarchy is:

1. An actual run's `run_config`, resolved priors, data identity, and saved Stan code.
2. Its frozen experiment's `cases.tsv`, model constructor, data preparation, and environment.
3. The delivered package and its PLAN.
4. A dated analytical review.
5. Historical notebook comments or conversational proposals.

A planned case is not proof of an accepted Slurm submission. A completed fitting job is not proof of completed checking. A completed check is not proof of diagnostic validity. A low RMSE is not proof of good onset/uncertainty/generalization.

## Where we left off

The latest verified experiment is:

`/ibex/user/meleg/DCM_MODEL/experiments/shape150_20261022T102851Z_3426481`

Preparation job 54250664 completed. Initial dispatcher job 54250665 started. At the archived dispatcher time **2026-10-25T11:14:41Z**, 79 cases were complete, 48 fitting, and 155 not submitted. Requested CPUs for active model stages were 768. These are historical observations, not live queue claims.

The 79 completed cases comprise 67 posterior and 12 prior-only fits. All new 40%, 50%, 75%, and 100% training cases were still incomplete. The best-supported next action is to collect those existing results before designing another large campaign.

Regional MLD-varying models, especially shape cases 34, 129, and prior sensitivity 249, are promising. Retain no-MLD, half-shift and full-shift controls. Region-specific tensors and unseen-date results remain important unanswered comparisons.

The verified queue archive from **27 October 2026 at 10:47 UTC / 13:47 Saudi time** contains **19 running DCM fits, one other running job, and 315 pending jobs without DCM manifest matches**. The DCM fits use 304 CPUs and all belong to the shape campaign. See [CURRENT_JOBS.md](CURRENT_JOBS.md) for every job, linked model settings, source evidence and the separate October 25 historical registrations. **No DCM controller or check job appears in this snapshot.** Inspect the current controller ledger/accounting and saved dispatch timestamp before expecting further cases to launch. The missing controller is an observation, not a diagnosed failure; the uploaded queue supplies no new scientific results.

## First working session

1. Obtain authorized Ibex access and the exact raw/prepared data assets listed in [EXTERNAL_ASSETS.md](EXTERNAL_ASSETS.md). Do not substitute a new preparation into an existing frozen campaign.
2. Read the latest pointer; record the resolved experiment path and current UTC time.
3. Inspect queue/accounting and dispatcher status. Do not rerun a historical installer merely because a model is absent from a finished-results table.
4. Check `FIT_COMPLETE.txt`, `CHECKS_COMPLETE.txt`, failures, and the manifest. Fitting and checking are separate stages.
5. Refresh compact results, export a timestamped archive, and retain the previous snapshot.
6. Compare **every** completed case on common cohorts. Keep invalid, failed, prior-only, and incomplete entries labelled in the register.
7. Update this repository's dated state and decision log, including changes that overturn earlier conclusions.

## Essential conceptual distinctions

- Enhancement is not absolute chlorophyll, a ratio without subtraction, or percent unless multiplied by 100.
- A 20% sample is 20% of eligible **casts**, not random rows.
- `iter` includes warmup; `samples` is retained draws per chain.
- A numeric `by` smooth is a varying coefficient; a factor `by` smooth creates regional curves.
- A half/full MLD shift moves the profile along the depth coordinate; a regional MLD coefficient or tensor can also change shape.
- Regional `sigma` changes observation variability. Date/ENSO random intercepts shift the mean.
- `re_formula=NA` sets group effects to zero. It is not an average over new-date uncertainty.
- Full-data results are reconstruction. Former holdouts cease to be held out when included in training.
- Prior-only simulations are plausibility checks, not posterior accuracy results.
- Physical within-cast residual correlation and MCMC chain autocorrelation are different problems.

## What cannot yet be reproduced from GitHub alone

The compact evidence supports reading configurations, reviewing findings, and reconstructing many tables. It does not contain the CTD/MLD data, all posterior draws, or all archived cast-level prediction files. Some historical analysis scripts retain original absolute paths. The current review script accepts an extracted shape experiment; the source archive is required for its full numerical audit.

Package installation records are available, but there is no complete historical package lockfile or container. Repeating old fits exactly also needs data ordering, frozen scaling, contrasts/factor levels, CmdStan/compiler versions, seeds, and hardware context. Record any unavoidable deviation.

## Research decisions requiring a new result

Do not call regional MLD a definitive winner; the validation gain is partly an offset improvement and influenced by extreme casts. Do not assume more k, more tree depth, or higher adapt_delta always helps. Do not infer that less training data is better because a different training selection happened to score better. Do not delete difficult Red Sea casts to improve a reported primary score.

The user's standing priorities are whole-profile fidelity, explicit onset/sharpness evaluation, regional biology, computational practicality, and comparison of all models rather than an error-only shortlist.
