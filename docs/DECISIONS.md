# Decision and rationale record

| Decision/constraint | Reason | Status |
|---|---|---|
| Predict fractional enhancement relative to cast surface | Scientific target supplied by project owner | Current |
| Fit raw unsmoothed, unbinned observations | Explicit user requirement; retain vertical information | Current for reference fits; 5 m bin medians allowed for screening since 28 Oct |
| Main current domain 10–150 m | Reduce cost and focus DCM profile; original fit extended to 250 m | Current, retain depth-domain bridge |
| Include Region depth structure | MLD-only models lost important shape | Supported |
| Nested date and ENSO intercepts | Control shared contextual offsets | Current; not causal adjustment guarantee |
| Train on at least 10% whole casts | User rejected splitting 10% into 5% training | Current |
| Keep common validation/reserve cohorts | Controlled comparisons | Current development sets |
| Keep every run in comparison | User explicitly requested all models, not only best RMSE | Required |
| Regional residual sigma baseline | Strong calibration improvement | Supported, recheck new sizes |
| Continue shifts and shape-changing MLD models | Shift alone cannot express the user's sharpening hypothesis | Shape-changing MLD now tested through the peak profile (MLD on amplitude and width) |
| Test k7–15, larger k in Central/South RS | Regional complexity concern plus early results | Current |
| Test 10/20/40/50/75/100% and repeated 20% selections | Sample-size and cohort sensitivity | Already designed/partly complete |
| No automatic single-score winner | Shape, onset, calibration, validity and cost differ | Required |
| Repair checking from saved fits | Reporting/support bugs do not justify refitting | Applied |
| Preserve frozen source and registration IDs | Avoid duplicate fits and loss of provenance | Required |
| Refresh existing campaign before another broad grid | Many directly relevant tests were pending October 25 | Next action |
| No final production equation yet | Larger-data/tensor/new-date evidence incomplete | Open |
| **Choose structure on unseen dates (date-blocked split)** | Region-specific terms win seen-date tests and fail on new dates (shape150 270/279; next150 9/10); literature review of 27 Oct | Adopted 27 Oct 2026 |
| Score every case by joint per-date elpd as well as RMSE/onset | Log score is the general criterion; RMSE alone hides interval quality | Adopted for next150 onward |
| **Lognormal likelihood on ratio = chl_norm + 1** | next150: +5,600 / +4,500 elpd vs Gaussian; no impossible negative Chl | Adopted as default, 9 Nov 2026 |
| Half-shift `depth − 0.5·MLD` as the depth coordinate | Beat no-MLD in all 8 date-blocked comparisons; the peak model finds no extra MLD displacement | Current |
| Stop regional tensors and pooled regional MLD | 72–124 h per fit; worse on unseen dates | Decided Oct–Nov 2026 |
| Allow 5 m bin medians for fitting/screening | Gage approved the binned test pair on 28 Oct; about 35× faster at equal accuracy; scored on raw rows | Current for screening |
| Test depth as independent bins ("histogram" models) at 5 m and 1 m | Gage's request, 9 Nov | Done; not better than smooths on unseen dates |
| Gaussian versions of the lognormal models | Gage's request, 9 Nov | Done (round3); Gaussian trails lognormal everywhere |
| Wait for running tests before queueing more | Gage, 12 Nov | Current |

This record reflects decisions and constraints evidenced in the project history. Proposed future work is not a claim that new jobs were submitted.
