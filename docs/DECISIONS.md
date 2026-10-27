# Decision and rationale record

| Decision/constraint | Reason | Status |
|---|---|---|
| Predict fractional enhancement relative to cast surface | Scientific target supplied by project owner | Current |
| Fit raw unsmoothed, unbinned observations | Explicit user requirement; retain vertical information | Current |
| Main current domain 10–150 m | Reduce cost and focus DCM profile; original fit extended to 250 m | Current, retain depth-domain bridge |
| Include Region depth structure | MLD-only models lost important shape | Supported |
| Nested date and ENSO intercepts | Control shared contextual offsets | Current; not causal adjustment guarantee |
| Train on at least 10% whole casts | User rejected splitting 10% into 5% training | Current |
| Keep common validation/reserve cohorts | Controlled comparisons | Current development sets |
| Keep every run in comparison | User explicitly requested all models, not only best RMSE | Required |
| Regional residual sigma baseline | Strong calibration improvement | Supported, recheck new sizes |
| Continue shifts and shape-changing MLD models | Shift alone cannot express the user's sharpening hypothesis | Current |
| Test k7–15, larger k in Central/South RS | Regional complexity concern plus early results | Current |
| Test 10/20/40/50/75/100% and repeated 20% selections | Sample-size and cohort sensitivity | Already designed/partly complete |
| No automatic single-score winner | Shape, onset, calibration, validity and cost differ | Required |
| Repair checking from saved fits | Reporting/support bugs do not justify refitting | Applied |
| Preserve frozen source and registration IDs | Avoid duplicate fits and loss of provenance | Required |
| Refresh existing campaign before another broad grid | Many directly relevant tests were pending October 25 | Next action |
| No final production equation yet | Larger-data/tensor/new-date evidence incomplete | Open |

This record reflects decisions and constraints evidenced in the project history. Proposed future work is not a claim that new jobs were submitted.
