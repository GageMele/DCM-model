# Glossary

| Term | Project meaning |
|---|---|
| Cast | A CTD vertical profile; identity includes source metadata, not merely a row number |
| DCM | Deep chlorophyll maximum; operational detection definitions differ across helpers |
| chl_norm | Fractional enhancement (Chl−surface_chl)/surface_chl |
| MLD | Mixed-layer depth from MLD_temp_MIX |
| Region | Fixed mean contrasts and usually distinct depth curves |
| ENSO | Historical categorical grouping encoded from year |
| Numeric by smooth | Predictor multiplied by a depth smooth, e.g. M g(z) |
| Factor by smooth | Separate smooth for each Region |
| Tensor/t2 | Smooth surface over multiple coordinates |
| Half/full shift | Evaluate the depth curve at z−0.5m or z−m |
| Warp | Prescribed coordinate scaling changing physical width |
| k | Basis dimension; not fitted effective degrees of freedom |
| sds | Standard deviation governing smooth coefficients/penalties |
| sd | Group-level effect scale |
| sigma | Observation residual SD; can have its own regression |
| R-hat | Between/within-chain diagnostic; use rank-normalized version |
| ESS | Effective sample size, distinct from retained draw count |
| E-BFMI | Energy-transition diagnostic |
| Divergence | HMC integration/geometry warning |
| Tree hit | Sampler reaches its configured maximum trajectory depth |
| Warmup | Adaptation phase, not retained posterior draws |
| adapt_delta | Target acceptance control affecting adapted step size |
| Conditional prediction | Includes group effects under recorded handling of new levels |
| Population prediction here | Group effects set to zero, not marginal integration |
| Reconstruction | Prediction on data included in fitting |
| Cast-balanced RMSE | Square root of average cast MSE |
| Shape RMSE | Error after per-cast demeaning |
| Reserve | Separate development cohort repeatedly inspected |
| Left censoring | Enhancement already above threshold at earliest supported search point |
| Inadequate support | Data do not permit locating a profile feature |
| Manifest/ledger | Case/job registrations and accepted submission state |
| Controller/dispatcher | Lightweight scheduler that submits fits/checks as capacity permits |
| Frozen campaign | Copied code, cases, data selection and identities for a named experiment |
