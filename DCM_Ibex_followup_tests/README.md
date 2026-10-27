# DCM follow-up package

Read PLAN.md for the34-case design. Run `bash queue_followup_dcm.sh` from the project root on Ibex after uploading the single-file installer. Repeating it resumes the same experiment's job registrations; do not set FOLLOWUP_NEW unless intentionally starting another full campaign.

The installer verifies its embedded payload and package checksums. Before scheduling, it runs the base-R scheduling/reporting tests on Ibex. Then a preparation job validates the real model specifications, copied data, new holdout and priors before any fit can start. Failed preparation prevents fitting.

The monitor includes all previously recorded campaigns. Results include every registered fit rather than only low-error or diagnostic-passing models. Old fits are not refitted; completed older common-validation evaluations are merged with their source fits.

Local R scheduling/reporting validation was blocked by automatic approval review after an attempted connection to an external Visual Studio telemetry endpoint. Those checks have been moved to the authorized Linux Ibex environment. The package includes existing completed data-only audit counts; model construction is still verified by the Ibex preparation job. See VALIDATION.md for actual local checks and limits.
