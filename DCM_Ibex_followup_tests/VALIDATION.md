# Validation status

The installer contains34 cases:30 posterior and4 prior-only. Python checks passed for case IDs, unique labels, resource arithmetic, block counts, sampler iterations, matched variance controls, prior controls, and lack of exact posterior duplicates among231 earlier registered cases.

Existing completed data-only audit files were independently inspected: original108 and alternate108 training casts form the216-cast pool; reserve108 is outside primary540, date-blocked training108, and both previous holdouts108/117; date groups remain disjoint for date-blocked training/validation. All copied source files retain identical MD5 checksums. The completed data-only R audit compared raw response/predictor rows before writing its completion marker; it did not fit models.

R delimiters/strings and installer payload/checksums were checked without executing R. These are static checks, not a complete R parser or a substitute for running the included tests.

Automatic approval review blocked local R scheduling/reporting validation after detecting an attempted HTTPS connection to an external Visual Studio telemetry endpoint. It was not retried or bypassed. Windows also prevented local Bash startup with an OS named-object permission error. Local R scheduler tests and full Bash syntax validation are therefore not claimed as passed.

On Ibex, the installer runs bash -n and all base-R scheduling/reporting recovery tests before submitting jobs. The preparation job then constructs every actual formula/prior/data combination before releasing fits through afterok dependencies. Failure prevents fitting. These runtime gates are required, not optional.

Local data audit datasets: {"date_blocked_p10_d150": 108, "date_blocked_reserve_d150": 108, "date_blocked_validation_d150": 117, "primary_p10_alt_d150": 108, "primary_p10_d150": 108, "primary_p20_d150": 216, "primary_p50_d150": 540, "primary_reserve_d150": 108, "primary_validation_d150": 108}

Aggregate requested capacity if every fitting job ran simultaneously: 296 CPUs and1440GB. Slurm controls actual concurrency. No runtime/total-cost forecast is implied.
