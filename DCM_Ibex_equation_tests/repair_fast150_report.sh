#!/bin/bash
# Repair only the failed old report; do not touch fits, data, or Slurm jobs.
set -euo pipefail
cd "${DCM_ROOT:-/ibex/user/meleg/DCM_MODEL}"
export REPORT_EXPERIMENT="${1:-$(cat logs/latest_fast150_experiment.txt)}"
source "$REPORT_EXPERIMENT/ibex_environment.sh"
Rscript --vanilla - <<'RSCRIPT'
experiment <- Sys.getenv("REPORT_EXPERIMENT")
original <- file.path(experiment, "04_collect_results.R")
code <- paste(readLines(original), collapse = "\n")
old <- c("hardware$same_cpu_model <- FALSE", "hardware$sampling_speedup <- NA_real_")
new <- c("hardware$same_cpu_model <- rep(FALSE, nrow(hardware))",
         "hardware$sampling_speedup <- rep(NA_real_, nrow(hardware))")
for (i in seq_along(old)) {
  if (grepl(old[i], code, fixed = TRUE)) code <- sub(old[i], new[i], code, fixed = TRUE) else
    if (!grepl(new[i], code, fixed = TRUE)) stop("Unexpected collector version; original left unchanged.")
}
invisible(parse(text = code))  # Validate before writing a separate repaired script.
repaired <- file.path(experiment, "04_collect_results_repaired.R")
writeLines(code, repaired)
cat("Original preserved. Repaired report script:", repaired, "\n")
status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(repaired), shQuote(experiment)))
if (status != 0L) stop("Repaired report returned status ", status)
RSCRIPT
