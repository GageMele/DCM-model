# Regression check for the confirmed prior-only/zero-row collector failure.
root <- file.path(getwd(), "validation_equation150", "report_repair")
dir.create(file.path(root, "runs"), recursive = TRUE, showWarnings = FALSE)
cases <- read.delim("DCM_Ibex_fast_tests/cases.tsv", stringsAsFactors = FALSE)
write.table(cases, file.path(root, "cases.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
manifest <- data.frame(case_id = cases$case_id, label = cases$label,
  fit_job = 1000 + cases$case_id, check_job = 2000 + cases$case_id,
  cpus = cases$cpus, run_dir = file.path(root, "runs", sprintf("case_%02d", cases$case_id)))
write.table(manifest, file.path(root, "manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
file.copy("DCM_Ibex_fast_tests/04_collect_results.R", root, overwrite = TRUE)
for (i in 23:24) {
  dir.create(manifest$run_dir[i], showWarnings = FALSE)
  row <- data.frame(cases[i, ], casts = 108, observations = 15228,
    fit_seconds = 200, sampling_seconds = 2,
    target_ess_per_sampling_second = 300, max_rhat = 1.05, min_bulk_ess = 160,
    divergences = 0, treedepth_hits = 0, diagnostic_screen_pass = FALSE,
    cpu_model = "mock CPU", data_md5 = "same", source_md5 = "same")
  write.csv(row, file.path(manifest$run_dir[i], "test_summary.csv"), row.names = FALSE)
}
original <- file.path(root, "04_collect_results.R")
original_hash <- tools::md5sum(original)
env <- new.env(parent = globalenv())
env$commandArgs <- function(trailingOnly = FALSE) root
failure <- try(sys.source(original, envir = env), silent = TRUE)
stopifnot(inherits(failure, "try-error"), grepl("replacement has 1 row", failure, fixed = TRUE))

# Execute the exact embedded R repair, substituting only the child Rscript launch.
lines <- readLines("DCM_Ibex_equation_tests/repair_fast150_report.sh")
start <- grep("^Rscript --vanilla - <<'RSCRIPT'$", lines) + 1L
end <- which(lines == "RSCRIPT") - 1L
repair_code <- paste(lines[start:end], collapse = "\n")
repair_env <- new.env(parent = globalenv())
repair_env$system2 <- function(command, args) {
  check_env <- new.env(parent = globalenv())
  check_env$commandArgs <- function(trailingOnly = FALSE) root
  sys.source(file.path(root, "04_collect_results_repaired.R"), envir = check_env)
  0L
}
Sys.setenv(REPORT_EXPERIMENT = root)
eval(parse(text = repair_code), envir = repair_env)
stopifnot(identical(original_hash, tools::md5sum(original)),
  file.exists(file.path(root, "COLLECTION_COMPLETE.txt")),
  nrow(read.csv(file.path(root, "all_results.csv"))) == 2L,
  nrow(read.csv(file.path(root, "matched_threading_results.csv"))) == 0L)
cat("Confirmed zero-row failure reproduced; repaired report succeeds and original script is unchanged.\n")
