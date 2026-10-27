# Validate the real-input basis construction and all controlled comparisons.
package <- "DCM_Ibex_k_tests"
for (path in list.files(package, "\\.R$", full.names = TRUE)) parse(path)
source(file.path(package, "common.R"))
cases <- read.delim(file.path(package, "cases.tsv"), stringsAsFactors = FALSE)
stopifnot(nrow(cases) == 12, !anyDuplicated(cases$case_id),
  sum(cases$cpus) == 96, all(cases$cpus == cases$chains * cases$threads))
assert_matched <- function(ids, allowed) {
  x <- cases[cases$case_id %in% ids, , drop = FALSE]
  fixed <- setdiff(names(x), c("case_id", "label", allowed))
  stopifnot(all(vapply(x[fixed], function(z) length(unique(z)) == 1L, logical(1))))
}
for (id in c(2L, 4L, 6L, 8L)) assert_matched(c(id, id + 1L), "treedepth")
assert_matched(c(1, 2, 4, 6), "k_depth")
assert_matched(c(8, 11), c("warmup", "iter"))
stopifnot(cases$iter[8] - cases$warmup[8] == cases$iter[11] - cases$warmup[11])
assert_matched(c(8, 12), "adapt_delta")
stopifnot(inherits(try(case_k(data.frame(k_depth = 2, k_distance = 4, k_mld = 4)),
  silent = TRUE), "try-error"))

# Run the same preflight as Ibex, using the previously validated actual p10 data.
experiment <- file.path(getwd(), "validation_k150", "actual")
dir.create(file.path(experiment, "data"), recursive = TRUE, showWarnings = FALSE)
file.copy(list.files(package, "\\.R$", full.names = TRUE), experiment, overwrite = TRUE)
file.copy(file.path(package, "cases.tsv"), experiment, overwrite = TRUE)
file.copy("DCM_Ibex_workflow/02_fit_model.R", experiment, overwrite = TRUE)
source_data <- "validation_fast_tests/actual150_radix/data/p10.rds"
stopifnot(file.copy(source_data, file.path(experiment, "data", "p10.rds"), overwrite = TRUE))
stopifnot(unname(tools::md5sum(source_data)) ==
  unname(tools::md5sum(file.path(experiment, "data", "p10.rds"))))
Sys.setenv(EXPERIMENT_DIR = experiment)
source(file.path(experiment, "01_validate_campaign.R"))
report <- read.csv(file.path(experiment, "basis_validation.csv"))
stopifnot(nrow(report) == 6L, all(report$observations == 15228L),
  all(report$casts == 108L), all(report$stan_penalized_spline_columns > 0L),
  all(report$stan_fixed_effect_columns > 0L))
print(report)
cat("Actual data, complete brms tensor construction, and matched comparisons passed.\n")

# Exercise the collector on complete, diagnostic-only and failed cases.
collector <- file.path(getwd(), "validation_k150", "collector")
dir.create(file.path(collector, "runs"), recursive = TRUE, showWarnings = FALSE)
write.table(cases, file.path(collector, "cases.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
manifest <- data.frame(case_id = cases$case_id, label = cases$label,
  fit_job = 1000L + cases$case_id, check_job = 2000L + cases$case_id,
  cpus = cases$cpus, run_dir = file.path(collector, "runs", sprintf("case_%02d", cases$case_id)))
write.table(manifest, file.path(collector, "manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
for (i in c(2, 3)) {
  run <- manifest$run_dir[i]
  dir.create(run, showWarnings = FALSE)
  row <- data.frame(cases[i, ], data_md5 = "frozen-input", source_md5 = "frozen-model",
    cpu_model = "mock CPU", fit_seconds = 1000 * i, sampling_seconds = 400 * i,
    max_rhat = 1, min_bulk_ess = 500, min_tail_ess = 500,
    divergences = 0, treedepth_hits = 0, diagnostic_screen_pass = TRUE)
  write.csv(row, file.path(run, if (i == 2) "test_summary.csv" else "diagnostic_summary.csv"), row.names = FALSE)
  if (i == 3) writeLines("mock plotting failure", file.path(run, "CHECKS_FAILED.txt"))
}
dir.create(manifest$run_dir[4], showWarnings = FALSE)
writeLines("mock fitting failure", file.path(manifest$run_dir[4], "FIT_FAILED.txt"))
path <- file.path(package, "04_collect_results.R")
env <- new.env(parent = globalenv())
env$commandArgs <- function(trailingOnly = FALSE) collector
sys.source(path, envir = env)
status <- read.csv(file.path(collector, "case_status.csv"))
pairs <- read.csv(file.path(collector, "matched_treedepth_pairs.csv"))
stopifnot(status$status[2] == "complete", status$status[3] == "checks_failed",
  status$status[4] == "fit_failed", status$status[1] == "incomplete",
  nrow(pairs) == 1L, pairs$matched, pairs$same_cpu_model)
cat("Collector preserves failures, partial diagnostics and matched tree-depth pairs.\n")
