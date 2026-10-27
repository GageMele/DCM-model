# Build deliberately incomplete synthetic results to test the collector's failure handling.
package <- "DCM_Ibex_fast_tests"
root <- "validation_fast_tests/collector"
dir.create(root, recursive = TRUE, showWarnings = FALSE)
cases <- read.delim(file.path(package, "cases.tsv"), stringsAsFactors = FALSE)
write.table(cases, file.path(root, "cases.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
manifest <- data.frame(case_id = c(3, 21), fit_job = c("synthetic_fit3", "synthetic_fit21"),
                      check_job = c("synthetic_check3", "synthetic_check21"))
write.table(manifest, file.path(root, "manifest.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
for (id in c(3L, 21L)) {
  run <- file.path(root, "runs", sprintf("case_%02d", id))
  dir.create(run, recursive = TRUE, showWarnings = FALSE)
  result <- data.frame(cases[cases$case_id == id, ],
    cpu_model = "synthetic_cpu", data_md5 = "synthetic_data", source_md5 = "synthetic_source",
    casts = 108, observations = 15228, fit_seconds = 10,
    sampling_seconds = if (id == 3) 5 else 4,
    target_ess_per_sampling_second = 20, max_rhat = 1.2, min_bulk_ess = 100,
    divergences = 5, treedepth_hits = 10, diagnostic_screen_pass = FALSE)
  if (id == 3) {
    result$training_rmse <- 1
    write.csv(result, file.path(run, "test_summary.csv"), row.names = FALSE)
  } else {
    write.csv(result, file.path(run, "diagnostic_summary.csv"), row.names = FALSE)
    writeLines("Synthetic plotting failure", file.path(run, "CHECKS_FAILED.txt"))
  }
}
cat(normalizePath(root, winslash = "/"), "\n")
