#!/usr/bin/env Rscript
# Collect successes and failures without confusing different models with CPU tests.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
experiment <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
cases <- read.delim(file.path(experiment, "cases.tsv"), stringsAsFactors = FALSE)
manifest <- read.delim(file.path(experiment, "manifest.tsv"), stringsAsFactors = FALSE)
summaries <- statuses <- list()
for (i in seq_len(nrow(cases))) {
  run <- file.path(experiment, "runs", sprintf("case_%02d", cases$case_id[i]))
  summary <- file.path(run, "test_summary.csv")
  diagnostic_summary <- file.path(run, "diagnostic_summary.csv")
  fit_error <- file.path(run, "FIT_FAILED.txt")
  check_error <- file.path(run, "CHECKS_FAILED.txt")
  status <- if (file.exists(summary)) "complete" else if (file.exists(fit_error))
    "fit_failed" else if (file.exists(check_error)) "checks_failed" else "incomplete"
  note <- if (file.exists(fit_error)) paste(readLines(fit_error), collapse = " ") else
    if (file.exists(check_error)) paste(readLines(check_error), collapse = " ") else ""
  jobs <- manifest[manifest$case_id == cases$case_id[i], ]
  statuses[[i]] <- data.frame(case_id = cases$case_id[i], label = cases$label[i],
    status = status,
    fit_job = if (nrow(jobs) == 1L) jobs$fit_job else NA,
    check_job = if (nrow(jobs) == 1L) jobs$check_job else NA, note = note)
  available <- if (file.exists(summary)) summary else diagnostic_summary
  if (file.exists(available)) summaries[[length(summaries) + 1L]] <-
    read.csv(available, stringsAsFactors = FALSE)
}
status_table <- do.call(rbind, statuses)
write.csv(status_table, file.path(experiment, "case_status.csv"), row.names = FALSE)
print(status_table, row.names = FALSE)
if (!length(summaries)) stop("No complete results. Inspect case_status.csv and Slurm logs.")
columns_union <- unique(unlist(lapply(summaries, names)))
summaries <- lapply(summaries, function(x) {
  for (name in setdiff(columns_union, names(x))) x[[name]] <- NA
  x[columns_union]
})
results <- do.call(rbind, summaries)
results <- results[order(results$mode, results$subset, results$fit_seconds), ]
write.csv(results, file.path(experiment, "all_results.csv"), row.names = FALSE)
eligible <- results[results$mode == "posterior" & results$diagnostic_screen_pass, ]
write.csv(eligible, file.path(experiment, "diagnostic_passes.csv"), row.names = FALSE)
columns <- c("case_id", "label", "casts", "observations", "fit_seconds", "sampling_seconds",
  "target_ess_per_sampling_second", "max_rhat", "min_bulk_ess", "divergences",
  "treedepth_hits", "diagnostic_screen_pass")
capture.output(print(results[, columns], row.names = FALSE),
  file = file.path(experiment, "comparison.txt"))
cat("\n", paste(readLines(file.path(experiment, "comparison.txt")), collapse = "\n"), "\n")
# Hardware scaling is evaluated only on identical data/model/sampler settings.
hardware <- results[results$mode == "posterior" & results$subset == "p10" &
  results$prior == "current" & results$chains == 2 & results$iter == 1000 &
  results$warmup == 500 & results$adapt_delta == .95 & results$treedepth == 10, ]
baseline <- hardware[hardware$threads == 0, ]
hardware$same_cpu_model <- FALSE
hardware$sampling_speedup <- NA_real_
if (nrow(baseline) == 1) {
  same <- hardware$cpu_model == baseline$cpu_model & hardware$cpu_model != "unknown" &
    hardware$data_md5 == baseline$data_md5 & hardware$source_md5 == baseline$source_md5
  hardware$same_cpu_model <- same
  hardware$sampling_speedup[same] <- baseline$sampling_seconds / hardware$sampling_seconds[same]
}
write.csv(hardware, file.path(experiment, "matched_threading_results.csv"), row.names = FALSE)
writeLines(c("Do not rank prior-only runs as faster posterior fits.",
  "Different priors define different models; sampler/subset settings change the experiment.",
  "Compare matched rows and inspect diagnostics, not runtime alone.",
  "Hardware speedups are calculated only for matching recorded CPU models.",
  "Empirical priors use full 150m response statistics; training PPC is not validation.",
  "A lower tree-depth cap can make a run faster while reducing effective exploration.",
  "Resolve divergences, tree hits, poor R-hat, and low ESS before production."),
  file.path(experiment, "READ_RESULTS_FIRST.txt"))
writeLines("Complete", file.path(experiment, "COLLECTION_COMPLETE.txt"))
