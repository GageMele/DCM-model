#!/usr/bin/env Rscript
# Keep completed and failed cases visible; never select a model on speed alone.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
experiment <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
cases <- read.delim(file.path(experiment, "cases.tsv"), stringsAsFactors = FALSE)
manifest <- read.delim(file.path(experiment, "manifest.tsv"), stringsAsFactors = FALSE)
summaries <- statuses <- list()
for (i in seq_len(nrow(cases))) {
  run <- file.path(experiment, "runs", sprintf("case_%02d", cases$case_id[i]))
  summary <- file.path(run, "test_summary.csv")
  diagnostics <- file.path(run, "diagnostic_summary.csv")
  fit_error <- file.path(run, "FIT_FAILED.txt")
  check_error <- file.path(run, "CHECKS_FAILED.txt")
  state <- if (file.exists(summary)) "complete" else if (file.exists(fit_error))
    "fit_failed" else if (file.exists(check_error)) "checks_failed" else "incomplete"
  note <- if (file.exists(fit_error)) paste(readLines(fit_error), collapse = " ") else
    if (file.exists(check_error)) paste(readLines(check_error), collapse = " ") else ""
  jobs <- manifest[manifest$case_id == cases$case_id[i], , drop = FALSE]
  statuses[[i]] <- data.frame(case_id = cases$case_id[i], label = cases$label[i],
    status = state, fit_job = if (nrow(jobs) == 1L) jobs$fit_job else NA,
    check_job = if (nrow(jobs) == 1L) jobs$check_job else NA, note = note)
  available <- if (file.exists(summary)) summary else diagnostics
  if (file.exists(available)) summaries[[length(summaries) + 1L]] <-
    read.csv(available, stringsAsFactors = FALSE)
}
status <- do.call(rbind, statuses)
write.csv(status, file.path(experiment, "case_status.csv"), row.names = FALSE)
print(status, row.names = FALSE)
if (!length(summaries)) stop("No completed diagnostic summaries; inspect Slurm accounting and logs.")
all_names <- unique(unlist(lapply(summaries, names)))
summaries <- lapply(summaries, function(x) {
  for (name in setdiff(all_names, names(x))) x[[name]] <- NA
  x[all_names]
})
results <- do.call(rbind, summaries)
results <- results[order(results$case_id), ]
write.csv(results, file.path(experiment, "all_results.csv"), row.names = FALSE)
columns <- c("case_id", "label", "k_depth", "k_distance", "k_mld", "treedepth",
  "warmup", "adapt_delta", "fit_seconds", "sampling_seconds", "max_rhat",
  "min_bulk_ess", "min_tail_ess", "divergences", "treedepth_hits",
  "diagnostic_screen_pass")
capture.output(print(results[columns], row.names = FALSE),
  file = file.path(experiment, "comparison.txt"))
cat(paste(readLines(file.path(experiment, "comparison.txt")), collapse = "\n"), "\n")
write.csv(results[results$diagnostic_screen_pass %in% TRUE, ],
  file.path(experiment, "diagnostic_passes.csv"), row.names = FALSE)
# Match the tree-depth pairs; they have identical data, k and sampler settings.
pairs <- list()
for (id in c(2L, 4L, 6L, 8L)) {
  a <- results[results$case_id == id, , drop = FALSE]
  b <- results[results$case_id == id + 1L, , drop = FALSE]
  if (nrow(a) != 1L || nrow(b) != 1L) next
  fields <- c("data_md5", "source_md5", "k_depth", "k_distance", "k_mld",
    "prior", "chains", "threads", "iter", "warmup", "adapt_delta")
  matched <- all(vapply(fields, function(field) identical(a[[field]], b[[field]]), logical(1)))
  pairs[[length(pairs) + 1L]] <- data.frame(case_td12 = id, case_td14 = id + 1L,
    matched = matched,
    same_cpu_model = identical(a$cpu_model, b$cpu_model) && a$cpu_model != "unknown",
    fit_seconds_td12 = a$fit_seconds, fit_seconds_td14 = b$fit_seconds,
    posterior_tree_hits_td12 = a$treedepth_hits, posterior_tree_hits_td14 = b$treedepth_hits,
    min_bulk_ess_td12 = a$min_bulk_ess, min_bulk_ess_td14 = b$min_bulk_ess)
}
if (length(pairs)) write.csv(do.call(rbind, pairs),
  file.path(experiment, "matched_treedepth_pairs.csv"), row.names = FALSE)
writeLines(c(
  "All cases use the same frozen 10% cast subset, 10-150 m, without smoothing or binning.",
  "Cases 2/3, 4/5, 6/7 and 8/9 isolate max_treedepth 12 versus 14.",
  "Cases 1/2/4/6 isolate depth k at tree depth 12; other marginal k values stay at 5.",
  "Cases 8/11 isolate warmup length while retaining 500 sampling iterations per chain.",
  "Cases 8/12 isolate adapt_delta .95 versus .99.",
  "Changing k changes the model and its induced function prior; numerical prior constants stay fixed.",
  "Warmup divergences and tree hits are not post-warmup diagnostic failures.",
  "Inspect posterior profiles, predictive checks, R-hat, ESS, divergences and tree hits together.",
  "Parameter ESS summaries across different k involve different parameter sets.",
  "In-sample predictive checks are not independent predictive validation.",
  "The shortest runtime does not establish that a smaller basis is scientifically adequate.",
  "All incomplete cases remain in case_status.csv; obtain timeout/OOM states from sacct."
), file.path(experiment, "READ_RESULTS_FIRST.txt"))
writeLines("Complete", file.path(experiment, "COLLECTION_COMPLETE.txt"))
