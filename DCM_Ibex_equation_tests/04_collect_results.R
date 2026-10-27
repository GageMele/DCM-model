#!/usr/bin/env Rscript
# Reporting must retain missing fits and must also work with zero completions.
collect_equation_results <- function(experiment) {
  cases <- read.delim(file.path(experiment, "cases.tsv"), stringsAsFactors = FALSE)
  manifest <- read.delim(file.path(experiment, "manifest.tsv"), stringsAsFactors = FALSE)
  summaries <- statuses <- list()
  # Accounting is optional locally; on Ibex preserve TIMEOUT/OOM/CANCELLED states.
  accounting <- data.frame()
  if (nzchar(Sys.which("sacct"))) {
    ids <- paste(unique(manifest$fit_job), collapse = ",")
    lines <- tryCatch(suppressWarnings(system2("sacct", c("-X", "-n", "-P", "-j", ids,
      "--format=JobIDRaw,State,ExitCode,Elapsed"), stdout = TRUE, stderr = FALSE)),
      error = function(e) character())
    if (length(lines) && is.null(attr(lines, "status"))) {
      accounting <- tryCatch(read.table(text = lines, sep = "|", fill = TRUE,
        col.names = c("job_id", "state", "exit_code", "elapsed", "unused"),
        stringsAsFactors = FALSE), error = function(e) data.frame())
      if (nrow(accounting)) write.csv(accounting, file.path(experiment, "fit_accounting.csv"), row.names = FALSE)
    }
  }
  for (i in seq_len(nrow(cases))) {
    run <- file.path(experiment, "runs", sprintf("case_%02d", cases$case_id[i]))
    summary <- file.path(run, "test_summary.csv")
    diagnostics <- file.path(run, "diagnostic_summary.csv")
    fit_error <- file.path(run, "FIT_FAILED.txt")
    check_error <- file.path(run, "CHECKS_FAILED.txt")
    jobs <- manifest[manifest$case_id == cases$case_id[i], , drop = FALSE]
    state <- "unknown"
    if (nrow(accounting) && nrow(jobs) == 1L) {
      hit <- accounting[as.character(accounting$job_id) == as.character(jobs$fit_job), , drop = FALSE]
      if (nrow(hit)) state <- hit$state[1]
    }
    status <- if (file.exists(summary)) "complete" else if (file.exists(fit_error))
      "fit_failed" else if (file.exists(check_error)) "checks_failed" else "incomplete"
    note <- if (file.exists(fit_error)) paste(readLines(fit_error), collapse = " ") else
      if (file.exists(check_error)) paste(readLines(check_error), collapse = " ") else ""
    statuses[[i]] <- data.frame(case_id = cases$case_id[i], label = cases$label[i],
      equation = cases$equation[i], status = status, slurm_state = state,
      fit_job = if (nrow(jobs) == 1L) jobs$fit_job else NA, note = note)
    available <- if (file.exists(summary)) summary else diagnostics
    if (file.exists(available)) summaries[[length(summaries) + 1L]] <-
      read.csv(available, stringsAsFactors = FALSE)
  }
  status <- do.call(rbind, statuses)
  write.csv(status, file.path(experiment, "case_status.csv"), row.names = FALSE)
  print(status, row.names = FALSE)
  if (!length(summaries)) {
    write.csv(cases[0, ], file.path(experiment, "all_results.csv"), row.names = FALSE)
    writeLines("No completed fit diagnostics. Inspect case_status.csv and Slurm states.",
      file.path(experiment, "comparison.txt"))
    writeLines("No completed fit summaries.", file.path(experiment, "COLLECTION_INCOMPLETE.txt"))
    return(invisible(FALSE))
  }
  all_names <- unique(unlist(lapply(summaries, names)))
  summaries <- lapply(summaries, function(x) {
    for (name in setdiff(all_names, names(x))) x[[name]] <- NA
    x[all_names]
  })
  results <- do.call(rbind, summaries)
  results <- results[order(results$case_id), ]
  write.csv(results, file.path(experiment, "all_results.csv"), row.names = FALSE)
  columns <- c("case_id", "equation", "label", "chains", "iter", "warmup",
    "adapt_delta", "treedepth", "fit_seconds", "sampling_seconds", "max_rhat",
    "min_bulk_ess", "min_tail_ess", "divergences", "treedepth_hits", "diagnostic_screen_pass")
  for (name in setdiff(columns, names(results))) results[[name]] <- NA
  capture.output(print(results[columns], row.names = FALSE), file = file.path(experiment, "comparison.txt"))
  cat(paste(readLines(file.path(experiment, "comparison.txt")), collapse = "\n"), "\n")
  write.csv(results[results$diagnostic_screen_pass %in% TRUE, ],
    file.path(experiment, "diagnostic_passes.csv"), row.names = FALSE)
  # These pairs differ in equation only; report matching hardware separately.
  pairs <- list()
  for (pair in list(c(0L, 6L), c(0L, 12L), c(6L, 12L))) for (i in 1:6) {
    a <- results[results$case_id == i + pair[1], , drop = FALSE]
    b <- results[results$case_id == i + pair[2], , drop = FALSE]
    if (nrow(a) != 1L || nrow(b) != 1L) next
    fields <- c("data_md5", "prior", "chains", "threads", "iter", "warmup", "adapt_delta", "treedepth")
    matched <- all(vapply(fields, function(name) identical(a[[name]], b[[name]]), logical(1)))
    pairs[[length(pairs) + 1L]] <- data.frame(case_a = a$case_id, case_b = b$case_id,
      equation_a = a$equation, equation_b = b$equation,
      matched_input_and_sampler = matched,
      same_cpu_model = identical(a$cpu_model, b$cpu_model) && a$cpu_model != "unknown",
      fit_seconds_a = a$fit_seconds, fit_seconds_b = b$fit_seconds,
      posterior_tree_hits_a = a$treedepth_hits, posterior_tree_hits_b = b$treedepth_hits)
  }
  if (length(pairs)) write.csv(do.call(rbind, pairs),
    file.path(experiment, "matched_equation_pairs.csv"), row.names = FALSE)
  if (file.exists(file.path(experiment, "COLLECTION_INCOMPLETE.txt")))
    file.remove(file.path(experiment, "COLLECTION_INCOMPLETE.txt"))
  writeLines(c("Compare equations using predictive behavior and sampling quality, not runtime alone.",
    "Training predictive checks are in-sample, not independent validation.",
    "Check coefficient R-hat and ESS, zero post-warmup divergences, tree hits and E-BFMI.",
    "The four-chain 2000/1000 cases have a larger retained-draw budget and are not pure speed comparisons.",
    "All formulas use raw numeric MLD and distance and exactly the requested k=15 depth terms."),
    file.path(experiment, "READ_RESULTS_FIRST.txt"))
  writeLines("Reporting finished; see case_status.csv for incomplete/failed fits.",
    file.path(experiment, "COLLECTION_COMPLETE.txt"))
  invisible(TRUE)
}
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  stopifnot(length(args) == 1L)
  collect_equation_results(normalizePath(args[1], winslash = "/", mustWork = TRUE))
}
