#!/usr/bin/env Rscript
# Report all cases, including missing fits, and compare matched sensitivity pairs.
collect_baseline_results <- function(experiment) {
  cases <- read.delim(file.path(experiment, "cases.tsv"), stringsAsFactors = FALSE)
  manifest <- read.delim(file.path(experiment, "manifest.tsv"), stringsAsFactors = FALSE)
  summaries <- statuses <- list()
  accounting <- data.frame()
  if (nzchar(Sys.which("sacct")) && nrow(manifest)) {
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
      subset = cases$subset[i], prior = cases$prior[i], variant = cases$variant[i],
      status = status, slurm_state = state,
      fit_job = if (nrow(jobs) == 1L) jobs$fit_job else NA, note = note)
    available <- if (file.exists(summary)) summary else diagnostics
    if (file.exists(available)) summaries[[length(summaries) + 1L]] <-
      read.csv(available, stringsAsFactors = FALSE)
  }
  status <- do.call(rbind, statuses)
  write.csv(status, file.path(experiment, "case_status.csv"), row.names = FALSE)
  print(status, row.names = FALSE)
  complete_marker <- file.path(experiment, "COLLECTION_COMPLETE.txt")
  if (file.exists(complete_marker)) file.remove(complete_marker)
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
  columns <- c("case_id", "label", "subset", "prior", "variant", "observations", "casts",
    "iter", "warmup", "adapt_delta", "treedepth", "fit_seconds", "sampling_seconds",
    "target_ess_per_sampling_second", "max_rhat", "min_bulk_ess", "min_tail_ess",
    "divergences", "treedepth_hits", "treedepth_hit_fraction", "mean_leapfrogs",
    "min_ebfmi", "diagnostic_screen_pass")
  for (name in setdiff(columns, names(results))) results[[name]] <- NA
  capture.output(print(results[columns], row.names = FALSE), file = file.path(experiment, "comparison.txt"))
  cat(paste(readLines(file.path(experiment, "comparison.txt")), collapse = "\n"), "\n")
  write.csv(results[results$diagnostic_screen_pass %in% TRUE, ],
    file.path(experiment, "diagnostic_passes.csv"), row.names = FALSE)
  pair_definitions <- data.frame(a = c(seq(1L, 11L, by = 2L), rep(1:2, 3L)),
    b = c(seq(2L, 12L, by = 2L), 7:12),
    comparison = c(rep("custom_vs_default", 6L), rep(c("tree12", "delta99", "warmup2500"), each = 2L)))
  pairs <- list()
  for (i in seq_len(nrow(pair_definitions))) {
    a <- results[results$case_id == pair_definitions$a[i], , drop = FALSE]
    b <- results[results$case_id == pair_definitions$b[i], , drop = FALSE]
    if (nrow(a) != 1L || nrow(b) != 1L) next
    fields <- c("data_md5", "source_md5", "subset", "prior", "chains", "cores", "threads",
                "iter", "warmup", "adapt_delta", "treedepth", "refresh", "k_depth")
    changed <- switch(pair_definitions$comparison[i],
      custom_vs_default = "prior", tree12 = "treedepth", delta99 = "adapt_delta",
      warmup2500 = c("warmup", "iter"))
    fields <- setdiff(fields, changed)
    matched <- all(vapply(fields, function(f)
      f %in% names(a) && f %in% names(b) && !anyNA(c(a[[f]], b[[f]])) &&
        identical(a[[f]], b[[f]]), logical(1)))
    same_cpu <- all(c("cpu_model") %in% names(a)) &&
      !anyNA(c(a$cpu_model, b$cpu_model)) && identical(a$cpu_model, b$cpu_model) &&
      a$cpu_model != "unknown"
    pairs[[length(pairs) + 1L]] <- data.frame(case_a = a$case_id, case_b = b$case_id,
      comparison = pair_definitions$comparison[i], subset = a$subset,
      matched_other_settings_and_data = matched, same_cpu_model = same_cpu,
      fit_seconds_a = a$fit_seconds, fit_seconds_b = b$fit_seconds,
      sampling_seconds_a = a$sampling_seconds, sampling_seconds_b = b$sampling_seconds,
      posterior_tree_hit_fraction_a = a$treedepth_hit_fraction,
      posterior_tree_hit_fraction_b = b$treedepth_hit_fraction,
      max_rhat_a = a$max_rhat, max_rhat_b = b$max_rhat,
      min_bulk_ess_a = a$min_bulk_ess, min_bulk_ess_b = b$min_bulk_ess)
  }
  if (length(pairs)) write.csv(do.call(rbind, pairs),
    file.path(experiment, "matched_baseline_pairs.csv"), row.names = FALSE)
  if (file.exists(file.path(experiment, "COLLECTION_INCOMPLETE.txt")))
    file.remove(file.path(experiment, "COLLECTION_INCOMPLETE.txt"))
  writeLines(c("All cases retain 4000 post-warmup draws per chain, 16000 total.",
    "Compare baseline and sensitivity cases within the same subset and prior choice.",
    "Defaults are resolved by the installed brms version and may depend on the data subset.",
    "Custom priors are the exact historical constants; both choices are saved before/after fitting.",
    "Tree saturation is an efficiency concern; examine R-hat, ESS, divergences and E-BFMI separately.",
    "A diagnostic-screen pass does not establish scientific validity.",
    "Training PPC checks are in-sample; use held-out casts for independent model comparison.",
    "Other campaigns have different retained-draw and threading budgets: raw fit time alone is not a fair comparison."),
    file.path(experiment, "READ_RESULTS_FIRST.txt"))
  writeLines("Reporting finished; see case_status.csv for incomplete/failed fits.", complete_marker)
  invisible(TRUE)
}
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  stopifnot(length(args) == 1L)
  collect_baseline_results(normalizePath(args[1], winslash = "/", mustWork = TRUE))
}
