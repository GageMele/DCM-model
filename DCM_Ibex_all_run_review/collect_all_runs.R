#!/usr/bin/env Rscript
# Read existing small results only. No refitting, prediction, or fit.rds loading.
collect_all <- function(root, output) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  if (dir.exists(output)) stop("Choose a new output directory: ", output)
  dir.create(output, recursive = TRUE)
  output <- normalizePath(output, winslash = "/", mustWork = TRUE)
  warnings <- character()
  note <- function(x) warnings <<- c(warnings, x)
  read_table <- function(path, sep = ",") {
    if (!file.exists(path)) return(data.frame())
    tryCatch(read.table(path, header = TRUE, sep = sep, quote = "\"", comment.char = "",
      stringsAsFactors = FALSE, check.names = FALSE, fill = FALSE), error = function(e) {
        note(paste("Unreadable table:", path, conditionMessage(e))); data.frame()
      })
  }
  bind <- function(xs) {
    xs <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, xs)
    if (!length(xs)) return(data.frame())
    cols <- unique(unlist(lapply(xs, names)))
    do.call(rbind, lapply(xs, function(x) {
      for (nm in setdiff(cols, names(x))) x[[nm]] <- NA
      x[cols]
    }))
  }
  csv <- function(x, name) write.csv(x, file.path(output, name), row.names = FALSE, na = "")
  rel <- function(path) {
    path <- normalizePath(path, winslash = "/", mustWork = FALSE)
    prefix <- paste0(root, "/")
    if (startsWith(path, prefix)) substring(path, nchar(prefix) + 1L) else path
  }
  resolve <- function(path) {
    if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) return(path)
    file.path(root, path)
  }
  scan <- function(base, pattern) {
    if (!dir.exists(base)) return(character())
    list.files(base, pattern = pattern, recursive = TRUE, full.names = TRUE)
  }
  bases <- file.path(root, c("experiments", "benchmarks", "runs"))
  manifests <- unique(unlist(lapply(bases, scan, pattern = "^manifest\\.tsv$")))
  registry <- list()
  for (m in sort(manifests)) {
    x <- read_table(m, "\t")
    if (!nrow(x)) next
    campaign <- rel(dirname(m))
    for (i in seq_len(nrow(x))) {
      rec <- as.list(x[i, , drop = FALSE])
      path <- rec$run_dir
      if (is.null(path) || is.na(path) || !nzchar(path)) {
        if (!is.null(rec$case_id) && grepl("^[0-9]+$", rec$case_id))
          path <- file.path(dirname(m), "runs", sprintf("case_%02d", as.integer(rec$case_id)))
        else {
          path <- file.path(dirname(m), paste0("UNRESOLVED_MANIFEST_ROW_", i))
          note(paste("No run_dir for", campaign, "row", i, "-- inspect copied manifest"))
        }
      } else path <- resolve(as.character(path))
      key <- normalizePath(path, winslash = "/", mustWork = FALSE)
      if (!is.null(registry[[key]])) note(paste("Duplicate manifest registration:", key))
      registry[[key]] <- list(run = key, campaign = campaign, manifest = rel(m), record = rec)
    }
  }
  configs <- unique(unlist(lapply(bases, scan, pattern = "^run_config\\.rds$")))
  summaries <- unique(unlist(lapply(bases, scan, pattern = "^test_summary\\.csv$")))
  markers <- unique(unlist(lapply(bases, scan, pattern = "^FIT_COMPLETE\\.txt$")))
  for (path in sort(unique(dirname(c(configs, summaries, markers))))) {
    key <- normalizePath(path, winslash = "/", mustWork = FALSE)
    if (is.null(registry[[key]])) {
      r <- rel(path)
      campaign <- if (grepl("^(experiments|benchmarks)/", r))
        paste(strsplit(r, "/", fixed = TRUE)[[1]][1:2], collapse = "/") else "original_or_unregistered"
      registry[[key]] <- list(run = key, campaign = campaign, manifest = NA_character_, record = list())
    }
  }
  if (length(registry)) registry <- registry[order(names(registry))]
  if (!length(registry)) note("No run manifests/configurations/results found under experiments, benchmarks, or runs.")
  inventory <- list(); long <- list(); copied <- list(); missing <- list()
  hash_cache <- new.env(parent = emptyenv())
  md5 <- function(path) {
    if (!file.exists(path)) return(NA_character_)
    if (!exists(path, hash_cache, inherits = FALSE)) assign(path, unname(tools::md5sum(path)), hash_cache)
    get(path, hash_cache, inherits = FALSE)
  }
  scalar <- function(x) {
    if (is.null(x) || !length(x)) return(NA)
    if (length(x) == 1L && !is.list(x)) return(unname(x))
    paste(unlist(x), collapse = " | ")
  }
  for (i in seq_along(registry)) {
    reg <- registry[[i]]; run <- reg$run; uid <- sprintf("run_%04d", i)
    config <- if (file.exists(file.path(run, "run_config.rds"))) tryCatch(
      readRDS(file.path(run, "run_config.rds")), error = function(e) {
        note(paste("Unreadable run_config:", run, conditionMessage(e))); list()
      }) else list()
    rec <- reg$record
    if (is.data.frame(config$case) && nrow(config$case)) {
      for (nm in names(config$case)) rec[[nm]] <- scalar(config$case[[nm]][1])
    }
    # Config stores actual values; use them when available.
    for (nm in c("chains", "cores", "threads", "iter", "warmup", "adapt_delta", "seed",
      "fit_seconds", "brm_seconds_including_compile", "hostname", "allocated_cpus", "formula",
      "data_md5", "source_md5", "original_source_md5", "input_md5", "cmdstan_version")) {
      if (!is.null(config[[nm]])) rec[[nm]] <- scalar(config[[nm]])
    }
    if (!is.null(config$max_treedepth)) rec$treedepth <- scalar(config$max_treedepth)
    if (!is.null(config$k)) rec$config_k <- scalar(config$k)
    diagnostic_summary <- read_table(file.path(run, "diagnostic_summary.csv"))
    if (nrow(diagnostic_summary)) for (nm in names(diagnostic_summary)) rec[[nm]] <- scalar(diagnostic_summary[[nm]][1])
    summary <- read_table(file.path(run, "test_summary.csv"))
    if (nrow(summary) > 1L) note(paste("Multiple test_summary rows:", run, "-- first used; full file copied"))
    if (nrow(summary)) for (nm in names(summary)) rec[[nm]] <- scalar(summary[[nm]][1])
    exists_run <- function(x) file.exists(file.path(run, x))
    complete <- exists_run("FIT_COMPLETE.txt")
    fit_exists <- exists_run("fit.rds")
    checks_failed <- exists_run("CHECKS_FAILED.txt")
    status <- if (complete && checks_failed) "fit_complete_checks_failed" else if (complete && exists_run("CHECKS_COMPLETE.txt"))
      "fit_and_checks_complete" else if (complete) "fit_complete_checks_unfinished" else if (exists_run("FIT_FAILED.txt"))
      "fit_failed" else if (fit_exists) "fit_file_without_completion_marker" else "incomplete_or_not_started"
    campaign_dir <- resolve(reg$campaign)
    prior <- (!is.null(rec$mode) && !is.na(rec$mode) && grepl("prior", rec$mode, ignore.case = TRUE)) ||
      (!is.null(config$mode) && identical(config$mode, "prior")) ||
      (!is.null(rec$label) && grepl("prior_predictive", rec$label, fixed = TRUE))
    has_validation <- any(grepl("^validation_(conditional|population)_", names(summary)))
    eval_type <- if (prior) "prior_predictive" else if (has_validation) "held_out_casts_plus_training" else
      if ("training_rmse" %in% names(summary)) "legacy_bounded_in_sample" else "prediction_metrics_missing"
    rec$review_uid <- uid; rec$campaign <- reg$campaign; rec$run_path <- rel(run)
    rec$source_manifest <- reg$manifest; rec$snapshot_status <- status
    rec$fit_complete_marker <- complete; rec$fit_file_present <- fit_exists
    rec$summary_available <- nrow(summary) > 0; rec$evaluation_type <- eval_type
    rec$diagnostic_summary_available <- nrow(diagnostic_summary) > 0
    rec$validation_set_md5 <- if (has_validation) md5(file.path(campaign_dir, "data", "validation_model.rds")) else NA_character_
    rec$training_evaluation_set_md5 <- if (has_validation) md5(file.path(campaign_dir, "data", "p10_model.rds")) else NA_character_
    rec$resolved_prior_file_md5 <- md5(file.path(run, "resolved_priors.csv"))
    rec$is_prior_predictive <- prior
    if (!is.null(rec$iter) && !is.null(rec$warmup)) {
      rec$retained_per_chain_configured <- suppressWarnings(as.numeric(rec$iter) - as.numeric(rec$warmup))
    }
    # Do not manufacture pass/fail flags for missing diagnostics or old definitions.
    inventory[[uid]] <- as.data.frame(rec, stringsAsFactors = FALSE, check.names = FALSE)
    run_files <- if (dir.exists(run)) list.files(run, recursive = TRUE, full.names = TRUE) else character()
    # Explicit compact evidence types. Raw chain draws, fit/data RDS and row predictions excluded.
    keep <- grepl("\\.(csv|tsv|txt|png|pdf|stan)$", run_files, ignore.case = TRUE) &
      !grepl("/stan_csv/", gsub("\\\\", "/", run_files)) &
      !grepl("(observations\\.csv$|sampler_iterations_including_warmup\\.csv$|sampler_diagnostics\\.csv$)", run_files)
    for (f in run_files[keep]) {
      bytes <- file.info(f)$size
      destination <- file.path(output, "evidence", uid, substring(gsub("\\\\", "/", f), nchar(gsub("\\\\", "/", run)) + 2L))
      reason <- if (!is.finite(bytes) || bytes > 20 * 1024^2) "over_20_MiB_or_unreadable" else "copied"
      if (reason == "copied") {
        dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
        if (!file.copy(f, destination, overwrite = FALSE)) { reason <- "copy_failed"; note(paste("Copy failed:", f)) }
      }
      copied[[length(copied) + 1L]] <- data.frame(review_uid = uid, source = rel(f), bytes = bytes,
        archive_path = if (reason == "copied") substring(destination, nchar(output) + 2L) else NA_character_, status = reason)
      nm <- basename(f)
      if (reason == "copied" && grepl("\\.csv$", nm) &&
          (grepl("(by_region|strata|cast_metrics|prediction_scores|sampler_summary_by_phase|chain_times|parameter_diagnostics|diagnostic_summary|prior_predictive_quantiles)\\.csv$", nm))) {
        tab <- read_table(f)
        if (nrow(tab)) {
          tab$review_uid <- uid; tab$campaign <- reg$campaign; tab$evidence_table <- nm
          long[[nm]] <- c(long[[nm]], list(tab))
        }
      }
    }
    gaps <- c(if (!complete) "fit_not_marked_complete", if (!nrow(summary)) "test_summary_missing",
      if (!prior && !has_validation) "common_heldout_metrics_unavailable",
      if (!exists_run("resolved_priors.csv")) "resolved_prior_table_unavailable",
      if (checks_failed) "checks_failed_metrics_may_be_partial")
    if (length(gaps)) missing[[uid]] <- data.frame(review_uid = uid, issue = gaps)
  }
  all <- bind(inventory)
  if (nrow(all)) {
    case_order <- if ("case_id" %in% names(all)) suppressWarnings(as.numeric(all$case_id)) else rep(NA_real_, nrow(all))
    all <- all[order(all$campaign, case_order, all$run_path, na.last = TRUE), , drop = FALSE]
    first <- c("review_uid", "campaign", "case_id", "label", "snapshot_status", "evaluation_type", "run_path")
    all <- all[c(intersect(first, names(all)), setdiff(names(all), first))]
  }
  csv(all, "all_registered_runs.csv")
  if (nrow(all)) {
    csv(all[all$fit_complete_marker & !all$is_prior_predictive, , drop = FALSE], "all_completed_posterior_runs.csv")
    csv(all[all$is_prior_predictive, , drop = FALSE], "all_prior_predictive_runs.csv")
    views <- list(
      model_and_compute = c("equation", "formula", "family", "subset", "observations", "casts", "maximum_depth", "prior",
        "k_depth", "k_mld", "k_distance", "config_k", "chains", "iter", "warmup", "retained_per_chain_configured",
        "adapt_delta", "treedepth", "threads", "cpus", "cpu_model", "fit_seconds", "sampling_seconds", "requested_core_hours",
        "data_md5", "source_md5", "resolved_prior_file_md5"),
      sampler_and_precision = c("max_rhat", "min_bulk_ess", "min_tail_ess", "max_mcse_over_sd", "divergences",
        "treedepth_hits", "treedepth_hit_fraction", "mean_leapfrogs", "min_ebfmi", "convergence_screen_pass",
        "diagnostic_screen_pass", "tree_efficiency_warning", "target_ess_per_sampling_second"),
      prediction_and_profiles = c("validation_set_md5", "training_evaluation_set_md5", "training_rmse", "training_ppc_95_coverage",
        grep("^(validation|training)_(conditional|population)_", names(all), value = TRUE)))
    id <- c("review_uid", "campaign", "case_id", "label", "snapshot_status", "evaluation_type", "is_prior_predictive")
    for (nm in names(views)) csv(all[unique(intersect(c(id, views[[nm]]), names(all)))], paste0(nm, ".csv"))
    count <- aggregate(list(runs = rep(1L, nrow(all))), all[c("campaign", "snapshot_status", "evaluation_type")], sum)
    csv(count, "coverage_counts.csv")
  }
  csv(bind(missing), "missing_or_incomplete_metrics.csv")
  csv(bind(copied), "evidence_index.csv")
  for (nm in names(long)) csv(bind(long[[nm]]), paste0("all_", nm))
  for (m in manifests) {
    dest <- file.path(output, "campaign_metadata", rel(m)); dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    file.copy(m, dest, overwrite = FALSE)
    for (name in c("cases.tsv", "PLAN.md", "README.md", "split_counts.csv", "basis_counts.csv")) {
      f <- file.path(dirname(m), name)
      if (file.exists(f)) file.copy(f, file.path(dirname(dest), name), overwrite = FALSE)
    }
  }
  writeLines(c("Snapshot UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE), "Root:", root,
    "No runs are ranked, selected by error, or excluded for poor diagnostics.",
    "All registered runs are retained; completed posterior and prior runs have separate tables.",
    "Fit/check completion is based on saved markers, not live Slurm state. Incomplete can include timed-out runs.",
    "Only existing results are read. No fits, predictions, or diagnostics are recomputed.",
    "Legacy bounded training RMSE is not comparable to extensive held-out RMSE.",
    "Compare held-out scores only on the same validation set and prediction mode.",
    "A validation hash identifies the saved evaluation file; different hashes require a data identity check.",
    "Per-region/stratum/cast tables retain every available entry, even from runs with poor diagnostics.",
    "Prior predictive checks concern plausibility and scale, not accuracy on the observed response.",
    "Missing columns/values mean unavailable; missing is never converted to zero or pass.",
    "Legacy diagnostic_screen_pass and newer convergence_screen_pass have different definitions.",
    "Small outputs and plots are copied. Fit/data RDS, raw chain draws, row predictions, and >20 MiB files are excluded.",
    "Fits outside experiments/, benchmarks/, or runs/ must be registered or collected separately.",
    "Registered manifests with missing run_dir and case_id are flagged as unresolved; inspect their copied manifests."),
    file.path(output, "READ_ME_FIRST.txt"))
  writeLines(unique(warnings), file.path(output, "collection_warnings.txt"))
  message("Registered runs: ", nrow(all), "; completed fit markers: ", if (nrow(all)) sum(all$fit_complete_marker) else 0L,
    "; warnings: ", length(unique(warnings)), ". Saved: ", output)
  invisible(all)
}
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L) stop("Usage: Rscript collect_all_runs.R PROJECT_ROOT NEW_OUTPUT_DIRECTORY")
  collect_all(args[1], args[2])
}
