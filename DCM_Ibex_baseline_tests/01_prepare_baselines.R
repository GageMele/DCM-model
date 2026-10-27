#!/usr/bin/env Rscript
# Reuse the frozen 10% and 20% files; extend the SAME cast selection to 50%.
prepare_baselines <- function(experiment) {
  source(file.path(experiment, "common.R"), local = TRUE)
  source(file.path(experiment, "baseline.R"), local = TRUE)
  if (file.exists(file.path(experiment, "PREPARATION_COMPLETE.txt")))
    stop("This experiment was already prepared.")
  for (p in c("brms", "cmdstanr", "mgcv", "posterior"))
    if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)
  cases <- read.delim(file.path(experiment, "cases.tsv"), stringsAsFactors = FALSE)
  validate_baseline_cases(cases)
  full_path <- file.path(experiment, "full_input_250m.rds")
  bundle <- readRDS(full_path)
  original <- validate_bundle(bundle, full = TRUE)
  d <- original[original$depth >= 10 & original$depth <= 150, , drop = FALSE]
  if (!nrow(d)) stop("No rows remain at 10-150 m.")
  check_baseline_support(d)
  metadata <- unique(d[c("cast_key", "Region", "log10_distance", "MLD_temp_MIX")])
  if (anyDuplicated(metadata$cast_key)) stop("Predictors vary within a cast.")
  metadata <- metadata[order(metadata$cast_key, method = "radix"), , drop = FALSE]
  regions <- sort(unique(as.character(d$Region)), method = "radix")
  # Reproduce the first campaign's supported, region-stratified cast order.
  orders <- lapply(seq_along(regions), function(i) supporting_cast_order(
    metadata[as.character(metadata$Region) == regions[i], ], 4100L + i))
  names(orders) <- regions
  sizes <- vapply(orders, function(x) length(x$ids), integer(1))
  counts <- vapply(orders, function(x) x$minimum, integer(1))
  reports <- selections <- list()
  previous_ids <- character()
  for (fraction in c(.05, .10, .20, .50)) {
    counts <- allocate_cast_counts(sizes, counts, fraction)
    ids <- unlist(lapply(seq_along(regions), function(i)
      head(orders[[i]]$ids, counts[i])), use.names = FALSE)
    stopifnot(all(previous_ids %in% ids))
    previous_ids <- ids
    if (fraction == .05) next  # Needed to reproduce the original nested selection.
    subset <- d[d$cast_key %in% ids, , drop = FALSE]
    set.seed(5000L + as.integer(100 * fraction))
    subset <- subset[sample.int(nrow(subset)), , drop = FALSE]
    subset$Region <- factor(as.character(subset$Region), levels = regions)
    subset$Ocean <- factor(as.character(subset$Ocean), levels = levels(d$Ocean))
    subset$enso <- factor(as.character(subset$enso), levels = levels(d$enso))
    subset$date <- factor(as.character(subset$date))
    rownames(subset) <- NULL
    label <- sprintf("p%02d", as.integer(fraction * 100))
    path <- file.path(experiment, "data", paste0(label, ".rds"))
    if (fraction < .5) {
      # Exact row/value/order identity also verifies that the copied parent is right.
      reused <- readRDS(path)
      stopifnot(identical(validate_bundle(reused), subset),
        isTRUE(all.equal(reused$config$requested_cast_fraction, fraction)),
        identical(reused$config$parent_bundle_md5, unname(tools::md5sum(full_path))))
    }
    report <- check_baseline_support(subset)
    report$subset <- label
    report$requested_fraction <- fraction
    report$full_region_casts <- sizes[match(report$Region, regions)]
    report$region_cast_fraction <- report$casts / report$full_region_casts
    report$global_cast_fraction <- length(ids) / nrow(metadata)
    report$global_observation_fraction <- nrow(subset) / nrow(d)
    reports[[label]] <- report
    selections[[label]] <- data.frame(subset = label,
      metadata[metadata$cast_key %in% ids, c("cast_key", "Region")])
    if (fraction == .5) {
      out <- bundle
      out$model_data <- subset
      out$region_summary <- report
      out$config$cast_fraction <- length(ids) / nrow(metadata)
      out$config$model_depth_m <- c(10, 150)
      out$config$exploratory_max_depth_m <- 150
      out$config$requested_cast_fraction <- fraction
      out$config$subset_seed <- 4100L
      out$config$row_order <- "fixed random order, observations unchanged"
      out$config$parent_bundle_md5 <- unname(tools::md5sum(full_path))
      out$config$parent_bundle <- normalizePath(full_path, winslash = "/")
      out$created_utc <- format(Sys.time(), tz = "UTC")
      saveRDS(out, path)
    }
    cat(label, ":", length(ids), "casts;", nrow(subset), "observations\n")
  }
  write.csv(do.call(rbind, reports), file.path(experiment, "subset_report.csv"), row.names = FALSE)
  write.csv(do.call(rbind, selections), file.path(experiment, "selected_casts.csv"), row.names = FALSE)
  write.csv(data.frame(observations = nrow(d), casts = nrow(metadata),
    minimum_depth = min(d$depth), maximum_depth = max(d$depth)),
    file.path(experiment, "full_150m_counts.csv"), row.names = FALSE)

  # Resolve custom AND actual brms-default priors separately on each dataset.
  validation <- list()
  for (subset in c("p10", "p20", "p50")) {
    data <- readRDS(file.path(experiment, "data", paste0(subset, ".rds")))$model_data
    for (prior in c("custom", "default")) {
      spec <- baseline_spec(prior)
      standata <- brms::make_standata(spec$formula, data = data, prior = spec$priors)
      label <- paste(subset, prior, sep = "_")
      code <- brms::make_stancode(spec$formula, data = data, prior = spec$priors,
        threads = brms::threading(2))
      writeLines(code, file.path(experiment, paste0("validated_", label, ".stan")))
      # validate_prior merges supplied overrides with version/data-specific defaults.
      resolved <- brms::validate_prior(spec$priors, spec$formula, data = data)
      write.csv(as.data.frame(resolved), file.path(experiment,
        paste0("validated_", label, "_priors.csv")), row.names = FALSE)
      capture.output(spec$formula, resolved, file = file.path(experiment,
        paste0("validated_", label, "_formula_priors.txt")))
      matrices <- standata[grepl("^Zs", names(standata))]
      validation[[label]] <- data.frame(subset = subset, prior = prior,
        observations = nrow(data), casts = length(unique(data$cast_key)),
        stan_fixed_effect_columns = standata$K,
        stan_penalized_spline_columns = sum(vapply(matrices, function(x)
          if (is.matrix(x)) ncol(x) else 0L, integer(1))))
      rm(standata); invisible(gc())
    }
  }
  write.csv(do.call(rbind, validation), file.path(experiment, "baseline_validation.csv"), row.names = FALSE)
  paths <- file.path(experiment, "data", paste0(c("p10", "p20", "p50"), ".rds"))
  write.csv(data.frame(subset = c("p10", "p20", "p50"), md5 = unname(tools::md5sum(paths))),
    file.path(experiment, "data_checksums.csv"), row.names = FALSE)
  capture.output(sessionInfo(), file = file.path(experiment, "validation_session_info.txt"))
  writeLines(c("All baseline cases, nested whole-cast subsets, formula and priors validated.",
    "10% and 20% files are unchanged copies from the original fast150 campaign.",
    "50% extends that same cast ordering; chlorophyll normalization is unchanged."),
    file.path(experiment, "PREPARATION_COMPLETE.txt"))
  invisible(TRUE)
}
if (sys.nframe() == 0L) prepare_baselines(Sys.getenv("EXPERIMENT_DIR"))
