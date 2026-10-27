#!/usr/bin/env Rscript
# Create nested whole-cast subsets from the FULL prepared dataframe.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
experiment <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
source(file.path(experiment, "common.R"))
full_path <- Sys.getenv("FULL_MODEL_DATA", "data/processed/full_unsmoothed_v2.rds")
if (file.exists(file.path(experiment, "PREPARATION_COMPLETE.txt")))
  stop("This experiment was already prepared.")
bundle <- readRDS(full_path)
d_original <- validate_bundle(bundle, full = TRUE)
# The existing 250 m jobs are untouched. New experiments stop at 150 m.
d <- d_original[d_original$depth <= 150, , drop = FALSE]
if (!nrow(d)) stop("No observations remain within 10-150 m.")
write.csv(data.frame(stage = c("original_prepared", "trimmed_150m"),
  observations = c(nrow(d_original), nrow(d)),
  casts = c(length(unique(d_original$cast_key)), length(unique(d$cast_key))),
  minimum_depth = c(min(d_original$depth), min(d$depth)),
  maximum_depth = c(max(d_original$depth), max(d$depth))),
  file.path(experiment, "depth_trim_report.csv"), row.names = FALSE)
write.csv(full_response_statistics(d_original),
  file.path(experiment, "original_250m_response_statistics.csv"), row.names = FALSE)
full_support <- check_support(d)
write.csv(full_support, file.path(experiment, "full_region_support.csv"), row.names = FALSE)
statistics <- full_response_statistics(d)
saveRDS(statistics, file.path(experiment, "full_response_statistics.rds"))
write.csv(statistics, file.path(experiment, "full_response_statistics.csv"), row.names = FALSE)

# Check that predictors used for cast selection are truly cast-level metadata.
metadata <- unique(d[c("cast_key", "Region", "log10_distance", "MLD_temp_MIX")])
if (anyDuplicated(metadata$cast_key))
  stop("Distance or MLD varies within a cast; review metadata before selecting casts.")
# Use bytewise ordering to make the seeded selection independent of OS locale.
metadata <- metadata[order(metadata$cast_key, method = "radix"), , drop = FALSE]
regions <- sort(unique(as.character(d$Region)), method = "radix")
orders <- lapply(seq_along(regions), function(i) {
  supporting_cast_order(metadata[as.character(metadata$Region) == regions[i], ], 4100L + i)
})
names(orders) <- regions
sizes <- vapply(orders, function(x) length(x$ids), integer(1))
counts <- vapply(orders, function(x) x$minimum, integer(1))
reports <- selections <- list()
dir.create(file.path(experiment, "data"), showWarnings = FALSE)
for (fraction in c(.05, .10, .20)) {
  counts <- allocate_cast_counts(sizes, counts, fraction)
  ids <- unlist(lapply(seq_along(regions), function(i) head(orders[[i]]$ids, counts[i])),
                use.names = FALSE)
  subset <- d[d$cast_key %in% ids, , drop = FALSE]
  # Shuffle rows once for threading balance; retain every observation unchanged.
  set.seed(5000L + as.integer(100 * fraction))
  subset <- subset[sample.int(nrow(subset)), , drop = FALSE]
  subset$Region <- factor(as.character(subset$Region), levels = regions)
  subset$Ocean <- factor(as.character(subset$Ocean), levels = levels(d$Ocean))
  subset$enso <- factor(as.character(subset$enso), levels = levels(d$enso))
  subset$date <- factor(as.character(subset$date))
  rownames(subset) <- NULL
  label <- sprintf("p%02d", as.integer(fraction * 100))
  report <- check_support(subset)
  report$subset <- label
  report$requested_fraction <- fraction
  report$full_region_casts <- sizes[match(report$Region, regions)]
  report$region_cast_fraction <- report$casts / report$full_region_casts
  report$global_cast_fraction <- length(ids) / nrow(metadata)
  report$global_observation_fraction <- nrow(subset) / nrow(d)
  reports[[label]] <- report
  selections[[label]] <- data.frame(subset = label,
    metadata[metadata$cast_key %in% ids, c("cast_key", "Region")])
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
  saveRDS(out, file.path(experiment, "data", paste0(label, ".rds")))
  cat(label, ":", length(ids), "casts;", nrow(subset), "observations;",
      round(100 * length(ids) / nrow(metadata), 2), "% of full casts\n")
}
write.csv(do.call(rbind, reports), file.path(experiment, "subset_report.csv"), row.names = FALSE)
write.csv(do.call(rbind, selections), file.path(experiment, "selected_casts.csv"), row.names = FALSE)
writeLines(c(paste("Full bundle:", normalizePath(full_path, winslash = "/")),
  paste("Full MD5:", unname(tools::md5sum(full_path))),
  "All subsets retain whole casts, with original chlorophyll and depths.",
  "All new test datasets are restricted to 10-150 m; surface baselines are unchanged.",
  "No smoothing or binning. Support constraints may increase the smallest subset.",
  "Full-data-derived priors are empirical Bayes exploratory sensitivity choices.",
  "They are not independent prior information or an independent validation."),
  file.path(experiment, "preparation_notes.txt"))
writeLines("Complete", file.path(experiment, "PREPARATION_COMPLETE.txt"))
