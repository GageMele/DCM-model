#!/usr/bin/env Rscript
# Create the model dataframe once. Run from the DCM_MODEL project directory.
# Inputs: data/raw/CTD_NOAA_WHOI_KAUST.csv and data/raw/mld_full.csv.
# Output: a bundle containing model_data, preprocessing choices, and QC tables.
# All settings can be changed using Sys.setenv(...) in R or export in Bash.

prepare_model_data <- function() {
  env_num <- function(name, default) {
    x <- suppressWarnings(as.numeric(Sys.getenv(name, as.character(default))))
    if (length(x) != 1L || !is.finite(x)) stop(name, " must be finite numeric.")
    x
  }
  ctd_file <- Sys.getenv("CTD_INPUT", "data/raw/CTD_NOAA_WHOI_KAUST.csv")
  mld_file <- Sys.getenv("MLD_INPUT", "data/raw/mld_full.csv")
  out_file <- Sys.getenv("MODEL_DATA", "data/processed/model_data.rds")
  fraction <- env_num("CAST_FRACTION", 1)
  seed <- env_num("DATA_SEED", 1234)
  # Retain each eligible input observation at its original depth and chlorophyll.
  # No running smoother, depth bins, or within-bin aggregation are applied.
  if (fraction <= 0 || fraction > 1) stop("CAST_FRACTION must be in (0, 1].")
  if (file.exists(out_file)) stop("Output exists. Choose a new MODEL_DATA path: ", out_file)
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  qc_dir <- paste0(tools::file_path_sans_ext(out_file), "_qc")
  dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
  csv <- function(x, name) write.csv(x, file.path(qc_dir, name), row.names = FALSE)
  needed <- function(d, names, label) {
    absent <- setdiff(names, colnames(d))
    if (length(absent)) stop(label, " lacks: ", paste(absent, collapse = ", "))
  }
  numeric_column <- function(x, name) {
    y <- suppressWarnings(as.numeric(as.character(x)))
    if (any(!is.na(x) & is.na(y))) stop("Non-numeric values in ", name)
    y
  }
  ocean_for <- function(region) {
    ifelse(region %in% c("North RS", "Central RS", "South RS", "KAUST"),
           "RedSea", ifelse(region %in% c("NWHI", "Samoa", "Mariana", "Hawiian", "Equatorial"),
                            "Pacific", "Other"))
  }
  canonical_region <- function(x) {
    x <- trimws(as.character(x))
    # Align the MLD spelling to the existing CTD/model region label.
    x[x == "Hawaiian" & !is.na(x)] <- "Hawiian"
    x
  }
  canonical_location <- function(x) {
    x <- trimws(as.character(x))
    # 'Not specified' represents unknown location metadata.
    missing <- is.na(x) | !nzchar(x) | tolower(x) == "not specified"
    x[missing] <- NA_character_
    x
  }
  # Use several identity fields to avoid mixing unrelated casts with the same ID.
  make_key <- function(d, columns) {
    needed(d, columns, "Key dataframe")
    # LocationCode may be unknown; other cast identifiers must be present.
    required <- setdiff(columns, "LocationCode")
    missing_counts <- vapply(d[required], function(x)
      sum(is.na(x) | !nzchar(trimws(as.character(x)))), integer(1))
    if (any(missing_counts > 0))
      stop("Missing or empty required identifiers: ", paste(
        paste0(names(missing_counts)[missing_counts > 0], " (", missing_counts[missing_counts > 0], " rows)"),
        collapse = ", "))
    if (any(vapply(d[columns], function(x) any(grepl("\034", as.character(x), fixed = TRUE)), logical(1))))
      stop("Reserved separator found in identifiers.")
    values <- lapply(d[columns], as.character)
    if ("LocationCode" %in% columns) {
      token <- "<missing:LocationCode>"
      if (any(values$LocationCode == token, na.rm = TRUE)) stop("Reserved missing-location token in input.")
      values$LocationCode[is.na(values$LocationCode)] <- token
    }
    do.call(paste, c(values, sep = "\034"))
  }
  raw <- read.csv(ctd_file, stringsAsFactors = FALSE, check.names = FALSE)
  # Prefer the notebook's original columns; accept standardized equivalents.
  aliases <- list(cast_id = c("DeepCTDID", "cast_id"), Region = c("Region.x", "Region"),
                  date_raw = c("UTCDate.x", "date_raw", "date"),
                  depth = c("Depth_m", "depth"), temperature = c("Temp_DegC", "temperature"))
  for (target in names(aliases)) {
    hit <- intersect(aliases[[target]], colnames(raw))
    if (!length(hit)) stop("No input column for ", target)
    raw[[target]] <- raw[[hit[1]]]
  }
  needed(raw, c("LocationCode", "cast_id", "Region", "date_raw", "depth", "temperature", "Chl", "distance"), "CTD")
  raw$Region <- canonical_region(raw$Region)
  raw$LocationCode <- canonical_location(raw$LocationCode)
  raw$Ocean <- ocean_for(raw$Region)
  # The notebook excludes KAUST, but retains Equatorial despite an older comment.
  raw <- raw[!is.na(raw$Region) & raw$Region != "KAUST", , drop = FALSE]
  for (v in c("depth", "temperature", "Chl", "distance")) raw[[v]] <- numeric_column(raw[[v]], v)
  date_chars <- substr(as.character(raw$date_raw), 1, 10)
  if (any(!grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", date_chars)))
    stop("CTD dates must be ISO yyyy-mm-dd. Convert ambiguous dates explicitly first.")
  raw$date_raw <- as.Date(date_chars)
  if (anyNA(raw$date_raw)) stop("Invalid CTD dates.")
  join_cols <- c("LocationCode", "Region", "cast_id", "Ocean")
  raw$cast_key <- make_key(raw, c(join_cols, "date_raw"))

  # Validate the right table before joining. Exact duplicate rows are collapsed;
  # conflicting MLDs for a key stop the workflow instead of multiplying CTD rows.
  mld <- read.csv(mld_file, stringsAsFactors = FALSE, check.names = FALSE)
  needed(mld, c(join_cols, "MLD_temp_MIX"), "MLD")
  # Standardize aliases before joining on the same four fields.
  mld$Region <- canonical_region(mld$Region)
  mld$Ocean <- trimws(as.character(mld$Ocean))
  mld$Ocean[mld$Ocean == "Red Sea" & !is.na(mld$Ocean)] <- "RedSea"
  mld$LocationCode <- canonical_location(mld$LocationCode)
  csv(data.frame(source = c("CTD", "MLD"), rows = c(nrow(raw), nrow(mld)),
    missing_location_rows = c(sum(is.na(raw$LocationCode)), sum(is.na(mld$LocationCode)))),
    "identifier_qc.csv")
  mld <- unique(mld[c(join_cols, "MLD_temp_MIX")])
  mld$MLD_temp_MIX <- numeric_column(mld$MLD_temp_MIX, "MLD_temp_MIX")
  mld_key <- make_key(mld, join_cols)
  conflict <- duplicated(mld_key) | duplicated(mld_key, fromLast = TRUE)
  if (any(conflict)) {
    csv(mld[conflict, ], "conflicting_mld_keys.csv")
    stop("Multiple MLD values per join key. See conflicting_mld_keys.csv.")
  }
  # A four-column join is only safe if the identifiers are unique across dates.
  identities <- unique(raw[c(join_cols, "date_raw")])
  identity_key <- make_key(identities, join_cols)
  if (anyDuplicated(identity_key))
    stop("A four-column MLD key spans multiple dates. Add date to BOTH input tables and join keys.")
  raw_join_key <- make_key(raw, join_cols)
  m <- match(raw_join_key, mld_key)
  raw$MLD_temp_MIX <- mld$MLD_temp_MIX[m]
  csv(unique(raw[is.na(m), c(join_cols, "date_raw")]), "unmatched_ctd_casts.csv")
  csv(mld[!mld_key %in% raw_join_key, ], "unused_mld_casts.csv")

  # Estimate each baseline BEFORE restricting profiles to the model depth range.
  groups <- split(seq_len(nrow(raw)), raw$cast_key)
  baseline <- vapply(groups, function(i) {
    ok <- raw$depth[i] >= 5 & raw$depth[i] <= 15 & is.finite(raw$depth[i]) & is.finite(raw$Chl[i])
    if (!any(ok)) NA_real_ else median(raw$Chl[i][ok])
  }, numeric(1))
  raw$surface_chl <- unname(baseline[raw$cast_key])
  raw$chl_norm <- (raw$Chl - raw$surface_chl) / raw$surface_chl
  raw$log10_distance <- suppressWarnings(log10(raw$distance))
  raw$year <- as.integer(format(raw$date_raw, "%Y"))
  raw$enso <- ifelse(raw$year == 2009, "el_nino", ifelse(raw$year %in% c(2010, 2011), "la_nina", "neutral"))
  # Preserve the source's year-based ENSO labels; these are not an updated ENSO index.
  raw$exclusion_reason <- ifelse(!is.finite(raw$surface_chl) | raw$surface_chl <= 0, "invalid_surface_baseline",
    ifelse(!is.finite(raw$depth) | raw$depth < 10 | raw$depth > 250, "outside_depth_range",
    ifelse(!is.finite(raw$chl_norm), "nonfinite_response",
    ifelse(!is.finite(raw$log10_distance), "invalid_distance",
    ifelse(!is.finite(raw$MLD_temp_MIX) | raw$MLD_temp_MIX < 0, "missing_or_invalid_mld", "retained")))))
  csv(as.data.frame(table(raw$Region, raw$exclusion_reason)), "row_exclusions_by_region.csv")
  csv(unique(raw[raw$exclusion_reason != "retained", c("cast_key", join_cols, "date_raw", "exclusion_reason")]),
      "excluded_cast_reasons.csv")
  d <- raw[raw$exclusion_reason == "retained", , drop = FALSE]
  if (!nrow(d)) stop("No eligible model observations.")
  d <- d[order(d$cast_key, d$depth), , drop = FALSE]
  # Draw complete casts ONCE, independently within regions. Full data is default.
  set.seed(seed)
  cast_rows <- unique(d[c("cast_key", "Region")])
  chosen <- unlist(lapply(split(cast_rows$cast_key, cast_rows$Region), function(ids) {
    if (fraction == 1) ids else sample(ids, max(1L, floor(length(ids) * fraction)))
  }), use.names = FALSE)
  d <- d[d$cast_key %in% chosen, , drop = FALSE]
  d$Region <- factor(d$Region)  # Same alphabetical reference-level convention.
  d$Ocean <- factor(d$Ocean)
  d$enso <- factor(d$enso)
  # 'date' is a grouping factor, not a continuous date effect. Scaling it in the
  # notebook changed labels, not group membership; retain readable ISO labels.
  d$date <- factor(as.character(d$date_raw))
  d$OceanRegion <- interaction(d$Ocean, d$Region, sep = "_", drop = TRUE)
  columns <- c("chl_norm", "Region", "depth", "date", "log10_distance", "enso", "Ocean", "OceanRegion",
               "MLD_temp_MIX", "cast_id", "cast_key", "LocationCode", "date_raw", "Chl", "surface_chl", "distance")
  model_data <- d[columns]
  rownames(model_data) <- NULL
  cast_summary <- do.call(rbind, lapply(split(model_data, model_data$Region), function(x) data.frame(
    Region = as.character(x$Region[1]), rows = nrow(x), casts = length(unique(x$cast_key)),
    depth_min = min(x$depth), depth_max = max(x$depth),
    distance_min = min(x$distance), distance_max = max(x$distance),
    mld_min = min(x$MLD_temp_MIX), mld_max = max(x$MLD_temp_MIX),
    distinct_depths = length(unique(x$depth)), distinct_distances = length(unique(x$log10_distance)),
    distinct_mlds = length(unique(x$MLD_temp_MIX)))))
  csv(cast_summary, "region_summary.csv")
  csv(unique(model_data[c("cast_key", "LocationCode", "cast_id", "Region", "date_raw", "MLD_temp_MIX")]), "selected_casts.csv")
  bundle <- list(model_data = model_data, region_summary = cast_summary,
    config = list(ctd_input = normalizePath(ctd_file), mld_input = normalizePath(mld_file),
      input_md5 = tools::md5sum(c(ctd_file, mld_file)), cast_fraction = fraction, seed = seed,
      smoothing = "none", binned = FALSE, baseline_depth_m = c(5, 15), model_depth_m = c(10, 250),
      identifier_aliases = c("Hawaiian -> Hawiian", "Red Sea -> RedSea", "Not specified/empty LocationCode -> NA"),
      response = "(Chl - surface_chl) / surface_chl", distance_units = "same as input distance; notebook treats these as km"),
    created_utc = format(Sys.time(), tz = "UTC"), session = capture.output(sessionInfo()))
  saveRDS(bundle, out_file, compress = "gzip")
  print(cast_summary)
  message("Saved model data: ", out_file, "\nQC directory: ", qc_dir)
  invisible(bundle)
}

# Command-line execution; in RStudio use source('01_prepare_data.R'); prepare_model_data().
if (sys.nframe() == 0L) prepare_model_data()
