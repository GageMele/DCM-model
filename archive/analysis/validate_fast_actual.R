# Confirm the actual 150 m subsets contain complete, unaltered eligible casts.
root <- commandArgs(trailingOnly = TRUE)[1]
full <- readRDS("validation_actual_inputs/full.rds")$model_data
full <- full[full$depth <= 150, ]
previous <- character()
for (label in c("p05", "p10", "p20")) {
  bundle <- readRDS(file.path(root, "data", paste0(label, ".rds")))
  d <- bundle$model_data
  ids <- unique(d$cast_key)
  expected <- full[full$cast_key %in% ids, ]
  stopifnot(nrow(d) == nrow(expected), all(previous %in% ids),
    all(d$depth >= 10 & d$depth <= 150),
    identical(bundle$config$smoothing, "none"), identical(bundle$config$binned, FALSE))
  a <- order(d$cast_key, d$depth, d$chl_norm, method = "radix")
  b <- order(expected$cast_key, expected$depth, expected$chl_norm, method = "radix")
  for (column in c("cast_key", "Region", "Ocean", "date", "enso"))
    stopifnot(identical(as.character(d[[column]][a]), as.character(expected[[column]][b])))
  for (column in c("depth", "chl_norm", "log10_distance", "MLD_temp_MIX", "Chl", "surface_chl"))
    stopifnot(identical(d[[column]][a], expected[[column]][b]))
  previous <- ids
  cat(label, ": all", nrow(d), "eligible observations preserved across", length(ids), "whole casts\n")
}
cat("Actual subsets are nested, bounded at 150 m, and contain unchanged observations.\n")
