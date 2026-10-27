# Validate actual-data subset identity, exact formula/priors, and reporting edges.
package <- "DCM_Ibex_baseline_tests"
for (path in list.files(package, "\\.R$", full.names = TRUE)) parse(path)
source(file.path(package, "common.R"))
source(file.path(package, "baseline.R"))
cases <- read.delim(file.path(package, "cases.tsv"), stringsAsFactors = FALSE)
validate_baseline_cases(cases)
expected <- chl_norm ~ 1 + Region + log10_distance +
  s(depth, by = Region, bs = "cr", k = 15) +
  (1 | Ocean:Region:date) + (1 | Ocean:Region:enso)
for (prior in c("custom", "default")) {
  spec <- baseline_spec(prior)
  stopifnot(identical(spec$formula$formula[[2]], expected[[2]]),
    identical(spec$formula$formula[[3]], expected[[3]]),
    spec$formula$family$family == "gaussian")
}
expected_priors <- c(
  brms::set_prior("student_t(3, 0.73, 5.448)", class = "Intercept"),
  brms::set_prior("normal(0, 2)", class = "b"),
  brms::set_prior("student_t(3, 0, 1.5)", class = "sd"),
  brms::set_prior("student_t(3, 0, 1.5)", class = "sds"),
  brms::set_prior("student_t(3, 0, 2.5)", class = "sigma"))
stopifnot(identical(baseline_spec("custom")$priors, expected_priors),
  is.null(baseline_spec("default")$priors), 2 * 2.724 == 5.448)

experiment <- file.path(getwd(), "validation_baseline150", "actual")
dir.create(file.path(experiment, "data"), recursive = TRUE, showWarnings = FALSE)
already_prepared <- file.exists(file.path(experiment, "PREPARATION_COMPLETE.txt"))
file.copy(list.files(package, "\\.R$", full.names = TRUE), experiment, overwrite = TRUE)
file.copy(file.path(package, "cases.tsv"), experiment, overwrite = TRUE)
stopifnot(file.copy("validation_actual_inputs/full.rds",
  file.path(experiment, "full_input_250m.rds"), overwrite = TRUE))
hashes <- character()
for (subset in c("p10", "p20")) {
  input <- file.path("validation_fast_tests/actual150_radix/data", paste0(subset, ".rds"))
  output <- file.path(experiment, "data", paste0(subset, ".rds"))
  stopifnot(file.copy(input, output, overwrite = TRUE))
  hashes[subset] <- unname(tools::md5sum(output))
}
source(file.path(experiment, "01_prepare_baselines.R"))
if (!already_prepared) prepare_baselines(experiment)
for (subset in names(hashes)) stopifnot(hashes[subset] ==
  unname(tools::md5sum(file.path(experiment, "data", paste0(subset, ".rds")))))
full <- readRDS(file.path(experiment, "full_input_250m.rds"))$model_data
full <- full[full$depth >= 10 & full$depth <= 150, ]
previous <- character()
for (label in c("p10", "p20", "p50")) {
  d <- readRDS(file.path(experiment, "data", paste0(label, ".rds")))$model_data
  ids <- unique(d$cast_key)
  expected <- full[full$cast_key %in% ids, ]
  stopifnot(all(previous %in% ids), nrow(d) == nrow(expected),
    length(ids) == c(p10 = 108L, p20 = 216L, p50 = 540L)[label])
  a <- order(d$cast_key, d$depth, d$chl_norm, method = "radix")
  b <- order(expected$cast_key, expected$depth, expected$chl_norm, method = "radix")
  for (column in c("cast_key", "Region", "Ocean", "date", "enso"))
    stopifnot(identical(as.character(d[[column]][a]), as.character(expected[[column]][b])))
  for (column in c("depth", "chl_norm", "log10_distance", "MLD_temp_MIX", "Chl", "surface_chl"))
    stopifnot(identical(d[[column]][a], expected[[column]][b]))
  previous <- ids
  cat(label, ":", length(ids), "whole casts and", nrow(d), "unaltered observations.\n")
}
report <- read.csv(file.path(experiment, "baseline_validation.csv"))
stopifnot(nrow(report) == 6L, all(report$stan_penalized_spline_columns == 104L),
  all(report$stan_fixed_effect_columns == 9L))
print(report)
for (subset in c("p10", "p20", "p50")) {
  for (prior in c("custom", "default")) {
    resolved <- read.csv(file.path(experiment, paste0("validated_", subset, "_", prior, "_priors.csv")),
      stringsAsFactors = FALSE)
    if (prior == "custom") stopifnot(resolved$prior[resolved$class == "b" & resolved$coef == ""] == "normal(0, 2)")
    if (prior == "default") stopifnot(resolved$prior[resolved$class == "b" & resolved$coef == ""] == "")
    code <- readLines(file.path(experiment, paste0("validated_", subset, "_", prior, ".stan")))
    stopifnot(any(grepl("reduce_sum", code, fixed = TRUE)))
  }
}

# Verify that actual fitting arguments preserve the requested historical settings.
parsed <- parse(file.path(package, "02_fit_case.R"))
find_calls <- function(x, name) {
  if (missing(x)) return(list())
  result <- list()
  if (is.call(x) && identical(x[[1]], parse(text = name)[[1]])) result <- list(x)
  if (is.call(x) || is.expression(x) || is.pairlist(x))
    for (element in as.list(x)) result <- c(result, find_calls(element, name))
  result
}
call <- find_calls(parsed, "brms::brm")
stopifnot(length(call) == 1L)
args <- as.list(call[[1]])[-1]
for (field in c("chains", "cores", "threads", "iter", "warmup", "refresh")) {
  expected_arg <- if (field == "threads") quote(thread_spec) else parse(text = paste0("case$", field))[[1]]
  stopifnot(identical(args[[field]], expected_arg))
}
stopifnot(identical(args$prior, quote(priors)), identical(args$save_warmup, TRUE),
  identical(args$backend, "cmdstanr"), identical(args$thin, 1L))

# Exercise no completions, partial diagnostics, prior pairs and sampler pairs.
source(file.path(package, "04_collect_results.R"))
collector <- file.path(getwd(), "validation_baseline150", "collector")
dir.create(file.path(collector, "runs"), recursive = TRUE, showWarnings = FALSE)
write.table(cases, file.path(collector, "cases.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
manifest <- data.frame(case_id = cases$case_id, label = cases$label,
  fit_job = 1000L + cases$case_id, check_job = 2000L + cases$case_id,
  cpus = cases$cpus, run_dir = file.path(collector, "runs", sprintf("case_%02d", cases$case_id)))
write.table(manifest, file.path(collector, "manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
stopifnot(identical(collect_baseline_results(collector), FALSE))
for (i in 1:12) {
  run <- manifest$run_dir[i]
  dir.create(run, showWarnings = FALSE)
  row <- data.frame(cases[i, ], data_md5 = cases$subset[i], source_md5 = "baseline-spec",
    cpu_model = "mock CPU", fit_seconds = 1000 + i, sampling_seconds = 400 + i,
    max_rhat = 1, min_bulk_ess = 500, min_tail_ess = 500,
    divergences = 0, treedepth_hits = 0, treedepth_hit_fraction = 0, diagnostic_screen_pass = TRUE)
  write.csv(row, file.path(run, if (i == 12) "diagnostic_summary.csv" else "test_summary.csv"), row.names = FALSE)
  if (i == 12) writeLines("mock plotting error", file.path(run, "CHECKS_FAILED.txt"))
}
stopifnot(identical(collect_baseline_results(collector), TRUE))
pairs <- read.csv(file.path(collector, "matched_baseline_pairs.csv"))
status <- read.csv(file.path(collector, "case_status.csv"))
stopifnot(nrow(pairs) == 12L, all(pairs$matched_other_settings_and_data),
  status$status[status$case_id == 12] == "checks_failed")
cat("Exact formula/priors, all six Stan constructions, data identity, fit arguments and reporting tests passed.\n")
