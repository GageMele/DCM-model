# Validate the exact equations, pairing and actual-input brms construction.
package <- "DCM_Ibex_equation_tests"
for (path in list.files(package, "\\.R$", full.names = TRUE)) parse(path)
source(file.path(package, "common.R"))
source(file.path(package, "equation.R"))
cases <- read.delim(file.path(package, "cases.tsv"), stringsAsFactors = FALSE)
stopifnot(nrow(cases) == 18L, !anyDuplicated(cases$case_id),
  sum(cases$cpus) == 168L, all(cases$k_depth == 15L))
for (i in 1:6) {
  matched <- cases[c(i, i + 6L, i + 12L), ]
  fixed <- setdiff(names(matched), c("case_id", "label", "equation"))
  stopifnot(all(vapply(matched[fixed], function(x) length(unique(x)) == 1L, logical(1))))
}
original <- "DCM_Ibex_workflow/02_fit_model.R"
expected <- list(
  additive = chl_norm ~ 1 + Region + log10_distance + MLD_temp_MIX +
    s(depth, by = Region, bs = "cr", k = 15) +
    (1 | Ocean:Region:date) + (1 | Ocean:Region:enso),
  mld_varying = chl_norm ~ 1 + Region + log10_distance +
    s(depth, by = MLD_temp_MIX, bs = "cr", k = 15) +
    s(depth, by = Region, bs = "cr", k = 15) +
    (1 | Ocean:Region:date) + (1 | Ocean:Region:enso),
  mld_only = chl_norm ~ 1 + Region + log10_distance +
    s(depth, by = MLD_temp_MIX, bs = "cr", k = 15) +
    (1 | Ocean:Region:date) + (1 | Ocean:Region:enso))
for (id in names(expected)) {
  spec <- load_equation_spec(id, original)
  stopifnot(identical(spec$formula$formula[[2L]], expected[[id]][[2L]]),
    identical(spec$formula$formula[[3L]], expected[[id]][[3L]]),
    spec$formula$family$family == "gaussian",
    identical(spec$current_priors, load_original_spec(original)$current_priors))
}
experiment <- file.path(getwd(), "validation_equation150", "actual")
dir.create(file.path(experiment, "data"), recursive = TRUE, showWarnings = FALSE)
file.copy(list.files(package, "\\.R$", full.names = TRUE), experiment, overwrite = TRUE)
file.copy(file.path(package, "cases.tsv"), experiment, overwrite = TRUE)
file.copy(original, experiment, overwrite = TRUE)
input <- "validation_fast_tests/actual150_radix/data/p10.rds"
stopifnot(file.copy(input, file.path(experiment, "data", "p10.rds"), overwrite = TRUE),
  unname(tools::md5sum(input)) == unname(tools::md5sum(file.path(experiment, "data", "p10.rds"))))
Sys.setenv(EXPERIMENT_DIR = experiment)
source(file.path(experiment, "01_validate_campaign.R"))
report <- read.csv(file.path(experiment, "equation_validation.csv"))
stopifnot(nrow(report) == 3L, all(report$observations == 15228L), all(report$casts == 108L))
print(report)
cat("All three exact formulas, unchanged priors, actual data and generated Stan code passed.\n")

# Exercise empty-report handling and complete/partial/failed paired cases.
source(file.path(package, "04_collect_results.R"))
collector <- file.path(getwd(), "validation_equation150", "collector")
dir.create(file.path(collector, "runs"), recursive = TRUE, showWarnings = FALSE)
write.table(cases, file.path(collector, "cases.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
manifest <- data.frame(case_id = cases$case_id, label = cases$label,
  fit_job = 1000L + cases$case_id, check_job = 2000L + cases$case_id,
  cpus = cases$cpus, run_dir = file.path(collector, "runs", sprintf("case_%02d", cases$case_id)))
write.table(manifest, file.path(collector, "manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
stopifnot(identical(collect_equation_results(collector), FALSE))
for (i in c(1L, 7L, 13L)) {
  run <- manifest$run_dir[i]
  dir.create(run, showWarnings = FALSE)
  row <- data.frame(cases[i, ], data_md5 = "frozen-input", source_md5 = "equation-spec",
    cpu_model = "mock CPU", fit_seconds = 1000 * i, sampling_seconds = 400 * i,
    max_rhat = 1, min_bulk_ess = 500, min_tail_ess = 500,
    divergences = 0, treedepth_hits = 0, diagnostic_screen_pass = TRUE)
  write.csv(row, file.path(run, if (i == 13) "diagnostic_summary.csv" else "test_summary.csv"), row.names = FALSE)
  if (i == 13) writeLines("mock plotting error", file.path(run, "CHECKS_FAILED.txt"))
}
stopifnot(identical(collect_equation_results(collector), TRUE))
pairs <- read.csv(file.path(collector, "matched_equation_pairs.csv"))
stopifnot(nrow(pairs) == 3L, all(pairs$matched_input_and_sampler))
cat("Collector handles zero completed fits, diagnostic-only rows and all three equation pairs.\n")
