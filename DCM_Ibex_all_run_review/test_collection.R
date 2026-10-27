#!/usr/bin/env Rscript
# Synthetic coverage/retention tests; does not touch project runs.
args <- commandArgs(trailingOnly = TRUE)
script_dir <- if (length(args)) args[1] else "DCM_Ibex_all_run_review"
source(file.path(script_dir, "collect_all_runs.R"))
base <- tempfile("dcm_all_review_test_"); dir.create(base)
root <- file.path(base, "project"); dir.create(root)
campaign <- file.path(root, "experiments", "extensive150_test")
dir.create(file.path(campaign, "data"), recursive = TRUE)
saveRDS(list(casts = c("heldout_a", "heldout_b")), file.path(campaign, "data", "validation_model.rds"))
saveRDS(list(casts = "train_a"), file.path(campaign, "data", "p10_model.rds"))
runs <- file.path(campaign, "runs", sprintf("case_%02d", 1:4))
for (r in runs) dir.create(file.path(r, "checks"), recursive = TRUE)
manifest <- data.frame(case_id = 1:4, label = c("good", "bad", "prior", "unfinished"),
  fit_job = 11:14, check_job = 21:24, cpus = 8,
  run_dir = paste0("experiments/extensive150_test/runs/case_0", 1:4))
write.table(manifest, file.path(campaign, "manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
for (i in 1:4) {
  cfg <- list(case = data.frame(case_id = i, mode = if (i == 3) "prior" else "posterior", equation = "test"),
    mode = "extensive", chains = 4, iter = 2000, warmup = 1000,
    max_treedepth = 14, formula = "chl_norm ~ Region + s(depth)")
  saveRDS(cfg, file.path(runs[i], "run_config.rds"))
  if (i < 4) writeLines("done", file.path(runs[i], "FIT_COMPLETE.txt"))
  if (i %in% c(1, 3)) writeLines("done", file.path(runs[i], "CHECKS_COMPLETE.txt"))
}
writeLines("plot failed", file.path(runs[2], "CHECKS_FAILED.txt"))
write.csv(data.frame(min_tail_ess = 12, min_ebfmi = 0.1),
  file.path(runs[2], "diagnostic_summary.csv"), row.names = FALSE)
write.csv(data.frame(min_bulk_ess = 21, divergences = 7),
  file.path(runs[4], "diagnostic_summary.csv"), row.names = FALSE)
for (i in 1:2) {
  write.csv(data.frame(validation_conditional_rmse = if (i == 1) 1.1 else 12,
    divergences = if (i == 1) 0 else 100, convergence_screen_pass = i == 1,
    max_rhat = if (i == 1) 1.001 else 1.5), file.path(runs[i], "test_summary.csv"), row.names = FALSE)
  write.csv(data.frame(cast = c("a", "b"), mse = c(i, 2*i)),
    file.path(runs[i], "checks", "validation_conditional_cast_metrics.csv"), row.names = FALSE)
  write.csv(data.frame(Region = "R1", rmse = i), file.path(runs[i], "checks", "validation_conditional_by_region.csv"), row.names = FALSE)
  write.csv(data.frame(x = 1), file.path(runs[i], "checks", "validation_conditional_observations.csv"), row.names = FALSE)
  dir.create(file.path(runs[i], "stan_csv"))
  write.csv(data.frame(lp__ = 1), file.path(runs[i], "stan_csv", "chain.csv"), row.names = FALSE)
}
# Legacy fit with a different summary schema, outside a manifest.
old <- file.path(root, "runs", "original_999"); dir.create(old, recursive = TRUE)
saveRDS(list(mode = "pilot", chains = 4), file.path(old, "run_config.rds"))
writeLines("done", file.path(old, "FIT_COMPLETE.txt"))
write.csv(data.frame(training_rmse = 0.01, min_tail_ess = 20), file.path(old, "test_summary.csv"), row.names = FALSE)
# Marker-only completion must be discovered, not silently omitted.
marker <- file.path(root, "runs", "marker_only"); dir.create(marker)
writeLines("done", file.path(marker, "FIT_COMPLETE.txt"))
# Historical benchmark manifest has no resolvable run path: keep and flag it.
bm <- file.path(root, "benchmarks", "threading_test"); dir.create(bm, recursive = TRUE)
write.table(data.frame(job_id = 123, threads = 2), file.path(bm, "manifest.tsv"), sep = "\t", row.names = FALSE)
out <- file.path(base, "report")
x <- collect_all(root, out)
stopifnot(nrow(x) == 7, sum(x$fit_complete_marker) == 5, sum(x$is_prior_predictive) == 1)
post <- read.csv(file.path(out, "all_completed_posterior_runs.csv"))
stopifnot(nrow(post) == 4, any(post$divergences == 100, na.rm = TRUE),
  any(post$snapshot_status == "fit_complete_checks_failed"), any(post$validation_conditional_rmse == 12, na.rm = TRUE))
stopifnot(x$min_tail_ess[x$label == "bad" & !is.na(x$label)] == 12,
  x$divergences[x$label == "unfinished" & !is.na(x$label)] == 7)
stopifnot(sum(!is.na(x$validation_set_md5)) == 2,
  length(unique(na.omit(x$validation_set_md5))) == 1,
  x$evaluation_type[grepl("original_999", x$run_path)] == "legacy_bounded_in_sample")
stopifnot(is.na(x$convergence_screen_pass[grepl("original_999", x$run_path)]))
casts <- read.csv(file.path(out, "all_validation_conditional_cast_metrics.csv"))
stopifnot(nrow(casts) == 4, length(unique(casts$review_uid)) == 2)
files <- list.files(file.path(out, "evidence"), recursive = TRUE)
stopifnot(!any(grepl("observations|stan_csv", files)))
stopifnot(length(readLines(file.path(out, "collection_warnings.txt"))) >= 1)
stopifnot(inherits(try(collect_all(root, out), silent = TRUE), "try-error"))
# A campaign with no completed results must not trigger the old zero-row assignment bug.
empty_root <- file.path(base, "empty"); dir.create(empty_root)
empty <- collect_all(empty_root, file.path(base, "empty_report"))
stopifnot(nrow(empty) == 0)
# Explicitly verify a full 102-case campaign is never truncated to a leader table.
many_root <- file.path(base, "many"); many_campaign <- file.path(many_root, "experiments", "many102")
dir.create(many_campaign, recursive = TRUE)
many <- data.frame(case_id = 1:102, label = paste0("case_", 1:102),
  run_dir = paste0("experiments/many102/runs/case_", sprintf("%02d", 1:102)))
write.table(many, file.path(many_campaign, "manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
for (i in 1:102) {
  r <- file.path(many_root, many$run_dir[i]); dir.create(r, recursive = TRUE)
  writeLines("done", file.path(r, "FIT_COMPLETE.txt"))
  write.csv(data.frame(training_rmse = 103 - i, divergences = i %% 3), file.path(r, "test_summary.csv"), row.names = FALSE)
}
many_out <- file.path(base, "many_report"); many_all <- collect_all(many_root, many_out)
many_post <- read.csv(file.path(many_out, "all_completed_posterior_runs.csv"))
stopifnot(nrow(many_all) == 102, nrow(many_post) == 102, setequal(many_post$case_id, 1:102),
  max(many_post$training_rmse) == 102, min(many_post$training_rmse) == 1)
cat("PASS: all-run retention, poor diagnostics, partial checks, prior separation, mixed schemas,\n",
    "legacy evaluation labels, shared validation identity, cast/Region evidence, raw-draw exclusions,\n",
    "unresolved benchmarks, empty campaigns, overwrite protection and all 102 cases without truncation.\n")
