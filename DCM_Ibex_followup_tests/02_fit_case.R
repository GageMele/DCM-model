#!/usr/bin/env Rscript
# Fit one case in a unique directory, using its requested equation exactly.
experiment <- Sys.getenv("EXPERIMENT_DIR")
source(file.path(experiment, "common.R"))
source(file.path(experiment, "campaign.R"))
if (!file.exists(file.path(experiment, "PREPARATION_COMPLETE.txt")))
  stop("Preparation has not completed successfully.")
case <- read_case(experiment, Sys.getenv("CASE_ID"))
for (package in c("brms", "cmdstanr", "mgcv", "posterior"))
  if (!requireNamespace(package, quietly = TRUE)) stop("Missing package: ", package)
run <- file.path(experiment, "runs", sprintf("case_%02d", case$case_id))
if (dir.exists(run)) stop("Run exists. Submit a fresh experiment to repeat it.")
dir.create(run, recursive = TRUE)
tryCatch({
  start <- proc.time()[["elapsed"]]
  k <- c(case$k_depth, case$k_distance, case$k_mld)
  data_path <- data_file(experiment,case)
  d <- readRDS(data_path)
  spec <- campaign_spec(case,readRDS(file.path(experiment,"region_levels.rds")))
  priors <- spec$priors
  if (case$mode == "prior" && case$prior == "mu_default")
    stop("Default flat coefficient priors do not define a prior-predictive distribution.")
  if (nzchar(Sys.getenv("CMDSTAN"))) cmdstanr::set_cmdstan_path(Sys.getenv("CMDSTAN"))
  cmdstanr::cmdstan_version()
  requested <- case$chains * max(1L, case$threads)
  allocated <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", as.character(requested)))
  if (allocated < requested) stop("Insufficient allocated CPUs.")
  Sys.setenv(STAN_NUM_THREADS = max(1L, case$threads))
  thread_spec <- thread_specification(case)
  config <- list(case = case, mode = "followup",
    chains = case$chains, cores = case$cores, threads = case$threads,
    iter = case$iter, warmup = case$warmup, adapt_delta = case$adapt_delta,
    max_treedepth = case$treedepth, refresh = case$refresh, k = k, seed = case$seed,
    job_id = Sys.getenv("SLURM_JOB_ID"), allocated_cpus = allocated,
    hostname = Sys.info()[["nodename"]],
    data_md5 = unname(tools::md5sum(data_path)),
    source_md5 = unname(tools::md5sum(file.path(experiment, "campaign.R"))),
    formula = paste(capture.output(spec$formula), collapse = " "),
    validation_md5 = unname(tools::md5sum(validation_file(experiment,case))),
    cmdstan_version = cmdstanr::cmdstan_version())
  saveRDS(config, file.path(run, "run_config.rds"))
  saveRDS(list(model_data=d), file.path(run, "model_data.rds"))
  saveRDS(priors, file.path(run, "requested_priors.rds"))
  capture.output(priors, file = file.path(run, "requested_priors.txt"))
  capture.output(config, file = file.path(run, "run_config.txt"))
  hardware <- system2("lscpu", stdout = TRUE)
  writeLines(hardware, file.path(run, "hardware.txt"))
  capture.output(sessionInfo(), file = file.path(run, "session_info.txt"))
  csv_dir <- file.path(run, "stan_csv")
  dir.create(csv_dir)
  writeLines(brms::make_stancode(spec$formula, data = d, prior = priors,
    threads = thread_spec, sample_prior=if(case$mode=="prior") "only" else "no"), file.path(run, "model.stan"))
  cat("Starting", case$label, "with", nrow(d), "rows and",
      length(unique(d$cast_key)), "whole casts\n")
  sampling_start <- proc.time()[["elapsed"]]
  fit <- brms::brm(spec$formula, data = d, prior = priors,
    backend = "cmdstanr", chains = case$chains, cores = case$cores,
    threads = thread_spec, iter = case$iter, warmup = case$warmup,
    seed = case$seed, refresh = case$refresh, save_warmup = TRUE, thin = 1L,
    control = list(adapt_delta = case$adapt_delta, max_treedepth = case$treedepth),
    save_pars = brms::save_pars(all = TRUE),
    sample_prior = if (case$mode == "prior") "only" else "no",
    output_dir = csv_dir, file = file.path(run, "fit"), file_refit = "always")
  brm_seconds <- proc.time()[["elapsed"]] - sampling_start
  # Hydrate draws before serialization so a downloaded fit remains usable.
  invisible(posterior::as_draws_array(fit))
  saveRDS(fit, file.path(run, "fit.rds"))
  config$brm_seconds_including_compile <- brm_seconds
  config$fit_seconds <- proc.time()[["elapsed"]] - start
  saveRDS(config, file.path(run, "run_config.rds"))
  resolved <- brms::prior_summary(fit)
  write.csv(as.data.frame(resolved), file.path(run, "resolved_priors.csv"), row.names = FALSE)
  capture.output(resolved, file = file.path(run, "resolved_priors.txt"))
  writeLines("Fit complete; inspect exploratory diagnostics.", file.path(run, "FIT_COMPLETE.txt"))
}, error = function(error) {
  writeLines(conditionMessage(error), file.path(run, "FIT_FAILED.txt"))
  stop(error)
})
