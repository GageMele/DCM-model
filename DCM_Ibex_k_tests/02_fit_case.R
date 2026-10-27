#!/usr/bin/env Rscript
# Fit one case in a unique directory, using the frozen original formula.
experiment <- Sys.getenv("EXPERIMENT_DIR")
source(file.path(experiment, "common.R"))
case <- read_case(experiment, Sys.getenv("CASE_ID"))
for (package in c("brms", "cmdstanr", "mgcv", "posterior"))
  if (!requireNamespace(package, quietly = TRUE)) stop("Missing package: ", package)
run <- file.path(experiment, "runs", sprintf("case_%02d", case$case_id))
if (dir.exists(run)) stop("Run exists. Submit a fresh experiment to repeat it.")
dir.create(run, recursive = TRUE)
tryCatch({
  start <- proc.time()[["elapsed"]]
  k <- case_k(case)
  data_path <- file.path(experiment, "data", paste0(case$subset, ".rds"))
  bundle <- readRDS(data_path)
  d <- validate_bundle(bundle)
  check_support(d, k)
  spec <- load_original_spec(file.path(experiment, "02_fit_model.R"), k)
  statistics <- readRDS(file.path(experiment, "full_response_statistics.rds"))
  priors <- selected_priors(case$prior, spec$current_priors, statistics)
  if (case$mode == "prior" && case$prior == "default")
    stop("Default flat coefficient priors do not define a prior-predictive distribution.")
  if (nzchar(Sys.getenv("CMDSTAN"))) cmdstanr::set_cmdstan_path(Sys.getenv("CMDSTAN"))
  cmdstanr::cmdstan_version()
  requested <- case$chains * max(1L, case$threads)
  allocated <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", as.character(requested)))
  if (allocated < requested) stop("Insufficient allocated CPUs.")
  Sys.setenv(STAN_NUM_THREADS = max(1L, case$threads))
  thread_spec <- if (case$threads == 0L) NULL else brms::threading(case$threads)
  config <- list(case = case, mode = if (case$mode == "prior") "prior" else "pilot",
    chains = case$chains, cores = case$chains, threads = case$threads,
    iter = case$iter, warmup = case$warmup, adapt_delta = case$adapt_delta,
    max_treedepth = case$treedepth, k = k, seed = 1234L,
    job_id = Sys.getenv("SLURM_JOB_ID"), allocated_cpus = allocated,
    hostname = Sys.info()[["nodename"]],
    data_md5 = unname(tools::md5sum(data_path)),
    source_md5 = unname(tools::md5sum(file.path(experiment, "02_fit_model.R"))),
    formula = paste(deparse(spec$formula$formula), collapse = " "),
    cmdstan_version = cmdstanr::cmdstan_version())
  saveRDS(config, file.path(run, "run_config.rds"))
  file.copy(data_path, file.path(run, "model_data.rds"))
  saveRDS(priors, file.path(run, "requested_priors.rds"))
  capture.output(priors, file = file.path(run, "requested_priors.txt"))
  capture.output(config, file = file.path(run, "run_config.txt"))
  hardware <- system2("lscpu", stdout = TRUE)
  writeLines(hardware, file.path(run, "hardware.txt"))
  capture.output(sessionInfo(), file = file.path(run, "session_info.txt"))
  csv_dir <- file.path(run, "stan_csv")
  dir.create(csv_dir)
  writeLines(brms::make_stancode(spec$formula, data = d, prior = priors,
    threads = thread_spec), file.path(run, "model.stan"))
  cat("Starting", case$label, "with", nrow(d), "rows and",
      length(unique(d$cast_key)), "whole casts\n")
  sampling_start <- proc.time()[["elapsed"]]
  fit <- brms::brm(spec$formula, data = d, prior = priors,
    backend = "cmdstanr", chains = case$chains, cores = case$chains,
    threads = thread_spec, iter = case$iter, warmup = case$warmup,
    seed = 1234L, refresh = 10L, save_warmup = TRUE,
    control = list(adapt_delta = case$adapt_delta, max_treedepth = case$treedepth),
    save_pars = brms::save_pars(all = TRUE),
    sample_prior = if (case$mode == "prior") "only" else "no",
    output_dir = csv_dir)
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
