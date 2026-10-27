#!/usr/bin/env Rscript
# Fit the current Gaussian model. This script does no wrangling or plotting.
# Run directories are unique; existing results are never silently reused.

fit_model <- function() {
  for (p in c("brms", "cmdstanr", "mgcv"))
    if (!requireNamespace(p, quietly = TRUE)) stop("Install required package: ", p)
  library(brms)
  env_num <- function(name, default, integer = FALSE) {
    x <- suppressWarnings(as.numeric(Sys.getenv(name, as.character(default))))
    if (!is.finite(x) || (integer && x != floor(x))) stop("Invalid ", name)
    if (integer) as.integer(x) else x
  }
  mode <- Sys.getenv("FIT_MODE", "pilot")
  if (!mode %in% c("pilot", "production", "prior")) stop("FIT_MODE: pilot, production, or prior.")
  data_file <- Sys.getenv("MODEL_DATA", "data/processed/model_data.rds")
  id <- Sys.getenv("SLURM_JOB_ID", format(Sys.time(), "%Y%m%dT%H%M%S"))
  run_dir <- Sys.getenv("RUN_DIR", file.path("runs", paste0(mode, "_", id)))
  if (dir.exists(run_dir)) stop("RUN_DIR already exists. Choose a new directory: ", run_dir)
  bundle <- readRDS(data_file)
  # Require freshly prepared observation-level data rather than an older bundle.
  if (!identical(bundle$config$smoothing, "none") || !identical(bundle$config$binned, FALSE))
    stop("Regenerate MODEL_DATA with the updated 01_prepare_data.R: unsmoothed, unbinned observations are required.")
  d <- bundle$model_data
  required <- c("chl_norm", "Region", "depth", "log10_distance", "MLD_temp_MIX", "date", "enso", "Ocean")
  if (!all(required %in% names(d))) stop("MODEL_DATA must be the output of 01_prepare_data.R.")
  if (!nrow(d) || anyNA(d[required]) || any(!is.finite(as.matrix(d[c("chl_norm", "depth", "log10_distance", "MLD_temp_MIX")]))))
    stop("Invalid model rows. Rerun preparation and inspect its QC tables.")
  chains <- env_num("MODEL_CHAINS", 4, TRUE)
  cores <- env_num("MODEL_CORES", chains, TRUE)
  threads <- env_num("MODEL_THREADS", 2, TRUE)
  iter <- env_num("MODEL_ITER", if (mode == "production") 5000 else 2000, TRUE)
  warmup <- env_num("MODEL_WARMUP", if (mode == "production") 2500 else 1000, TRUE)
  delta <- env_num("ADAPT_DELTA", 0.97)
  treedepth <- env_num("MAX_TREEDEPTH", 14, TRUE)
  seed <- env_num("MODEL_SEED", 1234, TRUE)
  if (chains < 2 || cores < 1 || cores > chains || threads < 1 || iter <= warmup || warmup < 1 ||
      delta <= 0 || delta >= 1 || treedepth < 1) stop("Invalid MCMC configuration.")
  allocated <- env_num("SLURM_CPUS_PER_TASK", cores * threads, TRUE)
  if (allocated < cores * threads) stop("Requested parallel chains x threads exceed allocated CPUs.")
  k <- c(env_num("K_DEPTH", 15, TRUE), env_num("K_DISTANCE", 5, TRUE), env_num("K_MLD", 5, TRUE))
  if (any(k < 3)) stop("Cubic regression spline marginal k must be at least 3.")
  # Conservative support check: do not fit a region-specific tensor to a
  # predictor with fewer distinct values than its requested marginal basis.
  support <- do.call(rbind, lapply(split(d, d$Region), function(x) {
    data.frame(Region = as.character(x$Region[1]),
      depth = length(unique(x$depth)), distance = length(unique(x$log10_distance)),
      mld = length(unique(x$MLD_temp_MIX)))
  }))
  if (any(as.matrix(support[-1]) < matrix(k, nrow(support), 3, byrow = TRUE))) {
    print(support)
    stop("Insufficient regional predictor support for k. Use more casts or an explicitly reviewed smaller k.")
  }
  stan_path <- Sys.getenv("CMDSTAN", "")
  if (nzchar(stan_path)) cmdstanr::set_cmdstan_path(stan_path)
  cmdstanr::cmdstan_version()  # Fail before creating outputs if CmdStan is absent.
  dir.create(run_dir, recursive = TRUE)
  run_dir <- normalizePath(run_dir, winslash = "/")
  stan_dir <- file.path(run_dir, "stan_csv")
  dir.create(stan_dir)
  # The k vector is inserted into the formula; all other model terms are unchanged.
  andrew_form <- bf(as.formula(sprintf(
    "chl_norm ~ 1 + Region + t2(depth, log10_distance, MLD_temp_MIX, by = Region, bs = c('cr', 'cr', 'cr'), k = c(%d, %d, %d)) + (1 | Ocean:Region:date) + (1 | Ocean:Region:enso)",
    k[1], k[2], k[3])), family = gaussian())
  # Preserve the LAST active priors in the supplied notebook. These historical
  # constants are not recalculated from this run or tightened to narrow intervals.
  priors <- c(
    prior(student_t(3, 0.73, 5.448), class = "Intercept"),
    prior(normal(0, 2), class = "b"),
    prior(student_t(3, 0, 1.5), class = "sd"),
    prior(student_t(3, 0, 1.5), class = "sds"),
    prior(student_t(3, 0, 2.5), class = "sigma")
  )
  # Validate the complete Stan data/basis before starting expensive sampling.
  stan_data <- make_standata(andrew_form, data = d, prior = priors)
  saveRDS(stan_data, file.path(run_dir, "stan_data.rds"), compress = "gzip")
  rm(stan_data); gc()
  writeLines(make_stancode(andrew_form, data = d, prior = priors,
    threads = threading(threads)), file.path(run_dir, "model.stan"))
  config <- list(mode = mode, chains = chains, cores = cores, threads = threads,
    iter = iter, warmup = warmup, adapt_delta = delta, max_treedepth = treedepth,
    k = k, seed = seed, input_md5 = tools::md5sum(data_file),
    formula = paste(deparse(andrew_form$formula), collapse = " "),
    cmdstan = cmdstanr::cmdstan_path(), cmdstan_version = cmdstanr::cmdstan_version())
  saveRDS(config, file.path(run_dir, "run_config.rds"))
  file.copy(data_file, file.path(run_dir, "model_data.rds"), overwrite = FALSE)
  for (script in c("01_prepare_data.R", "02_fit_model.R", "03_check_and_plot.R"))
    if (file.exists(script)) file.copy(script, file.path(run_dir, script))
  capture.output(print(config), file = file.path(run_dir, "run_config.txt"))
  capture.output(priors, file = file.path(run_dir, "priors.txt"))
  capture.output(sessionInfo(), file = file.path(run_dir, "session_info.txt"))
  # Persist CSVs on the HPC filesystem, never in a job's ephemeral /tmp.
  # Four parallel chains x two Stan threads = eight allocated CPUs.
  fit <- brm(andrew_form, data = d, prior = priors,
    chains = chains, cores = cores, threads = threading(threads),
    iter = iter, warmup = warmup, seed = seed, refresh = 50,
    backend = "cmdstanr", save_pars = save_pars(all = TRUE),
    sample_prior = if (mode == "prior") "only" else "no",
    control = list(adapt_delta = delta, max_treedepth = treedepth),
    output_dir = stan_dir, file = file.path(run_dir, "fit"), file_refit = "always")
  # Hydrate draws before serialization so plotting can load the RDS even if a
  # downloaded copy no longer has the original absolute CSV paths available.
  invisible(posterior::as_draws_array(fit))
  saveRDS(fit, file.path(run_dir, "fit.rds"), compress = "gzip")
  capture.output(summary(fit), file = file.path(run_dir, "model_summary.txt"))
  writeLines("Sampling completed; inspect diagnostics before scientific interpretation.",
             file.path(run_dir, "FIT_COMPLETE.txt"))
  message("Fit saved in ", run_dir)
  invisible(fit)
}

# In RStudio: source('02_fit_model.R'); fit_model().
if (sys.nframe() == 0L) fit_model()
