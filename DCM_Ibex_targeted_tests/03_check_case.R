#!/usr/bin/env Rscript
# Check one completed exploratory fit and create bounded predictive figures.
experiment <- Sys.getenv("EXPERIMENT_DIR")
source(file.path(experiment, "common.R"))
source(file.path(experiment, "campaign.R"))
case <- read_case(experiment, Sys.getenv("CASE_ID"))
run <- file.path(experiment, "runs", sprintf("case_%02d", case$case_id))
if (!file.exists(file.path(run, "FIT_COMPLETE.txt"))) {
  message("Fit did not complete; collector will retain its failure status.")
  quit(status = 0L)
}
tryCatch({
  for (package in c("brms", "posterior", "ggplot2", "bayesplot"))
    if (!requireNamespace(package, quietly = TRUE)) stop("Missing package: ", package)
  fit <- readRDS(file.path(run, "fit.rds"))
  d <- readRDS(file.path(run, "model_data.rds"))$model_data
  config <- readRDS(file.path(run, "run_config.rds"))
  checks <- file.path(run, "checks")
  dir.create(checks, showWarnings = FALSE)
  csv <- function(x, name) write.csv(x, file.path(checks, name), row.names = FALSE)
  draws <- posterior::as_draws_array(fit)
  diagnostics <- as.data.frame(posterior::summarise_draws(
    draws, "mean", "sd", "rhat", "ess_bulk", "ess_tail", "mcse_mean"))
  csv(diagnostics, "parameter_diagnostics.csv")
  varying <- diagnostics[is.finite(diagnostics$sd) & diagnostics$sd > 0, ]
  targets <- diagnostics[grepl("^(b_|sd_|sds_)|^sigma$", diagnostics$variable), ]
  target_ess <- if (nrow(targets) && all(is.finite(targets$ess_bulk)))
    median(targets$ess_bulk) else NA_real_
  rhat <- finite_max(varying$rhat)
  bulk <- finite_min(varying$ess_bulk)
  tail <- finite_min(varying$ess_tail)

  files <- list.files(file.path(run, "stan_csv"), "\\.csv$", full.names = TRUE)
  if (length(files) != case$chains) stop("Unexpected number of chain CSV files.")
  times <- do.call(rbind, lapply(files, chain_elapsed))
  csv(times, "chain_times.csv")
  sampler <- do.call(rbind, lapply(seq_along(files), function(chain) {
    x <- read_sampler_csv(files[chain], case$iter, case$warmup)
    x$chain <- chain
    x
  }))
  csv(sampler, "sampler_iterations_including_warmup.csv")
  phase_summary <- do.call(rbind, lapply(split(sampler, interaction(sampler$chain, sampler$phase)), function(x) {
    data.frame(chain = x$chain[1], phase = x$phase[1], iterations = nrow(x),
      divergences = sum(x$divergent__),
      treedepth_hits = sum(x$treedepth__ >= case$treedepth),
      median_leapfrogs = median(x$n_leapfrog__),
      maximum_leapfrogs = max(x$n_leapfrog__),
      ebfmi = if (nrow(x) > 2 && var(x$energy__) > 0)
        mean(diff(x$energy__)^2) / var(x$energy__) else NA_real_)
  }))
  csv(phase_summary, "sampler_summary_by_phase.csv")
  post <- phase_summary[phase_summary$phase == "sampling", ]
  divergences <- sum(post$divergences)
  tree_hits <- sum(post$treedepth_hits)
  tree_fraction <- tree_hits / sum(post$iterations)
  posterior_rows <- sampler[sampler$phase == "sampling", ]
  mean_leapfrogs <- mean(posterior_rows$n_leapfrog__)
  ebfmi <- finite_min(post$ebfmi)
  passed <- all(is.finite(c(rhat, bulk, tail, ebfmi))) &&
    rhat <= 1.01 && bulk >= 400 && tail >= 400 &&
    divergences == 0 && ebfmi >= .3
  sample_seconds <- finite_max(times$sampling_seconds)
  hardware <- readLines(file.path(run, "hardware.txt"))
  cpu <- grep("^Model name:", hardware, value = TRUE)
  cpu <- if (length(cpu)) trimws(sub("^Model name:", "", cpu[1])) else "unknown"
  core_hours <- config$allocated_cpus * config$fit_seconds / 3600
  result <- data.frame(case, job_id = config$job_id, cpu_model = cpu,
    hostname = config$hostname, data_md5 = config$data_md5,
    source_md5 = config$source_md5, observations = nrow(d),
    casts = length(unique(d$cast_key)), maximum_depth = max(d$depth),
    fit_seconds = config$fit_seconds,
    brm_seconds_including_compile = config$brm_seconds_including_compile,
    warmup_seconds = finite_max(times$warmup_seconds),
    sampling_seconds = sample_seconds, requested_core_hours = core_hours,
    retained_draws = posterior::ndraws(draws), median_target_bulk_ess = target_ess,
    target_ess_per_sampling_second = target_ess / sample_seconds,
    target_ess_per_requested_core_hour = target_ess / core_hours,
    max_rhat = rhat, min_bulk_ess = bulk, min_tail_ess = tail,
    divergences = divergences, treedepth_hits = tree_hits,
    treedepth_hit_fraction = tree_fraction, mean_leapfrogs = mean_leapfrogs, min_ebfmi = ebfmi,
    max_mcse_over_sd = finite_max(varying$mcse_mean / varying$sd),
    convergence_screen_pass = passed, tree_efficiency_warning = tree_hits > 0,
    diagnostic_screen_pass = passed && tree_hits == 0)
  # Retain useful runtime/diagnostic metrics even if a subsequent figure fails.
  write.csv(result, file.path(run, "diagnostic_summary.csv"), row.names = FALSE)

  trace_pars <- intersect(c("b_Intercept","sigma","nu"), posterior::variables(draws))
  if(length(trace_pars)) ggplot2::ggsave(file.path(checks,"trace_plots.png"),
    bayesplot::mcmc_trace(draws,pars=trace_pars),width=10,height=6,dpi=160,bg="white")

  source(file.path(experiment, "predictive.R"))
  result <- evaluate_predictions(fit, experiment, run, case, result)
  write.csv(result, file.path(run,"test_summary.csv"), row.names=FALSE)
  writeLines(c("Convergence and tree efficiency are reported separately.",
    "Convergence screen: Rhat <=1.01, bulk/tail ESS >=400, zero sampling divergences, E-BFMI >=0.3.",
    "Parameter MCSE/SD is also reported; inspect uncertainty of scientific predictions.",
    "Validation is on entire casts, excluded from all training subsets.",
    "Primary validation may share dates/ENSO levels; the date-blocked panel excludes nested date groups.",
    "Population predictive intervals condition on zero group effects; they do not integrate new-group variability.",
    "Raw observed peak-depth errors are noisy descriptive summaries, not the sole ranking criterion.",
    "Lag-1 residual correlation is a screening statistic; row spacing may be irregular.",
    "AR(1), where requested, describes adjacent observed rows and is not a continuous correlation length in metres.",
    "Depth-bin peak/shape/AUC metrics refer to that window; DCM classification is omitted in depth windows.",
    "Exclusive probability bins differ from exceedance probabilities; expected-curve and new-observation probabilities are labelled separately.",
    "Inspect matched repeat-seed profile predictions for Monte Carlo stability."),
    file.path(checks,"interpretation.txt"))
  writeLines("Complete",file.path(run,"CHECKS_COMPLETE.txt"))
}, error=function(error) {
  writeLines(conditionMessage(error),file.path(run,"CHECKS_FAILED.txt"))
  stop(error)
})
