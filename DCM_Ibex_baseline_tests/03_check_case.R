#!/usr/bin/env Rscript
# Check one completed exploratory fit and create bounded predictive figures.
experiment <- Sys.getenv("EXPERIMENT_DIR")
source(file.path(experiment, "common.R"))
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
    draws, "mean", "sd", "rhat", "ess_bulk", "ess_tail"))
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
    divergences == 0 && tree_hits == 0 && ebfmi >= .3
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
    diagnostic_screen_pass = passed)
  # Retain useful runtime/diagnostic metrics even if a subsequent figure fails.
  write.csv(result, file.path(run, "diagnostic_summary.csv"), row.names = FALSE)

  # Check a region-balanced set of existing observations: in-sample PPC only.
  set.seed(7201L)
  rows <- unlist(lapply(split(seq_len(nrow(d)), d$Region), function(ids) {
    ids[sample.int(length(ids), min(60L, length(ids)))]
  }), use.names = FALSE)
  check_data <- d[rows, , drop = FALSE]
  n_draws <- min(100L, posterior::ndraws(draws))
  draw_ids <- sort(sample.int(posterior::ndraws(draws), n_draws))
  mu <- brms::posterior_epred(fit, newdata = check_data,
                            re_formula = NULL, draw_ids = draw_ids)
  yrep <- brms::posterior_predict(fit, newdata = check_data,
                                 re_formula = NULL, draw_ids = draw_ids)
  q <- apply(yrep, 2L, quantile, probs = c(.025, .975))
  observed <- check_data$chl_norm
  result$training_check_observations <- length(observed)
  result$training_rmse <- sqrt(mean((observed - colMeans(mu))^2))
  result$training_ppc_95_coverage <- mean(observed >= q[1, ] & observed <= q[2, ])
  result$predictive_fraction_below_minus_one <- mean(yrep < -1)
  result$prediction_interpretation <- if (case$mode == "prior")
    "prior predictive; observed values shown for scale comparison" else
    "in-sample posterior predictive; not independent predictive accuracy"
  result$prior_interpretation <- if (case$prior == "empirical")
    "empirical Bayes from full 150m response; exploratory only" else case$prior
  csv(data.frame(check_data, fitted = colMeans(mu),
    predictive_lower = q[1, ], predictive_upper = q[2, ]), "check_observations.csv")
  trace_pars <- intersect(c("b_Intercept", "sigma"), posterior::variables(draws))
  if (length(trace_pars)) ggplot2::ggsave(file.path(checks, "trace_plots.png"),
    bayesplot::mcmc_trace(draws, pars = trace_pars), width = 10, height = 6, dpi = 160, bg = "white")
  ggplot2::ggsave(file.path(checks, "predictive_density.png"),
    bayesplot::ppc_dens_overlay(observed, head(yrep, 30L)),
    width = 10, height = 6, dpi = 160, bg = "white")
  # Profiles span only the tested depths and use real cast predictor pairs.
  casts <- unique(d[c("cast_key", "Region", "Ocean", "date", "enso",
                       "log10_distance", "MLD_temp_MIX")])
  examples <- do.call(rbind, lapply(split(casts, casts$Region), function(x) {
    x[order(x$MLD_temp_MIX, x$log10_distance)[ceiling(nrow(x) / 2)], , drop = FALSE]
  }))
  depth_grid <- seq(10, 150, by = 2)
  grid <- examples[rep(seq_len(nrow(examples)), each = length(depth_grid)), ]
  grid$depth <- rep(depth_grid, nrow(examples))
  profile <- brms::posterior_epred(fit, newdata = grid,
                                  re_formula = NA, draw_ids = draw_ids)
  grid$median <- apply(profile, 2L, median)
  grid$lower <- apply(profile, 2L, quantile, .025)
  grid$upper <- apply(profile, 2L, quantile, .975)
  csv(grid, "example_population_profiles.csv")
  p <- ggplot2::ggplot(grid, ggplot2::aes(x = depth, y = median)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = .25) +
    ggplot2::geom_line() + ggplot2::facet_wrap(~ Region, scales = "free_y") +
    ggplot2::labs(x = "Depth (m)", y = "Fractional chlorophyll enhancement",
      title = paste(case$label, "— exploratory population profiles")) +
    ggplot2::theme_minimal()
  ggplot2::ggsave(file.path(checks, "example_population_profiles.png"), p,
                  width = 12, height = 8, dpi = 160, bg = "white")
  write.csv(result, file.path(run, "test_summary.csv"), row.names = FALSE)
  writeLines(c(paste("Diagnostic screen:", passed),
    "R-hat <=1.01; bulk/tail ESS >=400; no divergences/tree hits; E-BFMI >=0.3.",
    "This conservative screen does not establish scientific validity.",
    "Tree saturation is an efficiency concern; inspect R-hat, ESS and predictive behavior.",
    "Training PPC values are not independent predictive accuracy.",
    "Custom priors are the historical constants; default priors are resolved by the installed brms version."),
    file.path(checks, "interpretation.txt"))
  writeLines("Complete", file.path(run, "CHECKS_COMPLETE.txt"))
}, error = function(error) {
  writeLines(conditionMessage(error), file.path(run, "CHECKS_FAILED.txt"))
  stop(error)
})
