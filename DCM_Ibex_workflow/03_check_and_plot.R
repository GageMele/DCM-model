#!/usr/bin/env Rscript
# Read an existing fit, check it, then create tables and figures without refitting.
# Required: RUN_DIR pointing at one completed run; no global CTD/fit objects needed.

check_and_plot <- function() {
  for (p in c("brms", "posterior", "ggplot2", "bayesplot"))
    if (!requireNamespace(p, quietly = TRUE)) stop("Install required package: ", p)
  library(ggplot2)
  run_dir <- Sys.getenv("RUN_DIR", "")
  if (!nzchar(run_dir)) stop("Set RUN_DIR to a completed model run.")
  out_dir <- Sys.getenv("PLOT_DIR", file.path(run_dir, "checks_and_plots"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fit <- readRDS(file.path(run_dir, "fit.rds"))
  d <- readRDS(file.path(run_dir, "model_data.rds"))$model_data
  config <- readRDS(file.path(run_dir, "run_config.rds"))
  draws_requested <- suppressWarnings(as.integer(Sys.getenv("PLOT_DRAWS", "1000")))
  if (is.na(draws_requested) || draws_requested < 50) stop("PLOT_DRAWS must be at least 50.")
  set.seed(42)
  csv <- function(x, name) write.csv(x, file.path(out_dir, name), row.names = FALSE)
  save_plot <- function(p, name, width = 12, height = 8.5) {
    for (ext in c("png", "pdf")) ggsave(file.path(out_dir, paste0(name, ".", ext)), p,
      width = width, height = height, units = "in", dpi = 300, bg = "white")
  }
  summarize_matrix <- function(x) {
    data.frame(mean = colMeans(x), median = apply(x, 2, median),
      lower = apply(x, 2, quantile, probs = 0.025), upper = apply(x, 2, quantile, probs = 0.975))
  }
  # Examine ALL sampled parameters, including spline and group coefficients.
  draws <- posterior::as_draws_array(fit)
  diag <- as.data.frame(posterior::summarise_draws(draws, "mean", "sd", "rhat", "ess_bulk", "ess_tail"))
  csv(diag, "parameter_diagnostics.csv")
  # Constants/quantities without a meaningful R-hat are excluded from these flags.
  bad <- diag[is.finite(diag$rhat) & (diag$rhat > 1.01 | diag$ess_bulk < 400 | diag$ess_tail < 400), ]
  csv(bad, "flagged_parameters.csv")
  # Select the SAME posterior draws across all prediction blocks; this preserves
  # posterior dependence for ocean averages and DCM metrics.
  total_draws <- posterior::ndraws(draws)
  draw_ids <- sort(sample.int(total_draws, min(draws_requested, total_draws)))
  trace_names <- intersect(c("b_Intercept", "sigma"), posterior::variables(draws))
  if (length(trace_names)) save_plot(bayesplot::mcmc_trace(draws, pars = trace_names), "trace_plots")
  rm(draws); gc()
  nuts <- brms::nuts_params(fit)
  csv(nuts, "sampler_diagnostics.csv")
  chain_sampler <- do.call(rbind, lapply(split(nuts, nuts$Chain), function(x) {
    energy <- x$Value[x$Parameter == "energy__"]
    data.frame(chain = x$Chain[1],
      divergences = sum(x$Value[x$Parameter == "divergent__"]),
      treedepth_hits = sum(x$Value[x$Parameter == "treedepth__"] >= config$max_treedepth),
      ebfmi = if (length(energy) > 2 && var(energy) > 0) mean(diff(energy)^2) / var(energy) else NA_real_)
  }))
  csv(chain_sampler, "sampler_summary.csv")
  diagnostic_pass <- !nrow(bad) && all(chain_sampler$divergences == 0) &&
    all(chain_sampler$treedepth_hits == 0) && all(is.finite(chain_sampler$ebfmi) & chain_sampler$ebfmi >= 0.3)
  writeLines(c(paste("Diagnostic screening passed:", diagnostic_pass),
    "Screen: R-hat <= 1.01; bulk/tail ESS >= 400; zero divergences/tree-depth hits; E-BFMI >= 0.3.",
    "These checks are necessary screening, not proof of model adequacy or predictive skill."),
    file.path(out_dir, "diagnostic_status.txt"))
  if (!diagnostic_pass) warning("Diagnostic flags found. Review tables before interpreting plots.")

  # Posterior predictive checks are bounded to 2,000 observations for memory.
  # Random effects are included for these existing observations (re_formula=NULL).
  check_rows <- sort(sample.int(nrow(d), min(2000L, nrow(d))))
  check_data <- d[check_rows, ]
  mu_check <- brms::posterior_epred(fit, newdata = check_data, re_formula = NULL, draw_ids = draw_ids)
  pp_ids <- head(draw_ids, min(100L, length(draw_ids)))
  yrep <- brms::posterior_predict(fit, newdata = check_data, re_formula = NULL, draw_ids = pp_ids)
  save_plot(bayesplot::ppc_dens_overlay(check_data$chl_norm, yrep[seq_len(min(50L, nrow(yrep))), , drop = FALSE]),
            if (config$mode == "prior") "prior_predictive_density" else "posterior_predictive_density")
  csv(data.frame(check_data, fitted = colMeans(mu_check), residual = check_data$chl_norm - colMeans(mu_check)),
      "check_observations.csv")
  if (config$mode == "prior") {
    message("Prior-only run: inspect prior predictive density; posterior interpretation skipped.")
    return(invisible(chain_sampler))
  }
  # In-sample Bayesian R2 on the documented checking subset, not full-data R2.
  r2 <- brms::bayes_R2(fit, newdata = check_data, re_formula = NULL,
                     ndraws = min(500L, length(draw_ids)))
  csv(as.data.frame(r2), "bayes_r2_check_subset.csv")
  residual_data <- data.frame(check_data, fitted = colMeans(mu_check),
                              residual = check_data$chl_norm - colMeans(mu_check))
  for (v in c("fitted", "depth", "log10_distance", "MLD_temp_MIX")) {
    p <- ggplot(residual_data, aes(x = .data[[v]], y = residual)) + geom_point(alpha = 0.3, size = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed") + facet_wrap(~ Region) + theme_minimal()
    save_plot(p, paste0("residuals_vs_", v))
  }
  save_plot(ggplot(residual_data, aes(sample = residual)) + stat_qq(alpha = 0.3) +
    stat_qq_line() + facet_wrap(~ Region) + theme_minimal(), "residual_qq")
  rm(mu_check, yrep); gc()

  # Use actual CAST-level distance/MLD pairs. This avoids constructing fictitious
  # combinations in regions where distance and MLD are strongly related.
  # These five scenarios change an observed pair; they do not isolate causality.
  profiles <- probabilities <- metrics <- grids <- list()
  reference_draws <- list()
  reference_grids <- list()
  thresholds_ratio <- c(1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
  # Metrics characterize positive fractional enhancement, not a fitted biomass DCM.
  # Half-maximum width is the contiguous component containing the maximum.
  metric_one <- function(v, z) {
    top <- max(v); peak <- which.max(v)
    ties <- which(abs(v - top) <= 1e-10)
    peak_depth <- mean(z[ties])
    has_peak <- top > 0 && peak_depth > 15
    edge <- any(ties %in% c(1L, length(z)))
    width <- NA_real_; clipped <- NA
    if (has_peak) {
      left <- right <- peak
      while (left > 1 && v[left - 1] >= top / 2) left <- left - 1L
      while (right < length(v) && v[right + 1] >= top / 2) right <- right + 1L
      low <- z[left]; high <- z[right]
      if (left > 1) low <- z[left - 1] + (top / 2 - v[left - 1]) * (z[left] - z[left - 1]) / (v[left] - v[left - 1])
      if (right < length(v)) high <- z[right] + (top / 2 - v[right]) * (z[right + 1] - z[right]) / (v[right + 1] - v[right])
      clipped <- left == 1L || right == length(v)
      if (!clipped) width <- high - low
    }
    w <- pmax(v, 0)
    wsd <- if (sum(w) > 0) sqrt(sum(w * (z - sum(w * z) / sum(w))^2) / sum(w)) else NA_real_
    c(peak_depth = peak_depth, max_fractional_enhancement = top, subsurface_enhancement_peak = has_peak,
      peak_at_depth_boundary = edge, enhancement_fwhm_m = width, fwhm_clipped = clipped, positive_enhancement_weighted_sd_m = wsd)
  }
  for (region in levels(d$Region)) {
    rd <- d[as.character(d$Region) == region, ]
    casts <- rd[!duplicated(rd$cast_key), ]
    sx <- sd(casts$log10_distance); sy <- sd(casts$MLD_temp_MIX)
    if (!is.finite(sx) || sx == 0) sx <- 1
    if (!is.finite(sy) || sy == 0) sy <- 1
    targets <- rbind(
      reference = c(median(casts$log10_distance), median(casts$MLD_temp_MIX)),
      near_shore = c(quantile(casts$log10_distance, 0.1), median(casts$MLD_temp_MIX)),
      offshore = c(quantile(casts$log10_distance, 0.9), median(casts$MLD_temp_MIX)),
      shallow_mld = c(median(casts$log10_distance), quantile(casts$MLD_temp_MIX, 0.1)),
      deep_mld = c(median(casts$log10_distance), quantile(casts$MLD_temp_MIX, 0.9)))
    for (scenario in rownames(targets)) {
      target <- targets[scenario, ]
      ix <- which.min(((casts$log10_distance - target[1]) / sx)^2 + ((casts$MLD_temp_MIX - target[2]) / sy)^2)
      representative <- casts[ix, ]
      # Same 1 m grid across regions, limited to observed depth support below.
      z <- seq(10, 250, by = 1)
      grid <- data.frame(Region = factor(region, levels = levels(d$Region)), depth = z,
        log10_distance = representative$log10_distance,
        MLD_temp_MIX = representative$MLD_temp_MIX,
        Ocean = representative$Ocean, date = representative$date, enso = representative$enso)
      grid$scenario <- scenario
      grid$reference_cast_key <- representative$cast_key
      grid$within_region_depth_range <- z >= min(rd$depth) & z <= max(rd$depth)
      # re_formula=NA sets group effects to zero. It does NOT integrate over the
      # random-effect distribution for a new date/cast/ENSO group.
      ep <- brms::posterior_epred(fit, newdata = grid, re_formula = NA, draw_ids = draw_ids)
      pp <- brms::posterior_predict(fit, newdata = grid, re_formula = NA, draw_ids = draw_ids)
      key <- paste(region, scenario, sep = "__")
      prof <- cbind(grid, summarize_matrix(ep))
      prof$ratio_median <- prof$median + 1
      prof$ratio_lower <- prof$lower + 1; prof$ratio_upper <- prof$upper + 1
      # Preserve unsupported-depth rows in CSV for audit; omit them from figures.
      profiles[[key]] <- prof
      grids[[key]] <- grid
      probabilities[[key]] <- do.call(rbind, lapply(thresholds_ratio, function(r) data.frame(
        grid, ratio_threshold = r, chl_norm_threshold = r - 1,
        probability_expected_profile = colMeans(ep > r - 1),
        probability_observation = colMeans(pp > r - 1))))
      probabilities[[key]]$probability_impossible_negative_ratio <- rep(colMeans(pp < -1), length(thresholds_ratio))
      if (scenario == "reference") {
        valid <- grid$within_region_depth_range
        metrics[[region]] <- data.frame(Region = region, draw_id = draw_ids,
          t(apply(ep[, valid, drop = FALSE], 1, metric_one, z = z[valid])))
        reference_draws[[region]] <- ep
        reference_grids[[region]] <- grid
      }
      rm(ep, pp); gc()
    }
  }
  profile <- do.call(rbind, profiles)
  prob <- do.call(rbind, probabilities)
  metric <- do.call(rbind, metrics)
  csv(profile, "posterior_profiles.csv"); csv(prob, "enhancement_exceedance_probabilities.csv")
  csv(do.call(rbind, grids), "prediction_scenarios.csv"); csv(metric, "reference_profile_peak_metrics_draws.csv")
  metric_summary <- do.call(rbind, lapply(split(metric, metric$Region), function(x) {
    measures <- setdiff(names(x), c("Region", "draw_id"))
    do.call(rbind, lapply(measures, function(v) {
      values <- x[[v]][is.finite(x[[v]])]
      data.frame(Region = as.character(x$Region[1]), metric = v, valid_draws = length(values),
        mean = if (length(values)) mean(values) else NA_real_,
        median = if (length(values)) median(values) else NA_real_,
        lower = if (length(values)) unname(quantile(values, 0.025)) else NA_real_,
        upper = if (length(values)) unname(quantile(values, 0.975)) else NA_real_)
    }))
  }))
  csv(metric_summary, "reference_profile_peak_metrics_summary.csv")
  ref <- profile[profile$scenario == "reference" & profile$within_region_depth_range, ]
  profile_plot <- function(x) ggplot(x, aes(depth, ratio_median, color = Region, fill = Region)) +
    geom_ribbon(aes(ymin = ratio_lower, ymax = ratio_upper), alpha = 0.2, color = NA) +
    geom_line(linewidth = 0.6) + geom_hline(yintercept = 1, linetype = "dashed") +
    scale_x_reverse() + coord_flip() + facet_wrap(~ Region, scales = "free_x") +
    labs(x = "Depth (m)", y = "Chlorophyll / surface baseline", subtitle = "Expected profile: median and 95% credible interval; group effects set to zero") +
    theme_minimal() + theme(legend.position = "none")
  save_plot(profile_plot(ref), "reference_profiles")
  for (ocean in unique(as.character(d$Ocean))) {
    regions <- unique(as.character(d$Region[as.character(d$Ocean) == ocean]))
    save_plot(profile_plot(ref[as.character(ref$Region) %in% regions, ]), paste0("reference_profiles_", ocean))
  }
  save_plot(ggplot(profile[profile$within_region_depth_range, ], aes(depth, ratio_median, color = scenario)) +
    geom_line(linewidth = 0.6) + scale_x_reverse() + coord_flip() + facet_wrap(~ Region, scales = "free_x") +
    labs(x = "Depth (m)", y = "Chlorophyll / surface baseline", subtitle = "Scenarios use observed distance/MLD pairs") +
    theme_minimal() + theme(legend.position = "bottom"), "observed_pair_scenarios")
  colors <- c("#20255c", "#4575b4", "#74add1", "#fdae61", "#d73027", "#a50026")
  for (ocean in c("Pacific", "RedSea")) {
    regions <- unique(as.character(d$Region[as.character(d$Ocean) == ocean]))
    limits <- if (ocean == "Pacific") c(1.5, 2, 3, 4, 5, 6) else c(7, 8, 9, 10, 11, 12)
    x <- prob[prob$scenario == "reference" & prob$within_region_depth_range &
      as.character(prob$Region) %in% regions & prob$ratio_threshold %in% limits, ]
    if (!nrow(x)) next
    x$threshold <- factor(paste0(">", x$ratio_threshold, "x"), levels = paste0(">", limits, "x"))
    for (type in c("observation", "expected_profile")) {
      field <- paste0("probability_", type)
      p <- ggplot(x, aes(depth, .data[[field]], color = threshold, fill = threshold)) +
        geom_ribbon(aes(ymin = 0, ymax = .data[[field]], group = threshold), alpha = 0.95, color = NA) +
        geom_line(color = "white", linewidth = 0.15) + scale_x_reverse() +
        scale_y_continuous(limits = c(0, 1)) + coord_flip() + facet_wrap(~ Region, ncol = 3) +
        scale_color_manual(values = colors, name = "Exceeds surface chlorophyll by") +
        scale_fill_manual(values = colors, name = "Exceeds surface chlorophyll by") +
        guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1)) +
        labs(x = "Depth (m)", y = "Exceedance probability", subtitle = paste("Prediction type:", type, "; group effects set to zero")) +
        theme_minimal() + theme(legend.position = "bottom")
      save_plot(p, paste0("probability_", ocean, "_", type))
    }
  }
  # Equal REGION weights within each ocean, with uncertainty combined DRAW BY DRAW.
  # This is a standardized descriptive average of regional reference profiles,
  # not an estimated fixed Ocean coefficient or an area/cast-weighted ocean mean.
  ocean_profiles <- list()
  for (ocean in unique(as.character(d$Ocean))) {
    regions <- unique(as.character(d$Region[as.character(d$Ocean) == ocean]))
    keep <- Reduce(`&`, lapply(reference_grids[regions], function(x) x$within_region_depth_range))
    avg <- Reduce(`+`, reference_draws[regions]) / length(regions)
    ocean_profiles[[ocean]] <- data.frame(Ocean = ocean, depth = seq(10, 250)[keep], summarize_matrix(avg[, keep, drop = FALSE]))
  }
  ocean <- do.call(rbind, ocean_profiles)
  csv(ocean, "equal_region_ocean_profiles.csv")
  save_plot(ggplot(ocean, aes(depth, mean + 1, color = Ocean, fill = Ocean)) +
    geom_ribbon(aes(ymin = lower + 1, ymax = upper + 1), alpha = 0.2, color = NA) + geom_line() +
    scale_x_reverse() + coord_flip() + theme_minimal() +
    labs(x = "Depth (m)", y = "Chlorophyll / surface baseline", subtitle = "Equal region weights at regional reference covariates"), "ocean_profiles")
  capture.output(sessionInfo(), file = file.path(out_dir, "plot_session_info.txt"))
  message("Checks and figures saved to ", out_dir)
  invisible(list(diagnostics = diag, profiles = profile, probabilities = prob))
}

# In RStudio: Sys.setenv(RUN_DIR='...'); source('03_check_and_plot.R'); check_and_plot().
if (sys.nframe() == 0L) check_and_plot()
