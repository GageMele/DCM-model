# Exact historical formula and explicitly supplied custom priors.
# No MLD term, tensor interaction, rescaling, smoothing or binning is added.
baseline_spec <- function(prior_mode) {
  formula <- brms::bf(chl_norm ~
    1 + Region + log10_distance +
    s(depth, by = Region, bs = "cr", k = 15) +
    (1 | Ocean:Region:date) + (1 | Ocean:Region:enso),
    family = gaussian())
  custom <- c(
    brms::set_prior("student_t(3, 0.73, 5.448)", class = "Intercept"),
    brms::set_prior("normal(0, 2)", class = "b"),
    brms::set_prior("student_t(3, 0, 1.5)", class = "sd"),
    brms::set_prior("student_t(3, 0, 1.5)", class = "sds"),
    brms::set_prior("student_t(3, 0, 2.5)", class = "sigma")
  )
  if (!prior_mode %in% c("custom", "default")) stop("Unknown prior mode.")
  # NULL invokes the installed brms version's defaults, resolved for this data.
  list(formula = formula, priors = if (prior_mode == "custom") custom else NULL)
}

validate_baseline_cases <- function(cases) {
  stopifnot(nrow(cases) == 12L, identical(cases$case_id, 1:12),
    all(cases$chains == 4L), all(cases$cores == 4L),
    all(cases$threads == 2L), all(cases$cpus == 8L),
    all(cases$iter - cases$warmup == 4000L), all(cases$refresh == 50L),
    all(cases$k_depth == 15L), all(cases$mode == "posterior"),
    all(cases$equation == "original_region_depth"),
    identical(cases$prior, rep(c("custom", "default"), 6L)),
    identical(cases$subset, c(rep(c("p10", "p20", "p50"), each = 2L), rep("p10", 6L))),
    all(cases$memory_gb %in% c(32L, 64L)))
  for (i in seq(1L, 11L, by = 2L)) {
    fields <- setdiff(names(cases), c("case_id", "label", "prior"))
    stopifnot(all(vapply(cases[i:(i + 1L), fields],
      function(x) length(unique(x)) == 1L, logical(1))))
  }
  baseline <- cases[1:6, ]
  stopifnot(all(baseline$iter == 5000L), all(baseline$warmup == 1000L),
    all(baseline$adapt_delta == .97), all(baseline$treedepth == 14L),
    all(baseline$variant == "baseline"))
  # Each extra pair changes exactly one sampler choice; retained draws stay fixed.
  for (i in 7:12) {
    reference <- cases[if (cases$prior[i] == "custom") 1L else 2L, ]
    fields <- c("chains", "cores", "threads", "refresh", "k_depth", "subset", "prior")
    stopifnot(all(vapply(fields, function(f) identical(cases[[f]][i], reference[[f]]), logical(1))))
    if (i %in% 7:8) stopifnot(cases$treedepth[i] == 12L,
      cases$adapt_delta[i] == .97, cases$warmup[i] == 1000L)
    if (i %in% 9:10) stopifnot(cases$treedepth[i] == 14L,
      cases$adapt_delta[i] == .99, cases$warmup[i] == 1000L)
    if (i %in% 11:12) stopifnot(cases$treedepth[i] == 14L,
      cases$adapt_delta[i] == .97, cases$warmup[i] == 2500L)
  }
  invisible(TRUE)
}

check_baseline_support <- function(d) {
  if (!is.factor(d$Region) || is.ordered(d$Region))
    stop("Region must be an unordered factor.")
  tab <- support_table(d)
  if (any(tab$depth < 15L)) stop("A region has fewer than 15 distinct depths.")
  if (length(unique(d$log10_distance)) < 2L) stop("Distance must vary.")
  invisible(tab)
}
