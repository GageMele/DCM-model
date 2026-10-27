# The three requested equations, written explicitly without rescaling predictors.
load_equation_spec <- function(equation, original_path) {
  original <- load_original_spec(original_path)
  if (equation == "additive") {
    andrew_form <- brms::bf(chl_norm ~
      1 + Region + log10_distance + MLD_temp_MIX +
      s(depth, by = Region, bs = "cr", k = 15) +
      (1 | Ocean:Region:date) + (1 | Ocean:Region:enso),
      family = gaussian())
  } else if (equation == "mld_varying") {
    andrew_form <- brms::bf(chl_norm ~
      1 + Region + log10_distance +
      s(depth, by = MLD_temp_MIX, bs = "cr", k = 15) +
      s(depth, by = Region, bs = "cr", k = 15) +
      (1 | Ocean:Region:date) + (1 | Ocean:Region:enso),
      family = gaussian())
  } else if (equation == "mld_only") {
    andrew_form <- brms::bf(chl_norm ~
      1 + Region + log10_distance +
      s(depth, by = MLD_temp_MIX, bs = "cr", k = 15) +
      (1 | Ocean:Region:date) + (1 | Ocean:Region:enso),
      family = gaussian())
  } else stop("Unknown equation: ", equation)
  list(formula = andrew_form, current_priors = original$current_priors)
}

check_equation_support <- function(d) {
  # Region must be an unordered factor; MLD must be a numeric multiplier.
  if (!is.factor(d$Region) || is.ordered(d$Region))
    stop("Region must be an unordered factor.")
  if (!is.numeric(d$MLD_temp_MIX)) stop("MLD_temp_MIX must remain numeric.")
  support <- support_table(d)
  if (any(support$depth < 15L)) stop("A region has fewer than 15 distinct depths.")
  if (length(unique(d$MLD_temp_MIX)) < 2L ||
      length(unique(d$log10_distance)) < 2L)
    stop("Distance and MLD must vary in the frozen data.")
  invisible(support)
}
