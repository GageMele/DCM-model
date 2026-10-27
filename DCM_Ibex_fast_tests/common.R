# Shared helpers for the isolated exploratory experiments.
read_case <- function(experiment, id) {
  cases <- read.delim(file.path(experiment, "cases.tsv"), stringsAsFactors = FALSE)
  hit <- cases[cases$case_id == as.integer(id), , drop = FALSE]
  if (nrow(hit) != 1L) stop("Unknown or duplicated CASE_ID: ", id)
  hit
}

validate_bundle <- function(bundle, full = FALSE) {
  if (!identical(bundle$config$smoothing, "none") ||
      !identical(bundle$config$binned, FALSE)) stop("Unsmoothed, unbinned data required.")
  if (full && !isTRUE(all.equal(bundle$config$cast_fraction, 1)))
    stop("FULL_MODEL_DATA must be the full prepared bundle, not the existing pilot.")
  d <- bundle$model_data
  needed <- c("cast_key", "chl_norm", "Region", "Ocean", "date", "enso",
              "depth", "log10_distance", "MLD_temp_MIX")
  if (!all(needed %in% names(d)) || !nrow(d) || anyNA(d[needed]))
    stop("Invalid prepared bundle.")
  if (any(!is.finite(as.matrix(d[c("chl_norm", "depth", "log10_distance", "MLD_temp_MIX")]))))
    stop("Nonfinite model values.")
  invisible(d)
}

support_table <- function(d) {
  do.call(rbind, lapply(split(d, droplevels(factor(d$Region))), function(x) {
    data.frame(Region = as.character(x$Region[1]), observations = nrow(x),
      casts = length(unique(x$cast_key)), depth = length(unique(x$depth)),
      distance = length(unique(x$log10_distance)), mld = length(unique(x$MLD_temp_MIX)))
  }))
}

check_support <- function(d, k = c(15L, 5L, 5L)) {
  tab <- support_table(d)
  if (any(as.matrix(tab[c("depth", "distance", "mld")]) <
          matrix(k, nrow(tab), 3L, byrow = TRUE))) {
    print(tab)
    stop("Insufficient predictor support; the spline dimensions were not changed.")
  }
  tab
}

weighted_quantile <- function(x, w, probabilities) {
  stopifnot(length(x) == length(w), all(is.finite(x)), all(w > 0))
  ord <- order(x)
  cumulative <- cumsum(w[ord]) / sum(w)
  vapply(probabilities, function(p) x[ord[which(cumulative >= p)[1L]]], numeric(1))
}

full_response_statistics <- function(d) {
  # Give every cast equal total weight, despite unequal profile lengths.
  index <- match(d$cast_key, unique(d$cast_key))
  weights <- 1 / tabulate(index)[index]
  y <- d$chl_norm
  q <- weighted_quantile(y, weights, c(.025, .25, .5, .75, .975))
  weighted_mad <- 1.4826 * weighted_quantile(abs(y - q[3]), weights, .5)
  scale <- max(.5, weighted_mad, (q[4] - q[2]) / 1.349)
  data.frame(observations = nrow(d), casts = length(unique(index)),
    observation_mean = mean(y), observation_sd = sd(y),
    observation_median = median(y), observation_mad = mad(y),
    cast_balanced_q025 = q[1], cast_balanced_q25 = q[2],
    cast_balanced_median = q[3], cast_balanced_q75 = q[4],
    cast_balanced_q975 = q[5], cast_balanced_mad = weighted_mad,
    empirical_location = q[3], empirical_scale = scale,
    minimum = min(y), maximum = max(y))
}

allocate_cast_counts <- function(region_sizes, lower, fraction) {
  # Constrained proportional allocation: preserve support and nested subsets.
  target <- max(as.integer(round(sum(region_sizes) * fraction)), sum(lower))
  counts <- lower
  quota <- fraction * region_sizes
  while (sum(counts) < target) {
    deficit <- quota - counts
    deficit[counts >= region_sizes] <- -Inf
    winner <- which.max(deficit)
    if (!is.finite(deficit[winner])) stop("Cannot allocate requested cast count.")
    counts[winner] <- counts[winner] + 1L
  }
  counts
}

supporting_cast_order <- function(metadata, seed, needed = c(5L, 5L)) {
  # Select real distance/MLD pairs; never manufacture predictor combinations.
  set.seed(seed)
  pool <- sample.int(nrow(metadata))
  chosen <- integer()
  distances <- mlds <- numeric()
  while (length(unique(distances)) < needed[1] || length(unique(mlds)) < needed[2]) {
    score <- (length(unique(distances)) < needed[1]) *
      as.integer(!(metadata$log10_distance[pool] %in% distances)) +
      (length(unique(mlds)) < needed[2]) *
      as.integer(!(metadata$MLD_temp_MIX[pool] %in% mlds))
    if (!length(score) || max(score) == 0) stop("Full region cannot support its basis.")
    pick <- which.max(score)
    selected <- pool[pick]
    chosen <- c(chosen, selected)
    distances <- c(distances, metadata$log10_distance[selected])
    mlds <- c(mlds, metadata$MLD_temp_MIX[selected])
    pool <- pool[-pick]
  }
  list(ids = metadata$cast_key[c(chosen, pool)], minimum = length(chosen))
}

load_original_spec <- function(path, k = c(15L, 5L, 5L)) {
  # Evaluate only the formula and prior assignments from the frozen fit script.
  # This retains the current scientific model without calling fit_model().
  original <- new.env(parent = asNamespace("brms"))
  sys.source(path, envir = original)
  if (!is.function(original$fit_model)) stop("Missing original fit_model().")
  expressions <- as.list(body(original$fit_model))[-1L]
  assignment <- function(name) {
    hits <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
      identical(x[[2]], as.name(name)), expressions)
    if (length(hits) != 1L) stop("Cannot locate original assignment: ", name)
    hits[[1L]]
  }
  spec <- new.env(parent = asNamespace("brms"))
  spec$k <- k
  eval(assignment("andrew_form"), envir = spec)
  eval(assignment("priors"), envir = spec)
  list(formula = spec$andrew_form, current_priors = spec$priors)
}

selected_priors <- function(mode, original, statistics) {
  if (mode == "current") return(original)
  if (mode == "default") return(NULL)  # Let brms determine its defaults.
  if (mode != "empirical") stop("Unknown prior mode.")
  m <- statistics$empirical_location
  s <- statistics$empirical_scale
  c(brms::set_prior(sprintf("student_t(3, %.12g, %.12g)", m, 2 * s), class = "Intercept"),
    brms::set_prior(sprintf("normal(0, %.12g)", s), class = "b"),
    brms::set_prior(sprintf("student_t(3, 0, %.12g)", s), class = "sd"),
    brms::set_prior(sprintf("student_t(3, 0, %.12g)", s), class = "sds"),
    brms::set_prior(sprintf("student_t(3, 0, %.12g)", s), class = "sigma"))
}

finite_min <- function(x) if (length(x) && all(is.finite(x))) min(x) else NA_real_
finite_max <- function(x) if (length(x) && all(is.finite(x))) max(x) else NA_real_

read_sampler_csv <- function(path, iter, warmup) {
  # Read just the seven sampler columns, including saved warmup diagnostics.
  connection <- file(path, "r")
  on.exit(close(connection))
  repeat {
    line <- readLines(connection, n = 1L, warn = FALSE)
    if (!length(line)) stop("CSV header absent: ", path)
    if (nzchar(line) && !startsWith(line, "#")) break
  }
  columns <- strsplit(line, ",", fixed = TRUE)[[1L]]
  keep <- c("lp__", "accept_stat__", "stepsize__", "treedepth__",
            "n_leapfrog__", "divergent__", "energy__")
  classes <- ifelse(columns %in% keep, "numeric", "NULL")
  data <- read.csv(path, comment.char = "#", colClasses = classes,
                   check.names = FALSE)
  if (!all(keep %in% names(data)) || nrow(data) != iter)
    stop("Unexpected sampler CSV rows; save_warmup=TRUE and thin=1 required.")
  data$phase <- c(rep("warmup", warmup), rep("sampling", iter - warmup))
  data
}

chain_elapsed <- function(path) {
  # Reading CSV lines does not parse or retain the large parameter table.
  lines <- readLines(path, warn = FALSE)
  seconds <- function(label) {
    line <- grep(label, lines, fixed = TRUE, value = TRUE)
    if (!length(line)) return(NA_real_)
    token <- regmatches(tail(line, 1L), regexpr(
      "[0-9]+(\\.[0-9]+)?([eE][+-]?[0-9]+)?", tail(line, 1L), perl = TRUE))
    as.numeric(token)
  }
  data.frame(csv = basename(path), warmup_seconds = seconds("(Warm-up)"),
             sampling_seconds = seconds("(Sampling)"))
}
