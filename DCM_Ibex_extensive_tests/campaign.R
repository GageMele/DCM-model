# The experiment matrix and scientific specifications. No model is fitted here.
campaign_cases <- function() {
  rows <- list()
  add <- function(block, equation, subset = "p10", kd = 7L, km = kd, kl = kd,
                  prior = "custom", chains = 4L, warmup = 1000L, samples = 1000L,
                  delta = .95, td = 12L, threads = 2L, seed = 1234L,
                  family = "gaussian", mode = "posterior") {
    id <- length(rows) + 1L
    rows[[id]] <<- data.frame(case_id = id,
      label = sprintf("c%03d_%s_%s_%s", id, block, equation, subset),
      block, equation, subset, prior, mode, chains, iter = warmup + samples,
      warmup, samples, adapt_delta = delta, treedepth = td, threads,
      cpus = chains * threads, k_depth = kd, k_mld = km, k_distance = kl,
      cores = chains, refresh = 50L, seed, family,
      memory_gb = if (subset == "p50") 64L else 32L)
  }
  equations <- c("original", "additive", "mld_vary", "distance_vary", "both_vary",
                 "global_both", "region_tensor_mld", "region_tensor_distance",
                 "global_tensor_mld", "global_tensor_distance")
  for (eq in equations) for (k in c(5L, 7L, 10L, 15L))
    add("equation_k", eq, kd = k,
        km = if (grepl("tensor", eq)) 3L else k,
        kl = if (grepl("tensor", eq)) 3L else k)
  for (k in list(c(5,3,3), c(7,3,3), c(5,4,3), c(5,3,4)))
    add("small_3d", "region_tensor_3d", kd=k[1], kl=k[2], km=k[3])
  for (k in list(c(10,5,5), c(15,5,5), c(10,10,5), c(10,5,10)))
    add("mixed_k", "both_vary", kd=k[1], km=k[2], kl=k[3])
  add("scaling", "both_scaled")
  anchors <- c("additive", "both_scaled")
  for (eq in anchors) {
    for (delta in c(.90, .97, .99)) add("delta", eq, delta=delta)
    for (td in c(10L, 14L)) add("tree", eq, td=td)
    # Hold retained draws fixed when changing warmup, and vice versa.
    for (warm in c(250L, 500L, 2000L)) add("warmup", eq, warmup=warm)
    for (n in c(500L, 2000L)) add("draws", eq, samples=n)
  }
  for (eq in anchors) {
    add("chains_same_total", eq, chains=2L, samples=2000L)
    add("chains_half_total", eq, chains=2L)
  }
  for (eq in anchors) for (p in c("default", "tighter", "half_normal"))
    add("priors", eq, prior=p)
  for (eq in anchors) for (threads in c(1L, 4L)) add("threads", eq, threads=threads)
  for (eq in c("additive", "mld_vary", "both_scaled", "region_tensor_mld"))
    for (subset in c("p20", "p50")) add("data_size", eq, subset=subset,
      km=if (eq == "region_tensor_mld") 3L else 7L,
      kl=if (eq == "region_tensor_mld") 3L else 7L)
  for (seed in c(20261016L, 20261017L)) add("repeat_seed", "both_scaled", seed=seed)
  for (eq in anchors) add("likelihood", eq, family="student")
  add("repeat_seed", "additive", seed=20261016L)
  for (eq in anchors) for (p in c("custom", "tighter", "half_normal"))
    add("prior_predictive", eq, prior=p, warmup=500L, samples=500L, mode="prior")
  ans <- do.call(rbind, rows)
  stopifnot(nrow(ans)==102L, sum(ans$mode=="posterior")==96L,
            !anyDuplicated(ans$label), all(ans$iter > ans$warmup))
  ans
}

campaign_spec <- function(case) {
  kd <- case$k_depth; km <- case$k_mld; kl <- case$k_distance
  region <- sprintf("s(depth, by=Region, bs='cr', k=%d)", kd)
  smooth <- function(by, k) sprintf("s(depth, by=%s, bs='cr', k=%d)", by, k)
  mld <- smooth("MLD_temp_MIX", km)
  distance <- smooth("log10_distance", kl)
  tensor <- function(x, k, by=FALSE) sprintf(
    "t2(depth, %s, bs=c('cr','cr'), k=c(%d,%d)%s)", x, kd, k,
    if (by) ", by=Region" else "")
  terms <- switch(case$equation,
    original = c("log10_distance", region),
    additive = c("log10_distance", "MLD_temp_MIX", region),
    mld_vary = c("log10_distance", mld, region),
    distance_vary = c("MLD_temp_MIX", distance, region),
    both_vary = c(mld, distance, region),
    both_scaled = c(smooth("mld_z", km), smooth("distance_z", kl), region),
    global_both = c(sprintf("s(depth, bs='cr', k=%d)", kd), mld, distance),
    region_tensor_mld = c("log10_distance", tensor("MLD_temp_MIX", km, TRUE)),
    region_tensor_distance = c("MLD_temp_MIX", tensor("log10_distance", kl, TRUE)),
    global_tensor_mld = c("log10_distance", tensor("MLD_temp_MIX", km)),
    global_tensor_distance = c("MLD_temp_MIX", tensor("log10_distance", kl)),
    region_tensor_3d = sprintf(paste0("t2(depth, log10_distance, MLD_temp_MIX, ",
      "by=Region, bs=c('cr','cr','cr'), k=c(%d,%d,%d))"), kd, kl, km),
    stop("Unknown equation: ", case$equation))
  # Numeric-by smooths include the associated predictor's constant-in-depth
  # effect. Do not also add that linear predictor and duplicate its null space.
  f <- as.formula(paste("chl_norm ~ 1 + Region +", paste(c(terms,
    "(1|Ocean:Region:date)", "(1|Ocean:Region:enso)"), collapse=" + ")))
  environment(f) <- asNamespace("brms")
  family <- if (case$family=="student") brms::student() else stats::gaussian()
  p <- if (case$prior=="default") NULL else {
    b <- if (case$prior=="tighter") "normal(0, 1)" else "normal(0, 2)"
    sd <- switch(case$prior, custom="student_t(3, 0, 1.5)",
                 tighter="normal(0, 0.75)", half_normal="normal(0, 1.5)",
                 stop("Unknown prior"))
    sigma <- if (case$prior=="half_normal") "normal(0, 2.5)" else "student_t(3, 0, 2.5)"
    c(brms::set_prior("student_t(3, 0.73, 5.448)", class="Intercept"),
      brms::set_prior(b, class="b"), brms::set_prior(sd, class="sd"),
      brms::set_prior(sd, class="sds"), brms::set_prior(sigma, class="sigma"))
  }
  list(formula=brms::bf(f, family=family), priors=p)
}

apply_scaling <- function(d, scaling) {
  d$mld_z <- (d$MLD_temp_MIX-scaling$mld_mean)/scaling$mld_sd
  d$distance_z <- (d$log10_distance-scaling$distance_mean)/scaling$distance_sd
  d
}

case_design_key <- function(cases) {
  cols <- cases[c("equation","subset","prior","mode","k_depth","k_mld",
                  "k_distance","family","threads")]
  # as.matrix.data.frame can pad numeric columns differently for a whole matrix
  # and for one row. Explicit character conversion keeps these keys identical.
  do.call(paste,c(lapply(cols,as.character),list(sep="_")))
}
