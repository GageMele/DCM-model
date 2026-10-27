#!/usr/bin/env Rscript
# Validate the actual frozen input and all bases before releasing any fit.
experiment <- Sys.getenv("EXPERIMENT_DIR")
source(file.path(experiment, "common.R"))
for (p in c("brms", "cmdstanr", "mgcv", "posterior"))
  if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)
cases <- read.delim(file.path(experiment, "cases.tsv"), stringsAsFactors = FALSE)
stopifnot(nrow(cases) == 12L, !anyDuplicated(cases$case_id),
  all(cases$chains == 2L), all(cases$threads == 4L),
  all(cases$cpus == cases$chains * cases$threads),
  all(cases$iter > cases$warmup), all(cases$subset == "p10"),
  all(cases$prior == "current"), all(cases$mode == "posterior"))
bundle <- readRDS(file.path(experiment, "data", "p10.rds"))
d <- validate_bundle(bundle)
stopifnot(all(d$depth >= 10 & d$depth <= 150),
  isTRUE(all.equal(bundle$config$requested_cast_fraction, .1)))
cat("Frozen data:", nrow(d), "observations;", length(unique(d$cast_key)), "casts\n")
report <- list()
unique_k <- unique(cases[c("k_depth", "k_distance", "k_mld")])
for (i in seq_len(nrow(unique_k))) {
  k <- case_k(unique_k[i, , drop = FALSE])
  check_support(d, k)
  spec <- load_original_spec(file.path(experiment, "02_fit_model.R"), k)
  # Construct the full brms/Stan data: catches unsupported bases or priors.
  standata <- brms::make_standata(spec$formula, data = d, prior = spec$current_priors)
  matrices <- standata[grepl("^Zs", names(standata))]
  spline_columns <- sum(vapply(matrices, function(x)
    if (is.matrix(x)) ncol(x) else 0L, integer(1)))
  report[[i]] <- data.frame(unique_k[i, , drop = FALSE],
    nominal_basis_product_per_region = prod(k),
    stan_fixed_effect_columns = standata$K,
    stan_penalized_spline_columns = spline_columns,
    observations = nrow(d), casts = length(unique(d$cast_key)))
  cat("Validated k =", paste(k, collapse = ","), "\n")
  rm(standata, matrices); gc()
}
write.csv(do.call(rbind, report), file.path(experiment, "basis_validation.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(experiment, "validation_session_info.txt"))
writeLines("All frozen data and proposed bases validated.",
  file.path(experiment, "PREPARATION_COMPLETE.txt"))
