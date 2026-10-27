#!/usr/bin/env Rscript
# Build all exact formulas on the frozen input before allowing any fit to run.
experiment <- Sys.getenv("EXPERIMENT_DIR")
source(file.path(experiment, "common.R"))
source(file.path(experiment, "equation.R"))
for (p in c("brms", "cmdstanr", "mgcv", "posterior"))
  if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)
cases <- read.delim(file.path(experiment, "cases.tsv"), stringsAsFactors = FALSE)
stopifnot(nrow(cases) == 18L, !anyDuplicated(cases$case_id),
  all(cases$chains %in% c(2L, 4L)), all(cases$threads == 4L),
  all(cases$cpus == cases$chains * cases$threads), all(cases$k_depth == 15L),
  all(cases$iter > cases$warmup), all(cases$subset == "p10"),
  all(cases$prior == "current"), all(cases$mode == "posterior"),
  setequal(cases$equation, c("additive", "mld_varying", "mld_only")))
bundle <- readRDS(file.path(experiment, "data", "p10.rds"))
d <- validate_bundle(bundle)
stopifnot(all(d$depth >= 10 & d$depth <= 150),
  isTRUE(all.equal(bundle$config$requested_cast_fraction, .1)))
check_equation_support(d)
cat("Frozen data:", nrow(d), "observations;", length(unique(d$cast_key)), "casts\n")
report <- list()
for (id in unique(cases$equation)) {
  spec <- load_equation_spec(id, file.path(experiment, "02_fit_model.R"))
  standata <- brms::make_standata(spec$formula, data = d, prior = spec$current_priors)
  matrices <- standata[grepl("^Zs", names(standata))]
  spline_columns <- sum(vapply(matrices, function(x)
    if (is.matrix(x)) ncol(x) else 0L, integer(1)))
  code <- brms::make_stancode(spec$formula, data = d, prior = spec$current_priors,
    threads = brms::threading(4))
  writeLines(code, file.path(experiment, paste0("validated_", id, ".stan")))
  capture.output(spec$formula, spec$current_priors,
    file = file.path(experiment, paste0("validated_", id, "_formula_priors.txt")))
  report[[id]] <- data.frame(equation = id, k_depth = 15L,
    observations = nrow(d), casts = length(unique(d$cast_key)),
    stan_fixed_effect_columns = standata$K,
    stan_penalized_spline_columns = spline_columns)
  cat("Validated", id, ":", spline_columns, "penalized spline columns\n")
}
write.csv(do.call(rbind, report), file.path(experiment, "equation_validation.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(experiment, "validation_session_info.txt"))
writeLines("All three exact equations and the frozen data validated.",
  file.path(experiment, "PREPARATION_COMPLETE.txt"))
