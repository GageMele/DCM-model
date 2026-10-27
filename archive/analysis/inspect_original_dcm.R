# Read saved objects without executing the original analysis or refitting.
out <- "original_model_review"
dir.create(out, showWarnings = FALSE)
fit <- readRDS("original_model_input/my_model_50_FINAL.rds")
cat("SAVED FIT CLASS:", class(fit), "\nNAMES:", names(fit), "\n")
cat("FORMULA:\n"); print(fit$formula)
cat("ENGINE CLASS:", class(fit$fit), "TYPE:", typeof(fit$fit), "\n")
cat("ENGINE ATTRIBUTES:", names(attributes(fit$fit)), "\n")
cat("BACKEND:", fit$backend, "\n")
d <- fit$data
cat("DATA DIM:", dim(d), "\nCOLUMNS:", names(d), "\n")
for (v in intersect(c("depth", "chl_norm", "log10_distance", "MLD_temp_MIX"), names(d))) {
  cat("VARIABLE", v, "RANGE", range(d[[v]], na.rm = TRUE), "\n")
  print(quantile(d[[v]], c(0,.01,.05,.5,.95,.99,1), na.rm = TRUE))
}
if ("Region" %in% names(d)) print(table(d$Region))
cat("PRIORS:\n"); print(fit$prior)
cat("VERSION:\n"); print(fit$version)
cat("TIME:\n"); print(fit$time)
if (isS4(fit$fit)) {
  for (nm in intersect(c("sim","stan_args","date"), names(attributes(fit$fit)))) {
    value <- attr(fit$fit,nm)
    if(nm=="sim") {
      cat("SIM FIELDS:",names(value),"\n")
      print(value[intersect(c("iter","warmup","chains","n_save","warmup2","thin"),names(value))])
      cat("SAMPLES CLASS/DIM:",class(value$samples),length(value$samples),"\n")
      cat("FIRST CHAIN ATTRIBUTES:",names(attributes(value$samples[[1]])),"\n")
      cat("FIRST CHAIN PARAMETERS:",head(names(value$samples[[1]]),20),"\n")
    } else if (nm=="stan_args") str(value,max.level=2) else print(value)
  }
}
e <- new.env(parent = emptyenv())
loaded <- load("original_model_input/fit_brms_model_50_FINAL.RData", envir=e)
cat("RDATA OBJECTS:",loaded,"\n")
for (nm in loaded) cat(nm, "CLASS",class(e[[nm]]),"SIZE",format(object.size(e[[nm]]),units="MB"),"\n")
saveRDS(list(formula=fit$formula, prior=fit$prior, version=fit$version,
             data_summary=summary(d), engine_class=class(fit$fit), backend=fit$backend),
        file.path(out,"original_metadata.rds"))
