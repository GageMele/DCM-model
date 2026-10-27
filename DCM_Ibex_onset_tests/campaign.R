# Six source controls are rescored; only execution="fit" cases start new models.
campaign_cases <- function() {
  rows <- list()
  add <- function(eq, subset="p10", split="primary", block="equation", reference=1L,
                  execution="fit", source_case=NA_integer_, seed=1234L,
                  delta=.95, warmup=1000L, samples=2000L, kd=7L, mode="posterior") {
    id <- length(rows)+1L
    threads <- if(execution=="reuse" || mode=="prior")2L else 4L
    rows[[id]] <<- data.frame(case_id=id,
      label=sprintf("o%03d_%s_%s_%s_%s",id,block,eq,split,subset), block,
      equation=eq,subset,split,prior="custom_z",mode,variance="region",cast="none",
      family="gaussian",autocor="none",depth_max=150L,chains=4L,
      iter=warmup+samples,warmup,samples,adapt_delta=delta,treedepth=14L,
      threads,cpus=4L*threads,k_depth=kd,k_mld=5L,k_distance=5L,
      sigma_smooth_sd=.5,cores=4L,refresh=50L,seed,
      memory_gb=if(subset=="p20")64L else 32L,
      reference_case=reference,execution,source_case)
    id
  }
  panels <- list(c("p10","primary"),c("p20","primary"),c("p10","date_blocked"))
  src <- c(1L,5L,13L,15L,22L,24L)
  for(j in seq_along(panels))for(eq in c("original_z","interaction_z")) {
    i <- length(rows)+1L
    add(eq,panels[[j]][1],panels[[j]][2],"existing_control",
        reference=if(i%%2L==1L)NA_integer_ else i-1L,
        execution="reuse",source_case=src[i])
  }
  equations <- c("mld_vary_z","mld_shift","mld_halfshift","mld_ratio",
                 "mld_hybrid","mld_shift_distance")
  for(j in seq_along(panels))for(eq in equations)
    add(eq,panels[[j]][1],panels[[j]][2],reference=2L*j-1L)
  # Isolate sampling stability, robust tuning and basis size for two new anchors.
  for(eq in c("mld_shift","mld_hybrid")) {
    ref <- if(eq=="mld_shift")8L else 11L
    add(eq,block="repeat_seed",reference=ref,seed=20261019L)
    add(eq,block="robust_sampler",reference=ref,delta=.99,warmup=2000L)
    add(eq,block="depth_basis",reference=ref,kd=10L)
  }
  for(eq in c("mld_shift","mld_hybrid"))for(seed in c(1234L,20261019L))
    add(eq,block="prior_predictive",reference=NA_integer_,seed=seed,
        mode="prior",warmup=500L,samples=500L)
  x <- do.call(rbind,rows)
  stopifnot(nrow(x)==34L,sum(x$execution=="fit" & x$mode=="posterior")==24L,
    sum(x$mode=="prior")==4L,sum(x$execution=="reuse")==6L,
    !anyDuplicated(x$label),all(x$subset %in% c("p10","p20")),
    all(x$cpus==x$chains*x$threads))
  x
}
.base_models <- new.env(parent=environment())
base_models <- function() {
  if(!exists("campaign_spec",envir=.base_models,inherits=FALSE))
    sys.source(file.path(Sys.getenv("EXPERIMENT_DIR",Sys.getenv("ONSET_PACKAGE")),
                         "base_models.R"),envir=.base_models)
  .base_models
}
apply_features <- function(d,scaling,region_levels) {
  d <- base_models()$apply_features(d,scaling,region_levels)
  if(any(!is.finite(d$MLD_temp_MIX)) || any(d$MLD_temp_MIX<=0))
    stop("MLD must be finite and positive for the predeclared depth/MLD tests; no rows dropped")
  # All transformations use predictors only and the observed MLD in metres.
  d$depth_mld <- d$depth-d$MLD_temp_MIX
  d$depth_mld_half <- d$depth-.5*d$MLD_temp_MIX
  d$depth_mld_ratio <- d$depth/d$MLD_temp_MIX
  d
}
campaign_spec <- function(case,region_levels) {
  if(case$equation %in% c("original_z","interaction_z","mld_vary_z"))
    return(base_models()$campaign_spec(case,region_levels))
  # Obtain exactly the same coefficient, smooth and random-effect priors.
  surrogate <- case;surrogate$equation <- "original_z"
  spec <- base_models()$campaign_spec(surrogate,region_levels)
  kd <- case$k_depth;ke <- case$k_mld
  region_s <- function(x)sprintf("s(%s,by=Region,bs='cr',k=%d)",x,kd)
  terms <- switch(case$equation,
    mld_shift=c("distance_z",region_s("depth_mld")),
    mld_halfshift=c("distance_z",region_s("depth_mld_half")),
    mld_ratio=c("distance_z",region_s("depth_mld_ratio")),
    mld_hybrid=c("distance_z",region_s("depth"),
      sprintf("s(depth_mld,bs='cr',k=%d)",ke)),
    mld_shift_distance=c(region_s("depth_mld"),
      sprintf("s(depth_mld,by=distance_z,bs='cr',k=%d)",ke)),
    stop("Unknown onset equation"))
  f <- as.formula(paste("chl_norm ~ 1 + Region +",paste(c(terms,
      "(1|Ocean:Region:date)","(1|Ocean:Region:enso)"),collapse=" + ")),
      env=asNamespace("brms"))
  spec$formula <- brms::bf(f,sigma~0+Region,family=gaussian())
  spec
}
data_file <- function(experiment,case)file.path(experiment,"data",
  paste0(case$split,"_",case$subset,"_d150.rds"))
training_eval_file <- function(experiment,case)file.path(experiment,"data",
  paste0(case$split,"_p10_d150.rds"))
validation_file <- function(experiment,case)file.path(experiment,"data",
  paste0(case$split,"_validation_d150.rds"))
reserve_validation_file <- function(experiment,case)file.path(experiment,"data",
  paste0(case$split,"_reserve_d150.rds"))
thread_specification <- function(case)brms::threading(case$threads)
