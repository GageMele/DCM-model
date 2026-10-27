# Focused follow-up to the 18 October all-run comparison.
campaign_cases <- function() {
  rows<-list()
  add<-function(block,eq="both_z",subset="p10",split="primary",seed=1234L,
    variance="region",warmup=1000L,samples=2000L,delta=.95,threads=2L,
    kd=7L,ke=5L,sigma_sd=.5,mode="posterior",reference=NA_integer_) {
    id<-length(rows)+1L
    rows[[id]]<<-data.frame(case_id=id,label=sprintf("f%03d_%s_%s_%s_s%d",id,block,eq,subset,seed),
      block,equation=eq,subset,split,prior="custom_z",mode,variance,cast="none",
      family="gaussian",autocor="none",depth_max=150L,chains=4L,iter=warmup+samples,
      warmup,samples,adapt_delta=delta,treedepth=14L,threads,cpus=4L*threads,
      k_depth=kd,k_mld=ke,k_distance=ke,sigma_smooth_sd=sigma_sd,
      cores=4L,refresh=50L,seed,memory_gb=if(subset=="p20" || variance=="region_depth")64L else 32L,
      reference_case=reference)
    id
  }
  equations<-c("original_z","both_z","interaction_z")
  # Equal settings across equations, plus an independent sampling seed.
  for(eq in equations)for(seed in c(1234L,20261018L))
    add("stability",eq,seed=seed,reference=if(eq=="original_z")if(seed==1234L)NA_integer_ else 1L else
      if(eq=="both_z")if(seed==1234L)1L else 3L else if(seed==1234L)3L else 5L)
  anchors<-c(original_z=1L,both_z=3L,interaction_z=5L)
  for(eq in equations)add("delta",eq,delta=.90,reference=anchors[[eq]])
  for(eq in equations)add("threads",eq,threads=4L,reference=anchors[[eq]])
  for(eq in equations)add("training_size",eq,subset="p20",reference=anchors[[eq]])
  # Matched robust-sampler controls isolate variance structure from tuning.
  for(eq in c("both_z","interaction_z")) {
    control<-add("sigma_control",eq,warmup=2000L,delta=.99,reference=anchors[[eq]])
    smooth<-add("sigma_depth",eq,variance="region_depth",warmup=2000L,delta=.99,reference=control)
    add("sigma_prior",eq,variance="region_depth",warmup=2000L,delta=.99,sigma_sd=.25,reference=smooth)
  }
  for(eq in equations)add("date_blocked",eq,split="date_blocked",reference=anchors[[eq]])
  # Alternative 108 training casts: the other half of the same 20% pool.
  for(eq in equations)add("training_composition",eq,subset="p10_alt",reference=anchors[[eq]])
  add("depth_basis","interaction_z",kd=10L,reference=5L)
  add("effect_basis_low","interaction_z",ke=3L,reference=5L)
  add("effect_basis_high","interaction_z",ke=7L,reference=5L)
  for(eq in c("both_z","interaction_z"))for(s in c(.5,.25))
    add("prior_sigma",eq,variance="region_depth",sigma_sd=s,mode="prior",warmup=500L,samples=500L)
  x<-do.call(rbind,rows)
  stopifnot(nrow(x)==34L,sum(x$mode=="posterior")==30L,!anyDuplicated(x$label),
    all(x$subset %in% c("p10","p10_alt","p20")),all(x$cpus==x$chains*x$threads))
  x
}

# Reuse the previously run mean equations and feature definitions exactly.
.base_models<-new.env(parent=environment())
base_models<-function() {
  if(!exists("campaign_spec",envir=.base_models,inherits=FALSE))
    sys.source(file.path(Sys.getenv("EXPERIMENT_DIR",Sys.getenv("FOLLOWUP_PACKAGE")),"base_models.R"),envir=.base_models)
  .base_models
}
campaign_spec<-function(case,region_levels) {
  spec<-base_models()$campaign_spec(case,region_levels)
  if(case$variance=="region_depth") {
    index<-spec$priors$class=="sds" & spec$priors$dpar=="sigma"
    stopifnot(sum(index)==1L)
    spec$priors$prior[index]<-sprintf("normal(0,%.12g)",case$sigma_smooth_sd)
  }
  spec
}
apply_features<-function(d,scaling,region_levels)base_models()$apply_features(d,scaling,region_levels)
data_file<-function(experiment,case)file.path(experiment,"data",paste0(case$split,"_",case$subset,"_d150.rds"))
training_eval_file<-function(experiment,case)file.path(experiment,"data",paste0(case$split,"_",if(case$subset=="p10_alt")"p10_alt"else"p10","_d150.rds"))
validation_file<-function(experiment,case)file.path(experiment,"data",paste0(case$split,"_validation_d150.rds"))
reserve_validation_file<-function(experiment,case)file.path(experiment,"data",paste0(case$split,"_reserve_d150.rds"))
thread_specification<-function(case)brms::threading(case$threads)
