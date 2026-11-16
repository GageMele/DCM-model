# Fit one frozen case. An accepted fitting job is never automatically repeated.
experiment<-Sys.getenv("EXPERIMENT_DIR")
source(file.path(experiment,"common.R"));source(file.path(experiment,"models.R"))
if(!file.exists(file.path(experiment,"PREPARATION_COMPLETE.txt")))stop("Preparation is incomplete")
case<-read_case(experiment,Sys.getenv("CASE_ID"));run<-run_directory(experiment,case)
if(dir.exists(run))stop("Run directory already exists; inspect it rather than overwriting a fit")
dir.create(run,recursive=TRUE)
if(case$mode=="reuse") {
  # Re-score an existing shape150 fit with this campaign's scoring (adds held-out
  # log predictive density). Nothing is refitted; the source run is only linked.
  tryCatch({
    source_exp<-readLines(file.path(experiment,"reuse_source.txt"))[1]
    src<-file.path(source_exp,"runs",case$source_run)
    if(!file.exists(file.path(src,"FIT_COMPLETE.txt")))stop("Source fit not complete: ",src)
    for(x in c("fit.rds","stan_csv","model_data.rds"))
      if(!file.symlink(file.path(src,x),file.path(run,x)))stop("Cannot link ",x)
    for(x in c("run_config.rds","run_config.txt","hardware.txt","session_info.txt","requested_priors.txt","resolved_priors.csv"))
      if(file.exists(file.path(src,x)))file.copy(file.path(src,x),file.path(run,x))
    writeLines(c(paste("Reused source fit",src),"Fit serialized with draws; diagnostics and predictions follow separately"),
      file.path(run,"FIT_COMPLETE.txt"))
  },error=function(e){writeLines(conditionMessage(e),file.path(run,"FIT_FAILED.txt"));stop(e)})
  quit(status=0L)
}
tryCatch({
  start<-proc.time()[["elapsed"]]
  for(p in c("brms","cmdstanr","mgcv","posterior"))if(!requireNamespace(p,quietly=TRUE))stop("Missing ",p)
  if(nzchar(Sys.getenv("CMDSTAN")))cmdstanr::set_cmdstan_path(Sys.getenv("CMDSTAN"))
  data_path<-data_file(experiment,case);d<-prepare_data(augment(readRDS(data_path)),case)
  spec<-campaign_spec(case,readRDS(file.path(experiment,"region_levels.rds")))
  allocated<-as.integer(Sys.getenv("SLURM_CPUS_PER_TASK",case$cpus))
  if(allocated<case$chains*case$threads)stop("Insufficient allocated CPUs")
  Sys.setenv(STAN_NUM_THREADS=case$threads)
  config<-list(case=case,job_id=Sys.getenv("SLURM_JOB_ID"),allocated_cpus=allocated,
    hostname=Sys.info()[["nodename"]],data_md5=unname(tools::md5sum(data_path)),
    source_md5=unname(tools::md5sum(file.path(experiment,"models.R"))),
    formula=paste(capture.output(spec$formula),collapse=" "),cmdstan_version=cmdstanr::cmdstan_version())
  saveRDS(config,file.path(run,"run_config.rds"));capture.output(config,file=file.path(run,"run_config.txt"))
  saveRDS(list(model_data=d),file.path(run,"model_data.rds"))
  saveRDS(spec$priors,file.path(run,"requested_priors.rds"))
  capture.output(spec$priors,file=file.path(run,"requested_priors.txt"))
  writeLines(system2("lscpu",stdout=TRUE),file.path(run,"hardware.txt"))
  capture.output(sessionInfo(),file=file.path(run,"session_info.txt"))
  csv_dir<-file.path(run,"stan_csv");dir.create(csv_dir)
  writeLines(brms::make_stancode(spec$formula,data=d,prior=spec$priors,
    threads=thread_specification(case),sample_prior=if(case$mode=="prior")"only"else"no"),file.path(run,"model.stan"))
  cat("Starting",case$label,"on",length(unique(d$cast_key)),"casts /",nrow(d),"rows\n")
  sampling_start<-proc.time()[["elapsed"]]
  fit<-brms::brm(spec$formula,data=d,prior=spec$priors,backend="cmdstanr",
    chains=case$chains,cores=case$cores,threads=thread_specification(case),
    iter=case$iter,warmup=case$warmup,seed=case$seed,refresh=case$refresh,
    save_warmup=TRUE,thin=1L,init=if(case$equation=="peak_nl")0 else NULL,
    control=list(adapt_delta=case$adapt_delta,max_treedepth=case$treedepth),
    save_pars=brms::save_pars(all=TRUE),sample_prior=if(case$mode=="prior")"only"else"no",
    output_dir=csv_dir,file=file.path(run,"fit"),file_refit="always")
  config$brm_seconds_including_compile<-proc.time()[["elapsed"]]-sampling_start
  invisible(posterior::as_draws_array(fit));saveRDS(fit,file.path(run,"fit.rds"))
  config$fit_seconds<-proc.time()[["elapsed"]]-start
  saveRDS(config,file.path(run,"run_config.rds"))
  write.csv(as.data.frame(brms::prior_summary(fit)),file.path(run,"resolved_priors.csv"),row.names=FALSE)
  writeLines("Fit serialized with draws; diagnostics and predictions follow separately",file.path(run,"FIT_COMPLETE.txt"))
},error=function(e){writeLines(conditionMessage(e),file.path(run,"FIT_FAILED.txt"));stop(e)})
