# Scheduler orchestration is separated from fitting and can be tested offline.
clean_job_ids<-function(x)unique(grep("^[0-9]+$",as.character(x),value=TRUE))

preparation_dependency<-function(experiment,job_id) {
  # Completed preparation is proved by its success marker. Its Slurm record
  # may already have expired, so do not keep adding a dependency on that ID.
  if(file.exists(file.path(experiment,"PREPARATION_COMPLETE.txt")))return(character())
  if(file.exists(file.path(experiment,"PREPARATION_FAILED.txt")))stop("Preparation failed; no new jobs submitted")
  c(paste0("--dependency=afterok:",job_id),"--kill-on-invalid-dep=yes")
}

active_scheduler_jobs<-function(known) {
  response<-suppressWarnings(system2("squeue",c("--me","--noheader","--format=%A"),stdout=TRUE,stderr=TRUE))
  status<-attr(response,"status")
  if(!is.null(status)&&status!=0)stop("Cannot verify live Slurm jobs: ",paste(response,collapse="\n"))
  intersect(clean_job_ids(known),clean_job_ids(trimws(response)))
}

tracked_work_jobs<-function(experiment) {
  ids<-character()
  manifest<-file.path(experiment,"manifest.tsv")
  if(file.exists(manifest)) {
    x<-read.delim(manifest,stringsAsFactors=FALSE)
    ids<-c(ids,x$fit_job,x$check_job)
  }
  evaluations<-file.path(experiment,"existing_evaluations.tsv")
  if(file.exists(evaluations))ids<-c(ids,read.delim(evaluations,stringsAsFactors=FALSE)$job_id)
  clean_job_ids(ids)
}

submit_collection<-function(experiment,root,submit_job,active_lookup=active_scheduler_jobs) {
  # Watch all recorded fits, checks and evaluations, not just newly added checks.
  known<-tracked_work_jobs(experiment)
  record<-file.path(experiment,"collection_job_id.txt")
  watched<-file.path(experiment,"collection_watch_job_ids.txt")
  if(file.exists(record)&&file.exists(watched)&&identical(sort(readLines(watched,warn=FALSE)),sort(known))) {
    previous<-clean_job_ids(readLines(record,warn=FALSE))
    if(length(intersect(previous,active_lookup(previous)))) {
      cat("Existing final collector still active:",previous,"\n")
      return(previous[1])
    }
  }
  Sys.setenv(DCM_ROOT=root,EXPERIMENT_DIR=experiment)
  for(attempt in 1:3) {
    live<-intersect(known,clean_job_ids(active_lookup(known)))
    dependency<-if(length(live))paste0("--dependency=afterany:",paste(live,collapse=":"))else character()
    args<-c("--parsable",paste0("--chdir=",root),"--export=ALL","--job-name=dcm-t150-results",
      dependency,file.path(experiment,"collect.slurm"))
    result<-tryCatch(submit_job(args),error=function(e)e)
    if(!inherits(result,"error")) {
      writeLines(result,record);writeLines(known,watched)
      cat("Final collector",result,"waits for",length(live),"currently active recorded jobs.\n")
      return(result)
    }
    if(!grepl("Job dependency problem|Invalid job id specified",conditionMessage(result),ignore.case=TRUE)||attempt==3L)stop(result)
    # A job can finish and disappear between squeue and sbatch. Refresh only
    # after a definite dependency rejection; never resubmit an accepted job.
  }
}

queue_frozen_experiment<-function(experiment,root,cases,submit_job,legacy=data.frame(),
  constraint="intel&cascadelake",fit_time="7-00:00:00",active_lookup=active_scheduler_jobs) {
  Sys.setenv(DCM_ROOT=root,EXPERIMENT_DIR=experiment)
  base<-c("--parsable",paste0("--chdir=",root),"--export=ALL")
  submit<-function(name,script,extra=character())submit_job(c(base,paste0("--job-name=",name),extra,file.path(experiment,script)))
  prepare<-submit("dcm-t150-prepare","prepare.slurm")
  writeLines(prepare,file.path(experiment,"preparation_job_id.txt"))
  # Preserve a monitorable pointer even if a later scheduler request fails.
  if(dir.exists(file.path(root,"logs")))writeLines(experiment,file.path(root,"logs","latest_targeted150_experiment.txt"))
  manifest<-data.frame(case_id=integer(),label=character(),fit_job=character(),check_job=character(),cpus=integer(),run_dir=character())
  write_manifest<-function()write.table(manifest,file.path(experiment,"manifest.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
  write_manifest();checks<-character()
  for(i in seq_len(nrow(cases))) {
    c<-cases[i,];Sys.setenv(CASE_ID=c$case_id)
    fit<-submit(paste0("dcm-t150-c",c$case_id),"fit.slurm",c(preparation_dependency(experiment,prepare),
      paste0("--cpus-per-task=",c$cpus),paste0("--mem=",c$memory_gb,"G"),paste0("--time=",fit_time),
      paste0("--constraint=",constraint),"--hint=nomultithread"))
    manifest<-rbind(manifest,data.frame(case_id=c$case_id,label=c$label,fit_job=fit,check_job="NA",cpus=c$cpus,
      run_dir=file.path(experiment,"runs",sprintf("case_%02d",c$case_id))))
    write_manifest()
    check<-submit(paste0("dcm-t150-check",c$case_id),"check.slurm",paste0("--dependency=afterany:",fit))
    manifest$check_job[nrow(manifest)]<-check;write_manifest();checks<-c(checks,check)
    cat("Queued case",c$case_id,c$label,"fit",fit,"check",check,"\n")
  }
  evaluations<-data.frame(external_id=character(),job_id=character(),source_run=character(),run_dir=character())
  for(i in seq_len(nrow(legacy))) {
    z<-legacy[i,];evaluation<-file.path(experiment,"evaluations",z$external_id)
    Sys.setenv(SOURCE_RUN=z$source_run,EVAL_DIR=evaluation,EXTERNAL_ID=z$external_id)
    id<-submit(paste0("dcm-t150-old",i),"score_existing.slurm",preparation_dependency(experiment,prepare))
    evaluations<-rbind(evaluations,data.frame(external_id=z$external_id,job_id=id,source_run=z$source_run,run_dir=evaluation))
    write.table(evaluations,file.path(experiment,"existing_evaluations.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
    checks<-c(checks,id)
  }
  if(!nrow(evaluations))write.table(evaluations,file.path(experiment,"existing_evaluations.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
  collector<-submit_collection(experiment,root,submit_job,active_lookup)
  list(prepare=prepare,manifest=manifest,evaluations=evaluations,collector=collector)
}

discover_existing<-function(root) {
  legacy<-list()
  for(prefix in c("baseline150","equation150","k150","fast150","extensive150"))for(m in Sys.glob(file.path(root,"experiments",paste0(prefix,"_*"),"manifest.tsv"))) {
    registry<-read.delim(m,stringsAsFactors=FALSE)
    old_cases<-read.delim(file.path(dirname(m),"cases.tsv"),stringsAsFactors=FALSE)
    for(i in seq_len(nrow(registry))) {
      z<-registry[i,];case<-old_cases[old_cases$case_id==z$case_id,,drop=FALSE]
      if(nrow(case)!=1L || case$mode!="posterior" || !case$subset %in% c("p10","p20","p50"))next
      if(!all(file.exists(file.path(z$run_dir,c("FIT_COMPLETE.txt","fit.rds","model_data.rds","run_config.rds")))))next
      legacy[[length(legacy)+1L]]<-data.frame(external_id=paste0(basename(dirname(m)),"_case_",z$case_id),source_run=z$run_dir)
    }
  }
  if(length(legacy))unique(do.call(rbind,legacy))else data.frame(external_id=character(),source_run=character())
}

scheduler_submit<-function(args) {
  response<-suppressWarnings(system2("sbatch",vapply(args,shQuote,character(1)),stdout=TRUE,stderr=TRUE))
  status<-attr(response,"status");if(!is.null(status)&&status!=0)stop(paste(response,collapse="\n"))
  ids<-grep("^[0-9]+(;[^[:space:]]+)?$",trimws(response),value=TRUE)
  if(length(ids)!=1L)stop("sbatch did not return exactly one job ID: ",paste(response,collapse="\n"))
  sub(";.*$","",ids)
}

submit_campaign<-function(root,package) {
  root<-normalizePath(root,mustWork=TRUE);package<-normalizePath(package,mustWork=TRUE)
  pointer<-file.path(root,"logs","latest_extensive150_experiment.txt")
  source_experiment<-Sys.getenv("SOURCE_EXPERIMENT",if(file.exists(pointer))readLines(pointer,n=1L)else "")
  if(!nzchar(source_experiment))stop("Need the frozen extensive150 source experiment")
  required<-c("PREPARATION_COMPLETE.txt","full_input_250m.rds","ibex_environment.sh",
    file.path("data",c("p10.rds","p20.rds","p50.rds","validation_model.rds")))
  if(!all(file.exists(file.path(source_experiment,required))))stop("Frozen source files are missing")
  source(file.path(package,"campaign.R"),local=TRUE)
  cases<-campaign_cases()
  checked<-read.delim(file.path(package,"cases.tsv"),stringsAsFactors=FALSE)
  stopifnot(isTRUE(all.equal(cases,checked,check.attributes=FALSE)))
  experiment<-file.path(root,"experiments",paste0("targeted150_",format(Sys.time(),"%Y%m%dT%H%M%SZ",tz="UTC"),"_",Sys.getpid()))
  if(dir.exists(experiment))stop("Experiment exists")
  dir.create(file.path(experiment,"data"),recursive=TRUE)
  files<-list.files(package,pattern="\\.(R|slurm|awk|sh|md|tsv)$",full.names=TRUE)
  if(!all(file.copy(files,experiment)))stop("Package copy failed")
  if(dir.exists(file.path(package,"original_review")) &&
    !file.copy(file.path(package,"original_review"),experiment,recursive=TRUE))stop("Original review copy failed")
  file.copy(file.path(source_experiment,"ibex_environment.sh"),file.path(experiment,"ibex_environment.sh"))
  file.copy(file.path(source_experiment,"full_input_250m.rds"),file.path(experiment,"full_input_250m.rds"))
  for(s in c("p10","p20","p50"))file.copy(file.path(source_experiment,"data",paste0(s,".rds")),file.path(experiment,"data",paste0(s,".rds")))
  file.copy(file.path(source_experiment,"data","validation_model.rds"),file.path(experiment,"source_validation.rds"))
  writeLines(source_experiment,file.path(experiment,"source_experiment.txt"))
  source_files<-file.path(source_experiment,c("full_input_250m.rds",file.path("data",c("p10.rds","p20.rds","p50.rds","validation_model.rds"))))
  copied_files<-c(file.path(experiment,"full_input_250m.rds"),file.path(experiment,"data",c("p10.rds","p20.rds","p50.rds")),file.path(experiment,"source_validation.rds"))
  stopifnot(identical(unname(tools::md5sum(source_files)),unname(tools::md5sum(copied_files))))
  for(f in list.files(experiment,pattern="\\.R$",full.names=TRUE))parse(f)
  legacy<-discover_existing(root)
  write.csv(legacy,file.path(experiment,"legacy_sources.csv"),row.names=FALSE)
  cat("Frozen targeted experiment:",experiment,"\nNew fits:",nrow(cases),"Legacy fits for prediction only:",nrow(legacy),"\n")
  if(Sys.getenv("TARGET_FREEZE_ONLY")=="1")return(invisible(experiment))
  # Validate the largest request before submitting the preparation or fit jobs.
  Sys.setenv(DCM_ROOT=root,EXPERIMENT_DIR=experiment,CASE_ID=1L)
  constraint<-Sys.getenv("TARGET_CONSTRAINT","intel&cascadelake")
  fit_time<-Sys.getenv("TARGET_TIME","7-00:00:00")
  test_args<-c("--test-only",paste0("--chdir=",root),"--cpus-per-task=16","--mem=64G",
    paste0("--constraint=",constraint),paste0("--time=",fit_time),"--hint=nomultithread","--export=ALL",file.path(experiment,"fit.slurm"))
  response<-system2("sbatch",vapply(test_args,shQuote,character(1)),stdout=TRUE,stderr=TRUE)
  if(!is.null(attr(response,"status"))&&attr(response,"status")!=0)stop(paste(response,collapse="\n"))
  result<-queue_frozen_experiment(experiment,root,cases,scheduler_submit,legacy,constraint,fit_time)
  writeLines(experiment,file.path(root,"logs","latest_targeted150_experiment.txt"))
  cat("Submitted",nrow(cases),"new fits plus",nrow(legacy),"legacy evaluations. Preparation job",result$prepare,
    "collector",result$collector,"\nMonitor: bash DCM_Ibex_targeted_tests/monitor_all_dcm_v5.sh --diagnostics\n")
  invisible(result)
}
if(sys.nframe()==0L) {
  a<-commandArgs(trailingOnly=TRUE)
  if(length(a)!=2L)stop("Usage: Rscript submit.R PROJECT_ROOT PACKAGE_DIRECTORY")
  submit_campaign(a[1],a[2])
}
