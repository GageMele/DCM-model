# Resumable Slurm submission. Accepted fits are never automatically replaced.
ids<-function(x)unique(grep("^[0-9]+$",trimws(as.character(x)),value=TRUE))
atomic_lines<-function(text,path) {
  temporary<-paste0(path,".tmp.",Sys.getpid());writeLines(text,temporary)
  if(!file.rename(temporary,path)) {
    if(.Platform$OS.type!="windows")stop("Cannot save scheduler record: ",path)
    stopifnot(file.copy(temporary,path,overwrite=TRUE));unlink(temporary)
  }
}
atomic_table<-function(x,path) {
  temporary<-paste0(path,".tmp.",Sys.getpid());write.table(x,temporary,sep="\t",row.names=FALSE,quote=FALSE)
  if(!file.rename(temporary,path)) {
    if(.Platform$OS.type!="windows")stop("Cannot save manifest")
    stopifnot(file.copy(temporary,path,overwrite=TRUE));unlink(temporary)
  }
}
hash_text<-function(text) {
  file<-tempfile();on.exit(unlink(file));writeLines(text,file);substr(unname(tools::md5sum(file)),1,10)
}
run_read<-function(command,args) {
  output<-suppressWarnings(system2(command,vapply(args,shQuote,character(1)),stdout=TRUE,stderr=TRUE))
  if(!is.null(attr(output,"status"))&&attr(output,"status")!=0)stop(command," failed: ",paste(output,collapse="\n"))
  output
}
active_jobs<-function(known)intersect(ids(known),ids(run_read("squeue",c("--me","--noheader","--format=%A"))))
find_named_job<-function(name) {
  live<-ids(run_read("squeue",c("--me","--noheader",paste0("--name=",name),"--format=%A")))
  archived<-run_read("sacct",c("-X","--noheader","--parsable2","--starttime=2026-10-18",
    paste0("--user=",Sys.getenv("USER")),paste0("--name=",name),"--format=JobIDRaw,JobName%100"))
  recovered<-vapply(strsplit(archived,"|",fixed=TRUE),function(z)if(length(z)>=2L && z[2]==name)z[1]else"",character(1))
  ids(c(live,recovered))
}
terminal_job<-function(id) {
  x<-run_read("sacct",c("-X","--noheader","--parsable2",paste0("--jobs=",id),"--format=JobIDRaw,State"))
  rows<-strsplit(x,"|",fixed=TRUE)
  states<-vapply(rows,function(z)if(length(z)>=2L&&z[1]==id)sub("[ +].*$","",z[2])else"",character(1))
  any(states %in% c("COMPLETED","FAILED","CANCELLED","TIMEOUT","OUT_OF_MEMORY","NODE_FAIL","PREEMPTED","BOOT_FAIL","DEADLINE","REVOKED"))
}
submit_job<-function(args) {
  output<-suppressWarnings(system2("sbatch",vapply(args,shQuote,character(1)),stdout=TRUE,stderr=TRUE))
  status<-attr(output,"status");accepted<-grep("^[0-9]+(;[^[:space:]]+)?$",trimws(output),value=TRUE)
  if((is.null(status)||status==0)&&length(accepted)==1L)return(sub(";.*$","",accepted))
  error<-simpleError(paste(output,collapse="\n"))
  if(!length(accepted)&&any(grepl("sbatch: error:",output,fixed=TRUE)))class(error)<-c("definite_rejection",class(error))
  stop(error)
}
ensure_job<-function(experiment,stage,name,args,submit=submit_job,lookup=find_named_job) {
  record<-file.path(experiment,"scheduler",paste0(stage,".job"))
  intent<-file.path(experiment,"scheduler",paste0(stage,".intent"))
  if(file.exists(record)) {
    id<-ids(readLines(record,warn=FALSE));if(length(id)!=1L)stop("Invalid saved job record")
    return(id)
  }
  found<-ids(lookup(name))
  if(length(found)>1L)stop("Multiple jobs share a unique stage name; refusing submission: ",name)
  if(length(found)==1L){atomic_lines(found,record);return(found)}
  if(file.exists(intent))stop("Earlier submission outcome is uncertain for ",name,". No duplicate submitted. Rerun after its Slurm record becomes visible; inspect scheduler/",stage,".intent if no job was accepted.")
  atomic_lines(c(name,format(Sys.time(),tz="UTC",usetz=TRUE)),intent)
  result<-tryCatch(submit(c("--parsable",paste0("--job-name=",name),args)),error=function(e)e)
  if(inherits(result,"error")) {
    if(inherits(result,"definite_rejection"))unlink(intent)
    stop(result)
  }
  stopifnot(length(ids(result))==1L);atomic_lines(result,record);result
}
all_work_jobs<-function(root) {
  paths<-c(Sys.glob(file.path(root,"experiments","*","manifest.tsv")),
    Sys.glob(file.path(root,"benchmarks","*","manifest.tsv")),
    Sys.glob(file.path(root,"experiments","*","existing_evaluations.tsv")))
  recorded<-unlist(lapply(paths,function(path){x<-read.delim(path,stringsAsFactors=FALSE)
    unlist(x[intersect(names(x),c("fit_job","check_job","job_id"))],use.names=FALSE)
  }),use.names=FALSE)
  pilot<-unlist(lapply(Sys.glob(file.path(root,"runs","*","run_config.rds")),function(path){
    tryCatch(readRDS(path)$job_id,error=function(e)character())
  }),use.names=FALSE)
  ids(c(recorded,pilot))
}
queue_cases<-function(experiment,root,cases,submit=submit_job,lookup=find_named_job,
  active=active_jobs,terminal=terminal_job) {
  Sys.setenv(DCM_ROOT=root,EXPERIMENT_DIR=experiment)
  token<-hash_text(experiment);base<-c(paste0("--chdir=",root),"--export=ALL")
  ensure<-function(stage,label,script,extra=character())ensure_job(experiment,stage,
    paste0("dcm-u150-",token,"-",label),c(base,extra,file.path(experiment,script)),submit,lookup)
  prep<-ensure("preparation","prep","prepare.slurm")
  atomic_lines(prep,file.path(experiment,"preparation_job_id.txt"))
  if(file.exists(file.path(experiment,"PREPARATION_FAILED.txt")))stop("Preparation failed; no replacement fits submitted. Inspect its error before starting a new campaign.")
  prep_dep<-function()if(file.exists(file.path(experiment,"PREPARATION_COMPLETE.txt")))character()else
    c(paste0("--dependency=afterok:",prep),"--kill-on-invalid-dep=yes")
  registry<-file.path(experiment,"manifest.tsv")
  manifest<-if(file.exists(registry))read.delim(registry,stringsAsFactors=FALSE)else
    data.frame(case_id=integer(),label=character(),fit_job=character(),check_job=character(),cpus=integer(),run_dir=character())
  if(anyDuplicated(manifest$case_id)||any(!manifest$case_id %in% cases$case_id))stop("Invalid saved manifest")
  for(i in seq_len(nrow(cases))) {
    case<-cases[i,];Sys.setenv(CASE_ID=case$case_id)
    fit<-ensure(sprintf("fit_%02d",case$case_id),sprintf("c%02d",case$case_id),"fit.slurm",
      c(prep_dep(),paste0("--cpus-per-task=",case$cpus),paste0("--mem=",case$memory_gb,"G"),
        paste0("--time=",if(case$mode=="prior")"02:00:00"else Sys.getenv("FOLLOWUP_TIME","7-00:00:00")),
        paste0("--constraint=",Sys.getenv("FOLLOWUP_CONSTRAINT","intel&cascadelake")),"--hint=nomultithread"))
    j<-match(case$case_id,manifest$case_id)
    if(is.na(j)) {
      manifest<-rbind(manifest,data.frame(case_id=case$case_id,label=case$label,fit_job=fit,check_job="NA",cpus=case$cpus,
        run_dir=file.path(experiment,"runs",sprintf("case_%02d",case$case_id))))
      j<-nrow(manifest)
    } else if(as.character(manifest$fit_job[j])!=fit)stop("Fit ledger disagrees with saved manifest")
    atomic_table(manifest,registry)
    # Never infer completion from a briefly absent squeue record.
    finished<-any(file.exists(file.path(manifest$run_dir[j],c("FIT_COMPLETE.txt","FIT_FAILED.txt"))))
    dep<-if(finished)character()else paste0("--dependency=afterany:",fit)
    check<-tryCatch(ensure(sprintf("check_%02d",case$case_id),sprintf("q%02d",case$case_id),"check.slurm",dep),error=function(e)e)
    if(inherits(check,"error")) {
      invalid<-inherits(check,"definite_rejection")&&grepl("Job dependency problem|Invalid job id",conditionMessage(check))
      if(!invalid||!terminal(fit))stop(check)
      check<-ensure(sprintf("check_%02d",case$case_id),sprintf("q%02d",case$case_id),"check.slurm")
    }
    manifest$check_job[j]<-check;atomic_table(manifest,registry)
    cat("Registered case",case$case_id,case$label,"fit",fit,"check",check,"\n")
  }
  known<-all_work_jobs(root);signature<-hash_text(sort(known))
  for(attempt in 1:3) {
    live<-active(known);dep<-if(length(live))paste0("--dependency=afterany:",paste(live,collapse=":"))else character()
    collector<-tryCatch(ensure(paste0("collect_",signature),paste0("results-",signature),"collect.slurm",dep),error=function(e)e)
    if(!inherits(collector,"error"))break
    invalid<-inherits(collector,"definite_rejection")&&grepl("Job dependency problem|Invalid job id",conditionMessage(collector))
    if(!invalid||attempt==3L)stop(collector)
  }
  atomic_lines(collector,file.path(experiment,"collection_job_id.txt"))
  atomic_lines(known,file.path(experiment,"collection_watch_job_ids.txt"))
  cat("All",nrow(cases),"fit/check pairs preserved. Collector",collector,"waits for",length(live),"active jobs across all campaigns.\n")
  invisible(manifest)
}
freeze_experiment<-function(root,package) {
  pointer<-file.path(root,"logs","latest_followup150_experiment.txt")
  if(file.exists(pointer)&&Sys.getenv("FOLLOWUP_NEW")!="1") {
    experiment<-readLines(pointer,n=1L)
    if(!file.exists(file.path(experiment,"CONFIG_FROZEN.txt")))stop("Initialization incomplete; inspect ",experiment)
    if(!identical(readLines(file.path(experiment,"package_md5.txt")),unname(tools::md5sum(file.path(package,"SHA256SUMS.txt")))))stop("Package differs from frozen experiment")
    return(experiment)
  }
  pointer_source<-file.path(root,"logs","latest_targeted150_experiment.txt")
  source<-Sys.getenv("FOLLOWUP_SOURCE",if(file.exists(pointer_source))readLines(pointer_source,n=1L)else"")
  required<-c("PREPARATION_COMPLETE.txt","ibex_environment.sh","scaling.rds","region_levels.rds","full_input_250m.rds",
    file.path("data",c("primary_p10_d150.rds","primary_p20_d150.rds","primary_p50_d150.rds","primary_validation_d150.rds","date_blocked_p10_d150.rds","date_blocked_validation_d150.rds")))
  if(!nzchar(source)||!all(file.exists(file.path(source,required))))stop("Need the successfully prepared targeted150 source campaign")
  experiment<-file.path(root,"experiments",paste0("followup150_",format(Sys.time(),"%Y%m%dT%H%M%SZ",tz="UTC"),"_",Sys.getpid()))
  dir.create(file.path(experiment,"scheduler"),recursive=TRUE)
  files<-list.files(package,pattern="\\.(R|slurm|sh|awk|md|tsv)$",full.names=TRUE)
  stopifnot(all(file.copy(files,experiment)),file.copy(file.path(source,"ibex_environment.sh"),experiment))
  writeLines(source,file.path(experiment,"source_experiment.txt"))
  writeLines(unname(tools::md5sum(file.path(package,"SHA256SUMS.txt"))),file.path(experiment,"package_md5.txt"))
  writeLines("Immutable case matrix and source pointer frozen",file.path(experiment,"CONFIG_FROZEN.txt"))
  atomic_lines(experiment,pointer);experiment
}
submit_campaign<-function(root,package) {
  root<-normalizePath(root,mustWork=TRUE);package<-normalizePath(package,mustWork=TRUE)
  source(file.path(package,"campaign.R"),local=TRUE);cases<-campaign_cases()
  stopifnot(isTRUE(all.equal(cases,read.delim(file.path(package,"cases.tsv")),check.attributes=FALSE)))
  experiment<-freeze_experiment(root,package)
  Sys.setenv(DCM_ROOT=root,EXPERIMENT_DIR=experiment,CASE_ID=1L)
  if(Sys.getenv("FOLLOWUP_FREEZE_ONLY")=="1")return(experiment)
  for(i in which(!duplicated(cases[c("cpus","memory_gb")])))run_read("sbatch",
    c("--test-only",paste0("--chdir=",root),"--export=ALL",paste0("--cpus-per-task=",cases$cpus[i]),
      paste0("--mem=",cases$memory_gb[i],"G"),paste0("--time=",Sys.getenv("FOLLOWUP_TIME","7-00:00:00")),
      paste0("--constraint=",Sys.getenv("FOLLOWUP_CONSTRAINT","intel&cascadelake")),"--hint=nomultithread",file.path(experiment,"fit.slurm")))
  result<-queue_cases(experiment,root,cases)
  cat("Experiment:",experiment,"\nMonitor: bash DCM_Ibex_followup_tests/monitor_all_dcm_v6.sh --diagnostics\n")
  invisible(result)
}
if(sys.nframe()==0L){a<-commandArgs(trailingOnly=TRUE);stopifnot(length(a)==2L);submit_campaign(a[1],a[2])}
