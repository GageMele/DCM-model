# The only allowed submission targets here are check.slurm and collect.slurm.
a<-commandArgs(trailingOnly=TRUE);root<-normalizePath(a[1]);package<-normalizePath(a[2])
experiment<-normalizePath(a[3]);repair<-normalizePath(Sys.getenv("ONSET_REPAIR_DIR"))
source(file.path(package,"ledger.R"))
Sys.setenv(DCM_ROOT=root,EXPERIMENT_DIR=experiment,ONSET_REPAIR_DIR=repair)
cases<-read.delim(file.path(experiment,"cases.tsv"),stringsAsFactors=FALSE)
original<-read.delim(file.path(experiment,"manifest.tsv"),stringsAsFactors=FALSE)
controls<-read.delim(file.path(experiment,"source_controls.tsv"),stringsAsFactors=FALSE)
stopifnot(nrow(cases)==34L,sum(cases$mode=="posterior")==30L,!anyDuplicated(original$case_id))
path<-file.path(repair,"repair_manifest.tsv")
manifest<-if(file.exists(path))read.delim(path,stringsAsFactors=FALSE)else
  data.frame(case_id=integer(),original_fit_job=character(),original_check_job=character(),check_job=character(),run_dir=character())
token<-hash_text(experiment);base<-c(paste0("--chdir=",root),"--export=ALL")
for(i in which(cases$mode=="posterior")) {
  case<-cases[i,];item<-original[original$case_id==case$case_id,,drop=FALSE];stopifnot(nrow(item)==1L)
  run<-as.character(item$run_dir)
  if(file.exists(file.path(run,"CHECKS_COMPLETE.txt"))) {cat("Already checked:",case$case_id,"\n");next}
  donor<-run;fit_id<-as.character(item$fit_job);work<-ids(c(item$fit_job,item$check_job))
  if(case$execution=="reuse") {
    old<-controls[controls$case_id==case$case_id,,drop=FALSE];stopifnot(nrow(old)==1L)
    donor<-as.character(old$source_run);fit_id<-as.character(old$fit_job);work<-ids(c(work,old$fit_job,old$check_job))
  }
  if(!file.exists(file.path(donor,"FIT_COMPLETE.txt")) && terminal_job(fit_id)) {
    cat("Source fit ended without completion; no replacement fit:",case$case_id,"\n");next
  }
  Sys.setenv(CASE_ID=case$case_id)
  for(attempt in 1:3) {
    waiting<-active_jobs(work)
    if(!file.exists(file.path(donor,"FIT_COMPLETE.txt")) && !terminal_job(fit_id))waiting<-ids(c(waiting,fit_id))
    if(!any(file.exists(file.path(run,c("CHECKS_COMPLETE.txt","CHECKS_FAILED.txt")))) &&
       !terminal_job(as.character(item$check_job)))waiting<-ids(c(waiting,item$check_job))
    dep<-if(length(waiting))paste0("--dependency=afterany:",paste(waiting,collapse=":"))else character()
    job<-tryCatch(ensure_job(repair,sprintf("check_%02d",case$case_id),sprintf("dcm-or1-%s-q%02d",token,case$case_id),
      c(base,dep,file.path(repair,"check.slurm"))),error=function(e)e)
    if(!inherits(job,"error"))break
    invalid<-inherits(job,"definite_rejection")&&grepl("Job dependency problem|Invalid job id",conditionMessage(job))
    if(!invalid||attempt==3L)stop(job)
  }
  j<-match(case$case_id,manifest$case_id)
  row<-data.frame(case_id=case$case_id,original_fit_job=as.character(item$fit_job),
    original_check_job=as.character(item$check_job),check_job=job,run_dir=run)
  if(is.na(j))manifest<-rbind(manifest,row)else if(as.character(manifest$check_job[j])!=job)stop("Repair ledger mismatch")
  atomic_table(manifest,path);cat("Registered evaluation only: case",case$case_id,"job",job,"\n")
}
known<-ids(manifest$check_job);signature<-hash_text(sort(known))
for(attempt in 1:3) {
  waiting<-active_jobs(known)
  dep<-if(length(waiting))paste0("--dependency=afterany:",paste(waiting,collapse=":"))else character()
  collector<-tryCatch(ensure_job(repair,paste0("collect_",signature),paste0("dcm-or1-",token,"-results-",signature),
    c(base,dep,file.path(repair,"collect.slurm"))),error=function(e)e)
  if(!inherits(collector,"error"))break
  invalid<-inherits(collector,"definite_rejection")&&grepl("Job dependency problem|Invalid job id",conditionMessage(collector))
  if(!invalid||attempt==3L)stop(collector)
}
atomic_lines(collector,file.path(repair,"collection_job_id.txt"))
cat("Collector:",collector,"; existing fit registrations unchanged; no fitting jobs submitted or cancelled.\n")
