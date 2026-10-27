# Finish scheduling an already-submitted targeted campaign. Never submit fits.
repair_existing_campaign<-function(experiment,root,submit_job=NULL,active_lookup=NULL) {
  experiment<-normalizePath(experiment,winslash="/",mustWork=TRUE)
  root<-normalizePath(root,winslash="/",mustWork=TRUE)
  if(!startsWith(experiment,paste0(root,"/experiments/targeted150_")))stop("Unexpected experiment path")
  needed<-c("PREPARATION_COMPLETE.txt","manifest.tsv","cases.tsv","existing_evaluations.tsv",
    "submit.R","submit_existing.R","score_existing.slurm","collect.slurm")
  if(!all(file.exists(file.path(experiment,needed))))stop("The submitted campaign is incomplete or preparation has not passed")
  cases<-read.delim(file.path(experiment,"cases.tsv"),stringsAsFactors=FALSE)
  manifest<-read.delim(file.path(experiment,"manifest.tsv"),stringsAsFactors=FALSE)
  stopifnot(nrow(cases)==63L,nrow(manifest)==63L,
    setequal(manifest$case_id,cases$case_id),!anyDuplicated(manifest$case_id),
    length(unique(c(manifest$fit_job,manifest$check_job)))==126L,
    all(grepl("^[0-9]+$",c(manifest$fit_job,manifest$check_job))))
  before<-nrow(read.delim(file.path(experiment,"existing_evaluations.tsv")))
  original_hash<-unname(tools::md5sum(file.path(experiment,"manifest.tsv")))
  source(file.path(experiment,"submit.R"),local=TRUE)
  source(file.path(experiment,"submit_existing.R"),local=TRUE)
  if(is.null(submit_job))submit_job<-scheduler_submit
  if(is.null(active_lookup))active_lookup<-active_scheduler_jobs
  recorded<-refresh_evaluations(experiment,root,submit_job,active_lookup)
  stopifnot(identical(original_hash,unname(tools::md5sum(file.path(experiment,"manifest.tsv")))))
  collector<-readLines(file.path(experiment,"collection_job_id.txt"),warn=FALSE)
  message<-c(paste("UTC:",format(Sys.time(),tz="UTC",usetz=TRUE)),
    paste("Experiment:",experiment),"All 63 original fit/check job registrations preserved.",
    paste("Previously registered old-fit evaluations:",before),
    paste("Additional evaluation jobs submitted:",nrow(recorded)-before),
    paste("Total registered old-fit evaluations:",nrow(recorded)),
    paste("Final collection job:",collector),"No model-fitting jobs submitted or cancelled.")
  writeLines(message,file.path(experiment,"DEPENDENCY_REPAIR_COMPLETE.txt"))
  cat(paste(message,collapse="\n"),"\n")
  invisible(recorded)
}
if(sys.nframe()==0L) {
  a<-commandArgs(trailingOnly=TRUE)
  if(length(a)!=2L)stop("Usage: Rscript repair_existing_campaign.R EXPERIMENT ROOT")
  repair_existing_campaign(a[1],a[2])
}
