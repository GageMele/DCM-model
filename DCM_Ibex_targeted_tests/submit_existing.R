# Append predictions for older fits that completed after initial submission.
# This entry point never submits a model-fitting job.
refresh_evaluations<-function(experiment,root,submit_job,active_lookup=active_scheduler_jobs) {
  if(!file.exists(file.path(experiment,"PREPARATION_COMPLETE.txt")))stop("Wait for targeted preparation to complete")
  legacy<-discover_existing(root)
  path<-file.path(experiment,"existing_evaluations.tsv")
  recorded<-read.delim(path,stringsAsFactors=FALSE)
  legacy<-legacy[!legacy$external_id %in% recorded$external_id,,drop=FALSE]
  if(!nrow(legacy))cat("All currently completed eligible fits already have evaluation jobs.\n")
  Sys.setenv(DCM_ROOT=root,EXPERIMENT_DIR=experiment)
  base<-c("--parsable",paste0("--chdir=",root),"--export=ALL")
  jobs<-character()
  for(i in seq_len(nrow(legacy))) {
    z<-legacy[i,];evaluation<-file.path(experiment,"evaluations",z$external_id)
    Sys.setenv(SOURCE_RUN=z$source_run,EVAL_DIR=evaluation,EXTERNAL_ID=z$external_id)
    id<-submit_job(c(base,paste0("--job-name=dcm-t150-old-extra",i),file.path(experiment,"score_existing.slurm")))
    recorded<-rbind(recorded,data.frame(external_id=z$external_id,job_id=id,source_run=z$source_run,run_dir=evaluation))
    write.table(recorded,path,sep="\t",row.names=FALSE,quote=FALSE)
    jobs<-c(jobs,id)
    cat("Queued predictions only:",z$external_id,"job",id,"\n")
  }
  submit_collection(experiment,root,submit_job,active_lookup)
  invisible(recorded)
}
if(sys.nframe()==0L) {
  a<-commandArgs(trailingOnly=TRUE)
  if(length(a)!=2L)stop("Usage: Rscript submit_existing.R EXPERIMENT PROJECT_ROOT")
  source(file.path(a[1],"submit.R"))
  refresh_evaluations(a[1],a[2],scheduler_submit)
}
