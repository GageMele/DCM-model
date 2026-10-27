# Offline integration checks for discovery and incremental prediction scheduling.
a<-commandArgs(trailingOnly=TRUE);package<-if(length(a))a[1]else"DCM_Ibex_targeted_tests"
source(file.path(package,"submit.R"));source(file.path(package,"submit_existing.R"))
root<-tempfile("targeted_discovery_");dir.create(root)
old<-file.path(root,"experiments","extensive150_fake");dir.create(old,recursive=TRUE)
cases<-data.frame(case_id=1:4,mode=c("posterior","posterior","prior","posterior"),subset=c("p10","p05","p10","p20"))
write.table(cases,file.path(old,"cases.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
registry<-data.frame(case_id=1:4,run_dir=file.path(old,paste0("run",1:4)))
write.table(registry,file.path(old,"manifest.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
for(run in registry$run_dir[1:3]) {
  dir.create(run)
  invisible(file.create(file.path(run,c("FIT_COMPLETE.txt","fit.rds","model_data.rds","run_config.rds"))))
}
found<-discover_existing(root)
stopifnot(nrow(found)==1,grepl("case_1$",found$external_id))
experiment<-file.path(root,"experiments","targeted150_fake");dir.create(experiment)
writeLines("prepared",file.path(experiment,"PREPARATION_COMPLETE.txt"))
write.table(data.frame(external_id=character(),job_id=character(),source_run=character(),run_dir=character()),
  file.path(experiment,"existing_evaluations.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
requests<-list();stub<-function(args){requests[[length(requests)+1L]]<<-args;as.character(800000L+length(requests))}
invisible(capture.output(recorded<-refresh_evaluations(experiment,root,stub,active_lookup=function(ids)ids)))
stopifnot(nrow(recorded)==1,length(requests)==2,
  any(grepl("score_existing.slurm$",requests[[1]])),!any(grepl("fit.slurm$",unlist(requests))),
  "--dependency=afterany:800001" %in% requests[[2]])
invisible(capture.output(refresh_evaluations(experiment,root,stub,active_lookup=function(ids)ids)))
stopifnot(length(requests)==2)
# A newly completed 20% model adds only its prediction and a collector job.
dir.create(registry$run_dir[4]);invisible(file.create(file.path(registry$run_dir[4],c("FIT_COMPLETE.txt","fit.rds","model_data.rds","run_config.rds"))))
invisible(capture.output(recorded<-refresh_evaluations(experiment,root,stub,active_lookup=function(ids)ids)))
stopifnot(nrow(recorded)==2,length(requests)==4,!anyDuplicated(recorded$external_id))
cat("PASS: all prior campaign discovery, >=10% posterior eligibility, incomplete exclusion, append-only predictions, duplicate prevention, zero refits.\n")
