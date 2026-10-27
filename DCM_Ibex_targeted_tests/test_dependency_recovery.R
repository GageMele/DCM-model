# Offline regression checks for the exact 63-fit / 67-of-117 recovery scenario.
# Only empty fixture files and scheduler stubs are used: no model is loaded.
a<-commandArgs(trailingOnly=TRUE);package<-if(length(a))a[1]else"DCM_Ibex_targeted_tests"
source(file.path(package,"campaign.R"));source(file.path(package,"submit.R"))
source(file.path(package,"submit_existing.R"))
source(file.path(package,"repair_existing_campaign.R"))
root<-tempfile("dcm_recovery_");dir.create(root)
old<-file.path(root,"experiments","extensive150_recovery");dir.create(old,recursive=TRUE)
old_cases<-data.frame(case_id=1:117,mode="posterior",subset="p10")
registry<-data.frame(case_id=1:117,run_dir=file.path(old,paste0("case_",1:117)))
write.table(old_cases,file.path(old,"cases.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(registry,file.path(old,"manifest.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
for(run in registry$run_dir) {
  dir.create(run)
  invisible(file.create(file.path(run,c("FIT_COMPLETE.txt","fit.rds","model_data.rds","run_config.rds"))))
}
legacy<-discover_existing(root);stopifnot(nrow(legacy)==117)
fixture<-function(name) {
  experiment<-file.path(root,"experiments",name);dir.create(experiment)
  writeLines("63 cases constructed",file.path(experiment,"PREPARATION_COMPLETE.txt"))
  x<-campaign_cases()
  write.table(x,file.path(experiment,"cases.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  manifest<-data.frame(case_id=x$case_id,label=x$label,fit_job=700000L+2*x$case_id-1L,
    check_job=700000L+2*x$case_id,cpus=x$cpus,run_dir=file.path(experiment,"runs",sprintf("case_%02d",x$case_id)))
  write.table(manifest,file.path(experiment,"manifest.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  recorded<-data.frame(legacy[1:67,],job_id=800000L+1:67,
    run_dir=file.path(experiment,"evaluations",legacy$external_id[1:67]))
  write.table(recorded,file.path(experiment,"existing_evaluations.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  experiment
}
experiment<-fixture("targeted150_recovery")
manifest_hash<-unname(tools::md5sum(file.path(experiment,"manifest.tsv")))
requests<-list();live<-c("700001","700002","800066","800067")
stub<-function(args) {
  requests[[length(requests)+1L]]<<-args
  id<-as.character(900000L+length(requests));live<<-c(live,id);id
}
lookup<-function(known)intersect(known,live)
invisible(capture.output(recorded<-refresh_evaluations(experiment,root,stub,lookup)))
stopifnot(nrow(recorded)==117,length(requests)==51,
  sum(vapply(requests,function(x)any(grepl("/score_existing.slurm$",x)),logical(1)))==50,
  !any(grepl("/fit.slurm$",unlist(requests))),!any(grepl("afterok",unlist(requests))),
  identical(as.character(recorded$job_id[1:67]),as.character(800000L+1:67)),
  identical(manifest_hash,unname(tools::md5sum(file.path(experiment,"manifest.tsv")))))
dependency<-requests[[51]][startsWith(requests[[51]],"--dependency=")]
watched<-strsplit(sub("--dependency=afterany:","",dependency),":",fixed=TRUE)[[1]]
stopifnot(setequal(watched,c("700001","700002","800066","800067",as.character(900001:900050))))
invisible(capture.output(refresh_evaluations(experiment,root,stub,lookup)))
stopifnot(length(requests)==51)

# Repeated recovery after an interrupted submission must preserve accepted jobs.
interrupted<-fixture("targeted150_interrupted")
accepted<-character();attempts<-0L
flaky<-function(args) {
  attempts<<-attempts+1L
  if(attempts==3L)stop("Submission rejected by scheduler")
  id<-as.character(910000L+attempts);accepted<<-c(accepted,id);id
}
invisible(capture.output(failed<-try(refresh_evaluations(interrupted,root,flaky,function(ids)ids),silent=TRUE)))
stopifnot(inherits(failed,"try-error"),nrow(read.delim(file.path(interrupted,"existing_evaluations.tsv")))==69)
resumed<-list()
resume_stub<-function(args){resumed[[length(resumed)+1L]]<<-args;as.character(920000L+length(resumed))}
invisible(capture.output(done<-refresh_evaluations(interrupted,root,resume_stub,function(ids)ids)))
stopifnot(nrow(done)==117,length(resumed)==49,all(accepted %in% done$job_id),!anyDuplicated(done$external_id))

# Refresh after a dependency disappears between squeue and sbatch.
race<-file.path(root,"race");dir.create(race)
write.table(data.frame(fit_job=101L,check_job=102L),file.path(race,"manifest.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
queries<-0L;race_requests<-list()
race_lookup<-function(ids){queries<<-queries+1L;intersect(ids,if(queries==1L)c("101","102")else"102")}
race_submit<-function(args) {
  race_requests[[length(race_requests)+1L]]<<-args
  if(any(grepl("afterany:101",args)))stop("sbatch: error: Batch job submission failed: Job dependency problem")
  "930001"
}
invisible(capture.output(submit_collection(race,root,race_submit,race_lookup)))
stopifnot(length(race_requests)==2,"--dependency=afterany:102" %in% race_requests[[2]])
instant<-list();instant_submit<-function(args){instant[[1]]<<-args;"930002"}
invisible(capture.output(submit_collection(race,root,instant_submit,function(ids)character())))
stopifnot(!any(startsWith(instant[[1]],"--dependency=")))

# Finished preparation must never generate a new dependency on an expired ID.
stopifnot(length(preparation_dependency(experiment,"53865144"))==0L)
pending<-file.path(root,"pending");dir.create(pending)
stopifnot("--dependency=afterok:53865144" %in% preparation_dependency(pending,"53865144"))
writeLines("failed",file.path(pending,"PREPARATION_FAILED.txt"))
stopifnot(inherits(try(preparation_dependency(pending,"53865144"),silent=TRUE),"try-error"))

# Exercise the delivered repair entry point with the same scheduler stubs.
entry<-fixture("targeted150_entry")
stopifnot(all(file.copy(file.path(package,c("submit.R","submit_existing.R")),entry)))
invisible(file.create(file.path(entry,c("score_existing.slurm","collect.slurm"))))
entry_requests<-list()
entry_submit<-function(args){entry_requests[[length(entry_requests)+1L]]<<-args;as.character(940000L+length(entry_requests))}
invisible(capture.output(entry_results<-repair_existing_campaign(entry,root,entry_submit,function(ids)ids)))
stopifnot(nrow(entry_results)==117,length(entry_requests)==51,
  file.exists(file.path(entry,"DEPENDENCY_REPAIR_COMPLETE.txt")),
  !any(grepl("/fit.slurm$",unlist(entry_requests))))
cat("PASS: exact 67/117 recovery adds 50 evaluations and one collector; zero refits; accepted IDs preserved; all active fits/checks/evaluations awaited; repeat/interruption recovery; expired-dependency race; completed/failed preparation handling.\n")
