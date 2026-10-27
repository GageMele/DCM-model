# Reuse the existing check scripts and saved fits; never enter the fitting script.
experiment<-normalizePath(Sys.getenv("EXPERIMENT_DIR"),mustWork=TRUE)
repair<-normalizePath(Sys.getenv("ONSET_REPAIR_DIR"),mustWork=TRUE)
source(file.path(experiment,"common.R"));source(file.path(experiment,"campaign.R"))
case<-read_case(experiment,Sys.getenv("CASE_ID"));stopifnot(case$mode=="posterior")
run<-file.path(experiment,"runs",sprintf("case_%02d",case$case_id))
dir.create(run,recursive=TRUE,showWarnings=FALSE)
stopifnot(startsWith(normalizePath(run),paste0(experiment,"/runs/")))
if(file.exists(file.path(run,"CHECKS_COMPLETE.txt")))quit(status=0L)
donor<-run
if(case$execution=="reuse") {
  c<-read.delim(file.path(experiment,"source_controls.tsv"),stringsAsFactors=FALSE)
  donor<-as.character(c$source_run[c$case_id==case$case_id]);stopifnot(length(donor)==1L)
}
if(!all(file.exists(file.path(donor,c("FIT_COMPLETE.txt","fit.rds"))))) {
  writeLines("Source fitting did not complete; no replacement fit submitted",
    file.path(run,"CHECKS_FAILED.txt"));stop("No completed source fit")
}
# Preserve earlier check outputs after any original checking job has ended.
backup<-file.path(repair,"case_backups",sprintf("case_%02d",case$case_id),
  paste0(format(Sys.time(),"%Y%m%dT%H%M%SZ",tz="UTC"),"_",Sys.getpid()))
dir.create(backup,recursive=TRUE)
checks<-file.path(run,"checks")
if(dir.exists(checks)) {
  if(!startsWith(normalizePath(checks),paste0(experiment,"/runs/")))stop("Unexpected checks path")
  if(!file.rename(checks,file.path(backup,"checks")))stop("Cannot preserve previous checks")
}
for(n in c("test_summary.csv","diagnostic_summary.csv","CHECKS_FAILED.txt","CHECK_STAGE.txt")) {
  p<-file.path(run,n)
  if(file.exists(p) && !file.copy(p,file.path(backup,n)))stop("Cannot preserve ",n)
}
failure<-file.path(run,"CHECKS_FAILED.txt")
if(file.exists(failure) && !file.rename(failure,file.path(backup,"previous_failure_marker.txt")))
  stop("Cannot archive previous failure marker")
writeLines(paste("Existing fit reused by evaluation job",Sys.getenv("SLURM_JOB_ID")),
  file.path(run,"ONSET_SUPPORT_REPAIR.txt"))
tryCatch({
  for(p in c("brms","posterior","mgcv","ggplot2","bayesplot"))
    if(!requireNamespace(p,quietly=TRUE))stop("Missing package: ",p)
  source(file.path(experiment,if(case$execution=="reuse")"reuse_check.R" else "03_check_case.R"))
  if(!file.exists(file.path(run,"CHECKS_COMPLETE.txt")))stop("Evaluation did not finish")
},error=function(e){writeLines(conditionMessage(e),file.path(run,"CHECKS_FAILED.txt"));stop(e)})
