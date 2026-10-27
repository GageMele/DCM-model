# Scheduling and reporting checks: base R only, no model fitting or networking.
a<-commandArgs(trailingOnly=TRUE);package<-if(length(a))a[1]else"DCM_Ibex_followup_tests"
for(f in list.files(package,pattern="\\.R$",full.names=TRUE))parse(f)
source(file.path(package,"campaign.R"));source(file.path(package,"submit.R"))
source(file.path(package,"04_collect_results.R"));source(file.path(package,"predictive.R"))
cases<-campaign_cases()
stopifnot(nrow(cases)==34L,sum(cases$mode=="posterior")==30L,sum(cases$mode=="prior")==4L,
  all(cases$subset %in% c("p10","p10_alt","p20")),all(cases$cpus==4L*cases$threads),
  all(cases$chains==4L),!anyDuplicated(cases$label))
for(pair in list(c(17,16),c(20,19))) {
  c<-cases[pair[1],];r<-cases[pair[2],]
  stopifnot(c$warmup==r$warmup,c$samples==r$samples,c$adapt_delta==r$adapt_delta,c$seed==r$seed,
    c$variance=="region_depth",r$variance=="region")
}
root<-tempfile("dcm_followup_test_");dir.create(file.path(root,"experiments","followup150_fixture","scheduler"),recursive=TRUE)
experiment<-file.path(root,"experiments","followup150_fixture")
write.table(cases,file.path(experiment,"cases.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
fake<-new.env(parent=emptyenv());fake$n<-1000L;fake$jobs<-list();fake$calls<-list()
lookup<-function(name)if(is.null(fake$jobs[[name]]))character()else fake$jobs[[name]]
submit<-function(args) {
  fake$n<-fake$n+1L;id<-as.character(fake$n)
  name<-sub("^--job-name=","",grep("^--job-name=",args,value=TRUE))
  fake$jobs[[name]]<-id;fake$calls[[length(fake$calls)+1L]]<-args;id
}
active<-function(known)ids(known)
invisible(capture.output(m<-queue_cases(experiment,root,cases,submit,lookup,active,function(id)TRUE)))
stopifnot(length(fake$calls)==70L,nrow(m)==34L,!anyDuplicated(c(m$fit_job,m$check_job)))
for(i in seq_len(nrow(cases))) {
  fit<-fake$calls[[2L*i]];check<-fake$calls[[2L*i+1L]]
  stopifnot(paste0("--cpus-per-task=",cases$cpus[i])%in%fit,paste0("--mem=",cases$memory_gb[i],"G")%in%fit,
    any(grepl("^--dependency=afterok:",fit)),paste0("--dependency=afterany:",m$fit_job[i])%in%check)
}
writeLines("done",file.path(experiment,"PREPARATION_COMPLETE.txt"));old_ids<-as.character(m$fit_job)
invisible(capture.output(repeat_m<-queue_cases(experiment,root,cases,submit,lookup,active,function(id)TRUE)))
stopifnot(length(fake$calls)==70L,identical(as.character(repeat_m$fit_job),old_ids))
# A returned job ID lost before local registration is recovered by its unique name.
unlink(file.path(experiment,"scheduler","fit_07.job"))
invisible(capture.output(recovered<-queue_cases(experiment,root,cases,submit,lookup,active,function(id)TRUE)))
stopifnot(length(fake$calls)==70L,identical(as.character(recovered$fit_job),old_ids))
uncertain<-try(ensure_job(experiment,"uncertain","uncertain_name",character(),function(args)stop("connection lost"),function(name)character()),silent=TRUE)
second<-try(ensure_job(experiment,"uncertain","uncertain_name",character(),function(args)stop("must not submit"),function(name)character()),silent=TRUE)
stopifnot(inherits(uncertain,"try-error"),inherits(second,"try-error"),grepl("uncertain",as.character(second)))
# Recover an expired dependency only after explicit rejection and terminal proof.
race_root<-tempfile("dcm_followup_race_");race_exp<-file.path(race_root,"experiments","followup150_fixture")
dir.create(file.path(race_exp,"scheduler"),recursive=TRUE)
writeLines("2000",file.path(race_exp,"scheduler","preparation.job"))
writeLines("2001",file.path(race_exp,"scheduler","fit_01.job"))
writeLines("done",file.path(race_exp,"PREPARATION_COMPLETE.txt"))
rejected<-FALSE
submit_race<-function(args){
  if(any(grepl("-q01$",args))&&any(grepl("^--dependency=",args))&&!rejected){
    rejected<<-TRUE;e<-simpleError("sbatch: error: Job dependency problem");class(e)<-c("definite_rejection",class(e));stop(e)
  };submit(args)
}
invisible(capture.output(race<-queue_cases(race_exp,race_root,cases[1,,drop=FALSE],submit_race,lookup,
  function(known)setdiff(ids(known),"2001"),function(id)TRUE)))
stopifnot(rejected,identical(as.character(race$fit_job),"2001"),length(fake$calls)==72L)
writeLines("fixture",file.path(experiment,"PREPARATION_FAILED.txt"));n<-length(fake$calls)
failure<-try(queue_cases(experiment,root,cases,submit,lookup,active,function(id)TRUE),silent=TRUE)
stopifnot(inherits(failure,"try-error"),length(fake$calls)==n)
unlink(file.path(experiment,"PREPARATION_FAILED.txt"))
cat("PASS: all34 requests/dependencies, repeats, completed preparation, orphan recovery, uncertainty stop and expired dependencies; zero replacement fits.\n")

invisible(capture.output(empty<-collect_results(experiment,root)));stopifnot(nrow(empty)==34L)
for(i in c(1L,3L,5L)) {
  run<-file.path(experiment,"runs",sprintf("case_%02d",i));dir.create(file.path(run,"checks"),recursive=TRUE)
  x<-data.frame(cases[i,],max_rhat=if(i==3)1.2 else 1.001,min_bulk_ess=800,min_tail_ess=800,
    divergences=if(i==3)2L else 0L,min_ebfmi=.7,validation_md5=if(i==5)"other"else"common",
    reserve_md5="reserve",validation_conditional_rmse=sqrt(mean(c(i,2*i))))
  write.csv(x,file.path(run,"test_summary.csv"),row.names=FALSE)
  writeLines("done",file.path(run,"FIT_COMPLETE.txt"));writeLines("done",file.path(run,"CHECKS_COMPLETE.txt"))
  cs<-data.frame(cast_key=c("A","B"),Region=c("R1","R2"),mse=c(i,2*i),shape_mse=c(i,2*i),
    mae=c(i,i),interval_score90=c(i,i),peak_depth_abs_error=c(i,i))
  write.csv(cs,file.path(run,"checks","validation_conditional_cast_metrics.csv"),row.names=FALSE)
  write.csv(cs,file.path(run,"checks","reserve_conditional_cast_metrics.csv"),row.names=FALSE)
}
write.csv(data.frame(cast_key=c("A","B"),Region=c("R1","R2"),Ocean="O",date=c("D1","D2")),file.path(experiment,"split_membership.csv"),row.names=FALSE)
invisible(capture.output(all<-collect_results(experiment,root)))
pairs<-read.csv(file.path(experiment,"results","all_matched_cast_comparisons.csv"))
stopifnot(nrow(all)==34L,nrow(pairs)==4L,sum(pairs$dataset=="validation")==1L,
  any(!pairs$a_convergence|!pairs$b_convergence),any(all$divergences==2L,na.rm=TRUE))
z<-data.frame(cast_key="A",Region="R",depth=c(10,20,30),chl_norm=c(0,2,0),predicted=c(0,1,0),lower=-1,upper=3,interval_score90=4)
pm<-profile_metrics(z);stopifnot(pm$peak_depth_abs_error==0,pm$positive_auc_abs_error==10,isTRUE(all.equal(pm$mse,1/3)))
cat("PASS: all34 statuses, diagnostic failures retained, separate primary/reserve cohorts, recomputed RMSE, date-bootstrap and profile/AUC summaries.\n")
