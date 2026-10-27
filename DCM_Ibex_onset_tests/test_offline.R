# Base-R gates run on Ibex BEFORE any sbatch submission; no fitting or networking.
a <- commandArgs(trailingOnly=TRUE);package <- a[1]
for(f in list.files(package,pattern="\\.R$",full.names=TRUE))parse(f)
Sys.setenv(ONSET_PACKAGE=package)
source(file.path(package,"campaign.R"));source(file.path(package,"submit.R"))
source(file.path(package,"onset_metrics.R"))
source(file.path(package,"04_collect_results.R"))
cases <- campaign_cases()
stopifnot(isTRUE(all.equal(cases,read.delim(file.path(package,"cases.tsv")),check.attributes=FALSE)),
  nrow(cases)==34L,sum(cases$execution=="fit")==28L,
  sum(cases$execution=="fit"&cases$mode=="posterior")==24L,
  sum(cases$execution=="reuse")==6L)
# Exact onset fixtures: known crossing, vertical displacement, absence, censoring,
# persistence, unsupported gaps, and a matching model with zero error.
z <- 10:150
curve <- onset_curve(z,pmax(0,(z-40)/20))
o <- onset_from_curve(curve,.5);stopifnot(o$status=="detected",o$depth==50)
shift <- onset_from_curve(onset_curve(z,pmax(0,(z-60)/20)),.5)
stopifnot(shift$depth-o$depth==20)
stopifnot(onset_from_curve(onset_curve(z,rep(0,length(z))),.5)$status=="not_reached_in_window")
stopifnot(onset_from_curve(onset_curve(z,rep(1,length(z))),.5)$status=="left_censored")
spike <- rep(0,length(z));spike[50:51] <- 2
stopifnot(onset_from_curve(onset_curve(z,spike),.5)$status=="not_reached_in_window")
gap <- onset_curve(c(10:30,50:150),c(rep(0,21),rep(1,101)))
stopifnot(onset_from_curve(gap,.5)$status=="gap_censored")
tab <- data.frame(cast_key="A",Region="R",Ocean="O",date="D",depth=z,
  MLD_temp_MIX=40,chl_norm=pmax(0,(z-40)/20),predicted=pmax(0,(z-40)/20))
m <- onset_profile_metrics(tab,rbind(tab$predicted,tab$predicted))
stopifnot(all(m$onset_abs_error==0),all(m$posterior_detected_draws==2),
          all(m$onset_interval_covers_observed_given_detected))
root <- tempfile("onset_scheduler_")
experiment <- file.path(root,"experiments","onset150_fixture")
dir.create(file.path(experiment,"scheduler"),recursive=TRUE)
write.table(cases,file.path(experiment,"cases.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
control <- data.frame(case_id=1:6,source_case=c(1,5,13,15,22,24),
  fit_job=as.character(101:106),check_job=as.character(201:206),source_run=file.path(root,paste0("old",1:6)))
write.table(control,file.path(experiment,"source_controls.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
fake <- new.env(parent=emptyenv());fake$n <- 1000L;fake$jobs <- list();fake$calls <- list()
lookup <- function(name)if(is.null(fake$jobs[[name]]))character()else fake$jobs[[name]]
submit <- function(args) {
  fake$n<-fake$n+1L;id<-as.character(fake$n)
  name<-sub("^--job-name=","",grep("^--job-name=",args,value=TRUE))
  fake$jobs[[name]]<-id;fake$calls[[length(fake$calls)+1L]]<-args;id
}
active <- function(known)ids(known)
invisible(capture.output(m <- queue_cases(experiment,root,cases,submit,lookup,active,function(id)TRUE)))
# 1 preparation + 28 new fits + 34 checks + 1 collector; six source IDs reused.
stopifnot(length(fake$calls)==64L,nrow(m)==34L,
  identical(as.character(m$fit_job[1:6]),as.character(101:106)))
fitcalls <- Filter(function(x)any(grepl("/fit.slurm$",x)),fake$calls)
reusecalls <- Filter(function(x)any(grepl("/reuse.slurm$",x)),fake$calls)
stopifnot(length(fitcalls)==28L,length(reusecalls)==6L,
  all(vapply(fitcalls,function(x)any(grepl("^--dependency=afterok:",x)),logical(1))),
  all(vapply(reusecalls,function(x)any(grepl("^--dependency=afterok:.*afterany:",x)),logical(1))))
for(j in seq_along(fitcalls)) {
  c <- cases[cases$execution=="fit",][j,]
  stopifnot(paste0("--cpus-per-task=",c$cpus)%in%fitcalls[[j]],
            paste0("--mem=",c$memory_gb,"G")%in%fitcalls[[j]])
}
writeLines("done",file.path(experiment,"PREPARATION_COMPLETE.txt"))
invisible(capture.output(again <- queue_cases(experiment,root,cases,submit,lookup,active,function(id)TRUE)))
same_registry <- function(a,b)identical(lapply(a,as.character),lapply(b,as.character))
stopifnot(length(fake$calls)==64L,same_registry(m,again))
unlink(file.path(experiment,"scheduler","fit_07.job"))
invisible(capture.output(recovered <- queue_cases(experiment,root,cases,submit,lookup,active,function(id)TRUE)))
stopifnot(length(fake$calls)==64L,same_registry(m,recovered))
uncertain <- try(ensure_job(experiment,"uncertain","uncertain_name",character(),
  function(args)stop("connection lost"),function(name)character()),silent=TRUE)
again_uncertain <- try(ensure_job(experiment,"uncertain","uncertain_name",character(),
  function(args)stop("must not submit"),function(name)character()),silent=TRUE)
stopifnot(inherits(uncertain,"try-error"),inherits(again_uncertain,"try-error"))
writeLines("fixture",file.path(experiment,"PREPARATION_FAILED.txt"))
blocked <- try(queue_cases(experiment,root,cases,submit,lookup,active,function(id)TRUE),silent=TRUE)
stopifnot(inherits(blocked,"try-error"),length(fake$calls)==64L)
unlink(file.path(experiment,"PREPARATION_FAILED.txt"))
# Keep incomplete and diagnostically poor fits in the comparison, with an onset
# contrast calculated from the same casts rather than a changing detectable subset.
for(i in c(1L,8L,10L)) {
  run <- file.path(experiment,"runs",sprintf("case_%02d",i))
  dir.create(file.path(run,"checks"),recursive=TRUE)
  fixture <- tab
  offset <- if(i==8L)0 else if(i==1L)5 else 10
  fixture$predicted <- pmax(0,(fixture$depth-40-offset)/20)
  pm <- onset_profile_metrics(fixture)
  write.csv(pm,file.path(run,"checks","validation_conditional_onset_cast_metrics.csv"),row.names=FALSE)
  summary <- data.frame(cases[i,],max_rhat=if(i==10L)1.2 else 1.001,
    min_bulk_ess=800,min_tail_ess=800,divergences=if(i==10L)2L else 0L,
    min_ebfmi=.7,validation_md5="same",reserve_md5="reserve")
  write.csv(summary,file.path(run,"test_summary.csv"),row.names=FALSE)
  writeLines("done",file.path(run,"FIT_COMPLETE.txt"))
  writeLines("done",file.path(run,"CHECKS_COMPLETE.txt"))
}
invisible(capture.output(all <- collect_results(experiment,root)))
invisible(capture.output(collect_onset(experiment)))
status <- read.csv(file.path(experiment,"results","onset_case_status.csv"))
scores <- read.csv(file.path(experiment,"results","onset_all_scores.csv"))
pairs <- read.csv(file.path(experiment,"results","onset_planned_paired_contrasts.csv"))
stopifnot(nrow(all)==34L,nrow(status)==34L,nrow(scores)==9L,
  any(all$divergences==2L,na.rm=TRUE),
  all(pairs$delta_onset_mae_common_casts[pairs$case_id==8L]== -5),
  all(pairs$jointly_detected_casts==1L))
cat("PASS: parse/case matrix, onset shifts/absence/censoring/gaps/persistence, posterior intervals;\n",
    "28 new fits, six source fits reused, dependencies/resources, repeat/orphan recovery, uncertain submission stop;\n",
    "all statuses, poor diagnostics retained and matched onset contrasts.\n")
