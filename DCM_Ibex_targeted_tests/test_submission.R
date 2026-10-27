# Offline scheduler tests. A stub records requests; no sbatch call is made.
args<-commandArgs(trailingOnly=TRUE)
package<-if(length(args))args[1]else"DCM_Ibex_targeted_tests"
source(file.path(package,"campaign.R"));source(file.path(package,"submit.R"))
x<-campaign_cases();stopifnot(nrow(x)==63,sum(x$mode=="posterior")==57,sum(x$mode=="prior")==6,
  min(as.numeric(sub("p","",x$subset)))>=10,all(x$chains==4L))
experiment<-tempfile("targeted_scheduler_");dir.create(experiment)
requests<-list();next_id<-700000L
stub<-function(args){next_id<<-next_id+1L;requests[[length(requests)+1L]]<<-args;as.character(next_id)}
legacy<-data.frame(external_id=c("old_a","old_b"),source_run=c("/fake/a","/fake/b"))
invisible(capture.output(result<-queue_frozen_experiment(experiment,"/fake/project",x,stub,legacy,active_lookup=function(ids)ids)))
stopifnot(length(requests)==130,nrow(result$manifest)==63,nrow(result$evaluations)==2)
for(i in 1:63) {
  fit_request<-requests[[2*i]];check_request<-requests[[2*i+1]]
  stopifnot(paste0("--dependency=afterok:",result$prepare)%in%fit_request,
    paste0("--cpus-per-task=",x$cpus[i])%in%fit_request,
    paste0("--mem=",x$memory_gb[i],"G")%in%fit_request,
    paste0("--dependency=afterany:",result$manifest$fit_job[i])%in%check_request)
}
stopifnot(length(unique(c(result$manifest$fit_job,result$manifest$check_job)))==126,
  any(grepl("afterany:",requests[[length(requests)]])),
  all(grepl("/fake/",result$evaluations$source_run)))
cat("PASS: 63 cases; minimum 10% training; 57 posterior/6 prior; all CPU/memory requests; preparation/fit/check dependencies; legacy evaluations; unique job IDs.\n")
