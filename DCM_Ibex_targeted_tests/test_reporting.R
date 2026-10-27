# Base-R synthetic checks: no Stan packages, network calls, or fitting.
a<-commandArgs(trailingOnly=TRUE);package<-if(length(a))a[1]else"DCM_Ibex_targeted_tests"
source(file.path(package,"campaign.R"));source(file.path(package,"04_collect_results.R"));source(file.path(package,"predictive.R"))
experiment<-tempfile("targeted_reporting_");dir.create(experiment)
cases<-campaign_cases()
write.table(cases,file.path(experiment,"cases.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
invisible(capture.output(empty<-collect_results(experiment)))
stopifnot(nrow(empty)==0,nrow(read.csv(file.path(experiment,"results","case_status.csv")))==63)
for(i in 1:3) {
  run<-file.path(experiment,"runs",sprintf("case_%02d",i));dir.create(file.path(run,"checks"),recursive=TRUE)
  result<-data.frame(cases[i,],max_rhat=if(i==2)2 else 1.001,min_bulk_ess=if(i==2)3 else 800,
    min_tail_ess=800,divergences=if(i==2)100 else 0,min_ebfmi=.7,
    validation_md5=if(i==3)"different_holdout"else"same_holdout")
  write.csv(result,file.path(run,"test_summary.csv"),row.names=FALSE)
  writeLines("done",file.path(run,"CHECKS_COMPLETE.txt"))
  cast<-data.frame(cast_key=c("A","B"),mse=c(i,2*i),shape_mse=c(i,2*i),mae=c(i,i),
    interval_score90=c(i,i),peak_depth_abs_error=c(i,i))
  write.csv(cast,file.path(run,"checks","validation_conditional_cast_metrics.csv"),row.names=FALSE)
}
# Reproduce the earlier zero-row report failure without excluding other results.
empty_run<-file.path(experiment,"runs","case_04");dir.create(empty_run,recursive=TRUE)
write.csv(data.frame(mode=character()),file.path(empty_run,"diagnostic_summary.csv"),row.names=FALSE)
invisible(capture.output(all<-collect_results(experiment)))
stopifnot(nrow(all)==3,any(all$max_rhat==2),any(all$divergences==100))
pairs<-read.csv(file.path(experiment,"results","all_paired_cast_comparisons.csv"))
stopifnot(nrow(pairs)==1,pairs$matched_casts==2,!pairs$b_convergence)
z<-data.frame(cast_key="A",Region="R",depth=c(10,20,30),chl_norm=c(0,2,0),predicted=c(0,1,0),
  lower=-1,upper=3,interval_score90=4)
m<-profile_metrics(z);stopifnot(m$observed_peak_depth==20,m$peak_depth_abs_error==0,m$positive_auc_abs_error==10,
  isTRUE(all.equal(m$mse,1/3)))
cat("PASS: zero-result and empty-file collection, all 63 statuses, diagnostic failures retained, all-pairs cohort matching, profile/DCM/AUC metrics.\n")
