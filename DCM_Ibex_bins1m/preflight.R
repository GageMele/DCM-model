# Preflight on the login node: build Stan code for every case without sampling.
# Optional second argument: comma-separated case ids, written to preflight_cases_<ids>.csv
# for merging (the scheduler only reads preflight_cases.csv).
args<-commandArgs(trailingOnly=TRUE);experiment<-args[1]
Sys.setenv(EXPERIMENT_DIR=experiment)
source(file.path(experiment,"common.R"));source(file.path(experiment,"models.R"))
cases<-read.delim(file.path(experiment,"cases.tsv"),stringsAsFactors=FALSE)
only<-if(length(args)>1)as.integer(strsplit(args[2],",")[[1]]) else integer(0)
if(length(only)){if(!all(only %in% cases$case_id))stop("Unknown case ids");cases<-cases[cases$case_id %in% only,]}
levels<-readRDS(file.path(experiment,"region_levels.rds"))
source_exp<-readLines(file.path(experiment,"reuse_source.txt"))[1]
out<-list()
for(i in seq_len(nrow(cases))) {
  case<-cases[i,];d<-prepare_data(augment(readRDS(data_file(experiment,case))),case)
  row<-data.frame(case_id=case$case_id,valid=TRUE,reason="",training_casts=length(unique(d$cast_key)),
    observations=nrow(d),formula="")
  r<-tryCatch({
    if(case$mode=="reuse") {
      src<-file.path(source_exp,"runs",case$source_run)
      if(!file.exists(file.path(src,"FIT_COMPLETE.txt")))stop("source fit not complete: ",src)
      row$formula<-paste("reuse",src)
    } else {
      spec<-campaign_spec(case,levels)
      code<-brms::make_stancode(spec$formula,data=d,prior=spec$priors,threads=thread_specification(case))
      if(!nzchar(paste(code,collapse="")))stop("empty Stan code")
      row$formula<-paste(capture.output(spec$formula),collapse=" ")
    }
    if(case$cast!="none")cat(sprintf("case %d: %s data, %d rows from %d casts\n",case$case_id,case$cast,
      nrow(d),length(unique(d$cast_key))))
    if(case$equation=="dcm_depth")cat("  DCM classes in training data:",tabulate(d$dcm_class,length(dcm_classes(case))),"\n")
    if(case$equation=="bin_class")cat("  enhancement classes in training data:",tabulate(d$enh_class,length(enh_labels)),"\n")
    if(case$family=="lognormal")cat(sprintf("case %d: %d of %d rows floored (Chl<=0)\n",case$case_id,
      sum(d$chl_norm+1<ratio_floor),nrow(d)))
    NULL},error=function(e)conditionMessage(e))
  if(!is.null(r)){row$valid<-FALSE;row$reason<-r}
  cat(sprintf("case %2d %-45s %s %s\n",case$case_id,substr(case$label,1,45),if(row$valid)"OK" else "FAIL",row$reason))
  out[[i]]<-row
}
p<-do.call(rbind,out)
write.csv(p,file.path(experiment,if(length(only))paste0("preflight_cases_",paste(only,collapse="_"),".csv") else "preflight_cases.csv"),row.names=FALSE)
if(!all(p$valid))quit(status=3L)
cat("All",nrow(p),"cases passed preflight\n")
