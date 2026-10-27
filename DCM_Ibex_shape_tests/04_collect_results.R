# Collect every case, including failed fits and full-data reconstruction.
args<-commandArgs(trailingOnly=TRUE);experiment<-normalizePath(args[1],mustWork=TRUE)
root<-dirname(dirname(experiment));output<-file.path(experiment,"results")
dir.create(output,showWarnings=FALSE)
read_table<-function(p,tsv=FALSE){if(!file.exists(p))return(data.frame())
  tryCatch(if(tsv)read.delim(p,stringsAsFactors=FALSE)else read.csv(p,stringsAsFactors=FALSE),error=function(e)data.frame())}
fill<-function(xs){xs<-Filter(function(x)is.data.frame(x)&&nrow(x)>0L,xs);if(!length(xs))return(data.frame())
  cols<-unique(unlist(lapply(xs,names)));do.call(rbind,lapply(xs,function(x){for(n in setdiff(cols,names(x)))x[[n]]<-rep(NA,nrow(x));x[cols]}))}
save_csv<-function(x,name){p<-file.path(output,name);tmp<-paste0(p,".tmp.",Sys.getpid());write.csv(x,tmp,row.names=FALSE);if(!file.rename(tmp,p))stop("Cannot publish ",p)}
status_of<-function(run){if(file.exists(file.path(run,"CHECKS_COMPLETE.txt")))"complete"else
  if(file.exists(file.path(run,"CHECKS_FAILED.txt")))"checks_failed"else
  if(file.exists(file.path(run,"FIT_FAILED.txt")))"fit_failed"else
  if(file.exists(file.path(run,"FIT_COMPLETE.txt")))"awaiting_checks"else"incomplete"}
safe_summary<-function(run){for(n in c("test_summary.csv","diagnostic_summary.csv")){
  x<-read_table(file.path(run,n));if(nrow(x)==1L)return(x)};data.frame()}
cases<-read_table(file.path(experiment,"cases.tsv"),TRUE);preflight<-read_table(file.path(experiment,"preflight_cases.csv"))
registry<-read_table(file.path(experiment,"manifest.tsv"),TRUE);records<-list();current<-list()

# Read current saved summaries from ALL older campaign manifests. Prediction-only
# evaluations are merged into their source fit rather than counted as refits.
evals<-fill(lapply(Sys.glob(file.path(root,"experiments","*","existing_evaluations.tsv")),read_table,tsv=TRUE))
manifests<-c(Sys.glob(file.path(root,"experiments","*","manifest.tsv")),Sys.glob(file.path(root,"benchmarks","*","manifest.tsv")))
for(path in manifests) {
  campaign<-dirname(path);m<-read_table(path,TRUE);cc<-read_table(file.path(campaign,"cases.tsv"),TRUE)
  if(!nrow(m)||!"run_dir"%in%names(m))next
  for(i in seq_len(nrow(m))) {
    id<-if("case_id"%in%names(m))m$case_id[i]else m$case[i]
    if(!length(id)||is.na(id))next
    c<-if("case_id"%in%names(cc))cc[cc$case_id==id,,drop=FALSE]else data.frame()
    row<-if(nrow(c)==1L)c else m[i,,drop=FALSE];run<-as.character(m$run_dir[i]);s<-safe_summary(run)
    if(nrow(s)==1L)for(n in names(s))row[[n]]<-s[[n]]
    if(nrow(evals)&&all(c("source_run","run_dir")%in%names(evals))) {
      e<-evals[evals$source_run==run,,drop=FALSE]
      for(j in seq_len(nrow(e))){p<-as.character(e$run_dir[j]);v<-safe_summary(p)
        if(file.exists(file.path(p,"CHECKS_COMPLETE.txt"))&&nrow(v)==1L&&
           "raw_training_identity_verified"%in%names(v)&&isTRUE(v$raw_training_identity_verified))
          for(n in grep("^(validation_|training_|reserve_|prediction_|raw_training_)",names(v),value=TRUE))row[[n]]<-v[[n]]}
    }
    key<-paste0(basename(campaign),":",sprintf("%03d",as.integer(id)))
    row$run_key<-key;row$campaign<-basename(campaign);row$run_dir<-run;row$status<-status_of(run)
    row$fit_completed<-file.exists(file.path(run,"FIT_COMPLETE.txt"))
    records[[key]]<-row
    if(normalizePath(campaign)==experiment)current[[as.character(id)]]<-row
  }
}
# Retain original pilot directories that predate campaign manifests.
for(p in Sys.glob(file.path(root,"runs","*","run_config.rds"))) {
  run<-dirname(p);key<-paste0("pilot:",basename(run))
  x<-safe_summary(run)
  if(!nrow(x))x<-data.frame(label=basename(run),mode=NA_character_)
  x$run_key<-key;x$campaign<-"pilot";x$run_dir<-run;x$status<-status_of(run)
  x$fit_completed<-file.exists(file.path(run,"FIT_COMPLETE.txt"));records[[key]]<-x
}
# Cases frozen before the first dispatch remain visible even without a manifest.
for(i in seq_len(nrow(cases))) {
  c<-cases[i,];id<-as.character(c$case_id);run<-file.path(experiment,"runs",sprintf("case_%03d",c$case_id))
  row<-if(!is.null(current[[id]]))current[[id]]else c
  row$run_key<-paste0(basename(experiment),":",sprintf("%03d",c$case_id));row$campaign<-basename(experiment)
  row$run_dir<-run;row$status<-status_of(run);row$heldout_eligible<-c$mode=="posterior"&&c$fraction<1
  if(nrow(preflight)&&any(preflight$case_id==c$case_id & !preflight$valid))row$status<-"preflight_failed"
  current[[id]]<-row;records[[row$run_key]]<-row
}
save_csv(fill(current),"shape_results.csv");save_csv(fill(records),"all_results.csv")

# Stream cast tables instead of retaining every model's draws or row predictions.
for(kind in c("cast_metrics","shape_cast_metrics","onset_cast_metrics","by_region","onset_scores")) {
  target<-file.path(output,paste0("all_shape_",kind,".csv"));temporary<-paste0(target,".tmp.",Sys.getpid())
  written<-FALSE;schema<-NULL
  for(i in seq_len(nrow(cases))) {
    c<-cases[i,];checks<-file.path(experiment,"runs",sprintf("case_%03d",c$case_id),"checks")
    paths<-list.files(checks,pattern=paste0("_",kind,"\\.csv$"),full.names=TRUE)
    if(kind=="cast_metrics")paths<-paths[!grepl("_(shape|onset)_cast_metrics\\.csv$",paths)]
    for(p in paths) {
      x<-read_table(p);if(!nrow(x))next
      tag<-sub(paste0("_",kind,"\\.csv$"),"",basename(p))
      mode<-if(endsWith(tag,"_conditional"))"conditional"else"population"
      dataset<-sub("_(conditional|population)$","",tag)
      x$case_id<-c$case_id;x$equation<-c$equation;x$subset<-c$subset;x$split<-c$split
      x$dataset<-dataset;x$prediction_mode<-mode;x$evaluation_role<-c$evaluation_role
      if(is.null(schema))schema<-names(x)
      if(!setequal(names(x),schema))stop("Inconsistent streamed metric schema: ",p)
      write.table(x[schema],temporary,sep=",",row.names=FALSE,col.names=!written,append=written)
      written<-TRUE
    }
  }
  if(written){if(!file.rename(temporary,target))stop("Cannot publish metric table")}
}

# Planned comparisons retain the same cast set. Full-data cases never receive
# held-out contrast labels, and onset/width errors use joint comparable casts.
contrasts<-list()
for(i in seq_len(nrow(cases))) {
  c<-cases[i,];if(c$mode!="posterior"||c$fraction==1)next
  for(type in c("equation_reference_case","settings_reference_case")) {
    ref<-c[[type]];if(is.na(ref)||ref==c$case_id)next
    a<-file.path(experiment,"runs",sprintf("case_%03d",c$case_id),"checks")
    b<-file.path(experiment,"runs",sprintf("case_%03d",ref),"checks")
    for(dataset in c("validation","reserve"))for(mode in c("conditional","population")) {
      tag<-paste(dataset,mode,sep="_")
      x<-read_table(file.path(a,paste0(tag,"_cast_metrics.csv")))
      y<-read_table(file.path(b,paste0(tag,"_cast_metrics.csv")))
      if(!nrow(x)||!setequal(x$cast_key,y$cast_key))next
      y<-y[match(x$cast_key,y$cast_key),,drop=FALSE]
      out<-data.frame(case_id=c$case_id,reference_case=ref,comparison=type,dataset,mode,matched_casts=nrow(x),
        delta_rmse=sqrt(mean(x$mse))-sqrt(mean(y$mse)),fraction_casts_lower_mse=mean(x$mse<y$mse),
        delta_shape_rmse=sqrt(mean(x$shape_mse))-sqrt(mean(y$shape_mse)),
        delta_peak_depth_mae=mean(x$peak_depth_abs_error-y$peak_depth_abs_error),
        delta_interval_score90=mean(x$interval_score90-y$interval_score90))
      for(kind in c("shape","onset")) {
        xx<-read_table(file.path(a,paste0(tag,"_",kind,"_cast_metrics.csv")))
        yy<-read_table(file.path(b,paste0(tag,"_",kind,"_cast_metrics.csv")))
        if(kind=="onset"){xx<-xx[xx$threshold==.5,,drop=FALSE];yy<-yy[yy$threshold==.5,,drop=FALSE]}
        if(!nrow(xx)||!setequal(xx$cast_key,yy$cast_key))next
        yy<-yy[match(xx$cast_key,yy$cast_key),,drop=FALSE]
        fields<-if(kind=="shape")c("rise_slope_abs_error","rise_width_abs_error")else"onset_abs_error"
        for(field in fields){good<-is.finite(xx[[field]])&is.finite(yy[[field]])
          out[[paste0(field,"_common_pairs")]]<-sum(good)
          out[[paste0("delta_",field,"_same_casts")]]<-if(any(good))mean(xx[[field]][good]-yy[[field]][good])else NA_real_}
      }
      contrasts[[length(contrasts)+1L]]<-out
    }
  }
}
save_csv(fill(contrasts),"planned_paired_contrasts.csv")
cat("Registered cases:",nrow(cases),"; all-campaign registrations:",length(records),"\n")
print(table(fill(current)$status))
