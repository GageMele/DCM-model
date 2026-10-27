# Keep all cases and all available pairwise cast comparisons, regardless of error.
collect_results<-function(experiment) {
  cases<-read.delim(file.path(experiment,"cases.tsv"),stringsAsFactors=FALSE)
  output<-file.path(experiment,"results");dir.create(output,showWarnings=FALSE)
  fill<-function(xs) {
    xs<-Filter(function(x)is.data.frame(x)&&nrow(x)>0L,xs)
    if(!length(xs))return(data.frame())
    columns<-unique(unlist(lapply(xs,names)))
    do.call(rbind,lapply(xs,function(x){for(n in setdiff(columns,names(x)))x[[n]]<-rep(NA,nrow(x));x[columns]}))
  }
  statuses<-rows<-list();casts<-list()
  read_one<-function(run,key,label,kind) {
    status<-if(file.exists(file.path(run,"CHECKS_COMPLETE.txt")))"complete" else
      if(file.exists(file.path(run,"CHECKS_FAILED.txt")))"checks_failed" else
      if(file.exists(file.path(run,"FIT_FAILED.txt")))"fit_failed" else
      if(file.exists(file.path(run,"FIT_COMPLETE.txt")))"awaiting_checks" else
      if(file.exists(file.path(experiment,"PREPARATION_FAILED.txt")))"blocked_preparation" else"incomplete"
    note_files<-file.path(run,c("CHECKS_FAILED.txt","FIT_FAILED.txt"));note_files<-note_files[file.exists(note_files)]
    note<-paste(unlist(lapply(note_files,readLines,warn=FALSE)),collapse="; ")
    if(status=="blocked_preparation")note<-paste(readLines(file.path(experiment,"PREPARATION_FAILED.txt"),warn=FALSE),collapse="; ")
    statuses[[length(statuses)+1L]]<<-data.frame(run_key=key,label,kind,status,note,run_dir=run)
    candidates<-file.path(run,c("test_summary.csv","diagnostic_summary.csv"));candidates<-candidates[file.exists(candidates)]
    if(length(candidates)) {
      loaded<-lapply(candidates,function(p)tryCatch(read.csv(p,stringsAsFactors=FALSE),error=function(e)NULL))
      loaded<-Filter(function(x)is.data.frame(x)&&nrow(x)==1L,loaded)
      if(!length(loaded)) {
        statuses[[length(statuses)]]$note<<-paste(note,"No readable one-row summary available yet.")
        return(invisible(NULL))
      }
      x<-loaded[[1]];x$run_key<-key;x$result_kind<-kind;x$status<-status
      required<-c("max_rhat","min_bulk_ess","min_tail_ess","divergences","min_ebfmi")
      x$harmonized_convergence_screen<-if(all(required %in% names(x)) && all(is.finite(unlist(x[required]))))
        x$max_rhat<=1.01 && x$min_bulk_ess>=400 && x$min_tail_ess>=400 && x$divergences==0 && x$min_ebfmi>=.3 else NA
      rows[[length(rows)+1L]]<<-x
      cast_file<-file.path(run,"checks","validation_conditional_cast_metrics.csv")
      if(file.exists(cast_file)&&"validation_md5"%in%names(x)&&length(x$validation_md5)==1L&&
        !is.na(x$validation_md5)&&nzchar(x$validation_md5))casts[[key]]<<-list(data=read.csv(cast_file,stringsAsFactors=FALSE),hash=x$validation_md5,
        convergence=x$harmonized_convergence_screen)
    }
  }
  for(i in seq_len(nrow(cases)))read_one(file.path(experiment,"runs",sprintf("case_%02d",cases$case_id[i])),
    sprintf("targeted:%03d",cases$case_id[i]),cases$label[i],"new_fit")
  external_file<-file.path(experiment,"existing_evaluations.tsv")
  if(file.exists(external_file)) {
    ext<-read.delim(external_file,stringsAsFactors=FALSE)
    for(i in seq_len(nrow(ext)))read_one(ext$run_dir[i],ext$external_id[i],ext$external_id[i],"legacy_evaluation")
  }
  status<-fill(statuses);all<-fill(rows)
  write.csv(status,file.path(output,"case_status.csv"),row.names=FALSE)
  write.csv(all,file.path(output,"all_results.csv"),row.names=FALSE)
  if(!nrow(all)){print(status,row.names=FALSE);return(invisible(all))}
  write.csv(all[all$mode=="prior",,drop=FALSE],file.path(output,"prior_checks.csv"),row.names=FALSE)
  write.csv(all[all$mode=="posterior",,drop=FALSE],file.path(output,"all_posterior_comparisons.csv"),row.names=FALSE)
  contrasts<-list();keys<-names(casts)
  if(length(keys)>1L)for(i in seq_len(length(keys)-1L))for(j in (i+1L):length(keys)) {
    a<-casts[[keys[i]]];b<-casts[[keys[j]]]
    if(!identical(a$hash,b$hash))next
    matched<-merge(a$data,b$data,by="cast_key",suffixes=c("_a","_b"))
    if(nrow(matched)!=nrow(a$data)||nrow(matched)!=nrow(b$data))next
    delta<-matched$mse_a-matched$mse_b
    contrasts[[length(contrasts)+1L]]<-data.frame(run_a=keys[i],run_b=keys[j],matched_casts=nrow(matched),
      a_convergence=a$convergence,b_convergence=b$convergence,delta_mean_cast_mse=mean(delta),
      descriptive_paired_se=sd(delta)/sqrt(length(delta)),fraction_casts_a_better=mean(delta<0),
      delta_rmse=sqrt(mean(matched$mse_a))-sqrt(mean(matched$mse_b)),
      delta_shape_rmse=sqrt(mean(matched$shape_mse_a))-sqrt(mean(matched$shape_mse_b)),
      delta_mae=mean(matched$mae_a-matched$mae_b),
      delta_interval_score90=mean(matched$interval_score90_a-matched$interval_score90_b),
      delta_peak_depth_mae=mean(matched$peak_depth_abs_error_a-matched$peak_depth_abs_error_b))
  }
  write.csv(fill(contrasts),file.path(output,"all_paired_cast_comparisons.csv"),row.names=FALSE)
  controls<-c("equation","variance","cast","family","autocor","prior","subset","split","depth_max",
    "k_depth","k_mld","k_distance","warmup","samples","adapt_delta","treedepth","threads","seed")
  controlled<-lapply(seq_len(nrow(cases)),function(i) {
    c<-cases[i,];ref<-cases[cases$case_id==c$reference_case & !is.na(c$reference_case),,drop=FALSE]
    if(nrow(ref)!=1L)return(NULL)
    changed<-controls[vapply(controls,function(n)!identical(c[[n]],ref[[n]]),logical(1))]
    data.frame(case_id=c$case_id,reference_case=ref$case_id,
      changed_settings=paste(vapply(changed,function(n)paste0(n,": ",ref[[n]]," -> ",c[[n]]),character(1)),collapse="; "),
      common_validation=c$split==ref$split,
      interpretation=if(c$split==ref$split)"Compare alongside diagnostics and hardware" else"Different training and validation cohorts; no paired inference")
  })
  write.csv(fill(controlled),file.path(output,"planned_controlled_contrasts.csv"),row.names=FALSE)
  cols<-intersect(c("run_key","equation","variance","cast","subset","status","max_rhat","min_bulk_ess",
    "divergences","treedepth_hits","validation_conditional_rmse","validation_conditional_interval_score90"),names(all))
  print(all[cols],row.names=FALSE)
  invisible(all)
}
if(sys.nframe()==0L) {
  a<-commandArgs(trailingOnly=TRUE)
  collect_results(if(length(a))a[1]else Sys.getenv("EXPERIMENT_DIR"))
}
