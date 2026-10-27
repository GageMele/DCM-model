# Compare every recorded fit; legacy prediction evaluations are merged, not refitted.
read_table<-function(path,tsv=FALSE) {
  if(!file.exists(path))return(data.frame())
  tryCatch(if(tsv)read.delim(path,stringsAsFactors=FALSE)else read.csv(path,stringsAsFactors=FALSE),error=function(e)data.frame())
}
fill_rows<-function(xs) {
  xs<-Filter(function(x)is.data.frame(x)&&nrow(x)>0L,xs)
  if(!length(xs))return(data.frame())
  cols<-unique(unlist(lapply(xs,names)))
  do.call(rbind,lapply(xs,function(x){for(n in setdiff(cols,names(x)))x[[n]]<-rep(NA,nrow(x));x[cols]}))
}
screen<-function(x) {
  ns<-c("max_rhat","min_bulk_ess","min_tail_ess","divergences","min_ebfmi")
  if(!all(ns %in% names(x)))return(NA)
  z<-suppressWarnings(as.numeric(unlist(x[ns])));if(any(!is.finite(z)))return(NA)
  z[1]<=1.01&&z[2]>=400&&z[3]>=400&&z[4]==0&&z[5]>=.3
}
safe_summary<-function(run) {
  for(n in c("test_summary.csv","diagnostic_summary.csv")){
    x<-read_table(file.path(run,n));if(nrow(x)==1L)return(x)
  };data.frame()
}
run_key<-function(campaign,id) {
  value<-suppressWarnings(as.numeric(id));suffix<-if(length(value)==1L&&is.finite(value))sprintf("%03d",as.integer(value))else as.character(id)
  paste0(basename(campaign),":",suffix)
}
status_of<-function(run) {
  if(!file.exists(file.path(run,"FIT_COMPLETE.txt"))&&file.exists(file.path(dirname(dirname(run)),"PREPARATION_FAILED.txt")))return("blocked_preparation")
  if(file.exists(file.path(run,"CHECKS_COMPLETE.txt")))"complete" else
  if(file.exists(file.path(run,"CHECKS_FAILED.txt")))"checks_failed" else
  if(file.exists(file.path(run,"FIT_FAILED.txt")))"fit_failed" else
  if(file.exists(file.path(run,"FIT_COMPLETE.txt")))"awaiting_checks" else"incomplete"
}
paired_metrics<-function(a,b) {
  delta<-a$mse-b$mse
  data.frame(matched_casts=nrow(a),delta_mean_cast_mse=mean(delta),
    descriptive_cast_se=sd(delta)/sqrt(length(delta)),fraction_casts_a_better=mean(delta<0),
    delta_rmse=sqrt(mean(a$mse))-sqrt(mean(b$mse)),
    delta_shape_rmse=sqrt(mean(a$shape_mse))-sqrt(mean(b$shape_mse)),
    delta_mae=mean(a$mae-b$mae),
    delta_peak_depth_mae=mean(a$peak_depth_abs_error-b$peak_depth_abs_error),
    delta_interval_score90=if(all(c("interval_score90") %in% names(a))&&"interval_score90"%in%names(b))mean(a$interval_score90-b$interval_score90)else NA_real_)
}
cluster_bootstrap<-function(a,b,metadata,B=500L) {
  meta<-metadata[match(a$cast_key,metadata$cast_key),,drop=FALSE]
  if(nrow(meta)!=nrow(a)||anyNA(meta[c("Ocean","Region","date")]))return(data.frame())
  groups<-split(seq_len(nrow(a)),paste(meta$Ocean,meta$Region,meta$date,sep="\r"))
  strata<-split(names(groups),vapply(groups,function(ii)as.character(meta$Region[ii[1]]),character(1)))
  delta<-a$mse-b$mse;set.seed(20261018L)
  boot<-replicate(B,{
    selected<-unlist(lapply(strata,function(keys)keys[sample.int(length(keys),length(keys),replace=TRUE)]),use.names=FALSE)
    mean(delta[unlist(groups[selected],use.names=FALSE)])
  })
  data.frame(nested_date_groups=length(groups),bootstrap_draws=B,
    bootstrap_low_mse_difference=unname(quantile(boot,.025)),bootstrap_high_mse_difference=unname(quantile(boot,.975)),
    interpretation="Exploratory region-stratified nested-date cluster bootstrap; reused development holdouts; not a final selection-adjusted interval")
}
collect_results<-function(experiment,root=dirname(dirname(experiment))) {
  output<-file.path(experiment,"results");dir.create(output,showWarnings=FALSE)
  manifests<-c(Sys.glob(file.path(root,"experiments","*","manifest.tsv")),Sys.glob(file.path(root,"benchmarks","*","manifest.tsv")))
  current_cases<-read_table(file.path(experiment,"cases.tsv"),TRUE)
  # Locate completed common evaluations from every older campaign.
  evaluations<-fill_rows(lapply(Sys.glob(file.path(root,"experiments","*","existing_evaluations.tsv")),function(p)read_table(p,TRUE)))
  records<-list();scores<-list();statuses<-list();tables<-list();verified<-0L;max_difference<-0
  for(manifest in manifests) {
    campaign<-dirname(manifest);registry<-read_table(manifest,TRUE)
    cases<-read_table(file.path(campaign,"cases.tsv"),TRUE)
    for(i in seq_len(nrow(registry))) {
      item<-registry[i,,drop=FALSE];id<-if("case_id"%in%names(item))item$case_id else item$case
      if(!length(id)||!"run_dir"%in%names(item))next
      run<-as.character(item$run_dir);key<-run_key(campaign,id)
      c<-if("case_id"%in%names(cases))cases[cases$case_id==id,,drop=FALSE]else data.frame()
      base<-if(nrow(c)==1L)c else item
      base$run_key<-key;base$campaign<-basename(campaign);base$run_dir<-run
      base$status<-status_of(run);base$fit_completed<-file.exists(file.path(run,"FIT_COMPLETE.txt"))
      statuses[[key]]<-base
      x<-safe_summary(run);metric_run<-run
      if(nrow(x)==1L)for(n in names(x))base[[n]]<-x[[n]]
      if(nrow(evaluations)&&all(c("source_run","run_dir")%in%names(evaluations))) {
        candidates<-evaluations[evaluations$source_run==run,,drop=FALSE]
        usable<-lapply(seq_len(nrow(candidates)),function(j){
          p<-as.character(candidates$run_dir[j]);z<-safe_summary(p)
          if(!file.exists(file.path(p,"CHECKS_COMPLETE.txt")))return(NULL)
          if(nrow(z)!=1L||!"raw_training_identity_verified"%in%names(z)||!isTRUE(z$raw_training_identity_verified))return(NULL)
          list(path=p,summary=z,time=as.numeric(file.info(file.path(p,"test_summary.csv"))$mtime))
        });usable<-Filter(Negate(is.null),usable)
        if(length(usable)) {
          chosen<-usable[[which.max(vapply(usable,function(z)z$time,numeric(1)))]]
          for(n in grep("^(validation_|training_|prediction_|raw_training_|evaluation_type)",names(chosen$summary),value=TRUE))base[[n]]<-chosen$summary[[n]]
          metric_run<-chosen$path;base$common_evaluation_run<-metric_run
        }
      }
      base$run_key<-key;base$status<-status_of(run);base$harmonized_convergence_screen<-screen(base)
      records[[key]]<-base
      for(dataset in c("training","validation","reserve"))for(mode in c("conditional","population")) {
        tag<-paste(dataset,mode,sep="_")
        cast<-read_table(file.path(metric_run,"checks",paste0(tag,"_cast_metrics.csv")))
        if(nrow(cast)) {
          if(anyDuplicated(cast$cast_key))stop("Duplicate cast in ",key," ",tag)
          reported<-base[[paste0(tag,"_rmse")]]
          if(length(reported)==1L&&is.finite(reported)) {
            difference<-abs(sqrt(mean(cast$mse))-reported)
            if(difference>1e-8)stop("RMSE audit failed: ",key," ",tag)
            verified<-verified+1L;max_difference<-max(max_difference,difference)
          }
          hash<-base[[if(dataset=="reserve")"reserve_md5"else"validation_md5"]]
          if(dataset!="training"&&mode=="conditional"&&length(hash)==1L&&!is.na(hash)&&nzchar(hash)&&
             "mode"%in%names(base)&&base$mode=="posterior") {
            cast<-cast[order(cast$cast_key,method="radix"),,drop=FALSE]
            scores[[paste(key,dataset,sep="|")]]<-list(key=key,dataset=dataset,hash=hash,data=cast,convergence=base$harmonized_convergence_screen)
          }
        }
        for(suffix in c("cast_metrics","by_region","strata")) {
          file<-paste0(tag,"_",suffix,".csv");tab<-if(suffix=="cast_metrics")cast else read_table(file.path(metric_run,"checks",file))
          if(nrow(tab)){tab$run_key<-key;tables[[file]]<-c(tables[[file]],list(tab))}
        }
      }
      for(file in c("predictive_batch_sensitivity.csv","prior_predictive_quantiles.csv")){
        tab<-read_table(file.path(metric_run,"checks",file));if(nrow(tab)){tab$run_key<-key;tables[[file]]<-c(tables[[file]],list(tab))}
      }
    }
  }
  # Retain the original unmanifested pilot as a status/configuration reference.
  for(p in Sys.glob(file.path(root,"runs","*","run_config.rds"))) {
    run<-dirname(p);key<-paste0("pilot:",basename(run))
    if(any(vapply(records,function(z)identical(z$run_dir,run),logical(1))))next
    z<-data.frame(run_key=key,campaign="pilot",run_dir=run,status=status_of(run),fit_completed=file.exists(file.path(run,"FIT_COMPLETE.txt")))
    x<-safe_summary(run);if(nrow(x)==1L)for(n in names(x))z[[n]]<-x[[n]]
    z$harmonized_convergence_screen<-screen(z);records[[key]]<-z;statuses[[key]]<-z
  }
  for(i in seq_len(nrow(current_cases))) {
    key<-run_key(experiment,current_cases$case_id[i])
    if(key %in% names(records))next
    z<-current_cases[i,,drop=FALSE];z$run_key<-key;z$campaign<-basename(experiment)
    z$run_dir<-file.path(experiment,"runs",sprintf("case_%02d",z$case_id));z$status<-"not_registered_yet"
    z$fit_completed<-FALSE;z$harmonized_convergence_screen<-NA
    records[[key]]<-z;statuses[[key]]<-z
  }
  all<-fill_rows(records);write.csv(all,file.path(output,"all_results.csv"),row.names=FALSE)
  write.csv(fill_rows(statuses),file.path(output,"case_status.csv"),row.names=FALSE)
  for(n in names(tables))write.csv(fill_rows(tables[[n]]),file.path(output,paste0("all_",n)),row.names=FALSE)
  current<-all[all$campaign==basename(experiment),,drop=FALSE]
  write.csv(current,file.path(output,"followup_results.csv"),row.names=FALSE)
  if("mode"%in%names(all))write.csv(all[!is.na(all$mode)&all$mode=="prior",,drop=FALSE],file.path(output,"prior_checks.csv"),row.names=FALSE)
  pairs<-list();keys<-names(scores)
  if(length(keys)>1L)for(i in seq_len(length(keys)-1L))for(j in seq.int(i+1L,length(keys))) {
    a<-scores[[i]];b<-scores[[j]]
    if(a$dataset!=b$dataset||a$hash!=b$hash||!identical(as.character(a$data$cast_key),as.character(b$data$cast_key)))next
    pairs[[length(pairs)+1L]]<-data.frame(run_a=a$key,run_b=b$key,dataset=a$dataset,
      a_convergence=a$convergence,b_convergence=b$convergence,paired_metrics(a$data,b$data))
  }
  write.csv(fill_rows(pairs),file.path(output,"all_matched_cast_comparisons.csv"),row.names=FALSE)
  # Planned comparisons add date-cluster uncertainty and regional differences.
  planned<-list();regional<-list();metadata<-read_table(file.path(experiment,"split_membership.csv"))
  if(nrow(metadata))metadata<-metadata[!duplicated(metadata$cast_key),,drop=FALSE]
  controls<-c("equation","subset","split","variance","sigma_smooth_sd","k_depth","k_mld","k_distance","seed","warmup","samples","adapt_delta","treedepth","threads")
  for(i in seq_len(nrow(current_cases))) {
    case<-current_cases[i,];reference<-case$reference_case
    if(is.na(reference)||case$mode!="posterior")next
    ref<-current_cases[current_cases$case_id==reference,,drop=FALSE]
    changed<-controls[vapply(controls,function(n)!identical(case[[n]],ref[[n]]),logical(1))]
    for(dataset in c("validation","reserve")) {
      a<-scores[[paste(run_key(experiment,case$case_id),dataset,sep="|")]]
      b<-scores[[paste(run_key(experiment,reference),dataset,sep="|")]]
      common<-!is.null(a)&&!is.null(b)&&a$hash==b$hash&&identical(as.character(a$data$cast_key),as.character(b$data$cast_key))
      z<-data.frame(case_id=case$case_id,reference_case=reference,dataset,
        changed_settings=paste(vapply(changed,function(n)paste0(n,": ",ref[[n]]," -> ",case[[n]]),character(1)),collapse="; "),common_cohort=common)
      if(common) {
        z<-cbind(z,paired_metrics(a$data,b$data));boot<-cluster_bootstrap(a$data,b$data,metadata)
        if(nrow(boot))z<-cbind(z,boot)
        for(region in unique(a$data$Region)){
          ii<-a$data$Region==region
          regional[[length(regional)+1L]]<-data.frame(case_id=case$case_id,reference_case=reference,dataset,Region=region,paired_metrics(a$data[ii,],b$data[ii,]))
        }
      }
      planned[[length(planned)+1L]]<-z
    }
  }
  write.csv(fill_rows(planned),file.path(output,"planned_controlled_contrasts.csv"),row.names=FALSE)
  write.csv(fill_rows(regional),file.path(output,"planned_regional_contrasts.csv"),row.names=FALSE)
  write.csv(data.frame(verified_rmse_aggregates=verified,maximum_difference=max_difference),file.path(output,"score_audit.csv"),row.names=FALSE)
  columns<-intersect(c("case_id","label","status","max_rhat","min_bulk_ess","divergences","treedepth_hits","validation_conditional_rmse","reserve_conditional_rmse"),names(current))
  print(current[columns],row.names=FALSE)
  cat("All registered runs:",nrow(all),"; recomputed RMSE aggregates:",verified,"; maximum discrepancy:",max_difference,"\n")
  invisible(all)
}
if(sys.nframe()==0L){a<-commandArgs(trailingOnly=TRUE);collect_results(a[1])}
