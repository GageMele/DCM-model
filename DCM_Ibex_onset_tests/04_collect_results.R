# Preserve every earlier registered run, plus focused onset comparisons for this campaign.
a <- commandArgs(trailingOnly=TRUE)
experiment <- if(length(a))a[1]else Sys.getenv("EXPERIMENT_DIR")
source(file.path(experiment,"collect_all_runs.R"))
source(file.path(experiment,"onset_metrics.R"))
collect_onset <- function(experiment) {
  cases <- read.delim(file.path(experiment,"cases.tsv"),stringsAsFactors=FALSE)
  output <- file.path(experiment,"results");dir.create(output,showWarnings=FALSE)
  casts <- list();aggregates <- list();local <- list();status <- list();contrasts <- list()
  for(i in seq_len(nrow(cases))) {
    c <- cases[i,];run <- file.path(experiment,"runs",sprintf("case_%02d",c$case_id))
    c$status <- status_of(run);status[[i]] <- c
    for(dataset in c("validation","reserve"))for(mode in c("conditional","population")) {
      tag <- paste(dataset,mode,sep="_")
      z <- read_table(file.path(run,"checks",paste0(tag,"_onset_cast_metrics.csv")))
      if(nrow(z)) {
        z$case_id <- c$case_id;z$equation <- c$equation;z$subset <- c$subset
        z$split <- c$split;z$dataset <- dataset;z$prediction_mode <- mode
        casts[[length(casts)+1L]] <- z
        for(threshold in unique(z$threshold)) {
          zz <- z[z$threshold==threshold,,drop=FALSE]
          aggregates[[length(aggregates)+1L]] <- data.frame(case_id=c$case_id,
            equation=c$equation,subset=c$subset,split=c$split,dataset,prediction_mode=mode,
            threshold,onset_scores(zz))
        }
      }
      w <- read_table(file.path(run,"checks",paste0(tag,"_mld_window_cast_metrics.csv")))
      if(nrow(w))local[[length(local)+1L]] <- data.frame(case_id=c$case_id,equation=c$equation,
        subset=c$subset,split=c$split,dataset,prediction_mode=mode,
        casts=nrow(w),mld_neighbourhood_rmse=sqrt(mean(w$mse)))
    }
  }
  cs <- fill_rows(casts)
  write.csv(fill_rows(status),file.path(output,"onset_case_status.csv"),row.names=FALSE)
  write.csv(cs,file.path(output,"onset_all_cast_metrics.csv"),row.names=FALSE)
  write.csv(fill_rows(aggregates),file.path(output,"onset_all_scores.csv"),row.names=FALSE)
  write.csv(fill_rows(local),file.path(output,"mld_neighbourhood_scores.csv"),row.names=FALSE)
  if(nrow(cs))for(i in which(cases$mode=="posterior")) {
    c <- cases[i,];references <- unique(na.omit(c(c$reference_case,
      if(c$execution=="fit" && c$block=="equation")if(c$split=="date_blocked")6L else if(c$subset=="p20")4L else 2L)))
    for(ref in references)for(dataset in c("validation","reserve"))
      for(mode in c("conditional","population"))for(t in c(.25,.5,1)) {
        aa <- cs[cs$case_id==c$case_id & cs$dataset==dataset & cs$prediction_mode==mode & cs$threshold==t,,drop=FALSE]
        bb <- cs[cs$case_id==ref & cs$dataset==dataset & cs$prediction_mode==mode & cs$threshold==t,,drop=FALSE]
        if(!nrow(aa)||!nrow(bb))next
        bb <- bb[match(aa$cast_key,bb$cast_key),,drop=FALSE]
        if(anyNA(bb$cast_key)||!identical(as.character(aa$cast_key),as.character(bb$cast_key)))next
        # Fair onset-MAE contrasts use the SAME casts detected by all three curves.
        eligible <- is.finite(aa$onset_abs_error)&is.finite(bb$onset_abs_error)
        avg <- function(v)if(length(v))mean(v)else NA_real_
        contrasts[[length(contrasts)+1L]] <- data.frame(case_id=c$case_id,reference_case=ref,
          dataset,prediction_mode=mode,threshold=t,matched_casts=nrow(aa),
          jointly_detected_casts=sum(eligible),
          delta_onset_mae_common_casts=avg(aa$onset_abs_error[eligible]-bb$onset_abs_error[eligible]),
          candidate_false_positives=sum(aa$false_positive,na.rm=TRUE),
          reference_false_positives=sum(bb$false_positive,na.rm=TRUE),
          candidate_false_negatives=sum(aa$false_negative,na.rm=TRUE),
          reference_false_negatives=sum(bb$false_negative,na.rm=TRUE))
      }
  }
  write.csv(fill_rows(contrasts),file.path(output,"onset_planned_paired_contrasts.csv"),row.names=FALSE)
  cat("Onset scores retain missing/censored detections explicitly; no automatic single-score winner.\n")
}
if(sys.nframe()==0L){collect_results(experiment);collect_onset(experiment)}
