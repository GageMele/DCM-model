#!/usr/bin/env Rscript
# Safe with zero completed fits and with fits whose predictive checks failed.
collect_results <- function(experiment) {
  cases <- read.delim(file.path(experiment,"cases.tsv"),stringsAsFactors=FALSE)
  output <- file.path(experiment,"results"); dir.create(output,showWarnings=FALSE)
  statuses <- summaries <- list()
  for(i in seq_len(nrow(cases))) {
    run <- file.path(experiment,"runs",sprintf("case_%02d",cases$case_id[i]))
    status <- if(file.exists(file.path(run,"CHECKS_COMPLETE.txt"))) "complete" else
      if(file.exists(file.path(run,"CHECKS_FAILED.txt"))) "checks_failed" else
      if(file.exists(file.path(run,"FIT_FAILED.txt"))) "fit_failed" else
      if(file.exists(file.path(run,"FIT_COMPLETE.txt"))) "awaiting_checks" else "incomplete"
    notes <- c(file.path(run,"FIT_FAILED.txt"),file.path(run,"CHECKS_FAILED.txt"))
    notes <- notes[file.exists(notes)]
    note <- if(length(notes)) paste(unlist(lapply(notes,readLines,warn=FALSE)),collapse="; ") else ""
    statuses[[i]] <- data.frame(case_id=cases$case_id[i],label=cases$label[i],status,note)
    candidates <- file.path(run,c("test_summary.csv","diagnostic_summary.csv"))
    candidates <- candidates[file.exists(candidates)]
    if(length(candidates)) {
      row <- read.csv(candidates[1],stringsAsFactors=FALSE)
      row$status <- status; summaries[[length(summaries)+1L]] <- row
    }
  }
  status <- do.call(rbind,statuses)
  write.csv(status,file.path(output,"case_status.csv"),row.names=FALSE)
  if(!length(summaries)) {print(status); message("No completed diagnostic summaries yet."); return(invisible(status))}
  columns <- unique(unlist(lapply(summaries,names)))
  summaries <- lapply(summaries,function(x) {
    for(n in setdiff(columns,names(x))) x[[n]] <- rep(NA,nrow(x)); x[columns]
  })
  all <- do.call(rbind,summaries)
  write.csv(all,file.path(output,"all_results.csv"),row.names=FALSE)
  write.csv(all[all$mode=="prior",,drop=FALSE],file.path(output,"prior_checks.csv"),row.names=FALSE)
  ranked <- all[all$mode=="posterior" & all$status=="complete",,drop=FALSE]
  if(nrow(ranked) && "validation_conditional_rmse" %in% names(ranked)) {
    ranked <- ranked[order(!ranked$convergence_screen_pass,ranked$validation_conditional_rmse,
                           ranked$validation_conditional_shape_rmse,ranked$requested_core_hours,na.last=TRUE),]
    write.csv(ranked,file.path(output,"predictive_ranking.csv"),row.names=FALSE)
  }
  # Every case is compared with the new original-formula k=15 case on exactly
  # the same validation casts. Paired SE is descriptive, not a selection-proof test.
  ref <- cases$case_id[cases$block=="equation_k" & cases$equation=="original" & cases$k_depth==15L]
  cast_file <- function(id) file.path(experiment,"runs",sprintf("case_%02d",id),
                                     "checks","validation_conditional_cast_metrics.csv")
  paired <- list()
  if(length(ref)==1L && file.exists(cast_file(ref))) {
    base <- read.csv(cast_file(ref),stringsAsFactors=FALSE)
    for(id in all$case_id[all$mode=="posterior"]) if(file.exists(cast_file(id))) {
      x <- read.csv(cast_file(id),stringsAsFactors=FALSE)
      matched <- merge(x,base,by="cast_key",suffixes=c("","_reference"))
      if(nrow(matched)!=nrow(base)) next
      delta <- matched$mse-matched$mse_reference
      paired[[length(paired)+1L]] <- data.frame(case_id=id,reference_case=ref,
        matched_casts=nrow(matched),delta_mean_cast_mse=mean(delta),
        descriptive_paired_se=sd(delta)/sqrt(length(delta)),
        delta_rmse=sqrt(mean(matched$mse))-sqrt(mean(matched$mse_reference)),
        delta_shape_rmse=sqrt(mean(matched$shape_mse))-sqrt(mean(matched$shape_mse_reference)))
    }
  }
  if(length(paired)) write.csv(do.call(rbind,paired),file.path(output,"paired_vs_original_k15.csv"),row.names=FALSE)
  # Factor-specific matched contrasts. Only one requested setting (or an explicitly
  # labelled package of settings) changes, including equal-total-draw chain tests.
  contrasts <- list()
  for(i in seq_len(nrow(all))) {
    row <- all[i,]; if(row$mode!="posterior") next
    pool <- cases[cases$block %in% c("equation_k","scaling") &
      cases$equation==row$equation & cases$k_depth==7L,,drop=FALSE]
    if(nrow(pool)!=1L || row$case_id==pool$case_id) next
    reference <- all[all$case_id==pool$case_id,,drop=FALSE]
    if(nrow(reference)!=1L) next
    same_cpu <- identical(row$cpu_model,reference$cpu_model) && row$cpu_model!="unknown"
    record <- data.frame(case_id=row$case_id,reference_case=reference$case_id,block=row$block,
      same_cpu_model=same_cpu,fit_time_ratio=row$fit_seconds/reference$fit_seconds,
      core_hour_ratio=row$requested_core_hours/reference$requested_core_hours,
      bulk_ess_ratio=row$min_bulk_ess/reference$min_bulk_ess,
      reference_converged=reference$convergence_screen_pass,case_converged=row$convergence_screen_pass)
    for(n in c("validation_conditional_rmse","validation_conditional_shape_rmse",
               "validation_population_rmse","validation_conditional_coverage90"))
      if(n %in% names(all)) record[[paste0("delta_",n)]] <- row[[n]]-reference[[n]]
    contrasts[[length(contrasts)+1L]] <- record
  }
  if(length(contrasts)) write.csv(do.call(rbind,contrasts),file.path(output,"matched_contrasts.csv"),row.names=FALSE)
  cols <- intersect(c("case_id","equation","subset","status","fit_seconds","max_rhat",
    "min_bulk_ess","divergences","treedepth_hit_fraction","convergence_screen_pass",
    "validation_conditional_rmse","validation_conditional_shape_rmse"),names(all))
  print(all[cols],row.names=FALSE)
  invisible(all)
}
if(sys.nframe()==0L) collect_results(Sys.getenv("EXPERIMENT_DIR"))
