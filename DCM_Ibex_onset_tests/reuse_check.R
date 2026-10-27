# Read an existing fit and write new predictions/metrics only to this campaign.
experiment <- Sys.getenv("EXPERIMENT_DIR")
source(file.path(experiment,"common.R"));source(file.path(experiment,"campaign.R"))
if(!file.exists(file.path(experiment,"PREPARATION_COMPLETE.txt")))stop("Preparation not complete")
case <- read_case(experiment,Sys.getenv("CASE_ID"))
stopifnot(case$execution=="reuse",case$mode=="posterior")
registry <- read.delim(file.path(experiment,"source_controls.tsv"),stringsAsFactors=FALSE)
item <- registry[registry$case_id==case$case_id,,drop=FALSE];stopifnot(nrow(item)==1L)
run <- file.path(experiment,"runs",sprintf("case_%02d",case$case_id))
dir.create(run,recursive=TRUE,showWarnings=FALSE)
if(file.exists(file.path(run,"CHECKS_COMPLETE.txt")))quit(status=0L)
tryCatch({
  if(!file.exists(file.path(item$source_run,"FIT_COMPLETE.txt")))
    stop("Source fitting job ended without a complete fit; it was not resubmitted")
  fit <- readRDS(file.path(item$source_run,"fit.rds"))
  old <- readRDS(file.path(item$source_run,"model_data.rds"))$model_data
  current <- readRDS(data_file(experiment,case))
  raw <- c("cast_key","Region","Ocean","date","enso","depth","log10_distance","MLD_temp_MIX","chl_norm")
  stopifnot(identical(old[raw],current[raw]),all(names(old) %in% names(current)),
            identical(old,current[names(old)]))
  # Hydrated saved fits are used directly; no brm/update/compilation call occurs here.
  result <- data.frame(case,source_run=item$source_run,source_fit_job=item$fit_job,
    evaluation_only=TRUE,raw_training_identity_verified=TRUE)
  for(name in c("test_summary.csv","diagnostic_summary.csv")) {
    path <- file.path(item$source_run,name)
    if(file.exists(path)) {
      x <- read.csv(path,stringsAsFactors=FALSE)
      keep <- intersect(names(x),c("max_rhat","min_bulk_ess","min_tail_ess","divergences",
        "treedepth_hits","min_ebfmi","fit_seconds","requested_core_hours",
        "diagnostic_screen_pass","convergence_screen_pass","cpu_model"))
      for(n in keep)result[[n]] <- x[[n]]
      break
    }
  }
  source(file.path(experiment,"predictive.R"))
  result <- evaluate_predictions(fit,experiment,run,case,result)
  write.csv(result,file.path(run,"test_summary.csv"),row.names=FALSE)
  writeLines(item$source_run,file.path(run,"REUSED_FIT.txt"))
  writeLines("Source fit reused without modification",file.path(run,"FIT_COMPLETE.txt"))
  writeLines("Complete",file.path(run,"CHECKS_COMPLETE.txt"))
},error=function(e){writeLines(conditionMessage(e),file.path(run,"CHECKS_FAILED.txt"));stop(e)})
