# Evaluate a completed legacy fit on the new common validation casts; never refit.
experiment<-Sys.getenv("EXPERIMENT_DIR")
source(file.path(experiment,"common.R"));source(file.path(experiment,"campaign.R"));source(file.path(experiment,"predictive.R"))
source_run<-Sys.getenv("SOURCE_RUN")
run<-Sys.getenv("EVAL_DIR")
if(!file.exists(file.path(source_run,"FIT_COMPLETE.txt")))stop("Source fit incomplete")
if(dir.exists(run))stop("Evaluation directory already exists")
dir.create(file.path(run,"checks"),recursive=TRUE)
tryCatch({
  config<-readRDS(file.path(source_run,"run_config.rds"))
  case<-config$case
  if(case$mode!="posterior")stop("Prior-only fit cannot be scored as posterior accuracy")
  case$split<-"primary";case$depth_max<-150L
  old<-readRDS(file.path(source_run,"model_data.rds"))$model_data
  expected<-readRDS(data_file(experiment,case))
  fields<-c("cast_key","chl_norm","Region","Ocean","date","enso","depth","log10_distance","MLD_temp_MIX")
  canon<-function(x) {
    x<-as.data.frame(lapply(x[fields],as.character),stringsAsFactors=FALSE)
    x<-x[do.call(order,c(unname(x),list(method="radix"))),,drop=FALSE];rownames(x)<-NULL;x
  }
  if(!identical(canon(old),canon(expected)))stop("Legacy training rows/values differ from the frozen subset; cannot claim a matched evaluation")
  validation<-readRDS(validation_file(experiment,case))
  if(any(old$cast_key %in% validation$cast_key))stop("Training/validation cast overlap")
  fit<-readRDS(file.path(source_run,"fit.rds"))
  diagnostic_file<-file.path(source_run,"diagnostic_summary.csv")
  result<-if(file.exists(diagnostic_file))read.csv(diagnostic_file,stringsAsFactors=FALSE)else case
  result$external_id<-Sys.getenv("EXTERNAL_ID");result$source_run<-source_run
  result$raw_training_identity_verified<-TRUE;result$evaluation_type<-"legacy_fit_common_validation"
  result<-evaluate_predictions(fit,experiment,run,case,result)
  write.csv(result,file.path(run,"test_summary.csv"),row.names=FALSE)
  writeLines("Evaluated without refitting",file.path(run,"CHECKS_COMPLETE.txt"))
},error=function(e){writeLines(conditionMessage(e),file.path(run,"CHECKS_FAILED.txt"));stop(e)})
