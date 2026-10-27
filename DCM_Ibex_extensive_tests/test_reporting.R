args <- commandArgs(TRUE)
package <- if(length(args)) args[1] else Sys.getenv("EXPERIMENT_DIR")
stopifnot(nzchar(package))
source(file.path(package,"04_collect_results.R"))
source(file.path(package,"campaign.R"))
cases <- campaign_cases()
fixture <- tempfile("report-fixture-",tmpdir=tempdir())
dir.create(fixture)
write.table(cases,file.path(fixture,"cases.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
add_result <- function(id,complete=TRUE,predictions=FALSE) {
  r <- cases[cases$case_id==id,]
  r$cpu_model <- "fixture_cpu"; r$fit_seconds <- id+100
  r$requested_core_hours <- 1; r$min_bulk_ess <- 700
  r$max_rhat <- 1.005; r$divergences <- 0L; r$treedepth_hit_fraction <- .001
  r$convergence_screen_pass <- TRUE
  if(predictions) {
    r$validation_conditional_rmse <- id/100+1
    r$validation_conditional_shape_rmse <- id/200+1
    r$validation_population_rmse <- 2
    r$validation_conditional_coverage90 <- .9
  }
  run <- file.path(fixture,"runs",sprintf("case_%02d",id))
  dir.create(file.path(run,"checks"),recursive=TRUE)
  write.csv(r,file.path(run,if(complete) "test_summary.csv" else "diagnostic_summary.csv"),row.names=FALSE)
  if(complete) writeLines("complete",file.path(run,"CHECKS_COMPLETE.txt")) else
    writeLines("fixture check failure",file.path(run,"CHECKS_FAILED.txt"))
  if(predictions) write.csv(data.frame(cast_key=c("a","b"),mse=c(1,4)+id/100,
    shape_mse=c(.2,.4)),file.path(run,"checks","validation_conditional_cast_metrics.csv"),row.names=FALSE)
}
add_result(97L)
invisible(capture.output(x <- collect_results(fixture)))
stopifnot(nrow(x)==1L,x$mode=="prior",!file.exists(file.path(fixture,"results","predictive_ranking.csv")))
add_result(4L, predictions=TRUE)
add_result(6L, predictions=TRUE)
add_result(50L, predictions=TRUE)
add_result(49L, complete=FALSE)
invisible(capture.output(x <- collect_results(fixture)))
ranked <- read.csv(file.path(fixture,"results","predictive_ranking.csv"))
contrasts <- read.csv(file.path(fixture,"results","matched_contrasts.csv"))
paired <- read.csv(file.path(fixture,"results","paired_vs_original_k15.csv"))
stopifnot(nrow(x)==5L,nrow(ranked)==3L,!any(ranked$mode=="prior"),
  !any(ranked$status!="complete"),50L %in% contrasts$case_id,
  all(contrasts$same_cpu_model),nrow(paired)==3L,
  abs(paired$delta_rmse[paired$case_id==4])<1e-12)
cat("Collector tests PASS: prior-only, mixed incomplete/complete, different summary columns, matched contrasts and paired cast results.\n")
