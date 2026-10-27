root <- normalizePath(".",winslash="/")
package <- file.path(root,"DCM_Ibex_extensive_tests")
for(f in list.files(package,"\\.R$",full.names=TRUE)) parse(f)
source(file.path(package,"campaign.R"))
source(file.path(package,"predictive.R"))
cases <- campaign_cases()
write.table(cases,file.path(package,"cases.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
stopifnot(sum(cases$cpus)==808L, sum(cases$subset=="p50")==4L,
  all(cases$cpus==cases$chains*cases$threads), all(cases$iter-cases$warmup==cases$samples),
  !anyDuplicated(cases[setdiff(names(cases),c("case_id","label","block"))]))
for(eq in c("additive","both_scaled")) {
  x <- cases[cases$equation==eq & cases$block=="warmup",]
  stopifnot(all(x$samples==1000L),identical(sort(x$warmup),c(250L,500L,2000L)))
  x <- cases[cases$equation==eq & cases$block=="draws",]
  stopifnot(all(x$warmup==1000L))
  x <- cases[cases$equation==eq & cases$block=="chains_same_total",]
  stopifnot(x$chains*x$samples==4000L)
}
# Metrics are checked against independently known values; unequal profile lengths
# must not give long casts more weight.
tab <- data.frame(cast_key=c("a","a","b","b","b"),Region="R",
  depth=c(10,20,10,20,30),chl_norm=c(1,3,2,2,2),predicted=c(2,4,0,0,0),
  lower=-10,upper=10)
metrics <- do.call(rbind,lapply(split(tab,tab$cast_key),profile_metrics))
scores <- score_casts(metrics)
stopifnot(abs(scores$rmse-sqrt(2.5))<1e-12,scores$mae==1.5,scores$shape_rmse==0,
          scores$coverage90==1)
experiment <- file.path(root,"validation_extensive150","actual")
dir.create(file.path(experiment,"data"),recursive=TRUE,showWarnings=FALSE)
file.copy(list.files(package,"\\.R$|cases.tsv$",full.names=TRUE),experiment,overwrite=TRUE)
parent <- file.path(root,"validation_baseline150","actual")
stopifnot(file.copy(file.path(parent,"full_input_250m.rds"),experiment,overwrite=TRUE))
for(s in c("p10","p20","p50")) stopifnot(file.copy(file.path(parent,"data",paste0(s,".rds")),
  file.path(experiment,"data"),overwrite=TRUE))
source(file.path(package,"01_prepare.R"))
if(!file.exists(file.path(experiment,"PREPARATION_COMPLETE.txt"))) print(prepare_campaign(experiment))
source(file.path(package,"common.R"))
for(s in c("p10","p20","p50")) stopifnot(
  unname(tools::md5sum(file.path(parent,"data",paste0(s,".rds"))))==
  unname(tools::md5sum(file.path(experiment,"data",paste0(s,".rds")))))
validation <- readRDS(file.path(experiment,"data","validation_model.rds"))
train <- readRDS(file.path(experiment,"data","p50_model.rds"))
stopifnot(!any(validation$cast_key %in% train$cast_key),length(unique(validation$cast_key))==108L)
# A collector with zero completed cases must succeed (the earlier report bug).
source(file.path(package,"04_collect_results.R"))
invisible(collect_results(experiment))
file.copy(file.path(experiment,c("basis_counts.csv","split_counts.csv","preflight_session.txt")),
          package,overwrite=TRUE)
cat("102 cases: R syntax, matched factors, actual-data Stan construction, unchanged source subsets, disjoint validation and empty collector PASS.\n")
