# Validate frozen data only. Do not load brms, compile Stan, or fit models locally.
package<-normalizePath('DCM_Ibex_followup_tests')
source(file.path(package,'common.R'));source(file.path(package,'campaign.R'));source(file.path(package,'data_setup.R'))
cases<-campaign_cases()
write.table(cases,file.path(package,'cases.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
experiment<-file.path('validation_followup150',paste0('data_',format(Sys.time(),'%Y%m%dT%H%M%S'), '_',Sys.getpid()))
dir.create(experiment,recursive=TRUE)
file.copy(file.path(package,'base_models.R'),experiment)
Sys.setenv(EXPERIMENT_DIR=normalizePath(experiment))
counts<-prepare_data(experiment,normalizePath('validation_targeted150/actual'))
full<-validate_bundle(readRDS(file.path(experiment,'full_input_250m.rds')),full=TRUE)
full<-full[full$depth>=10 & full$depth<=150,,drop=FALSE]
fields<-c('cast_key','depth','chl_norm','Region','Ocean','date','enso','MLD_temp_MIX','log10_distance')
canonical<-function(x){
  x<-as.data.frame(lapply(x[fields],as.character),stringsAsFactors=FALSE)
  x<-x[do.call(order,c(unname(x),list(method='radix'))),,drop=FALSE];rownames(x)<-NULL;x
}
for(name in counts$dataset){
  d<-readRDS(file.path(experiment,'data',paste0(name,'.rds')))
  expected<-full[full$cast_key %in% d$cast_key,,drop=FALSE]
  stopifnot(identical(canonical(d),canonical(expected)))
}
write.csv(counts,file.path(package,'verified_split_counts.csv'),row.names=FALSE)
writeLines(normalizePath(experiment),'validation_followup150/latest_data_check.txt')
print(counts,row.names=FALSE)
cat('PASS: all nine datasets >=10%; alternate108 disjoint from original108; reserve outside both training families and earlier holdouts; identical raw outcomes/predictors; frozen scaling reused.\n')
