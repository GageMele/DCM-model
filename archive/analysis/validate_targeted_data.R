root<-"validation_targeted150/actual"
dir.create(file.path(root,"data"),recursive=TRUE,showWarnings=FALSE)
for(s in c("p10","p20","p50"))file.copy(file.path("validation_baseline150/actual/data",paste0(s,".rds")),file.path(root,"data",paste0(s,".rds")),overwrite=TRUE)
file.copy("validation_baseline150/actual/full_input_250m.rds",file.path(root,"full_input_250m.rds"),overwrite=TRUE)
file.copy("validation_extensive150/actual/data/validation_model.rds",file.path(root,"source_validation.rds"),overwrite=TRUE)
for(f in c("common.R","campaign.R","cases.tsv"))file.copy(file.path("DCM_Ibex_targeted_tests",f),file.path(root,f),overwrite=TRUE)
source("DCM_Ibex_targeted_tests/01_prepare.R")
prepare_campaign(root)
