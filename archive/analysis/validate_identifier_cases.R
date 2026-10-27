# Meaningful checks for metadata aliases, unknown locations, and ambiguity guards.
source('DCM_Ibex_workflow/01_prepare_data.R')
root <- 'validation_metadata_cases'
dir.create(root,showWarnings=FALSE)
ctd <- do.call(rbind,lapply(c('Central RS','Hawiian'),function(r)
  data.frame(LocationCode=if(r=='Central RS') NA_character_ else 'HAW',
    DeepCTDID='cast1',Region.x=r,UTCDate.x='2011-01-01',Depth_m=c(5,10,15,20,25),
    Temp_DegC=27,Chl=c(1,2,3,7,1),distance=5)))
mld <- data.frame(LocationCode=c('Not specified','HAW'),Region=c('Central RS','Hawaiian'),
                  cast_id='cast1',Ocean=c('Red Sea','Pacific'),MLD_temp_MIX=c(20,30))
write.csv(ctd,file.path(root,'ctd.csv'),row.names=FALSE)
write.csv(mld,file.path(root,'mld.csv'),row.names=FALSE)
Sys.setenv(CTD_INPUT=file.path(root,'ctd.csv'),MLD_INPUT=file.path(root,'mld.csv'),
           MODEL_DATA=file.path(root,'pass.rds'),CAST_FRACTION='1')
b <- prepare_model_data()
stopifnot(nrow(b$model_data)==8,all(is.finite(b$model_data$MLD_temp_MIX)),
          all(is.na(b$model_data$LocationCode[b$model_data$Region=='Central RS'])),
          all(b$model_data$MLD_temp_MIX[b$model_data$Region=='Hawiian']==30))
# Conflicting MLDs remain prohibited even when location is unknown.
conflict <- mld[1,];conflict$MLD_temp_MIX <- 99
write.csv(rbind(mld,conflict),file.path(root,'mld.csv'),row.names=FALSE)
Sys.setenv(MODEL_DATA=file.path(root,'conflict.rds'))
err <- try(prepare_model_data(),silent=TRUE)
stopifnot(inherits(err,'try-error'),grepl('Multiple MLD',err))
# Four join fields spanning dates still cannot silently share one MLD.
write.csv(mld,file.path(root,'mld.csv'),row.names=FALSE)
extra <- ctd[1,];extra$UTCDate.x <- '2011-01-02'
write.csv(rbind(ctd,extra),file.path(root,'ctd.csv'),row.names=FALSE)
Sys.setenv(MODEL_DATA=file.path(root,'ambiguous.rds'))
err <- try(prepare_model_data(),silent=TRUE)
stopifnot(inherits(err,'try-error'),grepl('spans multiple dates',err))
# Missing cast IDs must produce a precise error and cannot act as matching keys.
ctd$DeepCTDID[1] <- NA_character_
write.csv(ctd,file.path(root,'ctd.csv'),row.names=FALSE)
Sys.setenv(MODEL_DATA=file.path(root,'missingcast.rds'))
err <- try(prepare_model_data(),silent=TRUE)
stopifnot(inherits(err,'try-error'),grepl('cast_id',err),!grepl('LocationCode',err))
cat('METADATA_TESTS_PASS: aliases; unknown-location matching; conflicting MLD rejection; date ambiguity rejection; required cast-ID validation.\n')
