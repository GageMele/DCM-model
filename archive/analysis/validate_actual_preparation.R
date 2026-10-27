# Validate preparation on the original CSVs without modifying those inputs.
source('DCM_Ibex_workflow/01_prepare_data.R')
dir.create('validation_actual_inputs',showWarnings=FALSE)
root <- 'K:/Gage Mele/chla/DCM distribution model Ibex'
Sys.setenv(CTD_INPUT=file.path(root,'CTD_NOAA_WHOI_KAUST.csv'),
           MLD_INPUT=file.path(root,'mld_full.csv'),CAST_FRACTION='1',
           MODEL_DATA='validation_actual_inputs/full.rds')
b <- prepare_model_data()
stopifnot(nlevels(b$model_data$Region)==8,
          identical(b$config$smoothing,'none'),identical(b$config$binned,FALSE),
          nrow(read.csv('validation_actual_inputs/full_qc/unmatched_ctd_casts.csv'))==0,
          nrow(read.csv('validation_actual_inputs/full_qc/unused_mld_casts.csv'))==0)
raw <- read.csv(file.path(root,'CTD_NOAA_WHOI_KAUST.csv'),stringsAsFactors=FALSE,check.names=FALSE)
keys <- function(region,cast,depth) paste(region,cast,depth,sep='|')
index <- match(keys(b$model_data$Region,b$model_data$cast_id,b$model_data$depth),
               keys(raw$Region.x,raw$DeepCTDID,raw$Depth_m))
stopifnot(!anyNA(index),all(b$model_data$Chl==raw$Chl[index]),
          all(b$model_data$depth==raw$Depth_m[index]))
Sys.setenv(CAST_FRACTION='0.25',MODEL_DATA='validation_actual_inputs/pilot.rds')
p <- prepare_model_data()
stopifnot(nlevels(p$model_data$Region)==8,
          all(p$region_summary$distinct_depths>=15),
          all(p$region_summary$distinct_distances>=5),
          all(p$region_summary$distinct_mlds>=5))
cat('ACTUAL_INPUTS_PASS: full and pilot preparation; all input cast identifiers match MLD; eight regions; raw Chl/depth retained; pilot marginal support passes.\n')
