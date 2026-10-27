# Local synthetic checks; no user data or model fitting is performed.
source('DCM_Ibex_workflow/01_prepare_data.R')
dir.create('validation_unsmoothed', showWarnings=FALSE)
ctd <- do.call(rbind, lapply(c('Central RS','Samoa'), function(region) {
  do.call(rbind,lapply(1:8,function(i) data.frame(
    LocationCode=paste0('L',i),DeepCTDID=paste0('C',i),Region.x=region,
    UTCDate.x=paste0('2011-01-',sprintf('%02d',i)),
    Depth_m=seq(5,100,5),Temp_DegC=27-seq(5,100,5)/100,
    Chl=(if(region=='Samoa') 1 else 2)*(1+seq(5,100,5)/100), distance=i)))
}))
mld <- unique(ctd[c('LocationCode','DeepCTDID','Region.x')])
names(mld) <- c('LocationCode','cast_id','Region')
mld$Ocean <- ifelse(mld$Region=='Samoa','Pacific','RedSea')
mld$MLD_temp_MIX <- seq_len(nrow(mld))*3
mld <- mld[-nrow(mld),]
write.csv(ctd,'validation_unsmoothed/ctd.csv',row.names=FALSE)
write.csv(rbind(mld,mld[1,]),'validation_unsmoothed/mld.csv',row.names=FALSE)
Sys.setenv(CTD_INPUT='validation_unsmoothed/ctd.csv',MLD_INPUT='validation_unsmoothed/mld.csv',
           MODEL_DATA='validation_unsmoothed/full.rds',SMOOTHING='notebook',CAST_FRACTION='1')
b <- prepare_model_data()
d <- b$model_data
stopifnot(nrow(d)==15*19, length(unique(d$cast_key))==15,
          all(is.finite(d$MLD_temp_MIX)),all(d$depth>=10),all(d$depth<=250),
          all(abs(d$surface_chl[d$Region=='Samoa']-1.1)<1e-10),
          all(abs(d$surface_chl[d$Region=='Central RS']-2.2)<1e-10),
          all(abs(d$chl_norm-(d$Chl/d$surface_chl-1))<1e-10))
stopifnot(identical(b$config$smoothing,'none'),identical(b$config$binned,FALSE))
# An obsolete SMOOTHING setting cannot activate removed processing code.
Sys.setenv(MODEL_DATA='validation_unsmoothed/subset.rds',CAST_FRACTION='0.5')
s <- prepare_model_data()
stopifnot(length(unique(s$model_data$cast_key))==7)
# Conflicting right-side join keys must fail, rather than expanding rows.
bad <- mld[1,];bad$MLD_temp_MIX <- bad$MLD_temp_MIX+1
write.csv(rbind(mld,bad),'validation_unsmoothed/mld.csv',row.names=FALSE)
Sys.setenv(MODEL_DATA='validation_unsmoothed/conflict.rds')
err <- try(prepare_model_data(),silent=TRUE)
stopifnot(inherits(err,'try-error'),grepl('Multiple MLD',err))
# Uneven depths in the same 5 m interval must remain separate observations,
# including a sharp chlorophyll excursion which would be removed by smoothing.
extra <- ctd[c(1,1),]
extra$Depth_m <- c(11.2,13.7)
extra$Chl <- c(0.31,4.8)
write.csv(rbind(ctd,extra),'validation_unsmoothed/ctd.csv',row.names=FALSE)
write.csv(mld,'validation_unsmoothed/mld.csv',row.names=FALSE)
Sys.setenv(MODEL_DATA='validation_unsmoothed/irregular.rds',CAST_FRACTION='1')
irregular <- prepare_model_data()$model_data
stopifnot(nrow(irregular)==287,
          identical(irregular$Chl[irregular$depth %in% c(11.2,13.7)],c(0.31,4.8)),
          identical(irregular$depth[irregular$depth %in% c(11.2,13.7)],c(11.2,13.7)))
key <- function(region,cast,depth) paste(region,cast,depth,sep='|')
index <- match(key(irregular$Region,irregular$cast_id,irregular$depth),
               key(ctd$Region.x,ctd$DeepCTDID,ctd$Depth_m))
stopifnot(all(irregular$Chl[!is.na(index)]==ctd$Chl[index[!is.na(index)]]),
          all(irregular$depth[!is.na(index)]==ctd$Depth_m[index[!is.na(index)]]))
# Check the actual mgcv three-way tensor basis with enough predictor support.
grid <- expand.grid(depth=seq(10,250,5),log10_distance=seq(0,2,length.out=5),MLD_temp_MIX=seq(10,100,length.out=5))
grid$Region <- factor('Samoa',levels=c('Central RS','Samoa'))
grid <- rbind(transform(grid,Region=factor('Central RS',levels=levels(grid$Region))),grid)
spec <- mgcv::interpret.gam(y ~ Region + t2(depth,log10_distance,MLD_temp_MIX,
  by=Region,bs=c('cr','cr','cr'),k=c(15,5,5)))$smooth.spec[[1]]
bases <- mgcv::smoothCon(spec,grid,absorb.cons=TRUE)
stopifnot(length(bases)==2,all(vapply(bases,function(x)all(is.finite(x$X)),logical(1))))
# Evaluate the metric helper from the plotting script without loading brms.
expr <- parse('DCM_Ibex_workflow/03_check_and_plot.R')
function_body <- expr[[1]][[3]][[3]]
metric_expr <- Filter(function(x)is.call(x)&&identical(x[[1]],as.name('<-'))&&identical(x[[2]],as.name('metric_one')),
                      as.list(function_body))[[1]]
eval(metric_expr)
v <- metric_one(c(0,1,2,1,0),c(10,20,30,40,50))
stopifnot(v['peak_depth']==30,v['enhancement_fwhm_m']==20,v['subsurface_enhancement_peak']==1)
v <- metric_one(rep(-1,5),seq(10,50,10))
stopifnot(v['subsurface_enhancement_peak']==0,is.na(v['enhancement_fwhm_m']))
cat('PASS: original Chl/depth retention, irregular observations without smoothing/binning, joins, cast identity, normalization, regional sampling, conflict rejection, tensor basis, profile metrics.\n')
