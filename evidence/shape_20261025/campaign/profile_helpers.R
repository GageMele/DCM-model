# Existing cast-weighted profile metrics, reused without changing observations.
trap_area <- function(depth,value)sum(diff(depth)*(head(value,-1)+tail(value,-1))/2)
profile_metrics <- function(x) {
  x<-x[order(x$depth),];e<-x$predicted-x$chl_norm
  shape<-(x$predicted-mean(x$predicted))-(x$chl_norm-mean(x$chl_norm))
  observed_peak<-x$depth[which.max(x$chl_norm)]
  predicted_peak<-x$depth[which.max(x$predicted)]
  data.frame(cast_key=as.character(x$cast_key[1]),Region=as.character(x$Region[1]),observations=nrow(x),
    mse=mean(e^2),mae=mean(abs(e)),bias=mean(e),shape_mse=mean(shape^2),
    peak_depth_abs_error=abs(observed_peak-predicted_peak),
    peak_height_abs_error=abs(max(x$chl_norm)-max(x$predicted)),
    coverage90=mean(x$chl_norm>=x$lower & x$chl_norm<=x$upper),
    interval_width90=mean(x$upper-x$lower),interval_score90=mean(x$interval_score90),
    positive_auc_abs_error=abs(trap_area(x$depth,pmax(x$predicted,0))-trap_area(x$depth,pmax(x$chl_norm,0))),
    observed_peak_depth=observed_peak,predicted_peak_depth=predicted_peak,
    observed_dcm_in_range=observed_peak>15,predicted_dcm_in_range=predicted_peak>15,
    peak_at_bottom=observed_peak==max(x$depth),
    residual_lag1=if(nrow(x)>3 && sd(e)>0)suppressWarnings(cor(head(e,-1),tail(e,-1)))else NA_real_)
}
score_casts <- function(x) {
  data.frame(casts=nrow(x),rmse=sqrt(mean(x$mse)),mae=mean(x$mae),bias=mean(x$bias),
    shape_rmse=sqrt(mean(x$shape_mse)),peak_depth_mae=mean(x$peak_depth_abs_error),
    peak_height_mae=mean(x$peak_height_abs_error),coverage90=mean(x$coverage90),
    interval_width90=mean(x$interval_width90),interval_score90=mean(x$interval_score90),
    positive_auc_mae=mean(x$positive_auc_abs_error),
    dcm_classification_agreement=mean(x$observed_dcm_in_range==x$predicted_dcm_in_range),
    median_residual_lag1=median(x$residual_lag1,na.rm=TRUE))
}

