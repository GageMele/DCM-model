# Whole-cast predictive checks. Every candidate uses identical evaluation casts.
# Conditional predictions retain learned date/ENSO effects; new levels draw from
# their estimated group distributions. Population predictions set those effects
# to zero and answer a different, explicitly labelled question.
profile_metrics <- function(x) {
  if(!nrow(x)) stop("Cannot score an empty cast")
  x <- x[order(x$depth),]
  e <- x$chl_norm-x$predicted
  shape <- (x$chl_norm-mean(x$chl_norm))-(x$predicted-mean(x$predicted))
  lag1 <- if(nrow(x)>3 && sd(e)>0) suppressWarnings(cor(head(e,-1),tail(e,-1))) else NA_real_
  data.frame(cast_key=as.character(x$cast_key[1]),Region=as.character(x$Region[1]),
    observations=nrow(x),mse=mean(e^2),mae=mean(abs(e)),bias=mean(x$predicted-x$chl_norm),
    shape_mse=mean(shape^2),
    peak_depth_abs_error=abs(x$depth[which.max(x$chl_norm)]-x$depth[which.max(x$predicted)]),
    peak_height_abs_error=abs(max(x$chl_norm)-max(x$predicted)),
    coverage90=mean(x$chl_norm>=x$lower & x$chl_norm<=x$upper),
    interval_width90=mean(x$upper-x$lower),residual_lag1=lag1)
}

score_casts <- function(tab) {
  data.frame(casts=nrow(tab),rmse=sqrt(mean(tab$mse)),mae=mean(tab$mae),bias=mean(tab$bias),
    shape_rmse=sqrt(mean(tab$shape_mse)),peak_depth_mae=mean(tab$peak_depth_abs_error),
    peak_height_mae=mean(tab$peak_height_abs_error),coverage90=mean(tab$coverage90),
    interval_width90=mean(tab$interval_width90),
    median_residual_lag1=if(any(is.finite(tab$residual_lag1))) median(tab$residual_lag1,na.rm=TRUE) else NA_real_)
}

evaluate_predictions <- function(fit, experiment, run, case, result, n_draws=400L) {
  checks <- file.path(run,"checks")
  train <- readRDS(file.path(experiment,"data","p10_model.rds"))
  fitted_data <- readRDS(file.path(experiment,"data",paste0(case$subset,"_model.rds")))
  validation <- readRDS(file.path(experiment,"data","validation_model.rds"))
  total <- posterior::ndraws(posterior::as_draws_array(fit))
  set.seed(8104L); draw_ids <- sort(sample.int(total,min(n_draws,total)))
  saveRDS(draw_ids,file.path(checks,"prediction_draw_ids.rds"))
  if(case$mode=="prior") {
    # Proper priors only. This is a scale check, never validation accuracy.
    yrep <- brms::posterior_predict(fit,newdata=train,draw_ids=draw_ids,re_formula=NULL)
    q <- apply(yrep,2L,quantile,probs=c(.05,.5,.95))
    write.csv(data.frame(train[c("cast_key","Region","depth","chl_norm")],
      prior_q05=q[1,],prior_median=q[2,],prior_q95=q[3,]),
      file.path(checks,"prior_predictive_observations.csv"),row.names=FALSE)
    write.csv(data.frame(t(quantile(as.vector(yrep),c(.001,.01,.05,.5,.95,.99,.999)))),
      file.path(checks,"prior_predictive_quantiles.csv"),row.names=FALSE)
    result$prediction_interpretation <- "prior predictive scale check; excluded from model ranking"
    return(result)
  }
  scores <- list(); all_rows <- list()
  for(dataset in c("training","validation")) {
    # The same 10% training casts are inspected even when a model fits 20%/50%.
    d <- if(dataset=="training") train else validation
    for(mode in c("conditional","population")) {
      cat("Prediction:",dataset,mode,nrow(d),"observations\n")
      set.seed(8201L)
      re <- if(mode=="conditional") NULL else NA
      mu <- brms::posterior_epred(fit,newdata=d,draw_ids=draw_ids,re_formula=re,
        allow_new_levels=TRUE,sample_new_levels="gaussian")
      set.seed(8201L)
      yrep <- brms::posterior_predict(fit,newdata=d,draw_ids=draw_ids,re_formula=re,
        allow_new_levels=TRUE,sample_new_levels="gaussian")
      q <- apply(yrep,2L,quantile,probs=c(.05,.95))
      tab <- data.frame(d[c("cast_key","Region","Ocean","date","enso","depth",
                            "log10_distance","MLD_temp_MIX","chl_norm")],
        predicted=colMeans(mu),lower=q[1,],upper=q[2,])
      group_key <- function(x, v) paste(x$Ocean,x$Region,x[[v]],sep="\r")
      tab$date_seen <- ifelse(group_key(d,"date") %in% group_key(fitted_data,"date"),"seen","new")
      tab$enso_seen <- ifelse(group_key(d,"enso") %in% group_key(fitted_data,"enso"),"seen","new")
      tab <- tab[order(tab$Region,tab$cast_key,tab$depth),]
      tag <- paste(dataset,mode,sep="_")
      cast <- do.call(rbind,lapply(split(tab,as.character(tab$cast_key)),profile_metrics))
      write.csv(cast,file.path(checks,paste0(tag,"_cast_metrics.csv")),row.names=FALSE)
      write.csv(tab,file.path(checks,paste0(tag,"_observations.csv")),row.names=FALSE)
      sc <- score_casts(cast); scores[[tag]] <- data.frame(dataset,mode,sc)
      for(name in names(sc)) result[[paste(tag,name,sep="_")]] <- sc[[name]]
      regional <- do.call(rbind,lapply(split(cast,cast$Region),function(x)
        data.frame(Region=x$Region[1],score_casts(x))))
      write.csv(regional,file.path(checks,paste0(tag,"_by_region.csv")),row.names=FALSE)
      # Quantile cuts depend only on training predictors, never validation y.
      breaks <- function(x) unique(c(-Inf,quantile(x,c(.25,.5,.75)),Inf))
      tab$depth_bin <- cut(tab$depth,c(9.999,30,60,100,150),include.lowest=TRUE)
      tab$mld_bin <- cut(tab$MLD_temp_MIX,breaks(train$MLD_temp_MIX),include.lowest=TRUE)
      tab$distance_bin <- cut(tab$log10_distance,breaks(train$log10_distance),include.lowest=TRUE)
      stratified <- do.call(rbind,lapply(c("depth_bin","mld_bin","distance_bin","date_seen","enso_seen"),function(v)
        do.call(rbind,lapply(split(tab,tab[[v]]),function(x) {
          if(!nrow(x)) return(NULL)
          cm <- do.call(rbind,lapply(split(x,as.character(x$cast_key)),profile_metrics))
          data.frame(stratum=v,bin=as.character(x[[v]][1]),score_casts(cm))
        }))))
      write.csv(stratified,file.path(checks,paste0(tag,"_strata.csv")),row.names=FALSE)
      if(mode=="conditional") all_rows[[dataset]] <- tab
      rm(mu,yrep); invisible(gc())
    }
  }
  write.csv(do.call(rbind,scores),file.path(checks,"prediction_scores.csv"),row.names=FALSE)
  result$prediction_interpretation <- paste("Held-out whole casts outside frozen 50%;",
    "conditional date/ENSO effects vs zero-effect population predictions; 400 draws;",
    "exploratory validation used for selection, not a final independent test.")
  # Save numerical results before optional plots, so rendering failure loses no scores.
  write.csv(result,file.path(run,"test_summary.csv"),row.names=FALSE)
  pdf(file.path(checks,"profile_checks.pdf"),width=11,height=8.5)
  on.exit(dev.off(),add=TRUE)
  par(mfrow=c(2,4),mar=c(4,4,3,1))
  for(dataset in names(all_rows)) for(region in unique(as.character(train$Region))) {
    tab <- all_rows[[dataset]]
    keys <- sort(unique(as.character(tab$cast_key[as.character(tab$Region)==region])))
    for(key in head(keys,2L)) {
      x <- tab[as.character(tab$cast_key)==key,]; x <- x[order(x$depth),]
      limits <- range(c(x$chl_norm,x$lower,x$upper),finite=TRUE)
      plot(x$chl_norm,x$depth,type="n",xlim=limits,ylim=c(150,10),
        xlab="chl_norm",ylab="Depth (m)",main=paste(dataset,region))
      polygon(c(x$lower,rev(x$upper)),c(x$depth,rev(x$depth)),
        border=NA,col=grDevices::adjustcolor("steelblue",alpha.f=.25))
      lines(x$predicted,x$depth,col="blue",lwd=2)
      points(x$chl_norm,x$depth,pch=16,cex=.3)
      mtext(substr(key,1,40),cex=.55)
    }
  }
  result
}
