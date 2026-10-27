# Common evaluation of full profiles; every primary run uses the same held-out casts.
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

plot_validation_examples<-function(tab,checks,mode,dataset="validation") {
  # Choose the same predictor-central validation cast per Region for every fit.
  # Selection uses neither response values nor model errors.
  metadata<-unique(tab[c("cast_key","Region","MLD_temp_MIX","log10_distance")])
  selected<-unlist(lapply(split(metadata,metadata$Region),function(x) {
    spread<-function(z)if(length(z)>1L && is.finite(sd(z)) && sd(z)>0)sd(z)else 1
    score<-abs(x$MLD_temp_MIX-median(x$MLD_temp_MIX))/spread(x$MLD_temp_MIX)+
      abs(x$log10_distance-median(x$log10_distance))/spread(x$log10_distance)
    x$cast_key[which.min(score)]
  }),use.names=FALSE)
  z<-tab[tab$cast_key %in% selected,,drop=FALSE]
  z$panel<-paste(z$Region,z$cast_key,sep=" | ")
  z<-z[order(z$panel,z$depth),,drop=FALSE]
  write.csv(z,file.path(checks,paste0(dataset,"_",mode,"_example_profiles.csv")),row.names=FALSE)
  p<-ggplot2::ggplot(z,ggplot2::aes(predicted,depth))+
    ggplot2::geom_ribbon(ggplot2::aes(xmin=lower,xmax=upper),orientation="y",fill="#94D2BD",alpha=.45)+
    ggplot2::geom_path(colour="#005F73",linewidth=.65)+
    ggplot2::geom_path(ggplot2::aes(x=chl_norm),colour="#222222",linewidth=.5)+
    ggplot2::scale_y_reverse()+ggplot2::facet_wrap(~panel,scales="free_x",ncol=3)+ggplot2::theme_bw()+
    ggplot2::labs(x="chl_norm",y="Depth (m)",title=paste("Independent",dataset,"casts:",mode),
      subtitle="Black: observed; blue: expected prediction; shading: 90% new-observation predictive interval")
  ggplot2::ggsave(file.path(checks,paste0(dataset,"_",mode,"_example_profiles.png")),p,width=13,height=10,dpi=130,bg="white")
}

make_original_style_plots <- function(fit,experiment,run,case,draw_ids) {
  reference<-readRDS(file.path(experiment,"profile_reference.rds"))
  grid<-reference[rep(seq_len(nrow(reference)),each=141L),,drop=FALSE]
  grid$depth<-rep(10:150,nrow(reference))
  levels<-readRDS(file.path(experiment,"region_levels.rds"))
  scaling<-readRDS(file.path(experiment,"scaling.rds"))[[case$split]]
  grid<-apply_features(grid,scaling,levels)
  # Group effects zero: intervals describe the expected population curve.
  mu<-brms::posterior_epred(fit,newdata=grid,re_formula=NA,draw_ids=draw_ids,allow_new_levels=TRUE)
  q<-apply(mu,2,quantile,c(.025,.25,.5,.75,.975))
  tab<-data.frame(grid[c("Region","Ocean","depth","log10_distance","MLD_temp_MIX")],
    mean=colMeans(mu),lower95=q[1,],lower50=q[2,],median=q[3,],upper50=q[4,],upper95=q[5,])
  write.csv(tab,file.path(run,"checks","population_reference_profiles.csv"),row.names=FALSE)
  peaks<-do.call(rbind,lapply(split(seq_len(nrow(grid)),grid$Region),function(ii) {
    j<-max.col(mu[,ii,drop=FALSE],ties.method="first")
    dep<-grid$depth[ii][j]
    height<-mu[cbind(seq_len(nrow(mu)),ii[j])]
    data.frame(Region=as.character(grid$Region[ii[1]]),
      peak_depth_q05=quantile(dep,.05),peak_depth_median=median(dep),peak_depth_q95=quantile(dep,.95),
      peak_height_q05=quantile(height,.05),peak_height_median=median(height),peak_height_q95=quantile(height,.95),
      probability_expected_peak_deeper_than_15m=mean(dep>15),probability_peak_at_150m=mean(dep==150))
  }))
  write.csv(peaks,file.path(run,"checks","population_expected_peak_summary.csv"),row.names=FALSE)
  # Average each joint draw first; never average confidence-interval endpoints.
  ocean<-do.call(rbind,lapply(unique(as.character(grid$Ocean)),function(o) {
    mm<-vapply(10:150,function(z)rowMeans(mu[,as.character(grid$Ocean)==o & grid$depth==z,drop=FALSE]),numeric(nrow(mu)))
    qq<-apply(mm,2,quantile,c(.025,.5,.975))
    data.frame(Ocean=o,depth=10:150,mean=colMeans(mm),lower95=qq[1,],median=qq[2,],upper95=qq[3,])
  }))
  write.csv(ocean,file.path(run,"checks","population_ocean_profiles_equal_region_weight.csv"),row.names=FALSE)
  library(ggplot2)
  p<-ggplot(tab,aes(mean,depth))+geom_ribbon(aes(xmin=lower95,xmax=upper95),orientation="y",fill="#94D2BD",alpha=.4)+
    geom_ribbon(aes(xmin=lower50,xmax=upper50),orientation="y",fill="#0A9396",alpha=.4)+geom_path(colour="#005F73")+
    scale_y_reverse()+facet_wrap(~Region,scales="free_x",ncol=3)+theme_bw()+
    labs(x="chl_norm",y="Depth (m)",title=case$label,
      subtitle="Expected population curves at fixed real training-cast predictor pairs; 50/95% posterior intervals")
  ggsave(file.path(run,"checks","population_profiles.png"),p,width=13,height=10,dpi=130,bg="white")
  p<-ggplot(ocean,aes(mean,depth))+geom_ribbon(aes(xmin=lower95,xmax=upper95),orientation="y",fill="#94D2BD",alpha=.5)+
    geom_path(colour="#005F73")+scale_y_reverse()+facet_wrap(~Ocean,scales="free_x")+theme_bw()+
    labs(x="chl_norm",y="Depth (m)",title="Expected ocean profiles with equal Region weights",
      subtitle="Joint posterior draws averaged across the fixed regional predictor pairs; 95% posterior intervals")
  ggsave(file.path(run,"checks","population_ocean_profiles.png"),p,width=11,height=7,dpi=130,bg="white")
  # Exceedance probabilities and exclusive bins are different quantities.
  cuts<-c(.5,1,2,3,4,5)
  exceed<-vapply(cuts,function(t)colMeans(mu>t),numeric(nrow(grid)))
  for(i in seq_along(cuts))tab[[paste0("prob_mu_gt_",cuts[i])]]<-exceed[,i]
  write.csv(tab,file.path(run,"checks","population_expected_exceedance.csv"),row.names=FALSE)
  probabilities<-cbind(1-exceed[,1],exceed[,-ncol(exceed),drop=FALSE]-exceed[,-1,drop=FALSE],exceed[,ncol(exceed)])
  bins<-c("<=1.5x","1.5-2x","2-3x","3-4x","4-5x","5-6x",">6x")
  long<-do.call(rbind,lapply(seq_along(bins),function(i)data.frame(tab[c("Region","depth")],bin=bins[i],probability=probabilities[,i])))
  long$bin<-factor(long$bin,levels=bins)
  write.csv(long,file.path(run,"checks","population_expected_probability_bins.csv"),row.names=FALSE)
  p<-ggplot(long,aes(depth,probability,fill=bin))+geom_area(position="stack")+coord_flip()+scale_x_reverse()+
    facet_wrap(~Region,ncol=3)+scale_fill_manual(values=c("#dddddd","#20255c","#4575b4","#74add1","#fdae61","#d73027","#a50026"))+
    theme_bw()+theme(legend.position="bottom")+
    labs(x="Depth (m)",y="Posterior probability",fill="Expected chlorophyll relative to surface",
      title="Exclusive enhancement bins for the expected population profile",
      subtitle="These are probabilities about the mean curve, not new-observation exceedance probabilities")
  ggsave(file.path(run,"checks","expected_enhancement_probability.png"),p,width=13,height=10,dpi=130,bg="white")
}

evaluate_predictions <- function(fit,experiment,run,case,result,n_draws=2000L) {
  checks<-file.path(run,"checks");dir.create(checks,showWarnings=FALSE,recursive=TRUE)
  train<-readRDS(training_eval_file(experiment,case))
  validation<-readRDS(validation_file(experiment,case))
  reserve<-readRDS(reserve_validation_file(experiment,case))
  fitted_data<-readRDS(data_file(experiment,case))
  stopifnot(!any(validation$cast_key %in% fitted_data$cast_key),
    !any(reserve$cast_key %in% fitted_data$cast_key),!any(reserve$cast_key %in% validation$cast_key))
  total<-posterior::ndraws(posterior::as_draws_array(fit))
  set.seed(8104L);draw_ids<-sort(sample.int(total,min(n_draws,total)))
  saveRDS(draw_ids,file.path(checks,"prediction_draw_ids.rds"))
  result$prediction_draws_used<-length(draw_ids)
  result$reserve_md5<-unname(tools::md5sum(reserve_validation_file(experiment,case)))
  result$validation_md5<-unname(tools::md5sum(validation_file(experiment,case)))
  result$training_evaluation_md5<-unname(tools::md5sum(training_eval_file(experiment,case)))
  if(case$mode=="prior") {
    set.seed(8201L)
    y<-brms::posterior_predict(fit,newdata=train[setdiff(names(train),"chl_norm")],draw_ids=draw_ids,re_formula=NULL)
    write.csv(data.frame(t(quantile(as.vector(y),c(.001,.01,.05,.5,.95,.99,.999)))),file.path(checks,"prior_predictive_quantiles.csv"),row.names=FALSE)
    result$prediction_interpretation<-"Prior-predictive plausibility only; no accuracy ranking"
    return(result)
  }
  scores<-list();batches<-list()
  for(dataset in c("training","validation","reserve"))for(mode in c("conditional","population")) {
    d<-switch(dataset,training=train,validation=validation,reserve=reserve)
    # Outcomes are retained only for scoring; never give validation y to prediction.
    predictors<-d[setdiff(names(d),"chl_norm")]
    re<-if(mode=="conditional")NULL else NA
    set.seed(8201L)
    mu<-brms::posterior_epred(fit,newdata=predictors,draw_ids=draw_ids,re_formula=re,allow_new_levels=TRUE,sample_new_levels="gaussian")
    set.seed(8201L)
    y<-brms::posterior_predict(fit,newdata=predictors,draw_ids=draw_ids,re_formula=re,allow_new_levels=TRUE,sample_new_levels="gaussian")
    q<-apply(y,2,quantile,c(.05,.95))
    tab<-data.frame(d[c("cast_key","Region","Ocean","date","enso","depth","log10_distance","MLD_temp_MIX","chl_norm")],
      predicted=colMeans(mu),lower=q[1,],upper=q[2,])
    tab$interval_score90<-tab$upper-tab$lower+20*pmax(tab$lower-tab$chl_norm,0)+20*pmax(tab$chl_norm-tab$upper,0)
    group_key<-function(x,v)paste(x$Ocean,x$Region,x[[v]],sep="\r")
    tab$date_seen<-ifelse(group_key(d,"date") %in% group_key(fitted_data,"date"),"seen","new")
    tab$enso_seen<-ifelse(group_key(d,"enso") %in% group_key(fitted_data,"enso"),"seen","new")
    tag<-paste(dataset,mode,sep="_")
    cast<-do.call(rbind,lapply(split(tab,as.character(tab$cast_key)),profile_metrics))
    write.csv(cast,file.path(checks,paste0(tag,"_cast_metrics.csv")),row.names=FALSE)
    # Per-row files retained locally for reconstruction and audit, not required for summaries.
    write.csv(tab,file.path(checks,paste0(tag,"_observations.csv")),row.names=FALSE)
    sc<-score_casts(cast);scores[[tag]]<-data.frame(dataset,mode,sc)
    for(n in names(sc))result[[paste(tag,n,sep="_")]]<-sc[[n]]
    regional<-do.call(rbind,lapply(split(cast,cast$Region),function(x)data.frame(Region=x$Region[1],score_casts(x))))
    write.csv(regional,file.path(checks,paste0(tag,"_by_region.csv")),row.names=FALSE)
    # Predictive batch sensitivity reuses existing draws; it is not a formal MCSE.
    result[[paste0(tag,"_fraction_observed_below_minus_one")]]<-mean(tab$chl_norm < -1)
    result[[paste0(tag,"_fraction_expected_below_minus_one")]]<-mean(tab$predicted < -1)
    result[[paste0(tag,"_fraction_predictive_below_minus_one")]]<-mean(y < -1)
    for(batch in 1:2) {
      ii<-if(batch==1L)seq_len(floor(nrow(mu)/2))else seq.int(floor(nrow(mu)/2)+1L,nrow(mu))
      qt<-apply(y[ii,,drop=FALSE],2,quantile,c(.05,.95))
      bt<-tab;bt$predicted<-colMeans(mu[ii,,drop=FALSE]);bt$lower<-qt[1,];bt$upper<-qt[2,]
      bt$interval_score90<-bt$upper-bt$lower+20*pmax(bt$lower-bt$chl_norm,0)+20*pmax(bt$chl_norm-bt$upper,0)
      bc<-do.call(rbind,lapply(split(bt,bt$cast_key),profile_metrics))
      batches[[length(batches)+1L]]<-data.frame(dataset,mode,batch,draws=length(ii),score_casts(bc))
    }
    # Store marginals for seen/new groups and fixed training-derived strata.
    metadata<-unique(train[c("cast_key","MLD_temp_MIX","log10_distance")])
    breaks<-function(v)unique(c(-Inf,as.numeric(quantile(metadata[[v]],c(.25,.5,.75))),Inf))
    tab$depth_bin<-cut(tab$depth,c(10,30,60,100,150),include.lowest=TRUE)
    tab$mld_bin<-cut(tab$MLD_temp_MIX,breaks("MLD_temp_MIX"),include.lowest=TRUE)
    tab$distance_bin<-cut(tab$log10_distance,breaks("log10_distance"),include.lowest=TRUE)
    strata<-list()
    for(v in c("depth_bin","mld_bin","distance_bin","date_seen","enso_seen"))for(b in unique(as.character(tab[[v]]))) {
      z<-tab[as.character(tab[[v]])==b,,drop=FALSE]
      cs<-do.call(rbind,lapply(split(z,z$cast_key),profile_metrics))
      score<-score_casts(cs)
      if(v=="depth_bin")score$dcm_classification_agreement<-NA_real_
      strata[[length(strata)+1L]]<-data.frame(stratum=v,bin=b,window_limited=v=="depth_bin",score)
    }
    write.csv(do.call(rbind,strata),file.path(checks,paste0(tag,"_strata.csv")),row.names=FALSE)
    # Predictive exceedance is explicitly labelled as a new observation probability.
    exceed<-data.frame(d[c("cast_key","Region","depth")])
    for(t in c(0,.5,1,2,3,4,5))exceed[[paste0("prob_new_observation_gt_",t)]]<-colMeans(y>t)
    write.csv(exceed,file.path(checks,paste0(tag,"_predictive_exceedance.csv")),row.names=FALSE)
    rm(mu,y);invisible(gc())
  }
  write.csv(do.call(rbind,scores),file.path(checks,"prediction_scores.csv"),row.names=FALSE)
  write.csv(do.call(rbind,batches),file.path(checks,"predictive_batch_sensitivity.csv"),row.names=FALSE)
  # Save scores before optional figure generation, retaining partial results if plotting fails.
  write.csv(result,file.path(run,"test_summary.csv"),row.names=FALSE)
  for(dataset in c("validation","reserve"))for(mode in c("conditional","population")) {
    tab<-read.csv(file.path(checks,paste0(dataset,"_",mode,"_observations.csv")),stringsAsFactors=FALSE)
    plot_validation_examples(tab,checks,mode,dataset)
  }
  make_original_style_plots(fit,experiment,run,case,draw_ids)
  result$prediction_interpretation<-paste("Whole-cast",case$split,"validation plus a disjoint reserve; DCM diagnostics confined to 10-150 m; population intervals set group effects to zero")
  result
}
