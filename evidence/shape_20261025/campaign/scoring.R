# Profile, onset and transition-shape scores use the same raw observations.
source(file.path(Sys.getenv("EXPERIMENT_DIR"),"profile_helpers.R"))
source(file.path(Sys.getenv("EXPERIMENT_DIR"),"onset_metrics.R"))
finite_mean<-function(x)if(any(is.finite(x)))mean(x,na.rm=TRUE)else NA_real_

curve_sharpness <- function(curve) {
  z<-curve$depth;y<-curve$value
  out<-list(slope=NA_real_,rise_depth=NA_real_,width=NA_real_,status="insufficient_support")
  if(length(z)<7L||sum(is.finite(y))<7L)return(out)
  rise<-rapid_rise(curve);out$slope<-unname(rise["slope"]);out$rise_depth<-unname(rise["depth"])
  peak<-which.max(replace(y,!is.finite(y),-Inf));height<-y[peak]
  if(!is.finite(height)||height<.5){out$status<-"insufficient_enhancement";return(out)}
  if(z[peak]<=15||peak>=length(z)-1L){out$status<-"boundary_peak";return(out)}
  rising<-list(depth=z[seq_len(peak)],value=y[seq_len(peak)])
  # Width from 25% to 75% of the profile's own positive enhancement peak.
  # These relative thresholds differ from the fixed +0.5 onset definition.
  lo<-onset_from_curve(rising,.25*height,persistence=3L)
  hi<-onset_from_curve(rising,.75*height,persistence=3L)
  if(lo$status!="detected"||hi$status!="detected") {
    out$status<-paste0("censored_",lo$status,"_",hi$status);return(out)
  }
  out$width<-hi$depth-lo$depth;out$status<-"detected";out
}
shape_metrics <- function(tab) {
  a<-curve_sharpness(onset_curve(tab$depth,tab$chl_norm))
  b<-curve_sharpness(onset_curve(tab$depth,tab$predicted))
  data.frame(cast_key=as.character(tab$cast_key[1]),Region=as.character(tab$Region[1]),
    MLD_temp_MIX=tab$MLD_temp_MIX[1],observed_max_rise_slope=a$slope,
    predicted_max_rise_slope=b$slope,rise_slope_abs_error=abs(a$slope-b$slope),
    observed_rise_depth=a$rise_depth,predicted_rise_depth=b$rise_depth,
    observed_rise_width_25_75=a$width,predicted_rise_width_25_75=b$width,
    rise_width_abs_error=abs(a$width-b$width),observed_width_status=a$status,predicted_width_status=b$status)
}

evaluation_sets <- function(experiment,case) {
  train<-readRDS(data_file(experiment,case))
  validation<-readRDS(file.path(experiment,"data",paste0(case$split,"_validation.rds")))
  reserve<-readRDS(file.path(experiment,"data",paste0(case$split,"_reserve.rds")))
  audit<-readRDS(file.path(experiment,"data",paste0(case$dataset_id,"_training_audit.rds")))
  if(case$fraction<1) {
    stopifnot(!any(validation$cast_key %in% train$cast_key),!any(reserve$cast_key %in% train$cast_key))
    if(case$split=="date_blocked")stopifnot(!any(group_key(validation) %in% group_key(train)))
    return(list(training=audit,validation=validation,reserve=reserve))
  }
  # Full-data fitting includes previous holdouts: never publish these as validation.
  stopifnot(all(validation$cast_key %in% train$cast_key),all(reserve$cast_key %in% train$cast_key))
  list(training_full=train,reconstruction_primary=validation,reconstruction_reserve=reserve)
}

score_predictions <- function(fit,experiment,run,case,result) {
  checks<-file.path(run,"checks");dir.create(checks,showWarnings=FALSE)
  total<-posterior::ndraws(posterior::as_draws_array(fit))
  set.seed(8104);ids<-sort(sample.int(total,min(as.integer(Sys.getenv("SHAPE_PRED_DRAWS","2000")),total)))
  saveRDS(ids,file.path(checks,"prediction_draw_ids.rds"));result$prediction_draws_used<-length(ids)
  sets<-evaluation_sets(experiment,case)
  if(case$mode=="prior") {
    d<-sets[[1]]
    # One actual cast per region gives a bounded, reproducible scale check.
    prior_keys<-vapply(split(as.character(d$cast_key),d$Region),function(x)sort(unique(x))[1],character(1))
    d<-d[d$cast_key %in% prior_keys,,drop=FALSE]
    set.seed(8201);y<-brms::posterior_predict(fit,newdata=d[setdiff(names(d),"chl_norm")],draw_ids=ids,re_formula=NULL)
    write.csv(data.frame(t(quantile(as.vector(y),c(.001,.01,.05,.5,.95,.99,.999)))),
      file.path(checks,"prior_predictive_quantiles.csv"),row.names=FALSE)
    result$prediction_interpretation<-"Prior-predictive scale check on one actual cast per region; not an accuracy comparison"
    return(result)
  }
  train<-readRDS(data_file(experiment,case));all_scores<-list()
  for(dataset in names(sets))for(mode in c("conditional","population")) {
    d<-sets[[dataset]];tag<-paste(dataset,mode,sep="_")
    groups<-split(seq_len(nrow(d)),group_key(d));casts<-shapes<-onsets<-examples<-list()
    negative_y<-negative_mu<-negative_observed<-n_y<-n_mu<-n_observed<-0
    for(gi in seq_along(groups)) {
      g<-d[groups[[gi]],,drop=FALSE];predictors<-g[setdiff(names(g),"chl_norm")]
      re<-if(mode=="conditional")NULL else NA
      # Keep all rows of a nested date group together so unseen group draws
      # are shared by its casts. Only marginal profile scores are aggregated.
      set.seed(8201+gi)
      mu<-brms::posterior_epred(fit,newdata=predictors,draw_ids=ids,re_formula=re,
        allow_new_levels=TRUE,sample_new_levels="gaussian")
      set.seed(8201+gi)
      y<-brms::posterior_predict(fit,newdata=predictors,draw_ids=ids,re_formula=re,
        allow_new_levels=TRUE,sample_new_levels="gaussian")
      q<-apply(y,2,quantile,c(.05,.95))
      tab<-data.frame(g[c("cast_key","Region","Ocean","date","enso","depth","log10_distance","MLD_temp_MIX","chl_norm")],
        predicted=colMeans(mu),lower=q[1,],upper=q[2,])
      tab$interval_score90<-tab$upper-tab$lower+20*pmax(tab$lower-tab$chl_norm,0)+20*pmax(tab$chl_norm-tab$upper,0)
      tab$date_seen<-group_key(g) %in% group_key(train)
      for(ii in split(seq_len(nrow(tab)),as.character(tab$cast_key))) {
        z<-tab[ii,,drop=FALSE];key<-as.character(z$cast_key[1])
        casts[[key]]<-profile_metrics(z);shapes[[key]]<-shape_metrics(z)
        onsets[[key]]<-onset_profile_metrics(z,mu[,ii,drop=FALSE],max_draws=400L)
        # Save one deterministic example per region, with observed raw values.
        region<-as.character(z$Region[1]);if(is.null(examples[[region]]))examples[[region]]<-z
      }
      negative_y<-negative_y+sum(y< -1);n_y<-n_y+length(y)
      negative_mu<-negative_mu+sum(colMeans(mu)< -1);n_mu<-n_mu+ncol(mu)
      negative_observed<-negative_observed+sum(g$chl_norm< -1);n_observed<-n_observed+nrow(g)
      rm(mu,y);invisible(gc())
    }
    cs<-do.call(rbind,casts);ss<-do.call(rbind,shapes);os<-do.call(rbind,onsets)
    write.csv(cs,file.path(checks,paste0(tag,"_cast_metrics.csv")),row.names=FALSE)
    write.csv(ss,file.path(checks,paste0(tag,"_shape_cast_metrics.csv")),row.names=FALSE)
    write.csv(os,file.path(checks,paste0(tag,"_onset_cast_metrics.csv")),row.names=FALSE)
    write.csv(do.call(rbind,examples),file.path(checks,paste0(tag,"_example_profiles.csv")),row.names=FALSE)
    score<-score_casts(cs);all_scores[[tag]]<-data.frame(dataset,mode,score)
    for(n in names(score))result[[paste(tag,n,sep="_")]]<-score[[n]]
    result[[paste0(tag,"_rise_slope_mae")]]<-finite_mean(ss$rise_slope_abs_error)
    result[[paste0(tag,"_rise_slope_pairs")]]<-sum(is.finite(ss$rise_slope_abs_error))
    result[[paste0(tag,"_rise_width_mae")]]<-finite_mean(ss$rise_width_abs_error)
    result[[paste0(tag,"_rise_width_pairs")]]<-sum(is.finite(ss$rise_width_abs_error))
    result[[paste0(tag,"_fraction_predictive_below_minus_one")]]<-negative_y/n_y
    result[[paste0(tag,"_fraction_expected_below_minus_one")]]<-negative_mu/n_mu
    result[[paste0(tag,"_fraction_observed_below_minus_one")]]<-negative_observed/n_observed
    rs<-do.call(rbind,lapply(split(cs,cs$Region),function(x)data.frame(Region=x$Region[1],score_casts(x))))
    write.csv(rs,file.path(checks,paste0(tag,"_by_region.csv")),row.names=FALSE)
    ons<-do.call(rbind,lapply(split(os,os$threshold),function(x)data.frame(threshold=x$threshold[1],onset_scores(x))))
    write.csv(ons,file.path(checks,paste0(tag,"_onset_scores.csv")),row.names=FALSE)
    primary<-ons[ons$threshold==.5,,drop=FALSE]
    for(n in setdiff(names(primary),"threshold"))result[[paste0(tag,"_onset50_",n)]]<-primary[[n]]
    # Checkpoint each dataset before any optional plotting operation.
    write.csv(result,file.path(run,"test_summary.csv"),row.names=FALSE)
    cat("Scored",tag,nrow(cs),"casts\n")
  }
  write.csv(do.call(rbind,all_scores),file.path(checks,"prediction_scores.csv"),row.names=FALSE)
  result$prediction_interpretation<-if(case$fraction==1)
    "100% training reconstruction; former holdout scores are apparent/in-training, not validation" else
    paste("Whole-cast",case$split,"validation and reserve; population predictions set date/ENSO effects to zero")
  result
}

inside_hull <- function(x,y,newx,newy) {
  h<-chull(x,y);if(length(h)<3L)return(rep(FALSE,length(newx)))
  px<-x[c(h,h[1])];py<-y[c(h,h[1])]
  vapply(seq_along(newx),function(j){cross<-(tail(px,-1)-head(px,-1))*(newy[j]-head(py,-1))-
    (tail(py,-1)-head(py,-1))*(newx[j]-head(px,-1));all(cross>= -1e-8)||all(cross<=1e-8)},logical(1))
}
plot_shape_sensitivity <- function(fit,experiment,run,case) {
  d<-readRDS(data_file(experiment,case));meta<-unique(d[c("cast_key","Region","Ocean","date","enso","MLD_temp_MIX","log10_distance")])
  levels<-readRDS(file.path(experiment,"region_levels.rds"));grids<-list()
  for(region in levels) {
    m<-meta[as.character(meta$Region)==region,,drop=FALSE]
    q<-quantile(m$MLD_temp_MIX,c(.1,.5,.9));base<-m[which.min(abs(m$log10_distance-median(m$log10_distance))),,drop=FALSE]
    g<-base[rep(1,3*141),,drop=FALSE];g$depth<-rep(10:150,3)
    g$MLD_temp_MIX<-rep(q,each=141);g$mld_level<-rep(c("Shallow MLD (10%)","Median MLD","Deep MLD (90%)"),each=141)
    g$inside_mld_distance_hull<-inside_hull(m$MLD_temp_MIX,m$log10_distance,g$MLD_temp_MIX,g$log10_distance)
    region_data<-d[as.character(d$Region)==region,,drop=FALSE]
    g$inside_depth_mld_hull<-inside_hull(region_data$depth,region_data$MLD_temp_MIX,g$depth,g$MLD_temp_MIX)
    grids[[region]]<-g
  }
  grid<-do.call(rbind,grids)
  grid<-apply_features(grid,readRDS(file.path(experiment,"scaling.rds"))[[case$split]],levels)
  total<-posterior::ndraws(posterior::as_draws_array(fit));set.seed(8104);ids<-sort(sample.int(total,min(total,1000L)))
  mu<-brms::posterior_epred(fit,newdata=grid,re_formula=NA,draw_ids=ids,allow_new_levels=TRUE)
  q<-apply(mu,2,quantile,c(.05,.5,.95));grid$mean<-colMeans(mu)
  grid$lower90<-q[1,];grid$median<-q[2,];grid$upper90<-q[3,]
  grid$joint_support<-grid$inside_mld_distance_hull & grid$inside_depth_mld_hull
  write.csv(grid,file.path(run,"checks","mld_shape_slices.csv"),row.names=FALSE)
  library(ggplot2)
  p<-ggplot(grid,aes(mean,depth,colour=mld_level,group=mld_level))+
    geom_path(linewidth=.7)+scale_y_reverse()+facet_wrap(~Region,scales="free_x",ncol=3)+
    theme_bw()+theme(legend.position="bottom")+
    labs(x="Expected chl_norm",y="Depth (m)",colour=NULL,title=case$label,
      subtitle="MLD quantile slices at fixed regional distance; date/ENSO effects zero; hypothetical associations, not causal effects")
  ggsave(file.path(run,"checks","mld_shape_slices.png"),p,width=13,height=10,dpi=140,bg="white")
}
