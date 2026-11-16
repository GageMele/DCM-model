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
  train<-augment(train);validation<-augment(validation);reserve<-augment(reserve);audit<-augment(audit)
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
  if(ordinal_case(case))return(score_ordinal(fit,experiment,run,case,result,sets,ids,checks))
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
    lpd_sum<-0;lpd_n<-0;date_elpd<-list()
    for(gi in seq_along(groups)) {
      g<-d[groups[[gi]],,drop=FALSE];predictors<-g[setdiff(names(g),"chl_norm")]
      predictors_ll<-g[setdiff(names(g),"chl_norm")]
      if(case$family!="lognormal")predictors_ll$chl_norm<-g$chl_norm
      if(case$family=="lognormal")predictors<-predictors[setdiff(names(predictors),"ratio")]
      re<-if(mode=="conditional")NULL else population_re(case)
      # Keep all rows of a nested date group together so unseen group draws
      # are shared by its casts. Only marginal profile scores are aggregated.
      set.seed(8201+gi)
      mu<-brms::posterior_epred(fit,newdata=predictors,draw_ids=ids,re_formula=re,
        allow_new_levels=TRUE,sample_new_levels="gaussian")
      set.seed(8201+gi)
      y<-brms::posterior_predict(fit,newdata=predictors,draw_ids=ids,re_formula=re,
        allow_new_levels=TRUE,sample_new_levels="gaussian")
      # Lognormal cases model ratio = chl_norm + 1; score on the chl_norm scale.
      if(case$family=="lognormal"){mu<-mu-1;y<-y-1}
      # Held-out log predictive density on the chl_norm scale (ratio is a pure
      # shift of chl_norm, so no Jacobian term). Per-row and joint per-date.
      set.seed(8201+gi)
      ll<-brms::log_lik(fit,newdata=predictors_ll,draw_ids=ids,re_formula=re,
        allow_new_levels=TRUE,sample_new_levels="gaussian")
      lme<-function(x){m<-max(x);m+log(mean(exp(x-m)))}
      lpd_sum<-lpd_sum+sum(apply(ll,2,lme));lpd_n<-lpd_n+ncol(ll)
      date_elpd[[length(date_elpd)+1L]]<-data.frame(group=group_key(g)[1],Region=as.character(g$Region[1]),
        rows=nrow(g),casts=length(unique(g$cast_key)),elpd_joint=lme(rowSums(ll)),lpd_rows=sum(apply(ll,2,lme)))
      rm(ll)
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
    de<-do.call(rbind,date_elpd);write.csv(de,file.path(checks,paste0(tag,"_date_elpd.csv")),row.names=FALSE)
    result[[paste0(tag,"_lpd_per_row")]]<-lpd_sum/lpd_n
    result[[paste0(tag,"_elpd_date_joint_sum")]]<-sum(de$elpd_joint)
    result[[paste0(tag,"_elpd_dates")]]<-nrow(de)
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
  if(ordinal_case(case))return(plot_bin_probabilities(fit,experiment,run,case))
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
  mu<-brms::posterior_epred(fit,newdata=grid,re_formula=population_re(case),draw_ids=ids,allow_new_levels=TRUE)
  if(case$family=="lognormal")mu<-mu-1
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

# Ordinal depth-bin models (bin_class, dcm_depth). Holdout casts are reduced exactly
# like the training data (5 m bin medians, or one DCM row per cast) and scored on the
# predicted class probabilities: held-out log predictive density (per row and joint per
# date group, unseen groups drawn from their fitted distribution), ranked probability
# score, Brier scores for each "at least this class" threshold, and calibration tables.
score_ordinal <- function(fit,experiment,run,case,result,sets,ids,checks) {
  y_name<-if(case$equation=="dcm_depth")"dcm_class" else "enh_class"
  labels<-if(case$equation=="dcm_depth")dcm_classes(case) else enh_labels
  bv<-bin_var(case)
  lme<-function(x){m<-max(x);m+log(mean(exp(x-m)))}
  for(dataset in names(sets))for(mode in c("conditional","population")) {
    d<-prepare_data(sets[[dataset]],case);tag<-paste(dataset,mode,sep="_")
    groups<-split(seq_len(nrow(d)),group_key(d));rows<-list();date_elpd<-list()
    re<-if(mode=="conditional")NULL else population_re(case)
    for(gi in seq_along(groups)) {
      g<-d[groups[[gi]],,drop=FALSE]
      set.seed(8201+gi)
      p<-brms::posterior_epred(fit,newdata=g,draw_ids=ids,re_formula=re,allow_new_levels=TRUE,sample_new_levels="gaussian")
      set.seed(8201+gi)
      ll<-brms::log_lik(fit,newdata=g,draw_ids=ids,re_formula=re,allow_new_levels=TRUE,sample_new_levels="gaussian")
      date_elpd[[gi]]<-data.frame(group=group_key(g)[1],Region=as.character(g$Region[1]),rows=nrow(g),
        casts=length(unique(g$cast_key)),elpd_joint=lme(rowSums(ll)),lpd_rows=sum(apply(ll,2,lme)))
      pbar<-apply(p,c(2,3),mean);K<-ncol(pbar);y<-g[[y_name]]
      cum_p<-t(apply(pbar,1,cumsum))[,seq_len(K-1),drop=FALSE]
      cum_o<-outer(y,seq_len(K-1),"<=")*1
      out<-data.frame(g[intersect(c("cast_key","Region","Ocean","date","enso","depth",bv,"MLD_temp_MIX",
        "log10_distance","chl_norm","dcm_depth","dcm_enhancement"),names(g))],observed=y,
        lpd=apply(ll,2,lme),rps=rowSums((cum_p-cum_o)^2)/(K-1),predicted_mode=max.col(pbar,"first"),
        predicted_median=1L+rowSums(cum_p<.5))
      # P(class >= k) for k = 2..K, i.e. 1 - P(class <= k-1).
      for(k in 2:K)out[[sprintf("p_at_least_%d",k)]]<-1-cum_p[,k-1]
      if(case$equation=="dcm_depth") {
        lo<-1L+rowSums(cum_p<.1);hi<-1L+rowSums(cum_p<.9)
        out$covered80<-y>=lo&y<=hi
      }
      rows[[gi]]<-out;rm(p,ll)
    }
    tab<-do.call(rbind,rows);de<-do.call(rbind,date_elpd);K<-length(labels)
    write.csv(tab,file.path(checks,paste0(tag,"_ordinal_rows.csv")),row.names=FALSE)
    write.csv(de,file.path(checks,paste0(tag,"_date_elpd.csv")),row.names=FALSE)
    result[[paste0(tag,"_rows")]]<-nrow(tab);result[[paste0(tag,"_casts")]]<-length(unique(tab$cast_key))
    result[[paste0(tag,"_lpd_per_row")]]<-mean(tab$lpd)
    result[[paste0(tag,"_elpd_date_joint_sum")]]<-sum(de$elpd_joint);result[[paste0(tag,"_elpd_dates")]]<-nrow(de)
    result[[paste0(tag,"_rps")]]<-mean(tab$rps)
    result[[paste0(tag,"_modal_accuracy")]]<-mean(tab$predicted_mode==tab$observed)
    result[[paste0(tag,"_median_class_mae")]]<-mean(abs(tab$predicted_median-tab$observed))
    for(k in 2:K)result[[paste0(tag,sprintf("_brier_at_least_%d",k))]]<-
      mean((tab[[sprintf("p_at_least_%d",k)]]-(tab$observed>=k))^2)
    if(case$equation=="dcm_depth") {
      # Class error times the class width is the DCM depth error in m.
      result[[paste0(tag,"_dcm_depth_mae_m")]]<-dcm_width(case)*mean(abs(tab$predicted_median-tab$observed))
      result[[paste0(tag,"_dcm_coverage80")]]<-mean(tab$covered80)
      by<-split(tab,tab$Region)
    } else by<-split(tab,list(tab$Region,tab[[bv]]),drop=TRUE)
    # Calibration: observed share vs mean predicted probability of each threshold.
    cal<-do.call(rbind,lapply(by,function(x){r<-data.frame(Region=x$Region[1],
      depth_bin=if(bv%in%names(x))as.character(x[[bv]][1]) else NA,rows=nrow(x))
      for(k in 2:K){r[[sprintf("observed_at_least_%d",k)]]<-mean(x$observed>=k)
        r[[sprintf("predicted_at_least_%d",k)]]<-mean(x[[sprintf("p_at_least_%d",k)]])};r}))
    write.csv(cal,file.path(checks,paste0(tag,"_calibration.csv")),row.names=FALSE)
    write.csv(result,file.path(run,"test_summary.csv"),row.names=FALSE)
    cat("Scored",tag,nrow(tab),"rows\n")
  }
  result$prediction_interpretation<-paste("Ordinal",case$equation,"model; whole-cast",case$split,
    "holdouts reduced like the training data; population predictions set date/ENSO effects to zero")
  result
}

# Population-level class probabilities by region and MLD (10/50/90% quantiles) at the
# region's median distance; date/ENSO effects zero. bin_class: probability that each
# 5 m bin reaches each enhancement class. dcm_depth: distribution of DCM depth.
plot_bin_probabilities <- function(fit,experiment,run,case) {
  d<-readRDS(data_file(experiment,case));meta<-unique(d[c("cast_key","Region","Ocean","date","enso","MLD_temp_MIX","log10_distance")])
  levels<-readRDS(file.path(experiment,"region_levels.rds"));grids<-list()
  bv<-bin_var(case)
  depths<-if(case$equation!="bin_class")50 else if(bv=="depth_m")10:150 else seq(12.5,147.5,5)
  for(region in levels) {
    m<-meta[as.character(meta$Region)==region,,drop=FALSE];if(!nrow(m))next
    q<-quantile(m$MLD_temp_MIX,c(.1,.5,.9));base<-m[which.min(abs(m$log10_distance-median(m$log10_distance))),,drop=FALSE]
    n<-length(depths);g<-base[rep(1,3*n),,drop=FALSE];g$depth<-rep(depths,3)
    g$MLD_temp_MIX<-rep(q,each=n);g$mld_level<-rep(c("Shallow MLD (10%)","Median MLD","Deep MLD (90%)"),each=n)
    grids[[region]]<-g
  }
  grid<-apply_features(do.call(rbind,grids),readRDS(file.path(experiment,"scaling.rds"))[[case$split]],levels)
  total<-posterior::ndraws(posterior::as_draws_array(fit));set.seed(8104);ids<-sort(sample.int(total,min(total,1000L)))
  p<-brms::posterior_epred(fit,newdata=grid,re_formula=population_re(case),draw_ids=ids,allow_new_levels=TRUE)
  K<-dim(p)[3];labels<-if(case$equation=="dcm_depth")dcm_classes(case) else enh_labels
  library(ggplot2)
  grid$mld_level<-factor(grid$mld_level,levels=c("Shallow MLD (10%)","Median MLD","Deep MLD (90%)"))
  if(case$equation=="bin_class") {
    long<-do.call(rbind,lapply(2:K,function(k){x<-apply(p[,,k:K,drop=FALSE],c(1,2),sum)
      q<-apply(x,2,quantile,c(.05,.5,.95))
      data.frame(grid[c("Region","depth",bv,"mld_level","MLD_temp_MIX")],
        threshold=sprintf("at least +%g%%",100*enh_breaks[k]),probability=colMeans(x),lower90=q[1,],upper90=q[3,])}))
    long$threshold<-factor(long$threshold,levels=unique(long$threshold))
    write.csv(long,file.path(run,"checks","bin_probability_slices.csv"),row.names=FALSE)
    pl<-ggplot(long[long$mld_level=="Median MLD",],aes(y=depth,x=probability,colour=threshold))+
      geom_errorbar(aes(xmin=lower90,xmax=upper90),width=0,orientation="y",alpha=.35)+geom_point(size=1)+
      scale_y_reverse()+facet_wrap(~Region,ncol=4)+theme_bw()+theme(legend.position="bottom")+
      labs(x=if(bv=="depth_m")"Probability the 1 m bin reaches the enhancement" else "Probability the 5 m bin's median reaches the enhancement",y="Depth bin centre (m)",colour=NULL,
        title=case$label,subtitle="Median MLD, median distance, date/ENSO effects zero; 90% intervals")
    ggsave(file.path(run,"checks","bin_probability_median_mld.png"),pl,width=13,height=8,dpi=140,bg="white")
    pl<-ggplot(long[long$threshold=="at least +50%",],aes(y=depth,x=probability,colour=mld_level))+
      geom_path()+geom_point(size=.8)+scale_y_reverse()+facet_wrap(~Region,ncol=4)+theme_bw()+theme(legend.position="bottom")+
      labs(x="Probability of at least +50% enhancement",y="Depth bin centre (m)",colour=NULL,title=case$label,
        subtitle="MLD quantile slices; hypothetical associations, not causal effects")
    ggsave(file.path(run,"checks","bin_probability_mld_slices.png"),pl,width=13,height=8,dpi=140,bg="white")
  } else {
    long<-do.call(rbind,lapply(seq_len(K),function(k){q<-apply(p[,,k],2,quantile,c(.05,.95))
      w<-dcm_width(case)
      data.frame(grid[c("Region","mld_level","MLD_temp_MIX")],dcm_bin=labels[k],dcm_mid=if(w==1)9+k else 10*k+5,
        probability=colMeans(p[,,k]),lower90=q[1,],upper90=q[2,])}))
    write.csv(long,file.path(run,"checks","dcm_depth_distribution.csv"),row.names=FALSE)
    # 10 m classes as grouped bars; 141 1 m classes as one line per MLD slice.
    pl<-if(dcm_width(case)==1)ggplot(long,aes(y=dcm_mid,x=probability,colour=mld_level))+geom_path() else
      ggplot(long,aes(y=dcm_mid,x=probability,fill=mld_level))+geom_col(position=position_dodge(width=8),width=8,orientation="y")
    pl<-pl+scale_y_reverse()+
      facet_wrap(~Region,ncol=4)+theme_bw()+theme(legend.position="bottom")+
      labs(x=sprintf("Probability the DCM is in this %g m bin",dcm_width(case)),y="Depth (m)",fill=NULL,colour=NULL,title=case$label,
        subtitle="Median distance, date/ENSO effects zero; hypothetical associations, not causal effects")
    ggsave(file.path(run,"checks","dcm_depth_distribution.png"),pl,width=13,height=8,dpi=140,bg="white")
  }
}
