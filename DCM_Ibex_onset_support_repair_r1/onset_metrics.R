# Onset is evaluated on 1-m grids within observed support, never outside a cast.
# Fit inputs remain raw. The 5-m running median is a scoring operation only.

# Only onset metrics need resolved vertical support; all observations stay in RMSE.
longest_supported_run <- function(good) {
  if(!length(good))return(0L)
  x<-rle(good);if(any(x$values))max(x$lengths[x$values])else 0L
}
onset_curve <- function(depth,value) {
  if(length(depth)!=length(value))stop("Depth/value lengths differ")
  empty<-function()list(depth=numeric(),value=numeric(),support_ok=FALSE)
  if(!length(depth))return(empty())
  if(any(!is.finite(depth))||any(!is.finite(value)))stop("Nonfinite input to onset scorer")
  z<-order(depth);depth<-depth[z];value<-value[z]
  if(anyDuplicated(depth))stop("Duplicate depths in onset evaluation")
  lo<-ceiling(max(10,min(depth)));hi<-floor(min(150,max(depth)))
  if(lo>hi)return(empty())
  grid<-seq.int(lo,hi)
  # approx needs two points; a singleton is retained and flagged below.
  y<-if(length(depth)==1L)rep(value,length(grid))else approx(depth,value,xout=grid,rule=1)$y
  j<-findInterval(grid,depth);between<-j>0L & j<length(depth)
  bad<-between & depth[pmin(j+1L,length(depth))]-depth[pmax(j,1L)]>3 & !grid %in% depth
  y[bad]<-NA_real_;sm<-rep(NA_real_,length(y));good<-which(is.finite(y))
  if(length(good))for(ii in split(good,cumsum(c(TRUE,diff(good)>1L))))
    sm[ii]<-if(length(ii)>=5L)stats::runmed(y[ii],k=5L,endrule="median")else y[ii]
  list(depth=grid,value=sm,support_ok=length(grid)>=7L && longest_supported_run(is.finite(sm))>=7L)
}
onset_from_curve <- function(curve,threshold=.5,persistence=5L) {
  z <- curve$depth;y <- curve$value
  eligible <- which(z>=15)

  if(length(z)<7L || longest_supported_run(is.finite(y))<7L ||
     longest_supported_run(is.finite(y)&z>=15)<persistence)
    return(list(depth=NA_real_,status="insufficient_support",boundary=NA_real_))
  good <- is.finite(y) & y>=threshold & z>=15
  # A run-length calculation keeps posterior onset scoring inexpensive.
  runs <- rle(good);ends <- cumsum(runs$lengths)
  starts <- (ends-runs$lengths+1L)[runs$values & runs$lengths>=persistence]
  if(length(starts)) {
    i <- starts[1]
    # A crossing already present at the upper evaluation boundary is censored.
    if(i==eligible[1])return(list(depth=NA_real_,status="left_censored",boundary=z[i]))
    # A first detection after a data gap cannot locate the onset precisely.
    if(any(!is.finite(y[seq.int(eligible[1],i)])))
      return(list(depth=NA_real_,status="gap_censored",boundary=z[i]))
    return(list(depth=z[i],status="detected",boundary=z[i]))
  }
  list(depth=NA_real_,status=if(any(!is.finite(y[eligible])))"gap_censored" else "not_reached_in_window",
       boundary=NA_real_)
}
rapid_rise <- function(curve) {
  z <- curve$depth;y <- curve$value;n <- length(z)
  ii <- which(z>=17 & z<=145 & seq_along(z)>2L & seq_along(z)<=n-2L)
  if(!length(ii) || sum(is.finite(y))<7L || diff(range(y,na.rm=TRUE))<.5)
    return(c(depth=NA_real_,slope=NA_real_))
  slope <- (y[ii+2L]-y[ii-2L])/4
  valid <- vapply(ii,function(i)all(is.finite(y[(i-2L):(i+2L)])),logical(1))
  slope[!valid] <- NA_real_
  if(!any(is.finite(slope)) || max(slope,na.rm=TRUE)<.02)
    return(c(depth=NA_real_,slope=NA_real_))
  j <- which.max(replace(slope,!is.finite(slope),-Inf))
  c(depth=z[ii[j]],slope=slope[j])
}
onset_profile_metrics <- function(tab,mu=NULL,max_draws=400L) {
  stopifnot(length(unique(tab$cast_key))==1L)
  observed <- onset_curve(tab$depth,tab$chl_norm)
  predicted <- onset_curve(tab$depth,tab$predicted)
  ro <- rapid_rise(observed);rp <- rapid_rise(predicted)
  # Spread the subsample over all supplied posterior draws, including all chains.
  draw_indices <- if(is.null(mu))integer()else unique(as.integer(round(
    seq(1,nrow(mu),length.out=min(max_draws,nrow(mu))))))
  draw_curves <- if(is.null(mu))NULL else lapply(draw_indices,
    function(j)onset_curve(tab$depth,mu[j,]))
  do.call(rbind,lapply(c(.25,.5,1),function(threshold) {
    a <- onset_from_curve(observed,threshold);b <- onset_from_curve(predicted,threshold)
    known <- c("detected","left_censored","not_reached_in_window")
    reached <- function(x)x %in% c("detected","left_censored")
    both <- a$status=="detected" && b$status=="detected"
    draw_onsets <- if(is.null(draw_curves))list()else
      lapply(draw_curves,onset_from_curve,threshold=threshold)
    statuses <- vapply(draw_onsets,function(c)c$status,character(1))
    draws <- vapply(draw_onsets,function(c)c$depth,numeric(1))
    ds <- draws[is.finite(draws)];interval <- if(length(ds))quantile(ds,c(.05,.5,.95))else rep(NA_real_,3)
    data.frame(cast_key=as.character(tab$cast_key[1]),Region=as.character(tab$Region[1]),
      Ocean=as.character(tab$Ocean[1]),date=as.character(tab$date[1]),threshold,

      mld=tab$MLD_temp_MIX[1],observed_rows=nrow(tab),
      minimum_depth=min(tab$depth),maximum_depth=max(tab$depth),
      grid_points=length(observed$depth),
      maximum_supported_run=longest_supported_run(is.finite(observed$value)),
      observed_status=a$status,predicted_status=b$status,
      observed_onset=a$depth,predicted_onset=b$depth,
      observed_onset_minus_mld=a$depth-tab$MLD_temp_MIX[1],
      predicted_onset_minus_mld=b$depth-tab$MLD_temp_MIX[1],
      onset_abs_error=if(both)abs(a$depth-b$depth)else NA_real_,
      onset_signed_error=if(both)b$depth-a$depth else NA_real_,
      reached_agreement=if(a$status %in% known && b$status %in% known)
        reached(a$status)==reached(b$status) else NA,
      false_positive=if(a$status %in% known && b$status %in% known)!reached(a$status)&&reached(b$status)else NA,
      false_negative=if(a$status %in% known && b$status %in% known)reached(a$status)&&!reached(b$status)else NA,
      observed_rapid_rise=unname(ro["depth"]),predicted_rapid_rise=unname(rp["depth"]),
      rapid_rise_abs_error=abs(unname(ro["depth"]-rp["depth"])),
      posterior_draws=length(statuses),posterior_detected_draws=length(ds),
      posterior_insufficient_support_draws=sum(statuses=="insufficient_support"),
      posterior_gap_censored_draws=sum(statuses=="gap_censored"),
      posterior_left_censored_draws=sum(statuses=="left_censored"),
      probability_reached_in_window=if(any(statuses %in% known))
        mean(reached(statuses[statuses %in% known]))else NA_real_,
      onset_q05_given_detected=interval[1],onset_median_given_detected=interval[2],
      onset_q95_given_detected=interval[3],
      onset_interval_covers_observed_given_detected=if(a$status=="detected"&&length(ds))
        a$depth>=interval[1]&&a$depth<=interval[3]else NA)
  }))
}
onset_scores <- function(x) {
  avg <- function(v)if(any(is.finite(v)))mean(v,na.rm=TRUE)else NA_real_
  data.frame(casts=nrow(x),
    observed_insufficient_support=sum(x$observed_status=="insufficient_support"),
    predicted_insufficient_support=sum(x$predicted_status=="insufficient_support"),
    reached_comparable=sum(!is.na(x$reached_agreement)),observed_detected=sum(x$observed_status=="detected"),
    predicted_detected=sum(x$predicted_status=="detected"),
    both_detected=sum(is.finite(x$onset_abs_error)),
    onset_mae_given_both_detected=avg(x$onset_abs_error),
    onset_bias_given_both_detected=avg(x$onset_signed_error),
    reached_agreement=avg(x$reached_agreement),false_positives=sum(x$false_positive,na.rm=TRUE),
    false_negatives=sum(x$false_negative,na.rm=TRUE),
    observed_left_censored=sum(x$observed_status=="left_censored"),
    predicted_left_censored=sum(x$predicted_status=="left_censored"),
    observed_not_reached=sum(x$observed_status=="not_reached_in_window"),
    predicted_not_reached=sum(x$predicted_status=="not_reached_in_window"),
    observed_gap_censored=sum(x$observed_status=="gap_censored"),
    predicted_gap_censored=sum(x$predicted_status=="gap_censored"),
    rapid_rise_comparable=sum(is.finite(x$rapid_rise_abs_error)),
    rapid_rise_mae=avg(x$rapid_rise_abs_error))
}
save_onset_metrics <- function(tab,mu,checks,tag) {
  indices <- split(seq_len(nrow(tab)),as.character(tab$cast_key))
  out <- do.call(rbind,lapply(indices,function(ii)
    onset_profile_metrics(tab[ii,,drop=FALSE],mu[,ii,drop=FALSE])))
  write.csv(out,file.path(checks,paste0(tag,"_onset_cast_metrics.csv")),row.names=FALSE)
  aggregate <- do.call(rbind,lapply(split(out,out$threshold),function(x)
    data.frame(threshold=x$threshold[1],onset_scores(x))))
  regional <- do.call(rbind,lapply(split(out,interaction(out$threshold,out$Region,drop=TRUE)),
    function(x)data.frame(threshold=x$threshold[1],Region=x$Region[1],onset_scores(x))))
  write.csv(aggregate,file.path(checks,paste0(tag,"_onset_scores.csv")),row.names=FALSE)
  write.csv(regional,file.path(checks,paste0(tag,"_onset_by_region.csv")),row.names=FALSE)
  # Explicitly score the MLD neighbourhood; equal total weight for each eligible cast.
  local <- tab[abs(tab$depth-tab$MLD_temp_MIX)<=20,,drop=FALSE]
  if(nrow(local)) {
    cs <- do.call(rbind,lapply(split(local,local$cast_key),function(x)data.frame(
      cast_key=x$cast_key[1],Region=x$Region[1],observations=nrow(x),
      mse=mean((x$predicted-x$chl_norm)^2))))
    write.csv(cs,file.path(checks,paste0(tag,"_mld_window_cast_metrics.csv")),row.names=FALSE)
  }
  # Use a predictor-central cast per Region, selected without response/error values.
  metadata <- unique(tab[c("cast_key","Region","MLD_temp_MIX","log10_distance")])
  eligible <- unique(out$cast_key[out$observed_status!="insufficient_support" &
                                  out$predicted_status!="insufficient_support"])
  metadata <- metadata[as.character(metadata$cast_key) %in% eligible,,drop=FALSE]
  if(!nrow(metadata)) {
    writeLines("No onset-resolved profiles available for this example plot",
      file.path(checks,paste0(tag,"_onset_plot_skipped.txt")))
    return(invisible(aggregate))
  }
  chosen <- unlist(lapply(split(metadata,metadata$Region,drop=TRUE),function(x) {
    spread <- function(v)if(length(v)>1L&&sd(v)>0)sd(v)else 1
    score <- abs(x$MLD_temp_MIX-median(x$MLD_temp_MIX))/spread(x$MLD_temp_MIX)+
      abs(x$log10_distance-median(x$log10_distance))/spread(x$log10_distance)
    as.character(x$cast_key[which.min(score)])
  }),use.names=FALSE)
  curves <- do.call(rbind,lapply(chosen,function(key) {
    x <- tab[as.character(tab$cast_key)==key,,drop=FALSE]
    obs <- onset_curve(x$depth,x$chl_norm);pred <- onset_curve(x$depth,x$predicted)
    data.frame(cast_key=key,Region=as.character(x$Region[1]),MLD=x$MLD_temp_MIX[1],
      depth=obs$depth,observed=obs$value,predicted=pred$value)
  }))
  curves$panel <- paste(curves$Region,curves$cast_key,sep=" | ")
  points <- out[out$cast_key %in% chosen & out$threshold==.5,,drop=FALSE]
  points$panel <- paste(points$Region,points$cast_key,sep=" | ")
  write.csv(curves,file.path(checks,paste0(tag,"_onset_example_profiles.csv")),row.names=FALSE)
  p <- ggplot2::ggplot(curves,ggplot2::aes(predicted,depth))+
    ggplot2::geom_path(colour="#005F73",linewidth=.7)+
    ggplot2::geom_path(ggplot2::aes(x=observed),colour="#222222",linewidth=.5)+
    ggplot2::geom_hline(data=unique(curves[c("panel","MLD")]),
      ggplot2::aes(yintercept=MLD),linetype="dashed",colour="#CA6702")+
    ggplot2::geom_point(data=points[points$observed_status=="detected",,drop=FALSE],
      ggplot2::aes(x=.5,y=observed_onset),inherit.aes=FALSE,colour="#222222",size=2)+
    ggplot2::geom_point(data=points[points$predicted_status=="detected",,drop=FALSE],
      ggplot2::aes(x=.5,y=predicted_onset),inherit.aes=FALSE,colour="#005F73",shape=17,size=2)+
    ggplot2::scale_y_reverse()+ggplot2::facet_wrap(~panel,scales="free_x",ncol=3)+
    ggplot2::theme_bw()+ggplot2::labs(x="chl_norm",y="Depth (m)",
      title=paste(tag,"enhancement onset"),
      subtitle="Black: observed; blue: expected; orange dashed: MLD; points: 50% onset; 5-m median for scoring")
  ggplot2::ggsave(file.path(checks,paste0(tag,"_onset_example_profiles.png")),p,
    width=13,height=10,dpi=130,bg="white")
  invisible(aggregate)
}

make_mld_sensitivity <- function(fit,experiment,run,case,draw_ids) {
  # Hypothetical population profiles within regional marginal predictor ranges.
  # These are model associations, not causal effects or independent validation.
  train <- readRDS(data_file(experiment,case))
  grid <- do.call(rbind,lapply(split(train,train$Region,drop=TRUE),function(d) {
    meta <- unique(d[c("cast_key","MLD_temp_MIX","log10_distance")])
    values <- as.numeric(quantile(meta$MLD_temp_MIX,c(.25,.5,.75)))
    do.call(rbind,lapply(seq_along(values),function(j) {
      x <- d[rep(1L,141L),,drop=FALSE];x$depth <- 10:150
      x$MLD_temp_MIX <- values[j];x$log10_distance <- median(meta$log10_distance)
      x$mld_scenario <- c("MLD Q25","MLD Q50","MLD Q75")[j];x
    }))
  }))
  grid <- apply_features(grid,readRDS(file.path(experiment,"scaling.rds"))[[case$split]],
    readRDS(file.path(experiment,"region_levels.rds")))
  # Keep a balanced spread over the saved prediction draw IDs.
  use <- draw_ids[unique(as.integer(round(seq(1,length(draw_ids),length.out=min(400L,length(draw_ids))))))]
  mu <- brms::posterior_epred(fit,newdata=grid[setdiff(names(grid),"chl_norm")],
    re_formula=NA,draw_ids=use,allow_new_levels=TRUE)
  q <- apply(mu,2,quantile,c(.05,.95))
  tab <- data.frame(grid[c("Region","depth","MLD_temp_MIX","log10_distance","mld_scenario")],
    predicted=colMeans(mu),lower90=q[1,],upper90=q[2,])
  write.csv(tab,file.path(run,"checks","population_mld_sensitivity.csv"),row.names=FALSE)
  p <- ggplot2::ggplot(tab,ggplot2::aes(predicted,depth,colour=mld_scenario,group=mld_scenario))+
    ggplot2::geom_path(linewidth=.7)+ggplot2::scale_y_reverse()+
    ggplot2::facet_wrap(~Region,scales="free_x",ncol=3)+ggplot2::theme_bw()+
    ggplot2::labs(x="Expected chl_norm",y="Depth (m)",colour=NULL,
      title="Population profile sensitivity to MLD",
      subtitle="MLD at regional training quartiles; distance fixed at regional training median; hypothetical combinations")
  ggplot2::ggsave(file.path(run,"checks","population_mld_sensitivity.png"),p,
    width=13,height=10,dpi=130,bg="white")
}
