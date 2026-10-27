"""Build the MLD-onset campaign without starting R, Stan, SSH, or network clients.

The existing frozen follow-up implementation supplies sampler checks and plotting.
All new source is written to a separate package; previous campaigns are untouched.
"""
from pathlib import Path, PurePosixPath
import base64, csv, hashlib, io, json, tarfile, textwrap

ROOT = Path(__file__).resolve().parent
OLD = ROOT / 'DCM_Ibex_followup_tests'
PKG = ROOT / 'DCM_Ibex_onset_tests'
PKG.mkdir(exist_ok=True)

def put(name, text):
    (PKG / name).write_text(text.strip() + '\n', encoding='utf-8', newline='\n')

for name in ['common.R', 'base_models.R', '02_fit_case.R', '03_check_case.R',
             'fit.slurm', 'check.slurm', 'prepare.slurm', 'collect.slurm', 'live_diagnostics.awk']:
    text = (OLD / name).read_text(encoding='utf-8-sig')
    text = text.replace('mode = "followup"', 'mode = "onset"')
    text = text.replace('Follow-up campaign', 'MLD-onset campaign')
    text = text.replace('all 34 requested models', 'all requested MLD-onset models')
    if name == 'check.slurm':
        text = text.replace('--time=08:00:00', '--time=16:00:00')
    put(name, text)

# Keep data transformations and priors identical to the existing standardized models.
put('campaign.R', r'''
# Six source controls are rescored; only execution="fit" cases start new models.
campaign_cases <- function() {
  rows <- list()
  add <- function(eq, subset="p10", split="primary", block="equation", reference=1L,
                  execution="fit", source_case=NA_integer_, seed=1234L,
                  delta=.95, warmup=1000L, samples=2000L, kd=7L, mode="posterior") {
    id <- length(rows)+1L
    threads <- if(execution=="reuse" || mode=="prior")2L else 4L
    rows[[id]] <<- data.frame(case_id=id,
      label=sprintf("o%03d_%s_%s_%s_%s",id,block,eq,split,subset), block,
      equation=eq,subset,split,prior="custom_z",mode,variance="region",cast="none",
      family="gaussian",autocor="none",depth_max=150L,chains=4L,
      iter=warmup+samples,warmup,samples,adapt_delta=delta,treedepth=14L,
      threads,cpus=4L*threads,k_depth=kd,k_mld=5L,k_distance=5L,
      sigma_smooth_sd=.5,cores=4L,refresh=50L,seed,
      memory_gb=if(subset=="p20")64L else 32L,
      reference_case=reference,execution,source_case)
    id
  }
  panels <- list(c("p10","primary"),c("p20","primary"),c("p10","date_blocked"))
  src <- c(1L,5L,13L,15L,22L,24L)
  for(j in seq_along(panels))for(eq in c("original_z","interaction_z")) {
    i <- length(rows)+1L
    add(eq,panels[[j]][1],panels[[j]][2],"existing_control",
        reference=if(i%%2L==1L)NA_integer_ else i-1L,
        execution="reuse",source_case=src[i])
  }
  equations <- c("mld_vary_z","mld_shift","mld_halfshift","mld_ratio",
                 "mld_hybrid","mld_shift_distance")
  for(j in seq_along(panels))for(eq in equations)
    add(eq,panels[[j]][1],panels[[j]][2],reference=2L*j-1L)
  # Isolate sampling stability, robust tuning and basis size for two new anchors.
  for(eq in c("mld_shift","mld_hybrid")) {
    ref <- if(eq=="mld_shift")8L else 11L
    add(eq,block="repeat_seed",reference=ref,seed=20261019L)
    add(eq,block="robust_sampler",reference=ref,delta=.99,warmup=2000L)
    add(eq,block="depth_basis",reference=ref,kd=10L)
  }
  for(eq in c("mld_shift","mld_hybrid"))for(seed in c(1234L,20261019L))
    add(eq,block="prior_predictive",reference=NA_integer_,seed=seed,
        mode="prior",warmup=500L,samples=500L)
  x <- do.call(rbind,rows)
  stopifnot(nrow(x)==34L,sum(x$execution=="fit" & x$mode=="posterior")==24L,
    sum(x$mode=="prior")==4L,sum(x$execution=="reuse")==6L,
    !anyDuplicated(x$label),all(x$subset %in% c("p10","p20")),
    all(x$cpus==x$chains*x$threads))
  x
}
.base_models <- new.env(parent=environment())
base_models <- function() {
  if(!exists("campaign_spec",envir=.base_models,inherits=FALSE))
    sys.source(file.path(Sys.getenv("EXPERIMENT_DIR",Sys.getenv("ONSET_PACKAGE")),
                         "base_models.R"),envir=.base_models)
  .base_models
}
apply_features <- function(d,scaling,region_levels) {
  d <- base_models()$apply_features(d,scaling,region_levels)
  if(any(!is.finite(d$MLD_temp_MIX)) || any(d$MLD_temp_MIX<=0))
    stop("MLD must be finite and positive for the predeclared depth/MLD tests; no rows dropped")
  # All transformations use predictors only and the observed MLD in metres.
  d$depth_mld <- d$depth-d$MLD_temp_MIX
  d$depth_mld_half <- d$depth-.5*d$MLD_temp_MIX
  d$depth_mld_ratio <- d$depth/d$MLD_temp_MIX
  d
}
campaign_spec <- function(case,region_levels) {
  if(case$equation %in% c("original_z","interaction_z","mld_vary_z"))
    return(base_models()$campaign_spec(case,region_levels))
  # Obtain exactly the same coefficient, smooth and random-effect priors.
  surrogate <- case;surrogate$equation <- "original_z"
  spec <- base_models()$campaign_spec(surrogate,region_levels)
  kd <- case$k_depth;ke <- case$k_mld
  region_s <- function(x)sprintf("s(%s,by=Region,bs='cr',k=%d)",x,kd)
  terms <- switch(case$equation,
    mld_shift=c("distance_z",region_s("depth_mld")),
    mld_halfshift=c("distance_z",region_s("depth_mld_half")),
    mld_ratio=c("distance_z",region_s("depth_mld_ratio")),
    mld_hybrid=c("distance_z",region_s("depth"),
      sprintf("s(depth_mld,bs='cr',k=%d)",ke)),
    mld_shift_distance=c(region_s("depth_mld"),
      sprintf("s(depth_mld,by=distance_z,bs='cr',k=%d)",ke)),
    stop("Unknown onset equation"))
  f <- as.formula(paste("chl_norm ~ 1 + Region +",paste(c(terms,
      "(1|Ocean:Region:date)","(1|Ocean:Region:enso)"),collapse=" + ")),
      env=asNamespace("brms"))
  spec$formula <- brms::bf(f,sigma~0+Region,family=gaussian())
  spec
}
data_file <- function(experiment,case)file.path(experiment,"data",
  paste0(case$split,"_",case$subset,"_d150.rds"))
training_eval_file <- function(experiment,case)file.path(experiment,"data",
  paste0(case$split,"_p10_d150.rds"))
validation_file <- function(experiment,case)file.path(experiment,"data",
  paste0(case$split,"_validation_d150.rds"))
reserve_validation_file <- function(experiment,case)file.path(experiment,"data",
  paste0(case$split,"_reserve_d150.rds"))
thread_specification <- function(case)brms::threading(case$threads)
''')

put('data_setup.R', r'''
# Copy previously frozen follow-up casts exactly; never choose casts from outcomes.
prepare_data <- function(experiment,source_experiment) {
  dir.create(file.path(experiment,"data"),recursive=TRUE,showWarnings=FALSE)
  required <- c("scaling.rds","region_levels.rds","profile_reference.rds",
    "full_cast_metadata.csv","profile_reference.csv","depth_spacing.csv","split_membership.csv")
  names <- c("primary_p10_d150","primary_p20_d150","primary_validation_d150",
    "primary_reserve_d150","date_blocked_p10_d150","date_blocked_validation_d150",
    "date_blocked_reserve_d150")
  stopifnot(all(file.exists(file.path(source_experiment,required))))
  stopifnot(all(file.copy(file.path(source_experiment,required),
                         file.path(experiment,required),overwrite=FALSE)))
  scaling <- readRDS(file.path(experiment,"scaling.rds"))
  levels <- readRDS(file.path(experiment,"region_levels.rds"))
  raw <- list();audits <- list()
  raw_cols <- c("cast_key","Region","Ocean","date","enso","depth",
                "log10_distance","MLD_temp_MIX","chl_norm")
  for(n in names) {
    path <- file.path(source_experiment,"data",paste0(n,".rds"))
    before <- readRDS(path);raw[[n]] <- before
    split <- if(startsWith(n,"date_blocked"))"date_blocked" else "primary"
    after <- apply_features(before,scaling[[split]],levels)
    # Row order, response, predictors and group membership must remain identical.
    stopifnot(identical(before[raw_cols],after[raw_cols]))
    saveRDS(after,file.path(experiment,"data",paste0(n,".rds")))
    audits[[n]] <- data.frame(dataset=n,source_file=path,
      source_md5=unname(tools::md5sum(path)),casts=length(unique(after$cast_key)),
      observations=nrow(after),minimum_mld=min(after$MLD_temp_MIX),
      source_raw_rows_verified=TRUE)
  }
  full <- read.csv(file.path(experiment,"full_cast_metadata.csv"),stringsAsFactors=FALSE)
  minimum <- ceiling(.1*length(unique(full$cast_key)))
  stopifnot(all(vapply(raw,function(d)length(unique(d$cast_key)),integer(1))>=minimum))
  for(split in c("primary","date_blocked")) {
    keys <- function(s)unique(raw[[paste0(split,"_",s,"_d150")]]$cast_key)
    training <- keys(if(split=="primary")"p20" else "p10")
    stopifnot(!length(intersect(training,keys("validation"))),
      !length(intersect(training,keys("reserve"))),
      !length(intersect(keys("validation"),keys("reserve"))))
  }
  date_key <- function(d)paste(d$Ocean,d$Region,d$date,sep="\r")
  stopifnot(!any(date_key(raw$date_blocked_p10_d150) %in%
                  date_key(raw$date_blocked_validation_d150)))
  write.csv(do.call(rbind,audits),file.path(experiment,"data_identity_audit.csv"),row.names=FALSE)
  invisible(do.call(rbind,audits))
}
''')

put('01_prepare.R', r'''
# Verify data identities, new model matrices and source controls before any fit.
experiment <- Sys.getenv("EXPERIMENT_DIR")
source(file.path(experiment,"common.R"));source(file.path(experiment,"campaign.R"))
source(file.path(experiment,"data_setup.R"))
tryCatch({
  if(file.exists(file.path(experiment,"PREPARATION_COMPLETE.txt")))stop("Already prepared")
  origin <- readLines(file.path(experiment,"source_experiment.txt"),n=1L)
  prepare_data(experiment,origin)
  cases <- read.delim(file.path(experiment,"cases.tsv"),stringsAsFactors=FALSE)
  stopifnot(isTRUE(all.equal(cases,campaign_cases(),check.attributes=FALSE)))
  for(p in c("brms","mgcv","cmdstanr","posterior"))
    if(!requireNamespace(p,quietly=TRUE))stop("Missing package: ",p)
  controls <- read.delim(file.path(experiment,"source_controls.tsv"),stringsAsFactors=FALSE)
  source_cases <- read.delim(file.path(origin,"cases.tsv"),stringsAsFactors=FALSE)
  settings <- c("equation","subset","split","prior","mode","variance","cast","family",
    "autocor","depth_max","chains","iter","warmup","adapt_delta","treedepth",
    "threads","k_depth","k_mld","k_distance","seed")
  for(i in which(cases$execution=="reuse")) {
    c <- cases[i,];old <- source_cases[source_cases$case_id==c$source_case,,drop=FALSE]
    stopifnot(nrow(old)==1L)
    for(n in settings)if(!isTRUE(all.equal(c[[n]],old[[n]],check.attributes=FALSE)))
      stop("Source control differs for ",c$case_id," ",n)
    stopifnot(sum(controls$case_id==c$case_id)==1L)
  }
  dir.create(file.path(experiment,"preflight"),showWarnings=FALSE)
  basis <- list();levels <- readRDS(file.path(experiment,"region_levels.rds"))
  for(i in seq_len(nrow(cases))) {
    c <- cases[i,];d <- readRDS(data_file(experiment,c));spec <- campaign_spec(c,levels)
    sd <- brms::make_standata(spec$formula,data=d,prior=spec$priors,
      threads=thread_specification(c),sample_prior=if(c$mode=="prior")"only" else "no")
    code <- brms::make_stancode(spec$formula,data=d,prior=spec$priors,
      threads=thread_specification(c),sample_prior=if(c$mode=="prior")"only" else "no")
    writeLines(code,file.path(experiment,"preflight",sprintf("case_%02d.stan",c$case_id)))
    write.csv(as.data.frame(brms::get_prior(spec$formula,data=d)),
      file.path(experiment,"preflight",sprintf("case_%02d_priors.csv",c$case_id)),row.names=FALSE)
    X <- sd$X
    if(is.matrix(X) && qr(X)$rank<ncol(X))stop("Rank-deficient fixed design in case ",c$case_id)
    basis[[i]] <- data.frame(case_id=c$case_id,fixed_columns=if(is.matrix(X))ncol(X)else NA_integer_,
      penalized_columns=sum(vapply(sd[grepl("^Zs",names(sd))],function(x)
        if(is.matrix(x))ncol(x)else 0L,integer(1))))
    cat("Validated onset case",c$case_id,c$label,"\n");rm(sd);invisible(gc())
  }
  write.csv(do.call(rbind,basis),file.path(experiment,"basis_counts.csv"),row.names=FALSE)
  files <- list.files(file.path(experiment,"data"),full.names=TRUE)
  write.csv(data.frame(file=basename(files),md5=unname(tools::md5sum(files))),
    file.path(experiment,"data_checksums.csv"),row.names=FALSE)
  writeLines("34 formula/data combinations validated; 28 new fits, 6 controls reused; no response changes",
             file.path(experiment,"PREPARATION_COMPLETE.txt"))
},error=function(e){writeLines(conditionMessage(e),file.path(experiment,"PREPARATION_FAILED.txt"));stop(e)})
''')

# The onset module is also usable on any older saved observation predictions.
put('onset_metrics.R', r'''
# Onset is evaluated on 1-m grids within observed support, never outside a cast.
# Fit inputs remain raw. The 5-m running median is a scoring operation only.
onset_curve <- function(depth,value) {
  z <- order(depth);depth <- depth[z];value <- value[z]
  if(anyDuplicated(depth))stop("Duplicate depths in onset evaluation")
  grid <- seq.int(ceiling(max(10,min(depth))),floor(min(150,max(depth))))
  if(length(grid)<7L)stop("Insufficient vertical support")
  y <- approx(depth,value,xout=grid,rule=1)$y
  # Do not interpolate across gaps larger than 3 metres.
  j <- findInterval(grid,depth)
  between <- j>0L & j<length(depth)
  bad <- between & depth[pmin(j+1L,length(depth))]-depth[pmax(j,1L)]>3 &
    !grid %in% depth
  y[bad] <- NA_real_
  sm <- rep(NA_real_,length(y))
  good <- which(is.finite(y))
  if(length(good))for(ii in split(good,cumsum(c(TRUE,diff(good)>1L)))) {
    sm[ii] <- if(length(ii)>=5L)stats::runmed(y[ii],k=5L,endrule="median") else y[ii]
  }
  list(depth=grid,value=sm)
}
onset_from_curve <- function(curve,threshold=.5,persistence=5L) {
  z <- curve$depth;y <- curve$value
  eligible <- which(z>=15)
  if(!length(eligible))return(list(depth=NA_real_,status="insufficient_support"))
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
      mld=tab$MLD_temp_MIX[1],observed_status=a$status,predicted_status=b$status,
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
  data.frame(casts=nrow(x),observed_detected=sum(x$observed_status=="detected"),
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
  chosen <- unlist(lapply(split(metadata,metadata$Region),function(x) {
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
  grid <- do.call(rbind,lapply(split(train,train$Region),function(d) {
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
''')

pred = (OLD / 'predictive.R').read_text(encoding='utf-8-sig')
pred = 'source(file.path(Sys.getenv("EXPERIMENT_DIR"),"onset_metrics.R"))\n' + pred
needle = '    cast<-do.call(rbind,lapply(split(tab,as.character(tab$cast_key)),profile_metrics))'
assert pred.count(needle) == 1
pred = pred.replace(needle, '''    # Evaluate threshold onset and the steepest sustained rise from expected curves.
    # Score held-out casts in both prediction modes; never change training responses.
    if(dataset!="training")save_onset_metrics(tab,mu,checks,tag)
''' + needle)
pred = pred.replace('  make_original_style_plots(fit,experiment,run,case,draw_ids)',
  '  make_original_style_plots(fit,experiment,run,case,draw_ids)\n  make_mld_sensitivity(fit,experiment,run,case,draw_ids)')
put('predictive.R', pred)

put('reuse_check.R', r'''
# Read an existing fit and write new predictions/metrics only to this campaign.
experiment <- Sys.getenv("EXPERIMENT_DIR")
source(file.path(experiment,"common.R"));source(file.path(experiment,"campaign.R"))
if(!file.exists(file.path(experiment,"PREPARATION_COMPLETE.txt")))stop("Preparation not complete")
case <- read_case(experiment,Sys.getenv("CASE_ID"))
stopifnot(case$execution=="reuse",case$mode=="posterior")
registry <- read.delim(file.path(experiment,"source_controls.tsv"),stringsAsFactors=FALSE)
item <- registry[registry$case_id==case$case_id,,drop=FALSE];stopifnot(nrow(item)==1L)
run <- file.path(experiment,"runs",sprintf("case_%02d",case$case_id))
dir.create(run,recursive=TRUE,showWarnings=FALSE)
if(file.exists(file.path(run,"CHECKS_COMPLETE.txt")))quit(status=0L)
tryCatch({
  if(!file.exists(file.path(item$source_run,"FIT_COMPLETE.txt")))
    stop("Source fitting job ended without a complete fit; it was not resubmitted")
  fit <- readRDS(file.path(item$source_run,"fit.rds"))
  old <- readRDS(file.path(item$source_run,"model_data.rds"))$model_data
  current <- readRDS(data_file(experiment,case))
  raw <- c("cast_key","Region","Ocean","date","enso","depth","log10_distance","MLD_temp_MIX","chl_norm")
  stopifnot(identical(old[raw],current[raw]),all(names(old) %in% names(current)),
            identical(old,current[names(old)]))
  # Hydrated saved fits are used directly; no brm/update/compilation call occurs here.
  result <- data.frame(case,source_run=item$source_run,source_fit_job=item$fit_job,
    evaluation_only=TRUE,raw_training_identity_verified=TRUE)
  for(name in c("test_summary.csv","diagnostic_summary.csv")) {
    path <- file.path(item$source_run,name)
    if(file.exists(path)) {
      x <- read.csv(path,stringsAsFactors=FALSE)
      keep <- intersect(names(x),c("max_rhat","min_bulk_ess","min_tail_ess","divergences",
        "treedepth_hits","min_ebfmi","fit_seconds","requested_core_hours",
        "diagnostic_screen_pass","convergence_screen_pass","cpu_model"))
      for(n in keep)result[[n]] <- x[[n]]
      break
    }
  }
  source(file.path(experiment,"predictive.R"))
  result <- evaluate_predictions(fit,experiment,run,case,result)
  write.csv(result,file.path(run,"test_summary.csv"),row.names=FALSE)
  writeLines(item$source_run,file.path(run,"REUSED_FIT.txt"))
  writeLines("Source fit reused without modification",file.path(run,"FIT_COMPLETE.txt"))
  writeLines("Complete",file.path(run,"CHECKS_COMPLETE.txt"))
},error=function(e){writeLines(conditionMessage(e),file.path(run,"CHECKS_FAILED.txt"));stop(e)})
''')
reuse_slurm = (PKG / 'check.slurm').read_text().replace('03_check_case.R', 'reuse_check.R')
put('reuse.slurm', reuse_slurm)

# Reuse the recovery ledger, dependency-race handling and all-campaign collector.
submit = (OLD / 'submit.R').read_text(encoding='utf-8-sig')
submit = submit.replace('FOLLOWUP', 'ONSET').replace('followup150', 'onset150')
submit = submit.replace('dcm-u150-', 'dcm-o150-').replace('2026-10-18', '2026-10-19')
submit = submit.replace('DCM_Ibex_followup_tests/monitor_all_dcm_v6.sh',
                        'DCM_Ibex_onset_tests/monitor_all_dcm_v7.sh')
submit = submit.replace('latest_targeted150_experiment.txt', 'latest_followup150_experiment.txt')
submit = submit.replace('"full_input_250m.rds",', '"split_membership.csv",')
submit = submit.replace('"primary_p50_d150.rds",', '"primary_reserve_d150.rds","date_blocked_reserve_d150.rds",')
submit = submit.replace('successfully prepared targeted150 source campaign',
                        'successfully prepared followup150 source campaign')
start = submit.index('    fit<-ensure(sprintf("fit_%02d",case$case_id)')
end = submit.index('    j<-match(case$case_id,manifest$case_id)', start)
submit = submit[:start] + r'''
    if(case$execution=="reuse") {
      controls<-read.delim(file.path(experiment,"source_controls.tsv"),stringsAsFactors=FALSE)
      source<-controls[controls$case_id==case$case_id,,drop=FALSE]
      stopifnot(nrow(source)==1L)
      fit<-as.character(source$fit_job)
    } else {
      fit<-ensure(sprintf("fit_%02d",case$case_id),sprintf("c%02d",case$case_id),"fit.slurm",
        c(prep_dep(),paste0("--cpus-per-task=",case$cpus),paste0("--mem=",case$memory_gb,"G"),
          paste0("--time=",if(case$mode=="prior")"02:00:00"else Sys.getenv("ONSET_TIME","7-00:00:00")),
          paste0("--constraint=",Sys.getenv("ONSET_CONSTRAINT","intel&cascadelake")),"--hint=nomultithread"))
    }
''' + submit[end:]
start = submit.index('    finished<-any(file.exists(')
end = submit.index('    manifest$check_job[j]<-check;', start)
submit = submit[:start] + r'''
    if(case$execution=="reuse") {
      # afterok preparation and afterany source checks, so failure is reported.
      work<-ids(c(source$fit_job,source$check_job))
      dep_reuse<-function() {
        # A missing squeue row alone does not prove the source fit has finished.
        active_work<-active(work)
        if(!file.exists(file.path(source$source_run,"FIT_COMPLETE.txt")) && !terminal(fit))
          active_work<-ids(c(active_work,fit))
        all<-c(if(!file.exists(file.path(experiment,"PREPARATION_COMPLETE.txt")))paste0("afterok:",prep),
          if(length(active_work))paste0("afterany:",paste(active_work,collapse=":")))
        if(length(all))paste0("--dependency=",paste(all,collapse=","))else character()
      }
      for(attempt in 1:3) {
        check<-tryCatch(ensure(sprintf("check_%02d",case$case_id),sprintf("q%02d",case$case_id),
          "reuse.slurm",dep_reuse()),error=function(e)e)
        if(!inherits(check,"error"))break
        invalid<-inherits(check,"definite_rejection")&&grepl("Job dependency problem|Invalid job id",conditionMessage(check))
        if(!invalid||attempt==3L)stop(check)
      }
    } else {
      finished<-any(file.exists(file.path(manifest$run_dir[j],c("FIT_COMPLETE.txt","FIT_FAILED.txt"))))
      dep<-if(finished)character()else paste0("--dependency=afterany:",fit)
      check<-tryCatch(ensure(sprintf("check_%02d",case$case_id),sprintf("q%02d",case$case_id),"check.slurm",dep),error=function(e)e)
      if(inherits(check,"error")) {
        invalid<-inherits(check,"definite_rejection")&&grepl("Job dependency problem|Invalid job id",conditionMessage(check))
        if(!invalid||!terminal(fit))stop(check)
        check<-ensure(sprintf("check_%02d",case$case_id),sprintf("q%02d",case$case_id),"check.slurm")
      }
    }
''' + submit[end:]
needle = '  writeLines(source,file.path(experiment,"source_experiment.txt"))'
assert submit.count(needle)==1
submit = submit.replace(needle, needle + r'''
  controls<-read.delim(file.path(package,"cases.tsv"),stringsAsFactors=FALSE)
  controls<-controls[controls$execution=="reuse",,drop=FALSE]
  origin<-read.delim(file.path(source,"manifest.tsv"),stringsAsFactors=FALSE)
  entries<-lapply(seq_len(nrow(controls)),function(i){
    old<-origin[origin$case_id==controls$source_case[i],,drop=FALSE]
    if(nrow(old)!=1L)stop("Missing unique registered source control ",controls$source_case[i])
    data.frame(case_id=controls$case_id[i],source_case=controls$source_case[i],
      fit_job=as.character(old$fit_job),check_job=as.character(old$check_job),source_run=as.character(old$run_dir))
  })
  atomic_table(do.call(rbind,entries),file.path(experiment,"source_controls.tsv"))
''')
put('submit.R', submit)

collect = (OLD / '04_collect_results.R').read_text(encoding='utf-8-sig')
collect = collect.replace('if(sys.nframe()==0L){a<-commandArgs(trailingOnly=TRUE);collect_results(a[1])}', '')
collect = collect.replace('"followup_results.csv"','"onset_results.csv"')
put('collect_all_runs.R', collect)
put('04_collect_results.R', r'''
# Preserve every earlier registered run, plus focused onset comparisons for this campaign.
a <- commandArgs(trailingOnly=TRUE)
experiment <- if(length(a))a[1]else Sys.getenv("EXPERIMENT_DIR")
source(file.path(experiment,"collect_all_runs.R"))
source(file.path(experiment,"onset_metrics.R"))
collect_onset <- function(experiment) {
  cases <- read.delim(file.path(experiment,"cases.tsv"),stringsAsFactors=FALSE)
  output <- file.path(experiment,"results");dir.create(output,showWarnings=FALSE)
  casts <- list();aggregates <- list();local <- list();status <- list();contrasts <- list()
  for(i in seq_len(nrow(cases))) {
    c <- cases[i,];run <- file.path(experiment,"runs",sprintf("case_%02d",c$case_id))
    c$status <- status_of(run);status[[i]] <- c
    for(dataset in c("validation","reserve"))for(mode in c("conditional","population")) {
      tag <- paste(dataset,mode,sep="_")
      z <- read_table(file.path(run,"checks",paste0(tag,"_onset_cast_metrics.csv")))
      if(nrow(z)) {
        z$case_id <- c$case_id;z$equation <- c$equation;z$subset <- c$subset
        z$split <- c$split;z$dataset <- dataset;z$prediction_mode <- mode
        casts[[length(casts)+1L]] <- z
        for(threshold in unique(z$threshold)) {
          zz <- z[z$threshold==threshold,,drop=FALSE]
          aggregates[[length(aggregates)+1L]] <- data.frame(case_id=c$case_id,
            equation=c$equation,subset=c$subset,split=c$split,dataset,prediction_mode=mode,
            threshold,onset_scores(zz))
        }
      }
      w <- read_table(file.path(run,"checks",paste0(tag,"_mld_window_cast_metrics.csv")))
      if(nrow(w))local[[length(local)+1L]] <- data.frame(case_id=c$case_id,equation=c$equation,
        subset=c$subset,split=c$split,dataset,prediction_mode=mode,
        casts=nrow(w),mld_neighbourhood_rmse=sqrt(mean(w$mse)))
    }
  }
  cs <- fill_rows(casts)
  write.csv(fill_rows(status),file.path(output,"onset_case_status.csv"),row.names=FALSE)
  write.csv(cs,file.path(output,"onset_all_cast_metrics.csv"),row.names=FALSE)
  write.csv(fill_rows(aggregates),file.path(output,"onset_all_scores.csv"),row.names=FALSE)
  write.csv(fill_rows(local),file.path(output,"mld_neighbourhood_scores.csv"),row.names=FALSE)
  if(nrow(cs))for(i in which(cases$mode=="posterior")) {
    c <- cases[i,];references <- unique(na.omit(c(c$reference_case,
      if(c$execution=="fit" && c$block=="equation")if(c$split=="date_blocked")6L else if(c$subset=="p20")4L else 2L)))
    for(ref in references)for(dataset in c("validation","reserve"))
      for(mode in c("conditional","population"))for(t in c(.25,.5,1)) {
        aa <- cs[cs$case_id==c$case_id & cs$dataset==dataset & cs$prediction_mode==mode & cs$threshold==t,,drop=FALSE]
        bb <- cs[cs$case_id==ref & cs$dataset==dataset & cs$prediction_mode==mode & cs$threshold==t,,drop=FALSE]
        if(!nrow(aa)||!nrow(bb))next
        bb <- bb[match(aa$cast_key,bb$cast_key),,drop=FALSE]
        if(anyNA(bb$cast_key)||!identical(as.character(aa$cast_key),as.character(bb$cast_key)))next
        # Fair onset-MAE contrasts use the SAME casts detected by all three curves.
        eligible <- is.finite(aa$onset_abs_error)&is.finite(bb$onset_abs_error)
        avg <- function(v)if(length(v))mean(v)else NA_real_
        contrasts[[length(contrasts)+1L]] <- data.frame(case_id=c$case_id,reference_case=ref,
          dataset,prediction_mode=mode,threshold=t,matched_casts=nrow(aa),
          jointly_detected_casts=sum(eligible),
          delta_onset_mae_common_casts=avg(aa$onset_abs_error[eligible]-bb$onset_abs_error[eligible]),
          candidate_false_positives=sum(aa$false_positive,na.rm=TRUE),
          reference_false_positives=sum(bb$false_positive,na.rm=TRUE),
          candidate_false_negatives=sum(aa$false_negative,na.rm=TRUE),
          reference_false_negatives=sum(bb$false_negative,na.rm=TRUE))
      }
  }
  write.csv(fill_rows(contrasts),file.path(output,"onset_planned_paired_contrasts.csv"),row.names=FALSE)
  cat("Onset scores retain missing/censored detections explicitly; no automatic single-score winner.\n")
}
if(sys.nframe()==0L){collect_results(experiment);collect_onset(experiment)}
''')

shell = (OLD/'submit_tests.sh').read_text(encoding='utf-8-sig')
shell = shell.replace('FOLLOWUP','ONSET').replace('followup150','onset150')
shell = shell.replace('follow-up','MLD-onset').replace('latest_targeted150','latest_followup150')
put('submit_tests.sh', shell)
put('collect_results.sh', r'''
#!/bin/bash
# Refresh completed results; this never starts a model fit.
set -euo pipefail
cd "${DCM_ROOT:-/ibex/user/meleg/DCM_MODEL}"
export EXPERIMENT_DIR="$(cat logs/latest_onset150_experiment.txt)"
source "$EXPERIMENT_DIR/ibex_environment.sh"
Rscript --vanilla "$EXPERIMENT_DIR/04_collect_results.R" "$EXPERIMENT_DIR"
''')
monitor = (OLD/'monitor_all_dcm_v6.sh').read_text(encoding='utf-8-sig')
monitor = monitor.replace('experiments/followup150_*/manifest.tsv',
                          'experiments/followup150_*/manifest.tsv experiments/onset150_*/manifest.tsv')
monitor = monitor.replace('logs/dcm-u150-*-c[0-9]*.out)',
                          'logs/dcm-u150-*-c[0-9]*.out logs/dcm-o150-*-c[0-9]*.out)')
put('monitor_all_dcm_v7.sh', monitor)
export = (OLD/'export_results.sh').read_text(encoding='utf-8-sig').replace('dcm_followup_comparison_','dcm_onset_comparison_')
put('export_results.sh', export)

# Cases.tsv is generated independently here and checked against campaign_cases() on Ibex.
rows=[]
def add(eq, subset='p10', split='primary', block='equation', reference=1,
        execution='fit', source_case='NA', seed=1234, delta=.95,
        warmup=1000, samples=2000, kd=7, mode='posterior'):
    i=len(rows)+1;threads=2 if execution=='reuse' or mode=='prior' else 4
    rows.append(dict(case_id=i,label=f'o{i:03}_{block}_{eq}_{split}_{subset}',block=block,
        equation=eq,subset=subset,split=split,prior='custom_z',mode=mode,variance='region',
        cast='none',family='gaussian',autocor='none',depth_max=150,chains=4,
        iter=warmup+samples,warmup=warmup,samples=samples,adapt_delta=delta,treedepth=14,
        threads=threads,cpus=4*threads,k_depth=kd,k_mld=5,k_distance=5,sigma_smooth_sd=.5,
        cores=4,refresh=50,seed=seed,memory_gb=64 if subset=='p20' else 32,
        reference_case=reference,execution=execution,source_case=source_case))
panels=[('p10','primary'),('p20','primary'),('p10','date_blocked')]
src=[1,5,13,15,22,24]
for j,(subset,split) in enumerate(panels):
    for eq in ['original_z','interaction_z']:
        i=len(rows)+1
        add(eq,subset,split,'existing_control','NA' if i%2 else i-1,'reuse',src[i-1])
equations=['mld_vary_z','mld_shift','mld_halfshift','mld_ratio','mld_hybrid','mld_shift_distance']
for j,(subset,split) in enumerate(panels):
    for eq in equations:add(eq,subset,split,reference=2*j+1)
for eq,ref in [('mld_shift',8),('mld_hybrid',11)]:
    add(eq,block='repeat_seed',reference=ref,seed=20261019)
    add(eq,block='robust_sampler',reference=ref,delta=.99,warmup=2000)
    add(eq,block='depth_basis',reference=ref,kd=10)
for eq in ['mld_shift','mld_hybrid']:
    for seed in [1234,20261019]:add(eq,block='prior_predictive',reference='NA',seed=seed,
                                  mode='prior',warmup=500,samples=500)
with (PKG/'cases.tsv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
    w.writeheader();w.writerows(rows)

put('test_offline.R', r'''
# Base-R gates run on Ibex BEFORE any sbatch submission; no fitting or networking.
a <- commandArgs(trailingOnly=TRUE);package <- a[1]
for(f in list.files(package,pattern="\\.R$",full.names=TRUE))parse(f)
Sys.setenv(ONSET_PACKAGE=package)
source(file.path(package,"campaign.R"));source(file.path(package,"submit.R"))
source(file.path(package,"onset_metrics.R"))
source(file.path(package,"04_collect_results.R"))
cases <- campaign_cases()
stopifnot(isTRUE(all.equal(cases,read.delim(file.path(package,"cases.tsv")),check.attributes=FALSE)),
  nrow(cases)==34L,sum(cases$execution=="fit")==28L,
  sum(cases$execution=="fit"&cases$mode=="posterior")==24L,
  sum(cases$execution=="reuse")==6L)
# Exact onset fixtures: known crossing, vertical displacement, absence, censoring,
# persistence, unsupported gaps, and a matching model with zero error.
z <- 10:150
curve <- onset_curve(z,pmax(0,(z-40)/20))
o <- onset_from_curve(curve,.5);stopifnot(o$status=="detected",o$depth==50)
shift <- onset_from_curve(onset_curve(z,pmax(0,(z-60)/20)),.5)
stopifnot(shift$depth-o$depth==20)
stopifnot(onset_from_curve(onset_curve(z,rep(0,length(z))),.5)$status=="not_reached_in_window")
stopifnot(onset_from_curve(onset_curve(z,rep(1,length(z))),.5)$status=="left_censored")
spike <- rep(0,length(z));spike[50:51] <- 2
stopifnot(onset_from_curve(onset_curve(z,spike),.5)$status=="not_reached_in_window")
gap <- onset_curve(c(10:30,50:150),c(rep(0,21),rep(1,101)))
stopifnot(onset_from_curve(gap,.5)$status=="gap_censored")
tab <- data.frame(cast_key="A",Region="R",Ocean="O",date="D",depth=z,
  MLD_temp_MIX=40,chl_norm=pmax(0,(z-40)/20),predicted=pmax(0,(z-40)/20))
m <- onset_profile_metrics(tab,rbind(tab$predicted,tab$predicted))
stopifnot(all(m$onset_abs_error==0),all(m$posterior_detected_draws==2),
          all(m$onset_interval_covers_observed_given_detected))
root <- tempfile("onset_scheduler_")
experiment <- file.path(root,"experiments","onset150_fixture")
dir.create(file.path(experiment,"scheduler"),recursive=TRUE)
write.table(cases,file.path(experiment,"cases.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
control <- data.frame(case_id=1:6,source_case=c(1,5,13,15,22,24),
  fit_job=as.character(101:106),check_job=as.character(201:206),source_run=file.path(root,paste0("old",1:6)))
write.table(control,file.path(experiment,"source_controls.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
fake <- new.env(parent=emptyenv());fake$n <- 1000L;fake$jobs <- list();fake$calls <- list()
lookup <- function(name)if(is.null(fake$jobs[[name]]))character()else fake$jobs[[name]]
submit <- function(args) {
  fake$n<-fake$n+1L;id<-as.character(fake$n)
  name<-sub("^--job-name=","",grep("^--job-name=",args,value=TRUE))
  fake$jobs[[name]]<-id;fake$calls[[length(fake$calls)+1L]]<-args;id
}
active <- function(known)ids(known)
invisible(capture.output(m <- queue_cases(experiment,root,cases,submit,lookup,active,function(id)TRUE)))
# 1 preparation + 28 new fits + 34 checks + 1 collector; six source IDs reused.
stopifnot(length(fake$calls)==64L,nrow(m)==34L,
  identical(as.character(m$fit_job[1:6]),as.character(101:106)))
fitcalls <- Filter(function(x)any(grepl("/fit.slurm$",x)),fake$calls)
reusecalls <- Filter(function(x)any(grepl("/reuse.slurm$",x)),fake$calls)
stopifnot(length(fitcalls)==28L,length(reusecalls)==6L,
  all(vapply(fitcalls,function(x)any(grepl("^--dependency=afterok:",x)),logical(1))),
  all(vapply(reusecalls,function(x)any(grepl("^--dependency=afterok:.*afterany:",x)),logical(1))))
for(j in seq_along(fitcalls)) {
  c <- cases[cases$execution=="fit",][j,]
  stopifnot(paste0("--cpus-per-task=",c$cpus)%in%fitcalls[[j]],
            paste0("--mem=",c$memory_gb,"G")%in%fitcalls[[j]])
}
writeLines("done",file.path(experiment,"PREPARATION_COMPLETE.txt"))
invisible(capture.output(again <- queue_cases(experiment,root,cases,submit,lookup,active,function(id)TRUE)))
same_registry <- function(a,b)identical(lapply(a,as.character),lapply(b,as.character))
stopifnot(length(fake$calls)==64L,same_registry(m,again))
unlink(file.path(experiment,"scheduler","fit_07.job"))
invisible(capture.output(recovered <- queue_cases(experiment,root,cases,submit,lookup,active,function(id)TRUE)))
stopifnot(length(fake$calls)==64L,same_registry(m,recovered))
uncertain <- try(ensure_job(experiment,"uncertain","uncertain_name",character(),
  function(args)stop("connection lost"),function(name)character()),silent=TRUE)
again_uncertain <- try(ensure_job(experiment,"uncertain","uncertain_name",character(),
  function(args)stop("must not submit"),function(name)character()),silent=TRUE)
stopifnot(inherits(uncertain,"try-error"),inherits(again_uncertain,"try-error"))
writeLines("fixture",file.path(experiment,"PREPARATION_FAILED.txt"))
blocked <- try(queue_cases(experiment,root,cases,submit,lookup,active,function(id)TRUE),silent=TRUE)
stopifnot(inherits(blocked,"try-error"),length(fake$calls)==64L)
unlink(file.path(experiment,"PREPARATION_FAILED.txt"))
# Keep incomplete and diagnostically poor fits in the comparison, with an onset
# contrast calculated from the same casts rather than a changing detectable subset.
for(i in c(1L,8L,10L)) {
  run <- file.path(experiment,"runs",sprintf("case_%02d",i))
  dir.create(file.path(run,"checks"),recursive=TRUE)
  fixture <- tab
  offset <- if(i==8L)0 else if(i==1L)5 else 10
  fixture$predicted <- pmax(0,(fixture$depth-40-offset)/20)
  pm <- onset_profile_metrics(fixture)
  write.csv(pm,file.path(run,"checks","validation_conditional_onset_cast_metrics.csv"),row.names=FALSE)
  summary <- data.frame(cases[i,],max_rhat=if(i==10L)1.2 else 1.001,
    min_bulk_ess=800,min_tail_ess=800,divergences=if(i==10L)2L else 0L,
    min_ebfmi=.7,validation_md5="same",reserve_md5="reserve")
  write.csv(summary,file.path(run,"test_summary.csv"),row.names=FALSE)
  writeLines("done",file.path(run,"FIT_COMPLETE.txt"))
  writeLines("done",file.path(run,"CHECKS_COMPLETE.txt"))
}
invisible(capture.output(all <- collect_results(experiment,root)))
invisible(capture.output(collect_onset(experiment)))
status <- read.csv(file.path(experiment,"results","onset_case_status.csv"))
scores <- read.csv(file.path(experiment,"results","onset_all_scores.csv"))
pairs <- read.csv(file.path(experiment,"results","onset_planned_paired_contrasts.csv"))
stopifnot(nrow(all)==34L,nrow(status)==34L,nrow(scores)==9L,
  any(all$divergences==2L,na.rm=TRUE),
  all(pairs$delta_onset_mae_common_casts[pairs$case_id==8L]== -5),
  all(pairs$jointly_detected_casts==1L))
cat("PASS: parse/case matrix, onset shifts/absence/censoring/gaps/persistence, posterior intervals;\n",
    "28 new fits, six source fits reused, dependencies/resources, repeat/orphan recovery, uncertain submission stop;\n",
    "all statuses, poor diagnostics retained and matched onset contrasts.\n")
''')

put('PLAN.md', '''
# MLD-onset tests — 19 October 2026

Question: does mixed-layer depth locate the beginning of subsurface chlorophyll
enhancement, rather than merely improving whole-profile error or peak depth?

This package registers 34 cases: **24 new posterior fits, 4 prior-predictive fits,
and 6 evaluations of existing fits**. No old model is refitted or cancelled.

## Equations

All retain Region, nested Ocean:Region:date and Ocean:Region:enso intercepts,
Gaussian observations and sigma ~ 0 + Region. Custom standardized priors are
copied from the existing follow-up. Distance is standardized log10_distance.

1. Absolute-depth varying MLD: existing s(depth,by=mld_z,k=5) alongside Region
   depth curves and a linear distance effect. This is a matched simpler MLD control.
2. MLD shift: replace depth in Region curves by depth - MLD, k=7. This tests
   translation with MLD; the enhancement can start above or below MLD.
3. Half shift: Region curves use depth - 0.5*MLD, k=7. This tests a weaker
   relationship instead of assuming a one-metre translation per metre of MLD.
4. MLD ratio: Region curves use depth / MLD, k=7. This tests stretching with MLD.
5. Hybrid: retain absolute-depth Region curves k=7 and add a shared smooth of
   depth - MLD, k=5. This allows regional absolute-depth structure and a separate
   MLD-relative contribution to coexist; check posterior dependence/diagnostics.
6. Shift plus distance: shifted Region curves k=7 plus a shared coefficient
   s(depth - MLD,by=distance_z,k=5), allowing distance to change shifted shape.

Every equation is run on the exact existing primary 108-cast training set, primary
216-cast training set, and 108-cast date-blocked training set: 18 posterior fits.
Six additional fits isolate repeat seed, robust sampler settings, and k=10 for
the shift and hybrid models. Four prior-only fits check those two new forms.
Controls reuse follow-up cases 1/5,13/15,22/24 (no-MLD and product interaction),
with raw training identity checks before scoring. Their original sampler settings
are preserved. New fits use four threads per chain; old controls use two, so
equation comparisons are statistical comparisons, not controlled speed benchmarks.

## Sampler settings

Default: four chains, 1000 warmup + 2000 retained draws per chain (iter=3000),
adapt_delta=.95, max_treedepth=14, four threads per chain (16 CPUs). Robust checks
use warmup=2000 and adapt_delta=.99, keeping 2000 retained draws. Prior-only fits
use 500 warmup + 500 retained and two threads (8 CPUs). Memory is 32 GB for 10%
fits and 64 GB for 20%. Wall limit 7 days is a scheduler limit, not a runtime estimate.
All 28 new fitting jobs together request 416 CPUs and 1088 GB if simultaneously
allocated; Slurm decides actual concurrency. Evaluation/preparation jobs are extra.

## Data and honest comparisons

No new split selection, response smoothing, binning, cast deletion or MLD imputation.
At least 10% (108 casts) is used to train every posterior model. Validation is
another >=10%, plus the existing disjoint 108-cast reserve. Primary validation may
share dates; date-blocked validation excludes training Ocean/Region/date groups.
Copied raw rows must be identical. New deterministic predictor columns change RDS
hashes, so a separate source-identity audit is retained. Reserve is now a development
cohort because we have already inspected it; it is not an untouched final test.
Nonpositive MLD stops preparation rather than dropping rows or silently changing
the ratio formula. Every actual formula, requested prior and design matrix is
constructed on Ibex before fits can start; rank deficiency fails preparation.

## New onset diagnostics

The response is (Chl-surface_chl)/surface_chl. Fixed thresholds .25,.5,1 therefore
mean 25%,50%,100% enhancement relative to the original 5-15 m surface baseline.
The primary onset threshold is .5; .25 and 1 are prespecified sensitivities.
Onset is the first depth at/after 15 m above threshold for five consecutive
one-metre grid points. For scoring only, both observed and expected profiles use
a five-metre running median. Fit inputs remain unsmoothed. Grid interpolation
stays inside observed 10-150 m support and never bridges gaps >3 m.

Threshold already reached at the upper boundary is left-censored. No threshold
crossing is labelled not reached in the observed window, never replaced by 150 m.
Gaps before first detection are flagged. Export detection agreement, false
positives/negatives, censor counts, conditional onset MAE, onset minus MLD,
and posterior onset probabilities/intervals from up to 400 expected-curve draws.
Paired onset-MAE contrasts use the same jointly detectable casts. Always inspect
false positives/negatives and the sample count alongside conditional MAE; a model
must not appear better just because it fails to predict difficult onsets.

Also export steepest-rise depth (positive centred four-metre slope >=.02 chl_norm
per metre, profile amplitude >=.5), error within +/-20 m of MLD, full-profile RMSE,
shape error, peak depth, interval calibration, regional metrics and MCMC diagnostics.
The steepest-rise location is a complementary operational metric, not a biological
definition of DCM onset. Do not force onset or peak depth to equal MLD. No physical
truncation or removal of large Red Sea errors is applied.

## Execution

Upload queue_onset_dcm.sh and run it from /ibex/user/meleg/DCM_MODEL. Repeating the
same installer resumes recorded submissions. ONSET_NEW=1 explicitly creates a new
campaign; leave it unset for recovery. Accepted old or new fit jobs are not replaced.
An ambiguous submission stops rather than risking duplicate jobs. Unique job names
permit recovery after interruption. All helper checks run before queueing; preparation
is an afterok dependency of new fits and reused evaluations. Source evaluations also
wait for active source fits/checks. The final collector includes all recorded earlier
campaigns and retains diagnostic failures and incomplete states.

Refresh: bash DCM_Ibex_onset_tests/collect_results.sh
Monitor all: bash DCM_Ibex_onset_tests/monitor_all_dcm_v7.sh --diagnostics
Export all: bash DCM_Ibex_onset_tests/export_results.sh

Scientific motivation: Mignot et al. 2014, Global Biogeochemical Cycles,
doi:10.1002/2013GB004781, describes mixed-layer and DCM structure while highlighting
light-driven DCM movements. It motivates testing the relationship, not fixing it.
Smooth syntax: https://paulbuerkner.com/brms/reference/s.html
''')
put('README.md', '''
# DCM MLD-onset campaign

Use the single-file installer queue_onset_dcm.sh. It installs into the new
DCM_Ibex_onset_tests directory and queues 28 new fits plus six existing-fit
evaluations. Prior campaigns and pointer files remain untouched. Read PLAN.md
for the scientific design and VALIDATION.md for checks and limits.

Run from the Ibex project directory:

    bash queue_onset_dcm.sh

To install for inspection without queueing:

    bash queue_onset_dcm.sh --extract-only

Repeating submission resumes the frozen campaign. No automatic refits or job
cancellations occur. Check results/onset_all_scores.csv,
results/onset_planned_paired_contrasts.csv and results/mld_neighbourhood_scores.csv
alongside all-run predictive/diagnostic tables. All completed or failed cases remain
visible. Lower conditional onset error alone is insufficient for model selection.
''')

def delimiters(path):
    text=path.read_text();stack=[];quote=None;escape=False;comment=False
    for ch in text:
        if comment:
            if ch=='\n':comment=False
            continue
        if quote:
            if escape:escape=False
            elif ch=='\\':escape=True
            elif ch==quote:quote=None
            continue
        if ch=='#':comment=True
        elif ch in ('"',"'",'`'):quote=ch
        elif ch in '({[':stack.append(ch)
        elif ch in ')}]':assert stack and stack.pop()=={')':'(',']':'[','}':'{'}[ch],path
    assert not quote and not stack,path

assert len(rows)==34 and sum(r['execution']=='fit' for r in rows)==28
assert sum(r['execution']=='fit' and r['mode']=='posterior' for r in rows)==24
assert len({r['label'] for r in rows})==34
assert all(r['iter']==r['warmup']+r['samples'] and r['cpus']==r['chains']*r['threads'] for r in rows)
assert sum(r['cpus'] for r in rows if r['execution']=='fit')==416
assert sum(r['memory_gb'] for r in rows if r['execution']=='fit')==1088
for r in rows:
    if isinstance(r['reference_case'],int):assert 1<=r['reference_case']<r['case_id']
# Check reuse identities against the source case matrix before shipping.
with (OLD/'cases.tsv').open(encoding='utf-8-sig') as f:source=list(csv.DictReader(f,delimiter='\t'))
check_fields=['equation','subset','split','prior','mode','variance','cast','family','autocor',
              'depth_max','chains','iter','warmup','adapt_delta','treedepth','threads',
              'k_depth','k_mld','k_distance','seed']
for r in rows[:6]:
    old=source[r['source_case']-1]
    for key in check_fields:
        try:equal=float(r[key])==float(old[key])
        except ValueError:equal=str(r[key])==old[key]
        assert equal,(r['case_id'],key,r[key],old[key])
for p in PKG.glob('*.R'):delimiters(p)
assert 'brms::brm(' not in (PKG/'reuse_check.R').read_text()
put('VALIDATION.md', '''
# Validation status

Passed locally without starting R/Stan: independent 34-row design generation;
24 new posterior/4 prior fits and six reused source fits; exact source-control
parameter matching; unique labels/IDs; CPU/memory/iteration arithmetic; reference
ordering; balanced R strings/delimiters; no fitting call in source reuse code;
Unix newlines; installer payload paths and SHA256 equality with source files.

R parser/runtime, model construction, Slurm and Bash are not claimed as locally
verified. Local R was previously blocked by automatic approval review because of
an attempted external telemetry connection; it is not retried here. This package
runs its base-R tests and Bash syntax checks on Ibex before submitting anything.
The tests cover onset shifts, absence, censoring, persistence, gaps, posterior
intervals, exact new-job counts, six reused fits, dependency/resources, repeats,
orphan recovery, uncertain-submission handling and failed-preparation blocking.

An Ibex preparation job then verifies the actual raw training/validation/reserve
rows, >=10% training, source-control settings, positive MLD, all requested brms
formulas/priors/Stan code and fixed-design rank before any model starts.
Failure stops fitting through afterok dependencies. New spline combinations
and posterior sampling therefore remain subject to real cluster preflight and
diagnostics; successful static checks do not guarantee model convergence.
''')

def bundle():
    files=sorted(p for p in PKG.iterdir() if p.is_file() and p.name!='SHA256SUMS.txt')
    put('SHA256SUMS.txt',''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.name+'\n' for p in files))
    files.append(PKG/'SHA256SUMS.txt')
    buffer=io.BytesIO()
    with tarfile.open(fileobj=buffer,mode='w:gz',format=tarfile.USTAR_FORMAT) as t:
        for p in files:
            data=p.read_bytes();info=tarfile.TarInfo(PKG.name+'/'+p.name)
            info.size=len(data);info.mode=0o644;info.mtime=1792378800
            t.addfile(info,io.BytesIO(data))
    payload=buffer.getvalue();digest=hashlib.sha256(payload).hexdigest()
    header='''#!/bin/bash
# Install and queue 24 posterior fits, 4 priors and 6 existing-fit evaluations.
# Repeating resumes registrations; --extract-only installs without queueing.
set -euo pipefail
bash -n "$0"
mode="${1:-submit}"
[[ "$mode" == submit || "$mode" == --extract-only ]] || { echo 'Usage: bash queue_onset_dcm.sh [--extract-only]' >&2; exit 2; }
root="${DCM_ROOT:-/ibex/user/meleg/DCM_MODEL}"
[[ -d "$root" ]] || { echo "Missing project directory: $root" >&2; exit 1; }
cd -- "$root"
root="$PWD"
package="$root/DCM_Ibex_onset_tests"
mkdir -p logs
for command in base64 tar sha256sum; do command -v "$command" >/dev/null || exit 1; done
if [[ "$mode" == submit ]]; then command -v sbatch >/dev/null || { echo 'Run this on Ibex.' >&2; exit 1; }; fi
payload=$(mktemp "$root/logs/onset_payload.XXXXXXXX.tar.gz")
trap 'rm -f -- "$payload"' EXIT
base64 --decode > "$payload" <<'__DCM_ONSET_PAYLOAD__'
'''
    footer='''
__DCM_ONSET_PAYLOAD__
printf '%s  %s\\n' '__DIGEST__' "$payload" | sha256sum --check --status
if [[ -L "$package" ]]; then echo 'Package path is a symlink; refusing extraction.' >&2; exit 1; fi
if [[ ! -e "$package" ]]; then
  staging=$(mktemp -d "$root/.dcm_onset_install.XXXXXXXX")
  tar --extract --gzip --file "$payload" --directory "$staging" --no-same-owner
  mv -- "$staging/DCM_Ibex_onset_tests" "$package"
  rmdir -- "$staging"
fi
[[ -f "$package/SHA256SUMS.txt" ]] || { echo 'Existing package incomplete; no overwrite.' >&2; exit 1; }
expected=$(tar -xOzf "$payload" DCM_Ibex_onset_tests/SHA256SUMS.txt | sha256sum | cut -d' ' -f1)
actual=$(sha256sum "$package/SHA256SUMS.txt" | cut -d' ' -f1)
[[ "$actual" == "$expected" ]] || { echo 'Existing package version differs; no overwrite or submissions.' >&2; exit 1; }
(cd -- "$package"; sha256sum --check --quiet SHA256SUMS.txt; for file in *.sh *.slurm; do bash -n "$file"; done)
printf '\\nInstalled: %s\\nPlan: %s/PLAN.md\\n' "$package" "$package"
if [[ "$mode" == --extract-only ]]; then echo 'No jobs submitted.'; exit 0; fi
bash "$package/submit_tests.sh"
'''.replace('__DIGEST__',digest)
    text=header+'\n'.join(textwrap.wrap(base64.b64encode(payload).decode(),76))+footer
    out=ROOT/'queue_onset_dcm.sh';out.write_text(text,encoding='ascii',newline='\n')
    start=text.index("<<'__DCM_ONSET_PAYLOAD__'\n")+len("<<'__DCM_ONSET_PAYLOAD__'\n")
    end=text.index('\n__DCM_ONSET_PAYLOAD__\n',start)
    saved=base64.b64decode(text[start:end]);assert hashlib.sha256(saved).hexdigest()==digest
    with tarfile.open(fileobj=io.BytesIO(saved),mode='r:gz') as t:
        for m in t.getmembers():
            path=PurePosixPath(m.name)
            assert not path.is_absolute() and '..' not in path.parts and path.parts[0]==PKG.name and m.isfile()
            assert t.extractfile(m).read()==(ROOT/m.name).read_bytes()
    print(json.dumps(dict(installer=str(out),bytes=out.stat().st_size,
          new_posterior_fits=24,prior_fits=4,reused_controls=6,payload_sha256=digest)))

if __name__=='__main__':bundle()
