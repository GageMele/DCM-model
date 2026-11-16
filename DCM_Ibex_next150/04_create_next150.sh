#!/bin/bash
# Create and launch the next150 campaign: 14 new literature-backed fits plus
# 7 re-scorings of existing shape150 anchor fits (adds held-out log predictive
# density). Reuses shape150's frozen data files (symlinked, read-only) and its
# scheduler/check code. Does not touch the running shape150 campaign.
set -euo pipefail
root=/ibex/user/meleg/DCM_MODEL
cd "$root"
shape=$(cat logs/latest_shape150_experiment.txt)
stamp=$(date -u +%Y%m%dT%H%M%SZ)
exp="$root/experiments/next150_${stamp}_$$"
echo "shape150 source: $shape"
echo "new experiment:  $exp"
mkdir -p "$exp/logs" "$exp/scheduler" "$exp/runs"
for f in common.R data.R 03_check_case.R 04_collect_results.R onset_metrics.R profile_helpers.R \
         scheduler.py worker.slurm ibex_environment.sh scaling.rds region_levels.rds; do
  cp "$shape/$f" "$exp/$f"
done
ln -s "$shape/data" "$exp/data"
echo "$shape" > "$exp/reuse_source.txt"
cat > "$exp/models.R" <<'EOF_NEXT150'
# Next150 equations: literature-backed follow-ups to the shape150 campaign.
# Inherits the shape150 feature construction; no fitting occurs when sourced.
run_directory <- function(experiment,case)file.path(experiment,"runs",sprintf("case_%03d",case$case_id))
data_file <- function(experiment,case)file.path(experiment,"data",paste0(case$dataset_id,".rds"))
thread_specification <- function(case)brms::threading(as.integer(case$threads))
group_key <- function(d)paste(d$Ocean,d$Region,d$date,sep="\034")

# Ratio of chlorophyll to the cast's own 5-15 m surface median (chl_norm + 1).
# The <0.2% of rows with Chl <= 0 are floored at 1% of surface for the lognormal
# likelihood only; scores are always computed against the unmodified chl_norm.
ratio_floor <- 0.01

apply_features <- function(d,scaling,region_levels) {
  d$Region<-factor(as.character(d$Region),levels=region_levels)
  d$Ocean<-factor(as.character(d$Ocean));d$date<-as.character(d$date);d$enso<-as.character(d$enso)
  d$mld_z<-(d$MLD_temp_MIX-scaling$mld_mean)/scaling$mld_sd
  d$distance_z<-(d$log10_distance-scaling$distance_mean)/scaling$distance_sd
  d$mld_distance_z<-d$mld_z*d$distance_z
  d$depth_mld<-d$depth-d$MLD_temp_MIX
  d$depth_mld_half<-d$depth-.5*d$MLD_temp_MIX
  for(i in seq_along(region_levels)) {
    yes<-as.character(d$Region)==region_levels[i]
    d[[paste0("region_only_",i)]]<-ordered(ifelse(yes,"on","off"),levels=c("off","on"))
    d[[paste0("region_mld_",i)]]<-as.numeric(yes)*d$mld_z
  }
  d<-augment(d)
  if(anyNA(d$Region)||any(!is.finite(d$mld_z)))stop("Invalid transformed predictors")
  d
}

# Columns added on top of the frozen shape150 data files (which already carry
# the shape150 features); applied to training, holdout and plotting data alike.
augment <- function(d) {
  d$depth_hm<-(d$depth-.5*d$MLD_temp_MIX)/100
  if("chl_norm" %in% names(d))d$ratio<-pmax(d$chl_norm+1,ratio_floor)
  d
}

groups_terms <- function(case)
  if(case$groups=="date_only")"(1|Ocean:Region:date)" else c("(1|Ocean:Region:date)","(1|Ocean:Region:enso)")

sigma_formula <- function(case)
  if(case$variance=="region_depth")
    sigma~0+Region+s(depth,by=Ocean,bs='cr',k=5) else sigma~0+Region

campaign_spec <- function(case,region_levels) {
  kd<-as.integer(case$k_depth)
  family<-if(case$family=="lognormal")brms::lognormal() else gaussian()
  response<-if(case$family=="lognormal")"ratio" else "chl_norm"
  u<-"depth_mld_half"
  regional<-function(v,k=kd,m=2)sprintf("s(%s,by=Region,bs='cr',k=%d%s)",v,k,if(m==1)",m=1" else "")
  pooled_mld<-vapply(seq_along(region_levels),function(i)
    sprintf("s(%s,by=region_mld_%d,bs='cr',k=5,m=1)",u,i),character(1))
  if(case$equation=="peak_nl") return(peak_spec(case,family))
  terms<-switch(case$equation,
    # Anchor equations re-expressed (Gaussian anchors live in shape150).
    half=c("distance_z",regional(u)),
    # Pedersen et al. 2019 'GI': shared Ocean-level curve + penalised regional deviations.
    half_GI=c("distance_z",sprintf("s(%s,by=Ocean,bs='cr',k=%d)",u,kd),regional(u,k=kd,m=1)),
    # Shared depth-varying MLD slope + regional deviations shrunk by m=1 penalty.
    half_mld_pooled=c("distance_z",regional(u),sprintf("s(%s,by=mld_z,bs='cr',k=%d)",u,kd),pooled_mld),
    # Island-mass effect: distance-from-shore effect allowed to vary with depth, by Ocean.
    half_dist_vary=c(regional(u),sprintf("s(%s,by=distance_z,bs='cr',k=5)",u)),
    stop("Unknown equation: ",case$equation))
  f<-as.formula(paste(response,"~",paste(c("1","Region",terms,groups_terms(case)),collapse=" + ")),
    env=asNamespace("brms"))
  form<-brms::bf(f,sigma_formula(case),family=family)
  scale_prior<-"student_t(3,0,1.5)"
  intercept<-if(case$family=="lognormal")"student_t(3,0,1)" else "student_t(3,0.73,5.448)"
  priors<-c(brms::set_prior(intercept,class="Intercept"),
    brms::set_prior(if(case$prior=="b_tight")"normal(0,1)" else "normal(0,2)",class="b"),
    brms::set_prior(scale_prior,class="sd"),brms::set_prior(scale_prior,class="sds"),
    brms::set_prior("normal(0,1)",class="b",dpar="sigma"))
  if(case$variance=="region_depth")priors<-c(priors,brms::set_prior("normal(0,0.5)",class="sds",dpar="sigma"))
  list(formula=form,priors=priors)
}

# Parametric Gaussian-peak-on-background profile (Lewis et al. 1983; Platt et al.
# 1988; Uitz et al. 2006; Mignot et al. 2011), on the log-ratio scale:
#   median ratio(z) = exp(lbg + slope*z) + exp(lamp) * exp(-0.5*((z - zm)/exp(lwid))^2)
# with z = (depth - 0.5*MLD)/100. MLD enters peak depth (shift), amplitude (erosion)
# and width separately, so 'MLD moves the DCM' and 'MLD weakens it' are distinguishable.
peak_spec <- function(case,family) {
  form<-brms::bf(ratio~log(exp(lbg+slope*depth_hm)+exp(lamp)*exp(-0.5*((depth_hm-zm)/exp(lwid))^2)),
    brms::lf(lbg~0+Region),brms::lf(slope~0+Region),
    brms::lf(lamp~0+Region+mld_z+distance_z+(1|g|Ocean:Region:date)),
    brms::lf(zm~0+Region+mld_z+(1|g|Ocean:Region:date)),
    brms::lf(lwid~0+Region+mld_z),
    sigma_formula(case),nl=TRUE,family=family)
  priors<-c(brms::set_prior("normal(0,0.5)",nlpar="lbg"),
    brms::set_prior("normal(-1,1)",nlpar="slope"),
    brms::set_prior("normal(0,1)",nlpar="lamp"),
    brms::set_prior("normal(0.6,0.3)",nlpar="zm"),
    brms::set_prior("normal(-1.4,0.5)",nlpar="lwid"),
    brms::set_prior("normal(0,0.5)",class="sd",nlpar="lamp"),
    brms::set_prior("normal(0,0.3)",class="sd",nlpar="zm"),
    brms::set_prior("normal(0,1)",class="b",dpar="sigma"))
  if(case$variance=="region_depth")priors<-c(priors,brms::set_prior("normal(0,0.5)",class="sds",dpar="sigma"))
  list(formula=form,priors=priors)
}
EOF_NEXT150
cat > "$exp/02_fit_case.R" <<'EOF_NEXT150'
# Fit one frozen case. An accepted fitting job is never automatically repeated.
experiment<-Sys.getenv("EXPERIMENT_DIR")
source(file.path(experiment,"common.R"));source(file.path(experiment,"models.R"))
if(!file.exists(file.path(experiment,"PREPARATION_COMPLETE.txt")))stop("Preparation is incomplete")
case<-read_case(experiment,Sys.getenv("CASE_ID"));run<-run_directory(experiment,case)
if(dir.exists(run))stop("Run directory already exists; inspect it rather than overwriting a fit")
dir.create(run,recursive=TRUE)
if(case$mode=="reuse") {
  # Re-score an existing shape150 fit with this campaign's scoring (adds held-out
  # log predictive density). Nothing is refitted; the source run is only linked.
  tryCatch({
    source_exp<-readLines(file.path(experiment,"reuse_source.txt"))[1]
    src<-file.path(source_exp,"runs",case$source_run)
    if(!file.exists(file.path(src,"FIT_COMPLETE.txt")))stop("Source fit not complete: ",src)
    for(x in c("fit.rds","stan_csv","model_data.rds"))
      if(!file.symlink(file.path(src,x),file.path(run,x)))stop("Cannot link ",x)
    for(x in c("run_config.rds","run_config.txt","hardware.txt","session_info.txt","requested_priors.txt","resolved_priors.csv"))
      if(file.exists(file.path(src,x)))file.copy(file.path(src,x),file.path(run,x))
    writeLines(c(paste("Reused source fit",src),"Fit serialized with draws; diagnostics and predictions follow separately"),
      file.path(run,"FIT_COMPLETE.txt"))
  },error=function(e){writeLines(conditionMessage(e),file.path(run,"FIT_FAILED.txt"));stop(e)})
  quit(status=0L)
}
tryCatch({
  start<-proc.time()[["elapsed"]]
  for(p in c("brms","cmdstanr","mgcv","posterior"))if(!requireNamespace(p,quietly=TRUE))stop("Missing ",p)
  if(nzchar(Sys.getenv("CMDSTAN")))cmdstanr::set_cmdstan_path(Sys.getenv("CMDSTAN"))
  data_path<-data_file(experiment,case);d<-augment(readRDS(data_path))
  spec<-campaign_spec(case,readRDS(file.path(experiment,"region_levels.rds")))
  allocated<-as.integer(Sys.getenv("SLURM_CPUS_PER_TASK",case$cpus))
  if(allocated<case$chains*case$threads)stop("Insufficient allocated CPUs")
  Sys.setenv(STAN_NUM_THREADS=case$threads)
  config<-list(case=case,job_id=Sys.getenv("SLURM_JOB_ID"),allocated_cpus=allocated,
    hostname=Sys.info()[["nodename"]],data_md5=unname(tools::md5sum(data_path)),
    source_md5=unname(tools::md5sum(file.path(experiment,"models.R"))),
    formula=paste(capture.output(spec$formula),collapse=" "),cmdstan_version=cmdstanr::cmdstan_version())
  saveRDS(config,file.path(run,"run_config.rds"));capture.output(config,file=file.path(run,"run_config.txt"))
  saveRDS(list(model_data=d),file.path(run,"model_data.rds"))
  saveRDS(spec$priors,file.path(run,"requested_priors.rds"))
  capture.output(spec$priors,file=file.path(run,"requested_priors.txt"))
  writeLines(system2("lscpu",stdout=TRUE),file.path(run,"hardware.txt"))
  capture.output(sessionInfo(),file=file.path(run,"session_info.txt"))
  csv_dir<-file.path(run,"stan_csv");dir.create(csv_dir)
  writeLines(brms::make_stancode(spec$formula,data=d,prior=spec$priors,
    threads=thread_specification(case),sample_prior=if(case$mode=="prior")"only"else"no"),file.path(run,"model.stan"))
  cat("Starting",case$label,"on",length(unique(d$cast_key)),"casts /",nrow(d),"rows\n")
  sampling_start<-proc.time()[["elapsed"]]
  fit<-brms::brm(spec$formula,data=d,prior=spec$priors,backend="cmdstanr",
    chains=case$chains,cores=case$cores,threads=thread_specification(case),
    iter=case$iter,warmup=case$warmup,seed=case$seed,refresh=case$refresh,
    save_warmup=TRUE,thin=1L,init=if(case$equation=="peak_nl")0 else NULL,
    control=list(adapt_delta=case$adapt_delta,max_treedepth=case$treedepth),
    save_pars=brms::save_pars(all=TRUE),sample_prior=if(case$mode=="prior")"only"else"no",
    output_dir=csv_dir,file=file.path(run,"fit"),file_refit="always")
  config$brm_seconds_including_compile<-proc.time()[["elapsed"]]-sampling_start
  invisible(posterior::as_draws_array(fit));saveRDS(fit,file.path(run,"fit.rds"))
  config$fit_seconds<-proc.time()[["elapsed"]]-start
  saveRDS(config,file.path(run,"run_config.rds"))
  write.csv(as.data.frame(brms::prior_summary(fit)),file.path(run,"resolved_priors.csv"),row.names=FALSE)
  writeLines("Fit serialized with draws; diagnostics and predictions follow separately",file.path(run,"FIT_COMPLETE.txt"))
},error=function(e){writeLines(conditionMessage(e),file.path(run,"FIT_FAILED.txt"));stop(e)})
EOF_NEXT150
cat > "$exp/scoring.R" <<'EOF_NEXT150'
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
      re<-if(mode=="conditional")NULL else NA
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
EOF_NEXT150
cat > "$exp/cases.tsv" <<'EOF_NEXT150'
case_id	label	block	equation	subset	fraction	split	dataset_id	selection	selection_seed	prior	mode	variance	family	cast	autocor	depth_max	k_depth	k_mld	k_distance	k_profile	chains	cores	threads	cpus	memory_gb	iter	warmup	samples	adapt_delta	treedepth	seed	refresh	evaluation_role	fit_days	priority	equation_reference_case	settings_reference_case	groups	source_run
1	n001_lognormal_likelihood_half_primary_p20	lognormal_likelihood	half	p20	0.2	primary	primary_p20_frozen	frozen	0	custom_z	posterior	region	lognormal	none	none	150	7	4	5	uniform	4	4	4	16	96	3000	1000	2000	0.95	14	1234	100	heldout	7	10	1	1	date_enso	NA
2	n002_lognormal_likelihood_half_date_blocked_p20	lognormal_likelihood	half	p20	0.2	date_blocked	date_blocked_p20_repeat_20261022	repeat	20261022	custom_z	posterior	region	lognormal	none	none	150	7	4	5	uniform	4	4	4	16	96	3000	1000	2000	0.95	14	1234	100	heldout	7	10	2	2	date_enso	NA
3	n003_lognormal_sigma_depth_half_primary_p20	lognormal_sigma_depth	half	p20	0.2	primary	primary_p20_frozen	frozen	0	custom_z	posterior	region_depth	lognormal	none	none	150	7	4	5	uniform	4	4	4	16	96	3000	1000	2000	0.95	14	1234	100	heldout	7	10	3	3	date_enso	NA
4	n004_lognormal_sigma_depth_half_date_blocked_p20	lognormal_sigma_depth	half	p20	0.2	date_blocked	date_blocked_p20_repeat_20261022	repeat	20261022	custom_z	posterior	region_depth	lognormal	none	none	150	7	4	5	uniform	4	4	4	16	96	3000	1000	2000	0.95	14	1234	100	heldout	7	10	4	4	date_enso	NA
5	n005_date_only_groups_half_primary_p20	date_only_groups	half	p20	0.2	primary	primary_p20_frozen	frozen	0	custom_z	posterior	region	gaussian	none	none	150	7	4	5	uniform	4	4	4	16	96	3000	1000	2000	0.95	14	1234	100	heldout	7	10	5	5	date_only	NA
6	n006_date_only_groups_half_date_blocked_p20	date_only_groups	half	p20	0.2	date_blocked	date_blocked_p20_repeat_20261022	repeat	20261022	custom_z	posterior	region	gaussian	none	none	150	7	4	5	uniform	4	4	4	16	96	3000	1000	2000	0.95	14	1234	100	heldout	7	10	6	6	date_only	NA
7	n007_pooled_shape_GI_half_GI_primary_p20	pooled_shape_GI	half_GI	p20	0.2	primary	primary_p20_frozen	frozen	0	custom_z	posterior	region_depth	lognormal	none	none	150	7	4	5	uniform	4	4	4	16	96	3000	1000	2000	0.95	14	1234	100	heldout	7	10	7	7	date_enso	NA
8	n008_pooled_shape_GI_half_GI_date_blocked_p20	pooled_shape_GI	half_GI	p20	0.2	date_blocked	date_blocked_p20_repeat_20261022	repeat	20261022	custom_z	posterior	region_depth	lognormal	none	none	150	7	4	5	uniform	4	4	4	16	96	3000	1000	2000	0.95	14	1234	100	heldout	7	10	8	8	date_enso	NA
9	n009_pooled_regional_mld_half_mld_pooled_primary_p20	pooled_regional_mld	half_mld_pooled	p20	0.2	primary	primary_p20_frozen	frozen	0	custom_z	posterior	region_depth	lognormal	none	none	150	7	4	5	uniform	4	4	4	16	96	3000	1000	2000	0.95	14	1234	100	heldout	7	10	9	9	date_enso	NA
10	n010_pooled_regional_mld_half_mld_pooled_date_blocked_p20	pooled_regional_mld	half_mld_pooled	p20	0.2	date_blocked	date_blocked_p20_repeat_20261022	repeat	20261022	custom_z	posterior	region_depth	lognormal	none	none	150	7	4	5	uniform	4	4	4	16	96	3000	1000	2000	0.95	14	1234	100	heldout	7	10	10	10	date_enso	NA
11	n011_distance_depth_vary_half_dist_vary_primary_p20	distance_depth_vary	half_dist_vary	p20	0.2	primary	primary_p20_frozen	frozen	0	custom_z	posterior	region_depth	lognormal	none	none	150	7	4	5	uniform	4	4	4	16	96	3000	1000	2000	0.95	14	1234	100	heldout	7	10	11	11	date_enso	NA
12	n012_distance_depth_vary_half_dist_vary_date_blocked_p20	distance_depth_vary	half_dist_vary	p20	0.2	date_blocked	date_blocked_p20_repeat_20261022	repeat	20261022	custom_z	posterior	region_depth	lognormal	none	none	150	7	4	5	uniform	4	4	4	16	96	3000	1000	2000	0.95	14	1234	100	heldout	7	10	12	12	date_enso	NA
13	n013_peak_profile_peak_nl_primary_p20	peak_profile	peak_nl	p20	0.2	primary	primary_p20_frozen	frozen	0	custom_z	posterior	region	lognormal	none	none	150	7	4	5	uniform	4	4	4	16	96	3000	1000	2000	0.95	14	1234	100	heldout	7	10	13	13	date_only	NA
14	n014_peak_profile_peak_nl_date_blocked_p20	peak_profile	peak_nl	p20	0.2	date_blocked	date_blocked_p20_repeat_20261022	repeat	20261022	custom_z	posterior	region	lognormal	none	none	150	7	4	5	uniform	4	4	4	16	96	3000	1000	2000	0.95	14	1234	100	heldout	7	10	14	14	date_only	NA
15	n015_rescore_shape150_case019_original_z_primary	rescore_anchor	original_z	p20	0.2	primary	primary_p20_frozen	frozen	0	custom_z	reuse	region	gaussian	none	none	150	7	4	5	uniform	4	1	4	1	8	3000	1000	2000	0.97	14	1234	100	heldout	1	0	15	15	date_enso	case_019
16	n016_rescore_shape150_case025_mld_halfshift_primary	rescore_anchor	mld_halfshift	p20	0.2	primary	primary_p20_frozen	frozen	0	custom_z	reuse	region	gaussian	none	none	150	7	4	5	uniform	4	1	4	1	8	3000	1000	2000	0.97	14	1234	100	heldout	1	0	16	16	date_enso	case_025
17	n017_rescore_shape150_case034_mld_vary_region_primary	rescore_anchor	mld_vary_region	p20	0.2	primary	primary_p20_frozen	frozen	0	custom_z	reuse	region	gaussian	none	none	150	7	4	5	uniform	4	1	4	1	8	3000	1000	2000	0.97	14	1234	100	heldout	1	0	17	17	date_enso	case_034
18	n018_rescore_shape150_case129_mld_vary_region_primary	rescore_anchor	mld_vary_region	p20	0.2	primary	primary_p20_frozen	frozen	0	custom_z	reuse	region	gaussian	none	none	150	7	4	5	rs15_other7	4	1	4	1	8	3000	1000	2000	0.97	14	1234	100	heldout	1	0	18	18	date_enso	case_129
19	n019_rescore_shape150_case265_original_z_date_blocked	rescore_anchor	original_z	p20	0.2	date_blocked	date_blocked_p20_repeat_20261022	repeat	20261022	custom_z	reuse	region	gaussian	none	none	150	7	4	5	uniform	4	1	4	1	8	3000	1000	2000	0.97	14	1234	100	heldout	1	0	19	19	date_enso	case_265
20	n020_rescore_shape150_case267_mld_halfshift_date_blocked	rescore_anchor	mld_halfshift	p20	0.2	date_blocked	date_blocked_p20_repeat_20261022	repeat	20261022	custom_z	reuse	region	gaussian	none	none	150	7	4	5	uniform	4	1	4	1	8	3000	1000	2000	0.97	14	1234	100	heldout	1	0	20	20	date_enso	case_267
21	n021_rescore_shape150_case270_mld_vary_region_date_blocked	rescore_anchor	mld_vary_region	p20	0.2	date_blocked	date_blocked_p20_repeat_20261022	repeat	20261022	custom_z	reuse	region	gaussian	none	none	150	7	4	5	uniform	4	1	4	1	8	3000	1000	2000	0.97	14	1234	100	heldout	1	0	21	21	date_enso	case_270
EOF_NEXT150
cat > "$exp/preflight.R" <<'EOF_NEXT150'
# Preflight on the login node: build Stan code for every case without sampling.
experiment<-commandArgs(trailingOnly=TRUE)[1]
Sys.setenv(EXPERIMENT_DIR=experiment)
source(file.path(experiment,"common.R"));source(file.path(experiment,"models.R"))
cases<-read.delim(file.path(experiment,"cases.tsv"),stringsAsFactors=FALSE)
levels<-readRDS(file.path(experiment,"region_levels.rds"))
source_exp<-readLines(file.path(experiment,"reuse_source.txt"))[1]
out<-list()
for(i in seq_len(nrow(cases))) {
  case<-cases[i,];d<-augment(readRDS(data_file(experiment,case)))
  row<-data.frame(case_id=case$case_id,valid=TRUE,reason="",training_casts=length(unique(d$cast_key)),
    observations=nrow(d),formula="")
  r<-tryCatch({
    if(case$mode=="reuse") {
      src<-file.path(source_exp,"runs",case$source_run)
      if(!file.exists(file.path(src,"FIT_COMPLETE.txt")))stop("source fit not complete: ",src)
      row$formula<-paste("reuse",src)
    } else {
      spec<-campaign_spec(case,levels)
      code<-brms::make_stancode(spec$formula,data=d,prior=spec$priors,threads=thread_specification(case))
      if(!nzchar(paste(code,collapse="")))stop("empty Stan code")
      row$formula<-paste(capture.output(spec$formula),collapse=" ")
    }
    if(case$family=="lognormal")cat(sprintf("case %d: %d of %d rows floored (Chl<=0)\n",case$case_id,
      sum(d$chl_norm+1<ratio_floor),nrow(d)))
    NULL},error=function(e)conditionMessage(e))
  if(!is.null(r)){row$valid<-FALSE;row$reason<-r}
  cat(sprintf("case %2d %-45s %s %s\n",case$case_id,substr(case$label,1,45),if(row$valid)"OK" else "FAIL",row$reason))
  out[[i]]<-row
}
p<-do.call(rbind,out);write.csv(p,file.path(experiment,"preflight_cases.csv"),row.names=FALSE)
if(!all(p$valid))quit(status=3L)
cat("All",nrow(p),"cases passed preflight\n")
EOF_NEXT150

echo '{"max_active": 14, "max_cpus": 240}' > "$exp/dispatch_limits.json"
python3 - "$exp" "$shape" <<'PY'
import json, sys, hashlib, datetime
from pathlib import Path
exp, shape = Path(sys.argv[1]), Path(sys.argv[2])
meta = json.loads((shape / "campaign.json").read_text())
digest = hashlib.sha256(b"".join((exp / f).read_bytes() for f in ["models.R", "02_fit_case.R", "scoring.R", "cases.tsv"])).hexdigest()
meta.update(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
            token=hashlib.sha256(str(exp).encode()).hexdigest()[:10], package_sha256="next150-" + digest)
(exp / "campaign.json").write_text(json.dumps(meta, indent=2) + "\n")
print("token", meta["token"])
PY

echo; echo "=== preflight (builds Stan code for every case; no sampling)"
source "$exp/ibex_environment.sh"
if ! Rscript --vanilla "$exp/preflight.R" "$exp"; then
  echo "PREFLIGHT FAILED: nothing submitted. Paste this output to Claude."; exit 3
fi
date -u +%FT%TZ > "$exp/PREPARATION_COMPLETE.txt"
echo "$exp" > logs/latest_next150_experiment.txt

echo; echo "=== submitting controller"
python3 - "$exp" <<'PY'
import sys
sys.path.insert(0, sys.argv[1])
import scheduler as s
ledger = s.Ledger(sys.argv[1])
print("controller job", ledger.ensure("controller_001", s.job_args(ledger, "controller")))
PY
sleep 5
squeue -u "$USER" -o "%.10i %.45j %.9T %.10M" | grep -E "JOBID|next|$(python3 -c "import json;print(json.load(open('$exp/campaign.json'))['token'])")" || true
