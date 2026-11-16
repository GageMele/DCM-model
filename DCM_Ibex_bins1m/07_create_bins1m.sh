#!/bin/bash
# bins1m: the depth-bin ("histogram") models at 1 m resolution. Every measured depth
# 10-150 m is its own bin (141 bins) and the raw rows are fitted directly (no 5 m
# medians); the DCM-depth model uses 141 1 m classes. Adds a lognormal bin-means pair.
# Creates and launches a new experiment. Reuses shape150's frozen data (symlinked,
# read-only) and the shared scheduler/check code. Does not touch running campaigns.
set -euo pipefail
root=/ibex/user/meleg/DCM_MODEL
cd "$root"
shape=$(cat logs/latest_shape150_experiment.txt)
r3=$(cat logs/latest_round3_experiment.txt)

stamp=$(date -u +%Y%m%dT%H%M%SZ)
exp="$root/experiments/bins1m_${stamp}_$$"
echo; echo "shape150 data:  $shape/data"
echo "new experiment: $exp"
mkdir -p "$exp/logs" "$exp/scheduler" "$exp/runs"
for f in common.R data.R 03_check_case.R 04_collect_results.R onset_metrics.R profile_helpers.R \
         scheduler.py worker.slurm ibex_environment.sh scaling.rds region_levels.rds; do
  cp "$r3/$f" "$exp/$f"
done
ln -s "$shape/data" "$exp/data"
echo "$shape" > "$exp/reuse_source.txt"
cat > "$exp/models.R" <<'EOF_BINS1M'
# Round3 equations: next150 (literature-backed follow-ups to shape150), plus
# Gaussian versions of the lognormal next150 models and depth-bin ("histogram")
# models at 5/10 m (round3) or 1 m (bins1m) resolution. Inherits the shape150
# feature construction; nothing is fitted when sourced.
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
  d$depth_bin<-factor(bin_labels[bin_index(d$depth)],levels=bin_labels)
  d$depth_m<-factor(m_labels[m_index(d$depth)],levels=m_labels)
  if("chl_norm" %in% names(d)) {
    d$ratio<-pmax(d$chl_norm+1,ratio_floor)
    d$enh_class<-as.integer(cut(d$chl_norm,enh_breaks,right=FALSE))
  }
  d
}

# Depth bins for the histogram models: 28 fixed 5 m bins from 10 m (10-15 ... 145-150 m,
# 150 m joins the last bin), labelled by their top depth.
bin_index <- function(depth)as.integer(pmin(pmax(floor((depth-10)/5),0),27))+1L
bin_labels <- sprintf("z%03d",seq(10,145,5))
# 1 m bins: every measured depth 10..150 m is its own bin (141 bins).
m_index <- function(depth)as.integer(pmin(pmax(round(depth),10),150))-9L
m_labels <- sprintf("m%03d",10:150)
# Enhancement classes relative to the 5-15 m surface median (chl_norm):
# 1 <+10%, 2 +10-50%, 3 +50-100%, 4 +100-200%, 5 >=+200% (>= 3x surface).
enh_breaks <- c(-Inf,.1,.5,1,2,Inf)
enh_labels <- c("<+10%","+10-50%","+50-100%","+100-200%",">=+200%")
# DCM depth classes (one per cast): 14 fixed 10 m bins, 10-20 ... 140-150 m (dcm10m),
# or 141 1 m bins (dcm1m).
dcm_index <- function(depth)as.integer(pmin(pmax(floor((depth-10)/10),0),13))+1L
dcm_labels <- sprintf("%d-%d m",seq(10,140,10),seq(20,150,10))
dcm_width <- function(case)if(identical(case$cast,"dcm1m"))1 else 10
dcm_classes <- function(case)if(dcm_width(case)==1)sprintf("%d m",10:150) else dcm_labels

# One row per cast: the bin with the highest enhancement is the DCM (5 m bin medians
# for dcm10m; the raw 1 m rows for dcm1m).
dcm_cast_data <- function(b,width=10) {
  b<-b[order(b$cast_key,-b$chl_norm,b$depth),,drop=FALSE];top<-b[!duplicated(b$cast_key),,drop=FALSE]
  top$dcm_depth<-top$depth;top$dcm_enhancement<-top$chl_norm
  top$dcm_class<-if(width==1)m_index(top$depth) else dcm_index(top$depth)
  rownames(top)<-NULL;top
}

# Data actually fitted for a case: raw rows, 5 m bin medians, or one DCM row per cast.
prepare_data <- function(d,case)switch(case$cast,none=,raw1m=d,bin5m=bin_depth(d,5),
  dcm10m=dcm_cast_data(bin_depth(d,5)),dcm1m=dcm_cast_data(d,1),
  stop("Unknown cast treatment: ",case$cast))
ordinal_case <- function(case)case$equation %in% c("bin_class","dcm_depth")
# The depth-bin factor a bin model uses: 1 m (depth_m) or 5 m (depth_bin).
bin_var <- function(case)if(case$cast %in% c("raw1m","dcm1m"))"depth_m" else "depth_bin"
bin_terms <- function(case)sprintf(c("(1+mld_z||%s)","(1+mld_z||Region:%s)"),bin_var(case))
# Population predictions drop date/ENSO effects only. In the bin models the depth-bin
# terms are group-level effects too, and they are the model's depth structure, so they stay.
# (The bin_means sigma uses fixed depth-bin effects because re_formula would drop them.)
population_re <- function(case)switch(case$equation,
  bin_means=,bin_class=as.formula(paste("~",paste(bin_terms(case),collapse="+"))),
  dcm_depth=~(0+mld_z|Region),NA)

# One-pair 5 m binning test (cast == "bin5m"): each cast's chl_norm is reduced to
# its median within 5 m depth bins, placed at the bin centre. Only the fitted data
# are binned; every holdout is still scored on the raw unbinned rows.
bin_depth <- function(d,width) {
  if(!width)return(d)
  derived<-c("depth","chl_norm","depth_mld","depth_mld_half","depth_hm","ratio","depth_bin","depth_m","enh_class",grep("^warp_",names(d),value=TRUE))
  varying<-vapply(setdiff(names(d),derived),function(v)
    any(tapply(as.character(d[[v]]),d$cast_key,function(x)length(unique(x))>1)),logical(1))
  # Other depth-varying columns are not used by the binned equation, so they are dropped.
  if(any(varying)){message("Binning drops depth-varying columns: ",paste(names(varying)[varying],collapse=", "))
    d<-d[setdiff(names(d),names(varying)[varying])]}
  lo<-floor(min(d$depth));top<-ceiling((max(d$depth)-lo)/width)-1
  b<-pmin(floor((d$depth-lo)/width),top);key<-paste(d$cast_key,b,sep="\034")
  first<-!duplicated(key);out<-d[first,,drop=FALSE]
  out$depth<-lo+width*(b[first]+.5)
  out$chl_norm<-as.numeric(tapply(d$chl_norm,key,median)[key[first]])
  out$rows_in_bin<-as.integer(table(key)[key[first]])
  out$depth_mld<-out$depth-out$MLD_temp_MIX;out$depth_mld_half<-out$depth-.5*out$MLD_temp_MIX
  for(v in grep("^warp_half_",names(out),value=TRUE)) {
    lambda<-as.numeric(sub("^warp_half_(pos|neg)","",v))/100*ifelse(grepl("_neg",v),-1,1)
    out[[v]]<-out$depth_mld_half*exp(lambda*out$mld_z)
  }
  out<-augment(out);rownames(out)<-NULL
  out[order(out$cast_key,out$depth),,drop=FALSE]
}

groups_terms <- function(case)
  if(case$groups=="date_only")"(1|Ocean:Region:date)" else c("(1|Ocean:Region:date)","(1|Ocean:Region:enso)")

sigma_formula <- function(case)
  if(case$variance=="region_depth")sigma~0+Region+s(depth,by=Ocean,bs='cr',k=5) else
  if(case$variance=="region_bin")as.formula(paste("sigma~0+Region+",bin_var(case))) else sigma~0+Region

campaign_spec <- function(case,region_levels) {
  kd<-as.integer(case$k_depth)
  family<-if(case$family=="lognormal")brms::lognormal() else gaussian()
  response<-if(case$family=="lognormal")"ratio" else "chl_norm"
  u<-"depth_mld_half"
  regional<-function(v,k=kd,m=2)sprintf("s(%s,by=Region,bs='cr',k=%d%s)",v,k,if(m==1)",m=1" else "")
  pooled_mld<-vapply(seq_along(region_levels),function(i)
    sprintf("s(%s,by=region_mld_%d,bs='cr',k=5,m=1)",u,i),character(1))
  if(case$equation=="peak_nl") return(peak_spec(case,family))
  if(case$equation %in% c("bin_means","bin_class","dcm_depth")) return(bin_spec(case))
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
# The Gaussian version models chl_norm = ratio - 1 directly on its own scale.
peak_spec <- function(case,family) {
  profile<-if(case$family=="lognormal")
    ratio~log(exp(lbg+slope*depth_hm)+exp(lamp)*exp(-0.5*((depth_hm-zm)/exp(lwid))^2)) else
    chl_norm~exp(lbg+slope*depth_hm)+exp(lamp)*exp(-0.5*((depth_hm-zm)/exp(lwid))^2)-1
  form<-brms::bf(profile,
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

# Depth-bin ("histogram") models: depth is a set of 28 separate 5 m bins rather than
# a smooth curve, so adjacent bins are not forced to be similar; bins and region x bin
# cells are partially pooled. MLD gets its own effect in every bin (and region x bin),
# the bin analogue of s(depth, by = MLD). Same distance, date and ENSO terms as the GAMMs.
#  bin_means: mean enhancement of each bin (a bar chart, not a curve); Gaussian, or
#             lognormal on ratio = chl_norm + 1.
#  bin_class: ordinal probability that each bin reaches each enhancement class.
#  dcm_depth: ordinal distribution of which 10 m bin holds each cast's DCM.
bin_spec <- function(case) {
  bins<-bin_terms(case)
  rhs<-switch(case$equation,
    bin_means=c("1","Region","distance_z","mld_z",bins),
    bin_class=c("1","Region","distance_z","mld_z",bins),
    dcm_depth=c("1","Region","distance_z","mld_z","(0+mld_z|Region)"))
  lhs<-switch(case$equation,bin_means=if(case$family=="lognormal")"ratio" else "chl_norm",
    bin_class=sprintf("enh_class|thres(%d)",length(enh_labels)-1),
    dcm_depth=sprintf("dcm_class|thres(%d)",length(dcm_classes(case))-1))
  f<-as.formula(paste(lhs,"~",paste(c(rhs,groups_terms(case)),collapse=" + ")),env=asNamespace("brms"))
  scale_prior<-"student_t(3,0,1.5)"
  if(case$equation=="bin_means") {
    lognormal<-case$family=="lognormal"
    form<-brms::bf(f,sigma_formula(case),family=if(lognormal)brms::lognormal() else gaussian())
    priors<-c(brms::set_prior(if(lognormal)"student_t(3,0,1)" else "student_t(3,0.73,5.448)",class="Intercept"),
      brms::set_prior("normal(0,2)",class="b"),brms::set_prior(scale_prior,class="sd"),
      brms::set_prior("normal(0,1)",class="b",dpar="sigma"))
  } else {
    form<-brms::bf(f,family=brms::cumulative("logit"))
    priors<-c(brms::set_prior("student_t(3,0,2.5)",class="Intercept"),
      brms::set_prior("normal(0,2)",class="b"),brms::set_prior(scale_prior,class="sd"))
  }
  list(formula=form,priors=priors)
}
EOF_BINS1M
cat > "$exp/02_fit_case.R" <<'EOF_BINS1M'
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
  data_path<-data_file(experiment,case);d<-prepare_data(augment(readRDS(data_path)),case)
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
EOF_BINS1M
cat > "$exp/scoring.R" <<'EOF_BINS1M'
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
EOF_BINS1M
cat > "$exp/preflight.R" <<'EOF_BINS1M'
# Preflight on the login node: build Stan code for every case without sampling.
# Optional second argument: comma-separated case ids, written to preflight_cases_<ids>.csv
# for merging (the scheduler only reads preflight_cases.csv).
args<-commandArgs(trailingOnly=TRUE);experiment<-args[1]
Sys.setenv(EXPERIMENT_DIR=experiment)
source(file.path(experiment,"common.R"));source(file.path(experiment,"models.R"))
cases<-read.delim(file.path(experiment,"cases.tsv"),stringsAsFactors=FALSE)
only<-if(length(args)>1)as.integer(strsplit(args[2],",")[[1]]) else integer(0)
if(length(only)){if(!all(only %in% cases$case_id))stop("Unknown case ids");cases<-cases[cases$case_id %in% only,]}
levels<-readRDS(file.path(experiment,"region_levels.rds"))
source_exp<-readLines(file.path(experiment,"reuse_source.txt"))[1]
out<-list()
for(i in seq_len(nrow(cases))) {
  case<-cases[i,];d<-prepare_data(augment(readRDS(data_file(experiment,case))),case)
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
    if(case$cast!="none")cat(sprintf("case %d: %s data, %d rows from %d casts\n",case$case_id,case$cast,
      nrow(d),length(unique(d$cast_key))))
    if(case$equation=="dcm_depth")cat("  DCM classes in training data:",tabulate(d$dcm_class,length(dcm_classes(case))),"\n")
    if(case$equation=="bin_class")cat("  enhancement classes in training data:",tabulate(d$enh_class,length(enh_labels)),"\n")
    if(case$family=="lognormal")cat(sprintf("case %d: %d of %d rows floored (Chl<=0)\n",case$case_id,
      sum(d$chl_norm+1<ratio_floor),nrow(d)))
    NULL},error=function(e)conditionMessage(e))
  if(!is.null(r)){row$valid<-FALSE;row$reason<-r}
  cat(sprintf("case %2d %-45s %s %s\n",case$case_id,substr(case$label,1,45),if(row$valid)"OK" else "FAIL",row$reason))
  out[[i]]<-row
}
p<-do.call(rbind,out)
write.csv(p,file.path(experiment,if(length(only))paste0("preflight_cases_",paste(only,collapse="_"),".csv") else "preflight_cases.csv"),row.names=FALSE)
if(!all(p$valid))quit(status=3L)
cat("All",nrow(p),"cases passed preflight\n")
EOF_BINS1M
cat > "$exp/cases.tsv" <<'EOF_BINS1M'
case_id	label	block	equation	subset	fraction	split	dataset_id	selection	selection_seed	prior	mode	variance	family	cast	autocor	depth_max	k_depth	k_mld	k_distance	k_profile	chains	cores	threads	cpus	memory_gb	iter	warmup	samples	adapt_delta	treedepth	seed	refresh	evaluation_role	fit_days	priority	equation_reference_case	settings_reference_case	groups	source_run
1	b001_hist1m_bin_means_primary_p20	hist1m_bin_means	bin_means	p20	0.2	primary	primary_p20_frozen	frozen	0	custom_z	posterior	region_bin	gaussian	raw1m	none	150	7	4	5	uniform	4	4	4	16	64	3000	1000	2000	0.95	14	1234	100	heldout	7	1	1	1	date_enso	NA
2	b002_hist1m_bin_means_date_blocked_p20	hist1m_bin_means	bin_means	p20	0.2	date_blocked	date_blocked_p20_repeat_20261022	repeat	20261022	custom_z	posterior	region_bin	gaussian	raw1m	none	150	7	4	5	uniform	4	4	4	16	64	3000	1000	2000	0.95	14	1234	100	heldout	7	2	2	2	date_enso	NA
3	b003_hist1m_bin_class_primary_p20	hist1m_bin_class	bin_class	p20	0.2	primary	primary_p20_frozen	frozen	0	custom_z	posterior	region	cumulative	raw1m	none	150	7	4	5	uniform	4	4	4	16	64	3000	1000	2000	0.95	14	1234	100	heldout	7	3	3	3	date_enso	NA
4	b004_hist1m_bin_class_date_blocked_p20	hist1m_bin_class	bin_class	p20	0.2	date_blocked	date_blocked_p20_repeat_20261022	repeat	20261022	custom_z	posterior	region	cumulative	raw1m	none	150	7	4	5	uniform	4	4	4	16	64	3000	1000	2000	0.95	14	1234	100	heldout	7	4	4	4	date_enso	NA
5	b005_hist1m_dcm_depth_primary_p20	hist1m_dcm_depth	dcm_depth	p20	0.2	primary	primary_p20_frozen	frozen	0	custom_z	posterior	region	cumulative	dcm1m	none	150	7	4	5	uniform	4	4	1	4	16	3000	1000	2000	0.95	14	1234	100	heldout	7	5	5	5	date_enso	NA
6	b006_hist1m_dcm_depth_date_blocked_p20	hist1m_dcm_depth	dcm_depth	p20	0.2	date_blocked	date_blocked_p20_repeat_20261022	repeat	20261022	custom_z	posterior	region	cumulative	dcm1m	none	150	7	4	5	uniform	4	4	1	4	16	3000	1000	2000	0.95	14	1234	100	heldout	7	6	6	6	date_enso	NA
7	b007_hist1m_bin_means_lognormal_primary_p20	hist1m_bin_means_lognormal	bin_means	p20	0.2	primary	primary_p20_frozen	frozen	0	custom_z	posterior	region_bin	lognormal	raw1m	none	150	7	4	5	uniform	4	4	4	16	64	3000	1000	2000	0.95	14	1234	100	heldout	7	7	7	7	date_enso	NA
8	b008_hist1m_bin_means_lognormal_date_blocked_p20	hist1m_bin_means_lognormal	bin_means	p20	0.2	date_blocked	date_blocked_p20_repeat_20261022	repeat	20261022	custom_z	posterior	region_bin	lognormal	raw1m	none	150	7	4	5	uniform	4	4	4	16	64	3000	1000	2000	0.95	14	1234	100	heldout	7	8	8	8	date_enso	NA
EOF_BINS1M
(cd "$exp" && md5sum -c --quiet) <<'SUMS' || { echo "embedded files corrupted (bad download/upload?); nothing submitted"; exit 1; }
3139ddd25d7bdc425a554b525601257b  models.R
a5daae56422eb8ab54a325b53592e42f  02_fit_case.R
efaa531eaf4533aae22e0602b7538b85  scoring.R
cdaa789771518be99b4a1dae0a9f20ce  preflight.R
1b2b506a135aff28cba83ff635efb1a4  cases.tsv
SUMS

echo '{"max_active": 8, "max_cpus": 120}' > "$exp/dispatch_limits.json"
python3 - "$exp" "$shape" <<'PY'
import json, sys, hashlib, datetime
from pathlib import Path
exp, shape = Path(sys.argv[1]), Path(sys.argv[2])
meta = json.loads((shape / "campaign.json").read_text())
digest = hashlib.sha256(b"".join((exp / f).read_bytes() for f in ["models.R", "02_fit_case.R", "scoring.R", "cases.tsv"])).hexdigest()
meta.update(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
            token=hashlib.sha256(str(exp).encode()).hexdigest()[:10], package_sha256="bins1m-" + digest)
(exp / "campaign.json").write_text(json.dumps(meta, indent=2) + "\n")
print("token", meta["token"])
PY

echo; echo "=== preflight (builds Stan code for every case; no sampling)"
source "$exp/ibex_environment.sh"
if ! Rscript --vanilla "$exp/preflight.R" "$exp"; then
  echo "PREFLIGHT FAILED: nothing submitted. Paste this output to Claude."; exit 3
fi
date -u +%FT%TZ > "$exp/PREPARATION_COMPLETE.txt"
echo "$exp" > logs/latest_bins1m_experiment.txt

echo; echo "=== submitting controller"
python3 - "$exp" <<'PY'
import sys
sys.path.insert(0, sys.argv[1])
import scheduler as s
ledger = s.Ledger(sys.argv[1])
print("controller job", ledger.ensure("controller_001", s.job_args(ledger, "controller")))
PY
sleep 5
tok=$(python3 -c "import json;print(json.load(open('$exp/campaign.json'))['token'])")
squeue -u "${USER:-$(whoami)}" -o "%.10i %.45j %.9T %.10M %.6C" | grep -E "JOBID|$tok" || true
