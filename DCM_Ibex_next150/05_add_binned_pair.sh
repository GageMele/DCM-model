#!/bin/bash
# Add the 5 m binned test pair (cases 22-23) to the running next150 campaign.
# Half-shift Gaussian (same equation as cases 25/267 in shape150) fitted to each
# cast's 5 m bin medians; every holdout is still scored on the raw 1 m rows.
# Safe while the controller runs: the two new fits are held until their
# preflight passes, and the files patched here are only read when a job starts.
set -euo pipefail
root=/ibex/user/meleg/DCM_MODEL
cd "$root"
exp=$(cat logs/latest_next150_experiment.txt)
echo "experiment: $exp"
[ -f "$exp/PREPARATION_COMPLETE.txt" ] || { echo "next150 is not prepared; nothing changed"; exit 1; }
if grep -q "^22	" "$exp/cases.tsv"; then echo "cases 22-23 already added; nothing changed"; exit 0; fi

# Only patch the exact files installed by 04_create_next150.sh.
cd "$exp"
md5sum -c --quiet <<'SUMS' || { echo "next150 files differ from the installed versions; nothing changed. Paste this output to Claude."; exit 1; }
0a5f75bb7201eeea175887e1defe2f9c  models.R
87025dd3fc36d82e74bac458a1a0a642  02_fit_case.R
6420edbc27b7be55290ba3bd93b92409  preflight.R
1ccd5d68eb8edc5a936307be3b214fc1  cases.tsv
SUMS
stamp=$(date -u +%Y%m%dT%H%M%SZ)
mkdir -p "backup_before_binned_$stamp" .binned_stage
cp -p models.R 02_fit_case.R preflight.R cases.tsv preflight_cases.csv "backup_before_binned_$stamp/"
echo "backup: $exp/backup_before_binned_$stamp"

cat > .binned_stage/models.R <<'EOF_BIN5'
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

# One-pair 5 m binning test (cast == "bin5m"): each cast's chl_norm is reduced to
# its median within 5 m depth bins, placed at the bin centre. Only the fitted data
# are binned; every holdout is still scored on the raw unbinned rows.
bin_width <- function(case)if(identical(case$cast,"bin5m"))5 else 0
bin_depth <- function(d,width) {
  if(!width)return(d)
  derived<-c("depth","chl_norm","depth_mld","depth_mld_half","depth_hm","ratio",grep("^warp_",names(d),value=TRUE))
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
EOF_BIN5
cat > .binned_stage/02_fit_case.R <<'EOF_BIN5'
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
  data_path<-data_file(experiment,case);d<-bin_depth(augment(readRDS(data_path)),bin_width(case))
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
EOF_BIN5
cat > .binned_stage/preflight.R <<'EOF_BIN5'
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
  case<-cases[i,];d<-bin_depth(augment(readRDS(data_file(experiment,case))),bin_width(case))
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
    if(bin_width(case))cat(sprintf("case %d: binned to %d rows (median %g raw rows per bin)\n",case$case_id,
      nrow(d),median(d$rows_in_bin)))
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
EOF_BIN5
cat > .binned_stage/cases.tsv <<'EOF_BIN5'
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
22	n022_binned_5m_half_primary_p20	binned_5m	half	p20	0.2	primary	primary_p20_frozen	frozen	0	custom_z	posterior	region	gaussian	bin5m	none	150	7	4	5	uniform	4	4	4	16	48	3000	1000	2000	0.95	14	1234	100	heldout	7	22	22	22	date_enso	NA
23	n023_binned_5m_half_date_blocked_p20	binned_5m	half	p20	0.2	date_blocked	date_blocked_p20_repeat_20261022	repeat	20261022	custom_z	posterior	region	gaussian	bin5m	none	150	7	4	5	uniform	4	4	4	16	48	3000	1000	2000	0.95	14	1234	100	heldout	7	23	23	23	date_enso	NA
EOF_BIN5
(cd .binned_stage && md5sum -c --quiet) <<'SUMS' || { echo "staged files corrupted (bad paste/upload?); nothing changed"; exit 1; }
b91afc5f5d01e933bd6992697066891f  models.R
2887ec2c374226090737427040325dc3  02_fit_case.R
7192e899dba9663df1a63fc6370ee822  preflight.R
7e6b601f26c30ddbb970e2b4d16e2a15  cases.tsv
SUMS

# Hold the two new fits so the controller cannot submit them before preflight.
for c in 022 023; do
  echo '{"utc": "'$stamp'", "reason": "HOLD until binned preflight passes (05_add_binned_pair.sh)", "args": []}' > "scheduler/fit_$c.rejected.json"
done
for f in models.R 02_fit_case.R preflight.R cases.tsv; do mv -f ".binned_stage/$f" "$f"; done
rmdir .binned_stage
echo "files patched; cases 22-23 added (held)"

echo; echo "=== preflight for cases 22-23 (builds Stan code; no sampling)"
source "$exp/ibex_environment.sh"
if ! Rscript --vanilla preflight.R "$exp" 22,23; then
  echo "PREFLIGHT FAILED: cases 22-23 stay held, nothing submitted. Paste this output to Claude."; exit 3
fi
python3 - "$exp" <<'PY'
import csv, json, os, sys
from pathlib import Path
exp = Path(sys.argv[1])
main = exp / "preflight_cases.csv"; add = exp / "preflight_cases_22_23.csv"
lines = main.read_text().splitlines(keepends=True); extra = add.read_text().splitlines(keepends=True)
assert lines[0] == extra[0] and len(extra) == 3
assert not any(l.split(",")[0] in ("22", "23") for l in lines[1:])
if not lines[-1].endswith("\n"): lines[-1] += "\n"
tmp = main.with_suffix(".csv.tmp"); tmp.write_text("".join(lines + extra[1:]))
os.replace(tmp, main)
# Room for the two extra fits alongside the 14 running ones and the re-score checks.
lim = exp / "dispatch_limits.json"; d = json.loads(lim.read_text())
d["max_active"] = max(int(d["max_active"]), 24); d["max_cpus"] = max(int(d["max_cpus"]), 300)
tmp = lim.with_suffix(".json.tmp"); tmp.write_text(json.dumps(d) + "\n"); os.replace(tmp, lim)
print("preflight merged; dispatch limits now", d)
PY
rm -f scheduler/fit_022.rejected.json scheduler/fit_023.rejected.json
echo "holds released: the controller submits cases 22-23 on its next pass"
echo; echo "=== queue"
tok=$(python3 -c "import json;print(json.load(open('campaign.json'))['token'])")
squeue -u "${USER:-$(whoami)}" -o "%.10i %.45j %.9T %.10M %.6C" | grep -E "JOBID|$tok" || true
