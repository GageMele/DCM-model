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
