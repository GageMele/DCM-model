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
