# MLD-shape equations. No fitting occurs when this file is sourced.
run_directory <- function(experiment,case)file.path(experiment,"runs",sprintf("case_%03d",case$case_id))
data_file <- function(experiment,case)file.path(experiment,"data",paste0(case$dataset_id,".rds"))
thread_specification <- function(case)brms::threading(as.integer(case$threads))
group_key <- function(d)paste(d$Ocean,d$Region,d$date,sep="\034")

apply_features <- function(d,scaling,region_levels) {
  d$Region<-factor(as.character(d$Region),levels=region_levels)
  d$Ocean<-factor(as.character(d$Ocean));d$date<-as.character(d$date);d$enso<-as.character(d$enso)
  d$mld_z<-(d$MLD_temp_MIX-scaling$mld_mean)/scaling$mld_sd
  d$distance_z<-(d$log10_distance-scaling$distance_mean)/scaling$distance_sd
  d$mld_distance_z<-d$mld_z*d$distance_z
  d$depth_mld<-d$depth-d$MLD_temp_MIX
  d$depth_mld_half<-d$depth-.5*d$MLD_temp_MIX
  # Positive lambda compresses the physical width when MLD is deeper.
  d$warp_half_pos025<-d$depth_mld_half*exp(.25*d$mld_z)
  d$warp_half_pos050<-d$depth_mld_half*exp(.50*d$mld_z)
  d$warp_half_neg025<-d$depth_mld_half*exp(-.25*d$mld_z)
  for(i in seq_along(region_levels)) {
    yes<-as.character(d$Region)==region_levels[i]
    # An ordered factor fits only the 'on' smooth. Its centred factor smooth
    # avoids duplicating the Region intercept when regional k values differ.
    d[[paste0("region_only_",i)]]<-ordered(ifelse(yes,"on","off"),levels=c("off","on"))
    d[[paste0("region_mld_",i)]]<-as.numeric(yes)*d$mld_z
  }
  if(anyNA(d$Region)||any(!is.finite(d$mld_z)))stop("Invalid transformed predictors")
  d
}

campaign_spec <- function(case,region_levels) {
  kd<-as.integer(case$k_depth);km<-as.integer(case$k_mld)
  kval<-function(i)if(case$k_profile=="rs15_other7")
    if(region_levels[i] %in% c("Central RS","South RS"))15L else 7L else kd
  regional<-function(variable) {
    if(case$k_profile=="uniform")return(sprintf("s(%s,by=Region,bs='cr',k=%d)",variable,kd))
    vapply(seq_along(region_levels),function(i)
      sprintf("s(%s,by=region_only_%d,bs='cr',k=%d)",variable,i,kval(i)),character(1))
  }
  regional_mld<-function(variable) vapply(seq_along(region_levels),function(i)
    sprintf("s(%s,by=region_mld_%d,bs='cr',k=%d)",variable,i,kval(i)),character(1))
  tensor<-function(variable)sprintf("t2(%s,mld_z,by=Region,bs=c('cr','cr'),k=c(%d,%d))",variable,kd,km)
  terms<-switch(case$equation,
    original_z=c("distance_z",regional("depth")),
    interaction_z=c(regional("depth"),"s(depth,by=mld_z,bs='cr',k=5)",
      "s(depth,by=distance_z,bs='cr',k=5)","s(depth,by=mld_distance_z,bs='cr',k=5)"),
    mld_halfshift=c("distance_z",regional("depth_mld_half")),
    mld_shift=c("distance_z",regional("depth_mld")),
    mld_vary_shared=c("distance_z",regional("depth"),sprintf("s(depth,by=mld_z,bs='cr',k=%d)",kd)),
    mld_vary_region=c("distance_z",regional("depth"),regional_mld("depth")),
    tensor_region=c("distance_z",tensor("depth")),
    tensor_shift_region=c("distance_z",tensor("depth_mld")),
    tensor_halfshift_region=c("distance_z",tensor("depth_mld_half")),
    warp_half_pos025=c("distance_z",regional("warp_half_pos025")),
    warp_half_pos050=c("distance_z",regional("warp_half_pos050")),
    warp_half_neg025=c("distance_z",regional("warp_half_neg025")),
    stop("Unknown equation: ",case$equation))
  f<-as.formula(paste("chl_norm ~",paste(c("1","Region",terms,
    "(1|Ocean:Region:date)","(1|Ocean:Region:enso)"),collapse=" + ")),env=asNamespace("brms"))
  form<-brms::bf(f,sigma~0+Region,family=gaussian())
  scale_prior<-if(case$prior=="scales_halfnormal")"normal(0,1.5)" else "student_t(3,0,1.5)"
  priors<-c(brms::set_prior("student_t(3,0.73,5.448)",class="Intercept"),
    brms::set_prior(if(case$prior=="b_tight")"normal(0,1)" else "normal(0,2)",class="b"),
    brms::set_prior(scale_prior,class="sd"),brms::set_prior(scale_prior,class="sds"),
    # sigma has a log link; these are proper priors on its regional coefficients.
    brms::set_prior("normal(0,1)",class="b",dpar="sigma"))
  list(formula=form,priors=priors)
}
