# A controlled follow-up campaign. No model is fitted in this file.
campaign_cases <- function() {
  rows <- list()
  add <- function(block,equation="additive_z",variance="region",cast="none",
    family="gaussian",prior="custom_z",kd=7L,ke=5L,subset="p10",split="primary",
    warmup=1000L,samples=1000L,delta=.95,td=12L,threads=2L,seed=1234L,
    mode="posterior",autocor="none",depth_max=150L,reference=NA_integer_) {
    id <- length(rows)+1L
    rows[[id]] <<- data.frame(case_id=id,label=sprintf("t%03d_%s_%s_%s",id,block,equation,subset),
      block,equation,subset,split,prior,mode,variance,cast,family,autocor,depth_max,
      chains=4L,iter=warmup+samples,warmup,samples,adapt_delta=delta,treedepth=td,
      threads,cpus=4L*threads,k_depth=kd,k_mld=ke,k_distance=ke,cores=4L,
      refresh=50L,seed,memory_gb=if(subset=="p50" || autocor!="none")64L else 32L,
      reference_case=reference)
    id
  }
  anchors <- c(additive_z=2L,both_z=8L)
  for(eq in names(anchors)) {
    constant <- add("structure",eq,variance="constant")
    region <- add("structure",eq,reference=constant)
    add("structure",eq,variance="region_depth",reference=region)
    add("structure",eq,cast="intercept",reference=region)
    add("structure",eq,cast="curve",reference=region)
    add("likelihood",eq,family="student4",reference=region)
  }
  for(eq in c("original_z","mld_vary_z","distance_vary_z","region_mld_z","region_both_z","interaction_z"))
    add("equation",eq,reference=if(eq %in% c("region_both_z","interaction_z"))8L else 2L)
  for(eq in names(anchors)) for(k in c(10L,15L)) add("depth_k",eq,kd=k,reference=anchors[[eq]])
  for(eq in names(anchors)) for(p in c("b_tight","scales_halfnormal","mu_default"))
    add("priors",eq,prior=p,reference=anchors[[eq]])
  for(eq in names(anchors)) {
    ref <- anchors[[eq]]
    for(delta in c(.90,.97)) add("delta",eq,delta=delta,reference=ref)
    add("tree",eq,td=14L,reference=ref)
    for(w in c(500L,2000L)) add("warmup",eq,warmup=w,reference=ref)
    add("draws",eq,samples=2000L,reference=ref)
    for(t in c(1L,4L)) add("threads",eq,threads=t,reference=ref)
    add("repeat_seed",eq,seed=20261017L,reference=ref)
  }
  for(eq in names(anchors)) for(s in c("p20","p50")) add("data_size",eq,subset=s,reference=anchors[[eq]])
  for(eq in names(anchors)) add("date_blocked",eq,split="date_blocked",reference=anchors[[eq]])
  for(eq in names(anchors)) for(p in c("custom_z","b_tight","scales_halfnormal"))
    add("prior_predictive",eq,prior=p,mode="prior",reference=anchors[[eq]])
  for(eq in names(anchors)) add("residual_ar",eq,variance="constant",autocor="ar1",threads=1L,
    reference=if(eq=="additive_z")62L else 63L)
  # Same historical settings as the existing 150 m p10 custom baseline.
  # This 250 m fit uses exactly the same training casts; score both at 10-150 m.
  add("historical_depth_bridge","original_raw",variance="constant",kd=15L,
    depth_max=250L,samples=4000L,delta=.97,td=14L)
  for(eq in names(anchors))add("ar_reference",eq,variance="constant",threads=1L,
    reference=if(eq=="additive_z")1L else 7L)
  x <- do.call(rbind,rows)
  stopifnot(nrow(x)==63L,sum(x$mode=="posterior")==57L,!anyDuplicated(x$label),
    all(x$subset %in% c("p10","p20","p50")),all(x$iter>x$warmup))
  x
}

campaign_spec <- function(case,region_levels) {
  kd <- case$k_depth; ke <- case$k_mld
  s <- function(by,k=ke) sprintf("s(depth,by=%s,bs='cr',k=%d)",by,k)
  region <- s("Region",kd)
  terms <- switch(case$equation,
    original_raw=c("log10_distance",region),
    original_z=c("distance_z",region),
    additive_z=c("distance_z","mld_z",region),
    mld_vary_z=c("distance_z",s("mld_z"),region),
    distance_vary_z=c("mld_z",s("distance_z"),region),
    both_z=c(s("mld_z"),s("distance_z"),region),
    region_mld_z=c("distance_z",vapply(seq_along(region_levels),function(i)s(paste0("region_mld_",i)),character(1)),region),
    region_both_z=c(vapply(seq_along(region_levels),function(i)s(paste0("region_mld_",i)),character(1)),
      vapply(seq_along(region_levels),function(i)s(paste0("region_distance_",i)),character(1)),region),
    interaction_z=c(s("mld_z"),s("distance_z"),s("mld_distance_z"),region),
    stop("Unknown equation"))
  terms <- c("1","Region",terms,"(1|Ocean:Region:date)","(1|Ocean:Region:enso)")
  if(case$cast=="intercept")terms<-c(terms,"(1|cast_key)")
  if(case$cast=="curve")terms<-c(terms,"(1+cast_b1+cast_b2+cast_b3+cast_b4||cast_key)")
  if(case$autocor=="ar1")terms<-c(terms,"ar(time=depth,gr=cast_key,p=1,cov=TRUE)")
  f<-as.formula(paste("chl_norm ~",paste(terms,collapse=" + ")),env=asNamespace("brms"))
  args<-list(f)
  if(case$variance=="region")args<-c(args,list(sigma~0+Region))
  if(case$variance=="region_depth")args<-c(args,list(sigma~0+Region+s(depth,by=Region,bs="cr",k=5)))
  if(case$family=="student4")args$nu<-4
  args$family<-if(case$family=="student4")brms::student() else gaussian()
  form<-do.call(brms::bf,args)
  priors<-NULL
  if(case$prior!="mu_default") {
    scale_prior<-if(case$prior=="scales_halfnormal")"normal(0,1.5)" else "student_t(3,0,1.5)"
    priors<-c(brms::set_prior("student_t(3,0.73,5.448)",class="Intercept"),
      brms::set_prior(if(case$prior=="b_tight")"normal(0,1)" else "normal(0,2)",class="b"),
      brms::set_prior(scale_prior,class="sd"),brms::set_prior(scale_prior,class="sds"))
  }
  if(case$variance=="constant") {
    if(case$prior!="mu_default")priors<-c(priors,brms::set_prior("student_t(3,0,2.5)",class="sigma"))
  } else {
    # Log-sigma coefficients; proper priors held constant in mean-prior contrasts.
    sigma_prior<-brms::set_prior("normal(0,1)",class="b",dpar="sigma")
    priors<-if(is.null(priors))sigma_prior else c(priors,sigma_prior)
    if(case$variance=="region_depth")priors<-c(priors,brms::set_prior("normal(0,0.5)",class="sds",dpar="sigma"))
  }
  if(case$cast!="none")priors<-c(priors,brms::set_prior("normal(0,1)",class="sd",group="cast_key"))
  if(case$autocor=="ar1")priors<-c(priors,brms::set_prior("normal(0,0.3)",class="ar"))
  list(formula=form,priors=priors)
}

apply_features <- function(d,scaling,region_levels) {
  d$mld_z<-(d$MLD_temp_MIX-scaling$mld_mean)/scaling$mld_sd
  d$distance_z<-(d$log10_distance-scaling$distance_mean)/scaling$distance_sd
  d$mld_distance_z<-d$mld_z*d$distance_z
  d$depth_z<-(d$depth-80)/70
  # Fixed knots avoid using validation data to select a cast-level basis.
  b<-splines::ns(d$depth,knots=c(45,80,115),Boundary.knots=c(10,150),intercept=FALSE)
  for(i in 1:4)d[[paste0("cast_b",i)]]<-b[,i]
  for(i in seq_along(region_levels)) {
    indicator<-as.numeric(as.character(d$Region)==region_levels[i])
    d[[paste0("region_mld_",i)]]<-d$mld_z*indicator
    d[[paste0("region_distance_",i)]]<-d$distance_z*indicator
  }
  d$Region<-factor(as.character(d$Region),levels=region_levels)
  d$cast_key<-as.character(d$cast_key)
  d<-d[order(d$cast_key,d$depth,method="radix"),,drop=FALSE]
  rownames(d)<-NULL
  d
}

data_file <- function(experiment,case)file.path(experiment,"data",paste0(case$split,"_",case$subset,"_d",case$depth_max,".rds"))
training_eval_file <- function(experiment,case)file.path(experiment,"data",paste0(case$split,"_p10_d150.rds"))
validation_file <- function(experiment,case)file.path(experiment,"data",paste0(case$split,"_validation_d150.rds"))
thread_specification <- function(case) {
  if(case$autocor=="ar1" || case$block=="ar_reference")NULL else brms::threading(case$threads)
}
