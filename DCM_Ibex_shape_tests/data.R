# Whole-cast sampling uses predictors and fixed seeds, never chlorophyll values.
read_source <- function(source,name) {
  p<-file.path(source,"data",paste0(name,".rds"))
  if(!file.exists(p))stop("Required frozen source data missing: ",p)
  readRDS(p)
}
cast_metadata <- function(d) {
  cols<-c("cast_key","Region","Ocean","date","enso","log10_distance","MLD_temp_MIX")
  m<-unique(d[cols]);m<-m[order(as.character(m$cast_key),method="radix"),,drop=FALSE]
  if(anyDuplicated(m$cast_key))stop("Predictors or groups vary within a cast")
  m
}
allocate_exact <- function(metadata,target,region_levels,lower=NULL,upper=NULL) {
  sizes<-as.integer(table(factor(metadata$Region,levels=region_levels)))
  if(is.null(lower))lower<-pmin(4L,sizes)
  if(is.null(upper))upper<-sizes
  if(sum(lower)>target||sum(upper)<target||any(lower>upper))stop("Cannot allocate requested whole-cast fraction")
  out<-lower;quota<-target*sizes/sum(sizes)
  while(sum(out)<target){gap<-quota-out;gap[out>=upper]<- -Inf;j<-which.max(gap);out[j]<-out[j]+1L}
  out
}
sample_keys <- function(full_meta,pool,target,region_levels,seed,mandatory=character(),allowed=NULL) {
  if(!is.null(allowed))pool<-pool[pool$cast_key %in% allowed,,drop=FALSE]
  must<-full_meta[full_meta$cast_key %in% mandatory,,drop=FALSE]
  if(length(setdiff(mandatory,pool$cast_key)))stop("Mandatory casts absent from sampling pool")
  lower<-as.integer(table(factor(must$Region,levels=region_levels)))
  counts<-allocate_exact(full_meta,target,region_levels,pmax(lower,pmin(4L,
    as.integer(table(factor(pool$Region,levels=region_levels))))),
    as.integer(table(factor(pool$Region,levels=region_levels))))
  out<-mandatory
  for(i in seq_along(region_levels)) {
    candidates<-sort(as.character(pool$cast_key[as.character(pool$Region)==region_levels[i] & !pool$cast_key %in% mandatory]))
    need<-counts[i]-lower[i];set.seed(as.integer(seed)+i)
    if(need>0L)out<-c(out,candidates[sample.int(length(candidates),need)])
  }
  stopifnot(length(unique(out))==target);out
}

prepare_data <- function(experiment,source) {
  dir.create(file.path(experiment,"data"),recursive=TRUE,showWarnings=FALSE)
  # The source bundle records the original quality filtering and preprocessing.
  full_path<-file.path(source,"full_input_250m.rds")
  full<-validate_bundle(readRDS(full_path),full=TRUE)
  full<-full[full$depth>=10 & full$depth<=150,,drop=FALSE]
  levels<-readRDS(file.path(source,"region_levels.rds"))
  scaling<-readRDS(file.path(source,"scaling.rds"))
  meta<-cast_metadata(full);total<-nrow(meta)
  validation<-read_source(source,"primary_validation_d150")
  reserve<-read_source(source,"primary_reserve_d150")
  blocked<-read_source(source,"date_blocked_validation_d150")
  anchors<-lapply(c("p10","p20","p50"),function(n)read_source(source,paste0("primary_",n,"_d150")))
  names(anchors)<-c("p10","p20","p50")
  keys<-function(d)unique(as.character(d$cast_key))
  stopifnot(length(keys(validation))>=ceiling(.1*total),length(keys(reserve))>=ceiling(.1*total),
    !any(keys(validation) %in% keys(reserve)),all(keys(anchors$p10) %in% keys(anchors$p20)),
    all(keys(anchors$p20) %in% keys(anchors$p50)))
  raw_cols<-c("cast_key","Region","Ocean","date","enso","depth","log10_distance","MLD_temp_MIX","chl_norm")
  canonical<-function(d){x<-as.data.frame(lapply(d[raw_cols],as.character),stringsAsFactors=FALSE)
    x<-x[order(x$cast_key,as.numeric(x$depth),method="radix"),,drop=FALSE];rownames(x)<-NULL;x}
  for(d in c(anchors,list(validation,reserve,blocked))) {
    expected<-full[full$cast_key %in% keys(d),,drop=FALSE]
    if(!identical(canonical(d),canonical(expected)))stop("Frozen source observations differ from full input")
  }
  pool<-meta[!meta$cast_key %in% c(keys(validation),keys(reserve)),,drop=FALSE]
  if(any(keys(anchors$p50) %in% c(keys(validation),keys(reserve))))stop("Source training/holdout overlap")
  selected<-list(p10=keys(anchors$p10),p20=keys(anchors$p20),p50=keys(anchors$p50))
  selected$p40<-sample_keys(meta,pool,round(.4*total),levels,20261022,
    mandatory=selected$p20,allowed=selected$p50)
  selected$p75<-sample_keys(meta,pool,round(.75*total),levels,20261022,mandatory=selected$p50)
  cases<-read.delim(file.path(experiment,"cases.tsv"),stringsAsFactors=FALSE)
  sets<-unique(cases[c("dataset_id","subset","fraction","selection","selection_seed","split")])
  membership<-list();split_counts<-list()
  for(i in seq_len(nrow(sets))) {
    s<-sets[i,];target<-if(s$fraction==1)total else max(ceiling(.1*total),as.integer(round(s$fraction*total)))
    if(s$selection=="all")chosen<-keys(full) else if(s$selection=="frozen")chosen<-selected[[s$subset]] else {
      eligible<-pool
      if(s$split=="date_blocked")eligible<-eligible[!group_key(eligible) %in% group_key(blocked),,drop=FALSE]
      chosen<-sample_keys(meta,eligible,target,levels,s$selection_seed)
    }
    if(length(chosen)!=target)stop("Requested training fraction was not preserved: ",s$dataset_id)
    train<-full[full$cast_key %in% chosen,,drop=FALSE]
    if(s$fraction<1) {
      stopifnot(!any(keys(train) %in% c(keys(validation),keys(reserve))))
      if(s$split=="date_blocked")stopifnot(!any(group_key(train) %in% group_key(blocked)))
    }
    d<-apply_features(train,scaling[[s$split]],levels)
    stopifnot(identical(canonical(train),canonical(d)))
    saveRDS(d,file.path(experiment,"data",paste0(s$dataset_id,".rds")))
    # A bounded training audit for subset fits; full-data fits score all casts.
    audit_keys<-if(s$fraction==1)chosen else sample_keys(meta,cast_metadata(train),
      ceiling(.1*total),levels,20261022)
    audit<-d[d$cast_key %in% audit_keys,,drop=FALSE]
    saveRDS(audit,file.path(experiment,"data",paste0(s$dataset_id,"_training_audit.rds")))
    z<-cast_metadata(d);z$dataset_id<-s$dataset_id;membership[[s$dataset_id]]<-z
    split_counts[[s$dataset_id]]<-data.frame(dataset_id=s$dataset_id,casts=length(chosen),
      observations=nrow(d),fraction=length(chosen)/total,split=s$split,
      heldout_eligible=s$fraction<1,selection_seed=s$selection_seed)
  }
  for(split in c("primary","date_blocked")) {
    for(name in c("validation","reserve")) {
      d<-if(name=="reserve")reserve else if(split=="primary")validation else blocked
      d<-apply_features(d,scaling[[split]],levels)
      saveRDS(d,file.path(experiment,"data",paste0(split,"_",name,".rds")))
      z<-cast_metadata(d);z$dataset_id<-paste0(split,"_",name);membership[[z$dataset_id[1]]]<-z
    }
  }
  saveRDS(scaling,file.path(experiment,"scaling.rds"));saveRDS(levels,file.path(experiment,"region_levels.rds"))
  write.csv(meta,file.path(experiment,"full_cast_metadata.csv"),row.names=FALSE)
  write.csv(do.call(rbind,membership),file.path(experiment,"split_membership.csv"),row.names=FALSE)
  write.csv(do.call(rbind,split_counts),file.path(experiment,"split_counts.csv"),row.names=FALSE)
  write.csv(data.frame(source=full_path,md5=unname(tools::md5sum(full_path)),casts=total,
    depth_min=min(full$depth),depth_max=max(full$depth),raw_identity_verified=TRUE),
    file.path(experiment,"source_identity.csv"),row.names=FALSE)
  invisible(list(casts=total,levels=levels,scaling=scaling))
}
