# Freeze valid whole-cast splits and construct every requested model before fitting.
prepare_campaign <- function(experiment) {
  source(file.path(experiment,"common.R"),local=TRUE)
  source(file.path(experiment,"campaign.R"),local=TRUE)
  if(file.exists(file.path(experiment,"PREPARATION_COMPLETE.txt")))stop("Already prepared")
  cases<-read.delim(file.path(experiment,"cases.tsv"),stringsAsFactors=FALSE)
  stopifnot(isTRUE(all.equal(cases,campaign_cases(),check.attributes=FALSE)))
  full250<-validate_bundle(readRDS(file.path(experiment,"full_input_250m.rds")),full=TRUE)
  full<-full250[full250$depth>=10 & full250$depth<=150,,drop=FALSE]
  sets<-setNames(lapply(c("p10","p20","p50"),function(s)
    validate_bundle(readRDS(file.path(experiment,"data",paste0(s,".rds"))))),c("p10","p20","p50"))
  validation<-readRDS(file.path(experiment,"source_validation.rds"))
  fields<-c("cast_key","chl_norm","Region","Ocean","date","enso","depth","log10_distance","MLD_temp_MIX")
  canonical<-function(x) {
    x<-as.data.frame(lapply(x[fields],as.character),stringsAsFactors=FALSE)
    x<-x[do.call(order,c(unname(x),list(method="radix"))),,drop=FALSE];rownames(x)<-NULL;x
  }
  for(d in c(sets,list(validation=validation))) {
    if(anyNA(d[fields]))stop("Missing source predictors")
    expected<-full[full$cast_key %in% d$cast_key,,drop=FALSE]
    stopifnot(identical(canonical(d),canonical(expected)))
  }
  meta<-unique(full[c("cast_key","Region","Ocean","date","enso","log10_distance","MLD_temp_MIX")])
  if(anyDuplicated(meta$cast_key))stop("Cast-level predictors or groups are not constant")
  meta<-meta[order(meta$cast_key,method="radix"),]
  total<-nrow(meta); minimum<-ceiling(.1*total)
  stopifnot(length(unique(sets$p10$cast_key))>=minimum,
    length(unique(validation$cast_key))>=minimum,
    !any(validation$cast_key %in% sets$p50$cast_key),
    all(sets$p10$cast_key %in% sets$p20$cast_key),all(sets$p20$cast_key %in% sets$p50$cast_key))
  region_levels<-levels(full$Region)
  if(is.null(region_levels))region_levels<-sort(unique(as.character(full$Region)))
  scales<-function(d) {
    x<-meta[meta$cast_key %in% d$cast_key,]
    s<-list(mld_mean=mean(x$MLD_temp_MIX),mld_sd=sd(x$MLD_temp_MIX),
      distance_mean=mean(x$log10_distance),distance_sd=sd(x$log10_distance))
    stopifnot(s$mld_sd>0,s$distance_sd>0);s
  }
  scaling<-list(primary=scales(sets$p10))
  # Block all casts belonging to the same nested Ocean/Region/date group.
  meta$date_block<-paste(meta$Ocean,meta$Region,meta$date,sep="\r")
  blocked_ids<-character()
  for(i in seq_along(region_levels)) {
    m<-meta[as.character(meta$Region)==region_levels[i],]
    groups<-split(as.character(m$cast_key),m$date_block)
    target<-ceiling(.1*nrow(m)); chosen<-character();n<-0L
    set.seed(20261016L+i); groups<-groups[sample.int(length(groups))]
    while(n<target) {
      sizes<-lengths(groups)
      eligible<-which(nrow(m)-n-sizes>=max(2L,ceiling(.1*nrow(m))))
      if(!length(eligible))stop("Cannot create disjoint date holdout for Region ",region_levels[i])
      j<-eligible[which.min(abs(n+sizes[eligible]-target))]
      chosen<-c(chosen,groups[[j]]);n<-length(chosen);groups<-groups[-j]
    }
    blocked_ids<-c(blocked_ids,chosen)
  }
  blocked_validation<-full[full$cast_key %in% blocked_ids,,drop=FALSE]
  pool<-meta[!meta$date_block %in% meta$date_block[meta$cast_key %in% blocked_ids],]
  sizes<-as.numeric(table(factor(meta$Region,levels=region_levels)))
  capacities<-as.numeric(table(factor(pool$Region,levels=region_levels)))
  counts<-pmin(5L,capacities)
  while(sum(counts)<minimum) {
    gap<-.1*sizes-counts;gap[counts>=capacities]<- -Inf
    j<-which.max(gap);if(!is.finite(gap[j]))stop("Not enough date-independent training casts")
    counts[j]<-counts[j]+1L
  }
  ids<-unlist(lapply(seq_along(region_levels),function(i) {
    z<-as.character(pool$cast_key[as.character(pool$Region)==region_levels[i]])
    set.seed(20261017L+i);z[sample.int(length(z),counts[i])]
  }),use.names=FALSE)
  blocked_train<-full[full$cast_key %in% ids,,drop=FALSE]
  stopifnot(length(unique(blocked_ids))>=minimum,length(unique(ids))>=minimum,
    !any(ids %in% blocked_ids),!any(meta$date_block[meta$cast_key %in% ids] %in%
      meta$date_block[meta$cast_key %in% blocked_ids]))
  scaling$date_blocked<-scales(blocked_train)
  saveRDS(scaling,file.path(experiment,"scaling.rds"))
  saveRDS(region_levels,file.path(experiment,"region_levels.rds"))
  datasets<-list()
  for(s in names(sets))datasets[[paste0("primary_",s,"_d150")]]<-apply_features(sets[[s]],scaling$primary,region_levels)
  datasets$primary_validation_d150<-apply_features(validation,scaling$primary,region_levels)
  datasets$primary_p10_d250<-apply_features(full250[full250$cast_key %in% sets$p10$cast_key,,drop=FALSE],scaling$primary,region_levels)
  datasets$date_blocked_p10_d150<-apply_features(blocked_train,scaling$date_blocked,region_levels)
  datasets$date_blocked_validation_d150<-apply_features(blocked_validation,scaling$date_blocked,region_levels)
  for(name in names(datasets))saveRDS(datasets[[name]],file.path(experiment,"data",paste0(name,".rds")))
  # All plots use fixed real cast predictor pairs selected from primary training.
  reference<-unique(datasets$primary_p10_d150[c("cast_key","Region","Ocean","date","enso","log10_distance","MLD_temp_MIX")])
  representative<-do.call(rbind,lapply(region_levels,function(region) {
    x<-reference[as.character(reference$Region)==region,]
    x[which.min(abs(x$log10_distance-median(x$log10_distance))/scaling$primary$distance_sd+
      abs(x$MLD_temp_MIX-median(x$MLD_temp_MIX))/scaling$primary$mld_sd),,drop=FALSE]
  }))
  saveRDS(representative,file.path(experiment,"profile_reference.rds"))
  write.csv(representative,file.path(experiment,"profile_reference.csv"),row.names=FALSE)
  split_counts<-do.call(rbind,lapply(names(datasets),function(name) {
    d<-datasets[[name]];data.frame(dataset=name,casts=length(unique(d$cast_key)),
      fraction_of_all_casts=length(unique(d$cast_key))/total,observations=nrow(d),minimum_depth=min(d$depth),maximum_depth=max(d$depth))
  }))
  write.csv(split_counts,file.path(experiment,"split_counts.csv"),row.names=FALSE)
  write.csv(meta,file.path(experiment,"full_cast_metadata.csv"),row.names=FALSE)
  write.csv(do.call(rbind,lapply(names(datasets),function(n) {
    x<-unique(datasets[[n]][fields[fields!="depth" & fields!="chl_norm"]]);x$dataset<-n;x
  })),file.path(experiment,"split_membership.csv"),row.names=FALSE)
  # AR uses adjacent observed rows, not a continuous correlation length in metres.
  spacing<-do.call(rbind,lapply(split(full,droplevels(factor(full$cast_key))),function(x) {
    z<-diff(sort(x$depth));data.frame(cast_key=x$cast_key[1],minimum_step=if(length(z))min(z)else NA_real_,maximum_step=if(length(z))max(z)else NA_real_,duplicate_depth=any(z==0))
  }))
  if(any(spacing$duplicate_depth))stop("Duplicate depth within cast: cannot fit requested AR model unchanged")
  write.csv(spacing,file.path(experiment,"depth_spacing.csv"),row.names=FALSE)
  dir.create(file.path(experiment,"preflight"),showWarnings=FALSE)
  signatures<-character();counts_out<-list()
  for(i in seq_len(nrow(cases))) {
    c<-cases[i,];d<-datasets[[paste0(c$split,"_",c$subset,"_d",c$depth_max)]]
    spec<-campaign_spec(c,region_levels)
    prior_only<-if(c$mode=="prior")"only" else "no"
    sd<-brms::make_standata(spec$formula,data=d,prior=spec$priors,sample_prior=prior_only)
    code<-brms::make_stancode(spec$formula,data=d,prior=spec$priors,threads=thread_specification(c),sample_prior=prior_only)
    prefix<-file.path(experiment,"preflight",sprintf("case_%03d",c$case_id))
    writeLines(code,paste0(prefix,".stan"))
    capture.output(spec$formula,file=paste0(prefix,"_formula.txt"))
    write.csv(as.data.frame(brms::validate_prior(spec$priors,spec$formula,data=d,sample_prior=prior_only)),paste0(prefix,"_priors.csv"),row.names=FALSE)
    zs<-sd[grepl("^Zs",names(sd))]
    counts_out[[i]]<-data.frame(case_id=c$case_id,observations=nrow(d),fixed_columns=sd$K,
      penalized_spline_columns=sum(vapply(zs,function(x)if(is.matrix(x))ncol(x)else 0L,integer(1))),
      stan_data_mb=as.numeric(object.size(sd))/1024^2)
    cat("Validated targeted case",c$case_id,c$label,"\n");rm(sd);invisible(gc())
  }
  write.csv(do.call(rbind,counts_out),file.path(experiment,"basis_counts.csv"),row.names=FALSE)
  files<-list.files(file.path(experiment,"data"),pattern="rds$",full.names=TRUE)
  write.csv(data.frame(file=basename(files),md5=unname(tools::md5sum(files))),file.path(experiment,"data_checksums.csv"),row.names=FALSE)
  writeLines("63 cases constructed; >=10% whole casts training plus >=10% independent validation; no response preprocessing changes.",file.path(experiment,"PREPARATION_COMPLETE.txt"))
  invisible(split_counts)
}
if(sys.nframe()==0L) {
  experiment<-Sys.getenv("EXPERIMENT_DIR")
  tryCatch(prepare_campaign(experiment),error=function(e) {
    writeLines(conditionMessage(e),file.path(experiment,"PREPARATION_FAILED.txt"))
    stop(e)
  })
}
