# Freeze new sets without using chlorophyll values for selection.
prepare_data<-function(experiment,source_experiment) {
  dir.create(file.path(experiment,"data"),showWarnings=FALSE,recursive=TRUE)
  required<-c("scaling.rds","region_levels.rds","profile_reference.rds","full_input_250m.rds",
    "full_cast_metadata.csv","profile_reference.csv","depth_spacing.csv")
  data_names<-c("primary_p10_d150.rds","primary_p20_d150.rds","primary_p50_d150.rds",
    "primary_validation_d150.rds","date_blocked_p10_d150.rds","date_blocked_validation_d150.rds")
  src<-c(file.path(source_experiment,required),file.path(source_experiment,"data",data_names))
  dst<-c(file.path(experiment,required),file.path(experiment,"data",data_names))
  stopifnot(all(file.exists(src)))
  if(!all(file.copy(src,dst,overwrite=FALSE)))stop("Frozen data copy failed or destination already exists")
  stopifnot(identical(unname(tools::md5sum(src)),unname(tools::md5sum(dst))))
  write.csv(data.frame(source=src,file=substring(dst,nchar(experiment)+2L),md5=unname(tools::md5sum(dst))),file.path(experiment,"source_checksums.csv"),row.names=FALSE)
  full<-validate_bundle(readRDS(file.path(experiment,"full_input_250m.rds")),full=TRUE)
  full<-full[full$depth>=10 & full$depth<=150,,drop=FALSE]
  read<-function(n)readRDS(file.path(experiment,"data",paste0(n,".rds")))
  p10<-read("primary_p10_d150");p20<-read("primary_p20_d150");p50<-read("primary_p50_d150")
  blocked<-read("date_blocked_p10_d150");bv<-read("date_blocked_validation_d150");v<-read("primary_validation_d150")
  stopifnot(all(p10$cast_key %in% p20$cast_key),all(p20$cast_key %in% p50$cast_key))
  alternate<-p20[!p20$cast_key %in% p10$cast_key,,drop=FALSE]
  saveRDS(alternate,file.path(experiment,"data","primary_p10_alt_d150.rds"))
  meta<-unique(full[c("cast_key","Region","Ocean","date","enso","log10_distance","MLD_temp_MIX")])
  meta<-meta[order(meta$cast_key,method="radix"),,drop=FALSE]
  stopifnot(!anyDuplicated(meta$cast_key));minimum<-ceiling(.1*nrow(meta))
  excluded<-unique(c(p50$cast_key,blocked$cast_key,v$cast_key,bv$cast_key))
  pool<-meta[!meta$cast_key %in% excluded,,drop=FALSE]
  levels<-readRDS(file.path(experiment,"region_levels.rds"))
  capacities<-as.integer(table(factor(pool$Region,levels=levels)))
  full_counts<-as.integer(table(factor(meta$Region,levels=levels)))
  counts<-pmin(5L,capacities)
  if(any(capacities<2L)||sum(capacities)<minimum)stop("Insufficient unused casts for a new >=10% reserve holdout")
  while(sum(counts)<minimum) {
    gap<-minimum*full_counts/sum(full_counts)-counts;gap[counts>=capacities]<- -Inf
    j<-which.max(gap);if(!is.finite(gap[j]))stop("Cannot fill reserve holdout")
    counts[j]<-counts[j]+1L
  }
  chosen<-unlist(lapply(seq_along(levels),function(i){
    candidates<-as.character(pool$cast_key[as.character(pool$Region)==levels[i]])
    set.seed(20261018L+i);candidates[sample.int(length(candidates),counts[i])]
  }),use.names=FALSE)
  reserve<-full[full$cast_key %in% chosen,,drop=FALSE]
  scaling<-readRDS(file.path(experiment,"scaling.rds"))
  for(split in c("primary","date_blocked"))saveRDS(apply_features(reserve,scaling[[split]],levels),
    file.path(experiment,"data",paste0(split,"_reserve_d150.rds")))
  key<-function(z)paste(z$Ocean,z$Region,z$date,sep="\r")
  stopifnot(length(unique(alternate$cast_key))>=minimum,!any(alternate$cast_key %in% p10$cast_key),
    length(unique(reserve$cast_key))>=minimum,!any(reserve$cast_key %in% excluded),
    !any(v$cast_key %in% p50$cast_key),!any(key(blocked) %in% key(bv)))
  sets<-list(primary_p10_d150=p10,primary_p10_alt_d150=alternate,primary_p20_d150=p20,
    primary_p50_d150=p50,primary_validation_d150=v,date_blocked_p10_d150=blocked,
    date_blocked_validation_d150=bv,primary_reserve_d150=read("primary_reserve_d150"),
    date_blocked_reserve_d150=read("date_blocked_reserve_d150"))
  counts_out<-do.call(rbind,lapply(names(sets),function(n){d<-sets[[n]]
    data.frame(dataset=n,casts=length(unique(d$cast_key)),fraction_of_all_casts=length(unique(d$cast_key))/nrow(meta),observations=nrow(d))}))
  stopifnot(all(counts_out$casts>=minimum))
  write.csv(counts_out,file.path(experiment,"split_counts.csv"),row.names=FALSE)
  write.csv(do.call(rbind,lapply(names(sets),function(n){
    d<-unique(sets[[n]][c("cast_key","Region","Ocean","date","enso","log10_distance","MLD_temp_MIX")]);d$dataset<-n;d
  })),file.path(experiment,"split_membership.csv"),row.names=FALSE)
  writeLines(c("Reserve casts excluded from recorded primary50% and date-blocked training and previous holdouts.",
    "Historical full-data fits cannot be given independent reserve-test claims.",
    "Reserve selection uses cast predictors and fixed RNG seeds only; outcomes are never inspected.",
    "This follow-up uses reserve results for development; keep a further untouched test for final claims."),file.path(experiment,"RESERVE_INTERPRETATION.txt"))
  invisible(counts_out)
}
