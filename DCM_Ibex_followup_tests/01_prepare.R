# Validate every model on Ibex before any fitting job is released.
experiment<-Sys.getenv("EXPERIMENT_DIR")
source(file.path(experiment,"common.R"));source(file.path(experiment,"campaign.R"))
source(file.path(experiment,"data_setup.R"))
tryCatch({
  if(file.exists(file.path(experiment,"PREPARATION_COMPLETE.txt")))stop("Already prepared")
  counts<-prepare_data(experiment,readLines(file.path(experiment,"source_experiment.txt"),n=1L))
  cases<-read.delim(file.path(experiment,"cases.tsv"),stringsAsFactors=FALSE)
  stopifnot(isTRUE(all.equal(cases,campaign_cases(),check.attributes=FALSE)))
  for(p in c("brms","mgcv","cmdstanr","posterior"))if(!requireNamespace(p,quietly=TRUE))stop("Missing package: ",p)
  output<-list();dir.create(file.path(experiment,"preflight"),showWarnings=FALSE)
  for(i in seq_len(nrow(cases))) {
    c<-cases[i,];d<-readRDS(data_file(experiment,c));spec<-campaign_spec(c,readRDS(file.path(experiment,"region_levels.rds")))
    stopifnot(length(unique(d$cast_key))>=108L)
    sd<-brms::make_standata(spec$formula,data=d,prior=spec$priors,threads=thread_specification(c),
      sample_prior=if(c$mode=="prior")"only"else"no")
    writeLines(brms::make_stancode(spec$formula,data=d,prior=spec$priors,threads=thread_specification(c),
      sample_prior=if(c$mode=="prior")"only"else"no"),file.path(experiment,"preflight",sprintf("case_%02d.stan",c$case_id)))
    write.csv(as.data.frame(brms::get_prior(spec$formula,data=d)),file.path(experiment,"preflight",sprintf("case_%02d_priors.csv",c$case_id)),row.names=FALSE)
    write.csv(as.data.frame(spec$priors),file.path(experiment,"preflight",sprintf("case_%02d_requested_priors.csv",c$case_id)),row.names=FALSE)
    zs<-sd[grepl("^Zs",names(sd))]
    output[[i]]<-data.frame(case_id=c$case_id,observations=nrow(d),fixed_columns=sd$K,
      penalized_spline_columns=sum(vapply(zs,function(x)if(is.matrix(x))ncol(x)else 0L,integer(1))))
    cat("Validated follow-up case",c$case_id,c$label,"\n");rm(sd);invisible(gc())
  }
  write.csv(do.call(rbind,output),file.path(experiment,"basis_counts.csv"),row.names=FALSE)
  files<-list.files(file.path(experiment,"data"),full.names=TRUE)
  write.csv(data.frame(file=basename(files),md5=unname(tools::md5sum(files))),file.path(experiment,"data_checksums.csv"),row.names=FALSE)
  writeLines("All 34 models constructed; training >=10%; old and reserve validation disjoint; raw responses preserved.",file.path(experiment,"PREPARATION_COMPLETE.txt"))
},error=function(e){writeLines(conditionMessage(e),file.path(experiment,"PREPARATION_FAILED.txt"));stop(e)})
