# Prepare immutable datasets and validate equation designs on a compute node.
experiment<-Sys.getenv("EXPERIMENT_DIR")
source(file.path(experiment,"common.R"));source(file.path(experiment,"models.R"))
source(file.path(experiment,"data.R"))
tryCatch({
  if(file.exists(file.path(experiment,"PREPARATION_COMPLETE.txt")))quit(status=0)
  for(p in c("brms","mgcv","cmdstanr","posterior","ggplot2","bayesplot"))
    if(!requireNamespace(p,quietly=TRUE))stop("Required package missing: ",p)
  origin<-readLines(file.path(experiment,"source_experiment.txt"),n=1L)
  info<-prepare_data(experiment,origin)
  cases<-read.delim(file.path(experiment,"cases.tsv"),stringsAsFactors=FALSE)
  dir.create(file.path(experiment,"preflight"),showWarnings=FALSE)
  # Use actual 20% training casts to check each distinct formula/prior design;
  # support is additionally checked on EVERY case's actual training dataset.
  template_case<-cases[cases$equation=="original_z" & cases$subset=="p20" &
    cases$selection=="frozen" & cases$split=="primary" & cases$mode=="posterior",][1,,drop=FALSE]
  template<-readRDS(data_file(experiment,template_case));tm<-cast_metadata(template)
  template_keys<-unlist(lapply(split(tm,tm$Region),function(m){
    m<-m[order(m$MLD_temp_MIX),,drop=FALSE]
    m$cast_key[unique(as.integer(round(seq(1,nrow(m),length.out=min(10,nrow(m))))))]
  }),use.names=FALSE)
  template<-template[template$cast_key %in% template_keys,,drop=FALSE]
  cache<-new.env(parent=emptyenv());prepared<-list()
  for(i in seq_len(nrow(cases))) {
    c<-cases[i,];reason<-NULL
    record<-tryCatch({
      d<-readRDS(data_file(experiment,c));support<-support_table(d)
      if(any(support$depth<c$k_depth))stop("Insufficient unique depths for requested k")
      if(startsWith(c$equation,"tensor")&&any(support$mld<c$k_mld))
        stop("Insufficient regional MLD support; requested k was preserved")
      spec<-campaign_spec(c,info$levels)
      signature<-paste(c$equation,c$k_depth,c$k_mld,c$k_profile,c$prior,c$mode,sep="|")
      if(!exists(signature,cache,inherits=FALSE)) {
        sd<-brms::make_standata(spec$formula,data=template,prior=spec$priors,
          threads=thread_specification(c),sample_prior=if(c$mode=="prior")"only" else "no")
        if(is.matrix(sd$X)&&qr(sd$X)$rank<ncol(sd$X))stop("Rank-deficient unpenalized design")
        code<-brms::make_stancode(spec$formula,data=template,prior=spec$priors,
          threads=thread_specification(c),sample_prior=if(c$mode=="prior")"only" else "no")
        assign(signature,code,cache);rm(sd);invisible(gc())
      }
      writeLines(get(signature,cache),file.path(experiment,"preflight",sprintf("case_%03d.stan",c$case_id)))
      write.csv(support,file.path(experiment,"preflight",sprintf("case_%03d_support.csv",c$case_id)),row.names=FALSE)
      data.frame(case_id=c$case_id,valid=TRUE,reason="",training_casts=length(unique(d$cast_key)),
        observations=nrow(d),formula=paste(capture.output(spec$formula),collapse=" "))
    },error=function(e)data.frame(case_id=c$case_id,valid=FALSE,reason=conditionMessage(e),
      training_casts=NA_integer_,observations=NA_integer_,formula=""))
    prepared[[i]]<-record
    cat("Preflight",c$case_id,c$label,if(record$valid)"PASS"else record$reason,"\n")
  }
  preflight<-do.call(rbind,prepared)
  write.csv(preflight,file.path(experiment,"preflight_cases.csv"),row.names=FALSE)
  if(!any(preflight$valid))stop("No valid cases; no model fitting will be submitted")
  files<-list.files(file.path(experiment,"data"),full.names=TRUE)
  write.csv(data.frame(file=basename(files),md5=unname(tools::md5sum(files))),
    file.path(experiment,"data_checksums.csv"),row.names=FALSE)
  writeLines(sprintf("Prepared %d cases; %d valid designs; %d full eligible casts; failures retained explicitly",
    nrow(cases),sum(preflight$valid),info$casts),file.path(experiment,"PREPARATION_COMPLETE.txt"))
},error=function(e){writeLines(conditionMessage(e),file.path(experiment,"PREPARATION_FAILED.txt"));stop(e)})
