# Verify data identities, new model matrices and source controls before any fit.
experiment <- Sys.getenv("EXPERIMENT_DIR")
source(file.path(experiment,"common.R"));source(file.path(experiment,"campaign.R"))
source(file.path(experiment,"data_setup.R"))
tryCatch({
  if(file.exists(file.path(experiment,"PREPARATION_COMPLETE.txt")))stop("Already prepared")
  origin <- readLines(file.path(experiment,"source_experiment.txt"),n=1L)
  prepare_data(experiment,origin)
  cases <- read.delim(file.path(experiment,"cases.tsv"),stringsAsFactors=FALSE)
  stopifnot(isTRUE(all.equal(cases,campaign_cases(),check.attributes=FALSE)))
  for(p in c("brms","mgcv","cmdstanr","posterior"))
    if(!requireNamespace(p,quietly=TRUE))stop("Missing package: ",p)
  controls <- read.delim(file.path(experiment,"source_controls.tsv"),stringsAsFactors=FALSE)
  source_cases <- read.delim(file.path(origin,"cases.tsv"),stringsAsFactors=FALSE)
  settings <- c("equation","subset","split","prior","mode","variance","cast","family",
    "autocor","depth_max","chains","iter","warmup","adapt_delta","treedepth",
    "threads","k_depth","k_mld","k_distance","seed")
  for(i in which(cases$execution=="reuse")) {
    c <- cases[i,];old <- source_cases[source_cases$case_id==c$source_case,,drop=FALSE]
    stopifnot(nrow(old)==1L)
    for(n in settings)if(!isTRUE(all.equal(c[[n]],old[[n]],check.attributes=FALSE)))
      stop("Source control differs for ",c$case_id," ",n)
    stopifnot(sum(controls$case_id==c$case_id)==1L)
  }
  dir.create(file.path(experiment,"preflight"),showWarnings=FALSE)
  basis <- list();levels <- readRDS(file.path(experiment,"region_levels.rds"))
  for(i in seq_len(nrow(cases))) {
    c <- cases[i,];d <- readRDS(data_file(experiment,c));spec <- campaign_spec(c,levels)
    sd <- brms::make_standata(spec$formula,data=d,prior=spec$priors,
      threads=thread_specification(c),sample_prior=if(c$mode=="prior")"only" else "no")
    code <- brms::make_stancode(spec$formula,data=d,prior=spec$priors,
      threads=thread_specification(c),sample_prior=if(c$mode=="prior")"only" else "no")
    writeLines(code,file.path(experiment,"preflight",sprintf("case_%02d.stan",c$case_id)))
    write.csv(as.data.frame(brms::get_prior(spec$formula,data=d)),
      file.path(experiment,"preflight",sprintf("case_%02d_priors.csv",c$case_id)),row.names=FALSE)
    X <- sd$X
    if(is.matrix(X) && qr(X)$rank<ncol(X))stop("Rank-deficient fixed design in case ",c$case_id)
    basis[[i]] <- data.frame(case_id=c$case_id,fixed_columns=if(is.matrix(X))ncol(X)else NA_integer_,
      penalized_columns=sum(vapply(sd[grepl("^Zs",names(sd))],function(x)
        if(is.matrix(x))ncol(x)else 0L,integer(1))))
    cat("Validated onset case",c$case_id,c$label,"\n");rm(sd);invisible(gc())
  }
  write.csv(do.call(rbind,basis),file.path(experiment,"basis_counts.csv"),row.names=FALSE)
  files <- list.files(file.path(experiment,"data"),full.names=TRUE)
  write.csv(data.frame(file=basename(files),md5=unname(tools::md5sum(files))),
    file.path(experiment,"data_checksums.csv"),row.names=FALSE)
  writeLines("34 formula/data combinations validated; 28 new fits, 6 controls reused; no response changes",
             file.path(experiment,"PREPARATION_COMPLETE.txt"))
},error=function(e){writeLines(conditionMessage(e),file.path(experiment,"PREPARATION_FAILED.txt"));stop(e)})
