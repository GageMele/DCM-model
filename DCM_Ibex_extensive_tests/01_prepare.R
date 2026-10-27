#!/usr/bin/env Rscript
prepare_campaign <- function(experiment) {
  source(file.path(experiment, "common.R"), local=TRUE)
  source(file.path(experiment, "campaign.R"), local=TRUE)
  if (file.exists(file.path(experiment,"PREPARATION_COMPLETE.txt"))) stop("Already prepared")
  cases <- read.delim(file.path(experiment,"cases.tsv"), stringsAsFactors=FALSE)
  stopifnot(isTRUE(all.equal(cases, campaign_cases(), check.attributes=FALSE)))
  full_bundle <- readRDS(file.path(experiment,"full_input_250m.rds"))
  full <- validate_bundle(full_bundle, full=TRUE)
  full <- full[full$depth>=10 & full$depth<=150,,drop=FALSE]
  sets <- lapply(c("p10","p20","p50"), function(s)
    validate_bundle(readRDS(file.path(experiment,"data",paste0(s,".rds")))))
  names(sets) <- c("p10","p20","p50")
  for (s in names(sets)) {
    d <- sets[[s]]
    stopifnot(all(d$depth>=10 & d$depth<=150),
      all(unique(d$cast_key) %in% full$cast_key))
    # Verify source rows/values despite the source subset's shuffled order.
    cols <- c("cast_key","chl_norm","Region","Ocean","date","enso",
              "depth","log10_distance","MLD_temp_MIX")
    canonical <- function(x) {
      x <- as.data.frame(lapply(x[cols], as.character), stringsAsFactors=FALSE)
      x <- x[do.call(order, c(unname(x), list(method="radix"))),,drop=FALSE]
      rownames(x) <- NULL; x
    }
    stopifnot(identical(canonical(d), canonical(full[full$cast_key %in% d$cast_key,])))
  }
  stopifnot(all(sets$p10$cast_key %in% sets$p20$cast_key),
            all(sets$p20$cast_key %in% sets$p50$cast_key))
  metadata <- unique(full[c("cast_key","Region","log10_distance","MLD_temp_MIX")])
  if (anyDuplicated(metadata$cast_key)) stop("Cast predictors are not constant")
  metadata <- metadata[order(metadata$cast_key, method="radix"),]
  regions <- sort(unique(as.character(metadata$Region)), method="radix")
  available <- metadata[!metadata$cast_key %in% sets$p50$cast_key,]
  quota <- .1 * as.numeric(table(factor(metadata$Region, levels=regions)))
  capacities <- as.numeric(table(factor(available$Region, levels=regions)))
  if (any(capacities<1)) stop("No unused validation cast in at least one Region")
  counts <- pmin(capacities, pmax(1L, floor(quota)))
  target <- round(.1*nrow(metadata))
  if (sum(counts)>target) stop("Validation minimum per Region exceeds target")
  while (sum(counts)<target) {
    deficit <- quota-counts; deficit[counts>=capacities] <- -Inf
    i <- which.max(deficit)
    if (!is.finite(deficit[i])) stop("Not enough unused validation casts")
    counts[i] <- counts[i]+1L
  }
  ids <- unlist(lapply(seq_along(regions), function(i) {
    pool <- as.character(available$cast_key[as.character(available$Region)==regions[i]])
    set.seed(20261015L+i); pool[sample.int(length(pool), counts[i])]
  }), use.names=FALSE)
  validation <- full[full$cast_key %in% ids,,drop=FALSE]
  stopifnot(!any(validation$cast_key %in% sets$p50$cast_key),
            length(unique(ids))==target)
  # Freeze these training-only scales across all data sizes. Response unchanged.
  meta10 <- metadata[metadata$cast_key %in% sets$p10$cast_key,]
  scaling <- list(mld_mean=mean(meta10$MLD_temp_MIX), mld_sd=sd(meta10$MLD_temp_MIX),
    distance_mean=mean(meta10$log10_distance), distance_sd=sd(meta10$log10_distance))
  stopifnot(scaling$mld_sd>0, scaling$distance_sd>0)
  saveRDS(scaling, file.path(experiment,"scaling.rds"))
  # Original bundle files remain byte-for-byte unchanged. Derived columns are separate.
  for (s in names(sets)) {
    sets[[s]] <- apply_scaling(sets[[s]], scaling)
    saveRDS(sets[[s]], file.path(experiment,"data",paste0(s,"_model.rds")))
  }
  validation <- apply_scaling(validation, scaling)
  # Match fixed-effect levels; unseen dates are allowed by prediction routines.
  for (s in c("Region","Ocean","enso")) {
    validation[[s]] <- factor(as.character(validation[[s]]), levels=levels(sets$p10[[s]]))
    if (anyNA(validation[[s]])) stop("Unrepresented validation fixed grouping level: ",s)
  }
  saveRDS(validation, file.path(experiment,"data","validation_model.rds"))
  write.csv(metadata[metadata$cast_key %in% ids,], file.path(experiment,"validation_casts.csv"),row.names=FALSE)
  write.csv(do.call(rbind,lapply(c(sets,list(validation=validation)), support_table)),
            file.path(experiment,"support_counts.csv"),row.names=TRUE)
  split_counts <- do.call(rbind,lapply(names(c(sets,list(validation=validation))),function(s) {
    x <- c(sets,list(validation=validation))[[s]]
    data.frame(dataset=s,casts=length(unique(x$cast_key)),observations=nrow(x))
  }))
  write.csv(split_counts,file.path(experiment,"split_counts.csv"),row.names=FALSE)
  dir.create(file.path(experiment,"preflight"),showWarnings=FALSE)
  representatives <- cases[!duplicated(case_design_key(cases)),,drop=FALSE]
  counts_out <- list()
  for (j in seq_len(nrow(representatives))) {
    c <- representatives[j,]; d <- sets[[c$subset]]; spec <- campaign_spec(c)
    cat("Preflight",j,"/",nrow(representatives),c$label,"\n")
    prior_only <- if(c$mode=="prior") "only" else "no"
    sd <- brms::make_standata(spec$formula,data=d,prior=spec$priors,sample_prior=prior_only)
    code <- brms::make_stancode(spec$formula,data=d,prior=spec$priors,
      threads=brms::threading(c$threads),sample_prior=prior_only)
    prefix <- file.path(experiment,"preflight",sprintf("case_%03d",c$case_id))
    writeLines(code,paste0(prefix,".stan"))
    resolved <- brms::validate_prior(spec$priors,spec$formula,data=d,sample_prior=prior_only)
    write.csv(as.data.frame(resolved),paste0(prefix,"_priors.csv"),row.names=FALSE)
    capture.output(spec$formula,file=paste0(prefix,"_formula.txt"))
    matrices <- sd[grepl("^Zs",names(sd))]
    counts_out[[j]] <- data.frame(case_id=c$case_id,key=case_design_key(c),
      observations=nrow(d),fixed_columns=sd$K,
      penalized_spline_columns=sum(vapply(matrices,function(x) if(is.matrix(x)) ncol(x) else 0L,integer(1))))
    rm(sd); invisible(gc())
  }
  write.csv(do.call(rbind,counts_out),file.path(experiment,"basis_counts.csv"),row.names=FALSE)
  paths <- c(file.path(experiment,"data",paste0(c("p10","p20","p50"),".rds")),
             file.path(experiment,"data","validation_model.rds"))
  write.csv(data.frame(file=basename(paths),md5=unname(tools::md5sum(paths))),
            file.path(experiment,"data_checksums.csv"),row.names=FALSE)
  capture.output(sessionInfo(),file=file.path(experiment,"preflight_session.txt"))
  writeLines("All cases passed Stan data/code and prior construction; held-out casts disjoint.",
             file.path(experiment,"PREPARATION_COMPLETE.txt"))
  invisible(split_counts)
}
if (sys.nframe()==0L) prepare_campaign(Sys.getenv("EXPERIMENT_DIR"))
