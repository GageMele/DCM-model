# Copy previously frozen follow-up casts exactly; never choose casts from outcomes.
prepare_data <- function(experiment,source_experiment) {
  dir.create(file.path(experiment,"data"),recursive=TRUE,showWarnings=FALSE)
  required <- c("scaling.rds","region_levels.rds","profile_reference.rds",
    "full_cast_metadata.csv","profile_reference.csv","depth_spacing.csv","split_membership.csv")
  names <- c("primary_p10_d150","primary_p20_d150","primary_validation_d150",
    "primary_reserve_d150","date_blocked_p10_d150","date_blocked_validation_d150",
    "date_blocked_reserve_d150")
  stopifnot(all(file.exists(file.path(source_experiment,required))))
  stopifnot(all(file.copy(file.path(source_experiment,required),
                         file.path(experiment,required),overwrite=FALSE)))
  scaling <- readRDS(file.path(experiment,"scaling.rds"))
  levels <- readRDS(file.path(experiment,"region_levels.rds"))
  raw <- list();audits <- list()
  raw_cols <- c("cast_key","Region","Ocean","date","enso","depth",
                "log10_distance","MLD_temp_MIX","chl_norm")
  for(n in names) {
    path <- file.path(source_experiment,"data",paste0(n,".rds"))
    before <- readRDS(path);raw[[n]] <- before
    split <- if(startsWith(n,"date_blocked"))"date_blocked" else "primary"
    after <- apply_features(before,scaling[[split]],levels)
    # Row order, response, predictors and group membership must remain identical.
    stopifnot(identical(before[raw_cols],after[raw_cols]))
    saveRDS(after,file.path(experiment,"data",paste0(n,".rds")))
    audits[[n]] <- data.frame(dataset=n,source_file=path,
      source_md5=unname(tools::md5sum(path)),casts=length(unique(after$cast_key)),
      observations=nrow(after),minimum_mld=min(after$MLD_temp_MIX),
      source_raw_rows_verified=TRUE)
  }
  full <- read.csv(file.path(experiment,"full_cast_metadata.csv"),stringsAsFactors=FALSE)
  minimum <- ceiling(.1*length(unique(full$cast_key)))
  stopifnot(all(vapply(raw,function(d)length(unique(d$cast_key)),integer(1))>=minimum))
  for(split in c("primary","date_blocked")) {
    keys <- function(s)unique(raw[[paste0(split,"_",s,"_d150")]]$cast_key)
    training <- keys(if(split=="primary")"p20" else "p10")
    stopifnot(!length(intersect(training,keys("validation"))),
      !length(intersect(training,keys("reserve"))),
      !length(intersect(keys("validation"),keys("reserve"))))
  }
  date_key <- function(d)paste(d$Ocean,d$Region,d$date,sep="\r")
  stopifnot(!any(date_key(raw$date_blocked_p10_d150) %in%
                  date_key(raw$date_blocked_validation_d150)))
  write.csv(do.call(rbind,audits),file.path(experiment,"data_identity_audit.csv"),row.names=FALSE)
  invisible(do.call(rbind,audits))
}
