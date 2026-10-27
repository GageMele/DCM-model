# Read configuration metadata only; no model fitting or package loading.
root <- 'dcm_snapshot_20261018_105731'
files <- list.files(root, pattern='^run_config\\.rds$', recursive=TRUE, full.names=TRUE)
rows <- lapply(files, function(path) {
  config <- readRDS(path)
  result <- list(run_path=gsub('\\\\','/',substring(dirname(path),nchar(root)+2L)))
  if (is.data.frame(config$case) && nrow(config$case)==1L)
    for (key in names(config$case)) result[[key]] <- as.character(config$case[[key]])
  for (key in names(config)) {
    value <- config[[key]]
    if (is.atomic(value) && length(value)==1L) result[[key]] <- as.character(value)
    if (key=='k') result$configured_k <- paste(value,collapse=',')
  }
  as.data.frame(result,stringsAsFactors=FALSE)
})
columns <- unique(unlist(lapply(rows,names)))
rows <- lapply(rows,function(x) {for (key in setdiff(columns,names(x)))x[[key]]<-NA_character_;x[columns]})
write.csv(do.call(rbind,rows),'DCM_review_20261018/config_metadata.csv',row.names=FALSE)
cat('Configuration records:',length(rows),'\n')
