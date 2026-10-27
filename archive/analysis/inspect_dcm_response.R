root <- "validation_extensive150/actual/data"
files <- file.path(root, c("p10_model.rds", "p20_model.rds", "p50_model.rds", "validation_model.rds"))
for (file in files[file.exists(files)]) {
  d <- readRDS(file)
  if (is.list(d) && !is.data.frame(d)) d <- d$model_data
  cat("FILE", basename(file), "MD5", tools::md5sum(file), "\n")
  cat("ROWS", nrow(d), "CASTS", length(unique(d$cast_key)), "DEPTH", range(d$depth), "\n")
  print(quantile(d$chl_norm, c(0,.001,.01,.05,.5,.95,.99,.999,1)))
  cat("mean", mean(d$chl_norm), "SD", sd(d$chl_norm), "\n")
  print(aggregate(chl_norm ~ Region, d, function(x) c(mean=mean(x),sd=sd(x),max=max(x))))
}
