# Construct the unchanged three-dimensional tensor on the smallest actual subset.
d <- readRDS("validation_fast_tests/actual150_radix/data/p05.rds")$model_data
smooths <- mgcv::smoothCon(mgcv::t2(depth, log10_distance, MLD_temp_MIX,
  by = Region, bs = c("cr", "cr", "cr"), k = c(15, 5, 5)),
  data = d, absorb.cons = TRUE)
stopifnot(length(smooths) == nlevels(d$Region))
report <- do.call(rbind, lapply(smooths, function(smooth) {
  stopifnot(nrow(smooth$X) == nrow(d), all(is.finite(smooth$X)))
  data.frame(region = smooth$by.level, basis_columns = ncol(smooth$X),
             null_space_dimension = smooth$null.space.dim,
             penalty_blocks = length(smooth$S))
}))
write.csv(report, "validation_fast_tests/basis_report.csv", row.names = FALSE)
print(report)
cat("Unchanged regional tensor bases constructed successfully on the actual 5% subset.\n")
