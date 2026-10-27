root <- "dcm_snapshot_20261016_090217"
files <- list.files(root, pattern = "^run_config\\.rds$", recursive = TRUE, full.names = TRUE)
configs <- setNames(lapply(files, readRDS), gsub("\\\\", "/", dirname(files)))
jsonlite::write_json(configs, "snapshot_configs.json", auto_unbox = TRUE, dataframe = "rows", null = "null", na = "null", digits = NA)
cat("Configs read:", length(configs), "\n")
for (p in c("ggplot2", "openxlsx", "writexl", "jsonlite", "patchwork", "ragg"))
  cat(p, requireNamespace(p, quietly = TRUE), "\n")
