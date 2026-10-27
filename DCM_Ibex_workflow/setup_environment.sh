#!/bin/bash -l
# Run inside an interactive COMPUTE allocation. Do not compile on a login node.
# This is a shell helper so the scientific workflow contains exactly three R files.
set -euo pipefail
: "${SLURM_JOB_ID:?Start an srun interactive compute allocation first.}"
cd "${DCM_ROOT:-/ibex/user/$USER/DCM_MODEL}"
source ./ibex_environment.sh
mkdir -p software/cmdstan logs data/raw data/processed runs
export CMDSTAN_VERSION="${CMDSTAN_VERSION:-2.37.0}"
# Package installation uses the compute allocation's CPUs, capped at four.
Rscript --vanilla - <<'RSCRIPT'
lib <- Sys.getenv("R_LIBS_USER")
.libPaths(c(lib, .libPaths()))
repos <- c(STAN = "https://stan-dev.r-universe.dev", CRAN = "https://cloud.r-project.org")
packages <- c("brms", "cmdstanr", "mgcv", "posterior", "bayesplot", "ggplot2")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
cores <- min(4L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))
if (length(missing)) install.packages(missing, repos = repos, Ncpus = cores)
for (p in packages) if (!requireNamespace(p, quietly = TRUE)) stop("Installation failed: ", p)
cmdstanr::check_cmdstan_toolchain(quiet = FALSE)
root <- Sys.getenv("DCM_ROOT")
version <- Sys.getenv("CMDSTAN_VERSION")
path <- file.path(root, "software", "cmdstan", paste0("cmdstan-", version))
if (!file.exists(file.path(path, "bin", "stanc"))) {
  cmdstanr::install_cmdstan(version = version, dir = dirname(path), cores = cores)
}
cmdstanr::set_cmdstan_path(path)
writeLines(cmdstanr::cmdstan_path(), file.path(root, "software", "cmdstan_path.txt"))
info <- data.frame(package = packages,
                   version = vapply(packages, function(p) as.character(packageVersion(p)), character(1)))
write.csv(info, file.path(root, "software", "package_versions.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(root, "software", "setup_session_info.txt"))
cat("CmdStan:", cmdstanr::cmdstan_path(), "\n")
RSCRIPT
module list
