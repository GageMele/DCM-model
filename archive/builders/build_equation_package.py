# Assemble a separate campaign from the validated shared fitting machinery.
from pathlib import Path

root = Path.cwd()
source = root / "DCM_Ibex_k_tests"
target = root / "DCM_Ibex_equation_tests"
target.mkdir(exist_ok=False)

for name in ("common.R", "02_fit_case.R", "03_check_case.R", "fit.slurm",
             "check.slurm", "prepare.slurm", "collect.slurm", "live_diagnostics.awk"):
    body = (source / name).read_text(encoding="utf-8")
    if name == "02_fit_case.R":
        body = body.replace("k <- case_k(case)", "k <- 15L")
        body = body.replace("check_support(d, k)", "check_additive_support(d)")
        body = body.replace(
            'spec <- load_original_spec(file.path(experiment, "02_fit_model.R"), k)',
            'source(file.path(experiment, "equation.R"))\n'
            '  spec <- additive_spec(file.path(experiment, "02_fit_model.R"))')
        body = body.replace(
            'source_md5 = unname(tools::md5sum(file.path(experiment, "02_fit_model.R"))),',
            'source_md5 = unname(tools::md5sum(file.path(experiment, "equation.R"))),\n'
            '    original_source_md5 = unname(tools::md5sum(file.path(experiment, "02_fit_model.R"))),')
        body = body.replace("using the frozen original formula", "using the requested additive formula")
    if name == "prepare.slurm":
        body = body.replace("tensor basis", "additive formula")
    if name == "fit.slurm":
        body = body.replace("Basis-size case", "Additive-equation case")
    (target / name).write_text(body.replace("\r\n", "\n"), encoding="utf-8", newline="\n")

# Extend the monitor by adding the new campaign; keep its existing behavior.
monitor = (root / "monitor_all_dcm.sh").read_text(encoding="utf-8")
monitor = monitor.replace("experiments/k150_*/manifest.tsv)",
    "experiments/k150_*/manifest.tsv experiments/additive150_*/manifest.tsv)")
monitor = monitor.replace("logs/dcm-k150-c[0-9]*.out)",
    "logs/dcm-k150-c[0-9]*.out logs/dcm-a150-c[0-9]*.out)")
monitor = monitor.replace(
    "for manifest in experiments/fast150_*/manifest.tsv experiments/k150_*/manifest.tsv; do",
    "for manifest in experiments/fast150_*/manifest.tsv experiments/k150_*/manifest.tsv experiments/additive150_*/manifest.tsv; do")
monitor = monitor.replace("both 150 m campaigns", "all 150 m campaigns")
monitor = monitor.replace("fast150 and k150", "fast150, k150 and additive150")
monitor = monitor.replace('diagnostics_awk="DCM_Ibex_k_tests/live_diagnostics.awk"',
    'diagnostics_awk="DCM_Ibex_equation_tests/live_diagnostics.awk"')
(target / "monitor_all_dcm_v2.sh").write_text(monitor, encoding="utf-8", newline="\n")
print("Created isolated additive-formula package.")
