# Build the isolated baseline campaign from the previously validated machinery.
from pathlib import Path
import csv

root = Path.cwd()
source = root / "DCM_Ibex_equation_tests"
target = root / "DCM_Ibex_baseline_tests"
target.mkdir(exist_ok=True)

for name in ("common.R", "03_check_case.R", "fit.slurm", "check.slurm",
             "prepare.slurm", "collect.slurm", "live_diagnostics.awk"):
    body = (source / name).read_text(encoding="utf-8")
    body = body.replace("Additive-equation case", "Historical baseline case")
    body = body.replace("# Check every proposed additive formula on the actual frozen 10% subset.",
                        "# Prepare nested subsets and validate both prior choices on every dataset.")
    body = body.replace("01_validate_campaign.R", "01_prepare_baselines.R")
    if name == "check.slurm":
        body = body.replace("--mem=8G", "--mem=32G").replace("--time=01:00:00", "--time=02:00:00")
    if name == "03_check_case.R":
        body = body.replace("  tree_hits <- sum(post$treedepth_hits)",
            "  tree_hits <- sum(post$treedepth_hits)\n"
            "  tree_fraction <- tree_hits / sum(post$iterations)\n"
            '  posterior_rows <- sampler[sampler$phase == "sampling", ]\n'
            "  mean_leapfrogs <- mean(posterior_rows$n_leapfrog__)")
        body = body.replace("divergences = divergences, treedepth_hits = tree_hits, min_ebfmi = ebfmi,",
            "divergences = divergences, treedepth_hits = tree_hits,\n"
            "    treedepth_hit_fraction = tree_fraction, mean_leapfrogs = mean_leapfrogs, min_ebfmi = ebfmi,")
        body = body.replace('"Short runs and these checks do not establish scientific validity.",',
            '"This conservative screen does not establish scientific validity.",\n'
            '    "Tree saturation is an efficiency concern; inspect R-hat, ESS and predictive behavior.",')
        body = body.replace('"Full-data-derived priors reuse the response and require sensitivity review."',
            '"Custom priors are the historical constants; default priors are resolved by the installed brms version."')
    (target / name).write_text(body, encoding="utf-8", newline="\n")

fit = (source / "02_fit_case.R").read_text(encoding="utf-8")
fit = fit.replace('"equation.R"', '"baseline.R"')
fit = fit.replace('case <- read_case(experiment, Sys.getenv("CASE_ID"))',
    'if (!file.exists(file.path(experiment, "PREPARATION_COMPLETE.txt")))\n'
    '  stop("Preparation has not completed successfully.")\n'
    'case <- read_case(experiment, Sys.getenv("CASE_ID"))')
fit = fit.replace("check_equation_support(d)", "check_baseline_support(d)")
fit = fit.replace('spec <- load_equation_spec(case$equation, file.path(experiment, "02_fit_model.R"))\n'
    '  statistics <- readRDS(file.path(experiment, "full_response_statistics.rds"))\n'
    '  priors <- selected_priors(case$prior, spec$current_priors, statistics)',
    'spec <- baseline_spec(case$prior)\n  priors <- spec$priors')
fit = fit.replace('mode = if (case$mode == "prior") "prior" else "pilot"', 'mode = "baseline"')
fit = fit.replace('cores = case$chains', 'cores = case$cores')
fit = fit.replace('    original_source_md5 = unname(tools::md5sum(file.path(experiment, "02_fit_model.R"))),\n', '')
fit = fit.replace('max_treedepth = case$treedepth, k = k, seed = 1234L,',
    'max_treedepth = case$treedepth, refresh = case$refresh, k = k, seed = 1234L,')
fit = fit.replace('seed = 1234L, refresh = 10L, save_warmup = TRUE,',
    'seed = 1234L, refresh = case$refresh, save_warmup = TRUE, thin = 1L,')
fit = fit.replace('output_dir = csv_dir)',
    'output_dir = csv_dir, file = file.path(run, "fit"), file_refit = "always")')
(target / "02_fit_case.R").write_text(fit, encoding="utf-8", newline="\n")

# Six baseline cases and three sampler sensitivities paired by prior choice.
rows = []
def add_pair(subset, variant, iterations=5000, warmup=1000, delta=.97, tree=14):
    for prior in ("custom", "default"):
        rows.append(dict(case_id=len(rows)+1, label=f"{subset}_{prior}_{variant}",
            subset=subset, prior=prior, mode="posterior", chains=4, iter=iterations,
            warmup=warmup, adapt_delta=delta, treedepth=tree, threads=2, cpus=8,
            k_depth=15, equation="original_region_depth", cores=4, refresh=50,
            variant=variant, memory_gb=64 if subset == "p50" else 32))
for subset in ("p10", "p20", "p50"):
    add_pair(subset, "baseline")
add_pair("p10", "tree12", tree=12)
add_pair("p10", "delta99", delta=.99)
add_pair("p10", "warmup2500", iterations=6500, warmup=2500)
with (target / "cases.tsv").open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)

monitor = (source / "monitor_all_dcm_v2.sh").read_text(encoding="utf-8")
monitor = monitor.replace("k150 and equation150.", "k150, equation150 and baseline150.")
monitor = monitor.replace("experiments/equation150_*/manifest.tsv)",
    "experiments/equation150_*/manifest.tsv experiments/baseline150_*/manifest.tsv)")
monitor = monitor.replace("logs/dcm-e150-c[0-9]*.out)",
    "logs/dcm-e150-c[0-9]*.out logs/dcm-base150-c[0-9]*.out)")
monitor = monitor.replace('diagnostics_awk="DCM_Ibex_equation_tests/live_diagnostics.awk"',
    'diagnostics_awk="DCM_Ibex_baseline_tests/live_diagnostics.awk"')
monitor = monitor.replace("fast150, k150 and equation150", "fast150, k150, equation150 and baseline150")
monitor = monitor.replace("experiments/equation150_*/manifest.tsv; do",
    "experiments/equation150_*/manifest.tsv experiments/baseline150_*/manifest.tsv; do")
(target / "monitor_all_dcm_v3.sh").write_text(monitor, encoding="utf-8", newline="\n")
print(f"Built {len(rows)} baseline cases; {sum(r['cpus'] for r in rows)} CPUs if all fit jobs run together.")
