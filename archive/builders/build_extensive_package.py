from pathlib import Path
import shutil

root = Path(__file__).resolve().parent
old = root / 'DCM_Ibex_baseline_tests'
out = root / 'DCM_Ibex_extensive_tests'
out.mkdir(exist_ok=True)
def read(name): return (old / name).read_text(encoding='utf-8')
def write(name, value): (out / name).write_text(value, encoding='utf-8', newline='\n')
for name in ('common.R', 'live_diagnostics.awk'):
    write(name, read(name))

fit = read('02_fit_case.R').replace('baseline.R','campaign.R')
fit = fit.replace('k <- 15L','k <- c(case$k_depth, case$k_distance, case$k_mld)')
fit = fit.replace('paste0(case$subset, ".rds")','paste0(case$subset, "_model.rds")')
fit = fit.replace('bundle <- readRDS(data_path)\n  d <- validate_bundle(bundle)\n  check_baseline_support(d)',
                  'd <- readRDS(data_path)')
fit = fit.replace('spec <- baseline_spec(case$prior)','spec <- campaign_spec(case)')
fit = fit.replace('mode = "baseline"','mode = "extensive"')
fit = fit.replace('seed = 1234L', 'seed = case$seed')
fit = fit.replace('file.copy(data_path, file.path(run, "model_data.rds"))',
                  'saveRDS(list(model_data=d), file.path(run, "model_data.rds"))')
fit = fit.replace('threads = thread_spec), file.path(run, "model.stan"))',
    'threads = thread_spec, sample_prior=if(case$mode=="prior") "only" else "no"), file.path(run, "model.stan"))')
write('02_fit_case.R',fit)

checks = read('03_check_case.R').split('  # Check a region-balanced')[0]
checks = checks.replace('"ess_bulk", "ess_tail"','"ess_bulk", "ess_tail", "mcse_mean"')
checks = checks.replace('divergences == 0 && tree_hits == 0 && ebfmi >= .3',
    'divergences == 0 && ebfmi >= .3')
checks = checks.replace('diagnostic_screen_pass = passed)',
    'max_mcse_over_sd = finite_max(varying$mcse_mean / varying$sd),\n'
    '    convergence_screen_pass = passed, tree_efficiency_warning = tree_hits > 0,\n'
    '    diagnostic_screen_pass = passed && tree_hits == 0)')
checks += '''  source(file.path(experiment, "predictive.R"))
  result <- evaluate_predictions(fit, experiment, run, case, result)
  write.csv(result, file.path(run,"test_summary.csv"), row.names=FALSE)
  writeLines(c("Convergence and tree efficiency are reported separately.",
    "Convergence screen: Rhat <=1.01, bulk/tail ESS >=400, zero sampling divergences, E-BFMI >=0.3.",
    "Parameter MCSE/SD is also reported; inspect uncertainty of scientific predictions.",
    "Validation is on entire casts, excluded from all training subsets.",
    "Validation may share dates and ENSO levels; it is not a new-cruise/date test.",
    "Population predictive intervals condition on zero group effects; they do not integrate new-group variability.",
    "Raw observed peak-depth errors are noisy descriptive summaries, not the sole ranking criterion.",
    "Lag-1 residual correlation is a screening statistic; row spacing may be irregular.",
    "Inspect matched repeat-seed profile predictions for Monte Carlo stability."),
    file.path(checks,"interpretation.txt"))
  writeLines("Complete",file.path(run,"CHECKS_COMPLETE.txt"))
}, error=function(error) {
  writeLines(conditionMessage(error),file.path(run,"CHECKS_FAILED.txt"))
  stop(error)
})
'''
write('03_check_case.R', checks)
for name in ('fit.slurm','check.slurm','prepare.slurm','collect.slurm'):
    content = read(name).replace('01_prepare_baselines.R','01_prepare.R')
    content = content.replace('Historical baseline', 'Extensive campaign')
    if name=='prepare.slurm':
        content=content.replace('--mem=32G','--mem=64G').replace('--time=00:30:00','--time=02:00:00')
    if name=='check.slurm': content=content.replace('--time=02:00:00','--time=04:00:00')
    write(name,content)
monitor = read('monitor_all_dcm_v3.sh')
monitor = monitor.replace('experiments/baseline150_*/manifest.tsv',
                         'experiments/baseline150_*/manifest.tsv experiments/extensive150_*/manifest.tsv')
monitor = monitor.replace('logs/dcm-base150-c[0-9]*.out',
                         'logs/dcm-base150-c[0-9]*.out logs/dcm-x150-c[0-9]*.out')
monitor = monitor.replace('DCM_Ibex_baseline_tests/live_diagnostics.awk',
                         'DCM_Ibex_extensive_tests/live_diagnostics.awk')
monitor = monitor.replace('and baseline150','baseline150 and extensive150')
write('monitor_all_dcm_v4.sh',monitor)

submit=read('submit_tests.sh')
submit=submit.replace('12 separately schedulable baseline fits','102 separately schedulable new fits')
submit=submit.replace('DCM_Ibex_baseline_tests','DCM_Ibex_extensive_tests')
submit=submit.replace('BASE_','EXT_').replace('latest_fast150','latest_baseline150')
submit=submit.replace('baseline150_','extensive150_').replace('latest_baseline150_experiment.txt"\nsubmission=',
                                                     'latest_extensive150_experiment.txt"\nsubmission=')
# The source pointer must remain baseline150; output pointer is extensive150.
submit=submit.replace('cat logs/latest_extensive150_experiment.txt','cat logs/latest_baseline150_experiment.txt')
submit=submit.replace('monitor_all_dcm_v3.sh','monitor_all_dcm_v4.sh').replace('baseline.R','campaign.R')
submit=submit.replace('"$source_experiment/data/p20.rds"',
                      '"$source_experiment/data/p20.rds" "$source_experiment/data/p50.rds"',1)
submit=submit.replace('for subset in p10 p20;','for subset in p10 p20 p50;')
submit=submit.replace('validate_baseline_cases(read.delim(file.path(a[1],"cases.tsv"))); cat("12 baseline cases passed\\n")',
  'stopifnot(isTRUE(all.equal(read.delim(file.path(a[1],"cases.tsv")),campaign_cases(),check.attributes=FALSE))); cat("102 cases passed\\n")')
submit=submit.replace('dcm-base150','dcm-x150').replace('--cpus-per-task=8 --constraint=', '--cpus-per-task=16 --constraint=')
submit=submit.replace('case_id label subset prior mode chains iterations warmup delta tree threads cpus kd equation cores refresh variant memory',
    'case_id label block equation subset prior mode chains iterations warmup samples delta tree threads cpus kd km kl cores refresh seed family memory')
submit=submit.replace('--dependency="afterok:$prepare_job"',
                      '--dependency="afterok:$prepare_job" --kill-on-invalid-dep=yes')
submit=submit.replace('12 fits submitted: 8 CPUs each, 32G (10%/20%) or 64G (50%), time limit $fit_time.',
    '102 fits submitted: 96 posterior + 6 prior predictive; 4/8/16 CPUs each; 32G or 64G; time limit $fit_time.')
write('submit_tests.sh',submit)

# Standalone mocked submission exercise, automatically run before real submission.
test=read('self_test_submission.sh').replace('DCM_Ibex_baseline_tests','DCM_Ibex_extensive_tests')
test=test.replace('dcm-baseline-selftest','dcm-extensive-selftest').replace('BASE_','EXT_')
test=test.replace('for subset in p10 p20;','for subset in p10 p20 p50;')
test=test.replace('"$cpus" == 8 && "$memory" == 64G','"$cpus" == 16 && "$memory" == 64G')
test=test.replace('latest_baseline150','latest_extensive150').replace('1026','1206')
test=test.replace('== 13 ]]','== 103 ]]').replace('== 25 ]]','== 205 ]]')
test=test.replace('&& "$cpus" == 8','&& "$cpus" -ge 4')
test=test.replace('-v id="$id"', '-v id="$id" -v cpus="$cpus"')
test=test.replace('$4!=8','$4!=cpus').replace('dcm-base150','dcm-x150')
test=test.replace('expected=($7==5 || $7==6)?"64G":"32G"',
                  'expected=($7==85 || $7==87 || $7==89 || $7==91)?"64G":"32G"')
test=test.replace('fits!=12','fits!=102').replace('12 fits, correct resources/dependencies, 26 jobs total.',
                                              '102 fits, correct resources/dependencies, 206 jobs total.')
write('self_test_submission.sh',test)
print('Generated fitting, checking, submission and monitoring scripts.')
