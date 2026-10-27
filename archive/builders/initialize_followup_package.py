"""Copy the tested execution helpers and add follow-up scoring and monitoring."""
from pathlib import Path
import shutil
ROOT=Path(__file__).resolve().parent
OLD=ROOT/'DCM_Ibex_targeted_tests'; NEW=ROOT/'DCM_Ibex_followup_tests'
NEW.mkdir(exist_ok=True)
for name in ('common.R','02_fit_case.R','03_check_case.R','fit.slurm','check.slurm','prepare.slurm','collect.slurm','live_diagnostics.awk'):
    text=(OLD/name).read_text(encoding='utf-8-sig')
    text=text.replace('config <- list(case = case, mode = "targeted"','config <- list(case = case, mode = "followup"')
    text=text.replace('Targeted campaign case','Follow-up campaign case').replace('all 63 requested models','all 34 requested models')
    if name=='check.slurm':text=text.replace('--time=04:00:00','--time=08:00:00')
    if name=='collect.slurm':text=text.replace('--mem=2G','--mem=4G').replace('--time=00:15:00','--time=01:00:00')
    (NEW/name).write_text(text,encoding='utf-8',newline='\n')
shutil.copyfile(OLD/'campaign.R',NEW/'base_models.R')

text=(OLD/'predictive.R').read_text(encoding='utf-8-sig')
text=text.replace('plot_validation_examples<-function(tab,checks,mode)', 'plot_validation_examples<-function(tab,checks,mode,dataset="validation")')
text=text.replace('paste0("validation_",mode,"_example_profiles', 'paste0(dataset,"_",mode,"_example_profiles')
text=text.replace('paste("Independent validation casts:",mode)', 'paste("Independent",dataset,"casts:",mode)')
text=text.replace('n_draws=1000L','n_draws=2000L')
text=text.replace('  fitted_data<-readRDS(data_file(experiment,case))', '  reserve<-readRDS(reserve_validation_file(experiment,case))\n  fitted_data<-readRDS(data_file(experiment,case))')
text=text.replace('  stopifnot(!any(validation$cast_key %in% fitted_data$cast_key))','  stopifnot(!any(validation$cast_key %in% fitted_data$cast_key),\n    !any(reserve$cast_key %in% fitted_data$cast_key),!any(reserve$cast_key %in% validation$cast_key))')
text=text.replace('  result$validation_md5<-', '  result$reserve_md5<-unname(tools::md5sum(reserve_validation_file(experiment,case)))\n  result$validation_md5<-')
text=text.replace('  scores<-list()', '  scores<-list();batches<-list()')
text=text.replace('for(dataset in c("training","validation"))','for(dataset in c("training","validation","reserve"))')
text=text.replace('    d<-if(dataset=="training")train else validation', '    d<-switch(dataset,training=train,validation=validation,reserve=reserve)')
needle='    # Store marginals for seen/new groups and fixed training-derived strata.'
insertion='''    # Predictive batch sensitivity reuses existing draws; it is not a formal MCSE.
    result[[paste0(tag,"_fraction_observed_below_minus_one")]]<-mean(tab$chl_norm < -1)
    result[[paste0(tag,"_fraction_expected_below_minus_one")]]<-mean(tab$predicted < -1)
    result[[paste0(tag,"_fraction_predictive_below_minus_one")]]<-mean(y < -1)
    for(batch in 1:2) {
      ii<-if(batch==1L)seq_len(floor(nrow(mu)/2))else seq.int(floor(nrow(mu)/2)+1L,nrow(mu))
      qt<-apply(y[ii,,drop=FALSE],2,quantile,c(.05,.95))
      bt<-tab;bt$predicted<-colMeans(mu[ii,,drop=FALSE]);bt$lower<-qt[1,];bt$upper<-qt[2,]
      bt$interval_score90<-bt$upper-bt$lower+20*pmax(bt$lower-bt$chl_norm,0)+20*pmax(bt$chl_norm-bt$upper,0)
      bc<-do.call(rbind,lapply(split(bt,bt$cast_key),profile_metrics))
      batches[[length(batches)+1L]]<-data.frame(dataset,mode,batch,draws=length(ii),score_casts(bc))
    }
'''
assert needle in text;text=text.replace(needle,insertion+needle)
text=text.replace('  write.csv(do.call(rbind,scores),file.path(checks,"prediction_scores.csv"),row.names=FALSE)', '  write.csv(do.call(rbind,scores),file.path(checks,"prediction_scores.csv"),row.names=FALSE)\n  write.csv(do.call(rbind,batches),file.path(checks,"predictive_batch_sensitivity.csv"),row.names=FALSE)')
text=text.replace('  for(mode in c("conditional","population")) {\n    tab<-read.csv(file.path(checks,paste0("validation_",mode,"_observations.csv")),stringsAsFactors=FALSE)\n    plot_validation_examples(tab,checks,mode)', '  for(dataset in c("validation","reserve"))for(mode in c("conditional","population")) {\n    tab<-read.csv(file.path(checks,paste0(dataset,"_",mode,"_observations.csv")),stringsAsFactors=FALSE)\n    plot_validation_examples(tab,checks,mode,dataset)')
text=text.replace('validation; DCM diagnostics confined', 'validation plus a disjoint reserve; DCM diagnostics confined')
(NEW/'predictive.R').write_text(text,encoding='utf-8',newline='\n')

monitor=(OLD/'monitor_all_dcm_v5.sh').read_text(encoding='utf-8-sig')
monitor=monitor.replace('experiments/targeted150_*/manifest.tsv', 'experiments/targeted150_*/manifest.tsv experiments/followup150_*/manifest.tsv')
monitor=monitor.replace('logs/dcm-t150-c[0-9]*.out)', 'logs/dcm-t150-c[0-9]*.out logs/dcm-u150-*-c[0-9]*.out)')
monitor=monitor.replace('all six 150 m campaigns','all recorded 150 m campaigns')
(NEW/'monitor_all_dcm_v6.sh').write_text(monitor,encoding='utf-8',newline='\n')
export=(OLD/'export_results.sh').read_text(encoding='utf-8-sig').replace('dcm_targeted_comparison_','dcm_followup_comparison_')
(NEW/'export_results.sh').write_text(export,encoding='utf-8',newline='\n')
print('Follow-up execution, predictive scoring, monitoring and export helpers initialized.')
