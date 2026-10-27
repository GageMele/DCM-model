"""Read saved scores only and print compact scientific comparison panels."""
from pathlib import Path
import csv,json,collections
root=Path(__file__).resolve().parent
res=next((root/'dcm_snapshot_20261021_124243').glob('ibex/user/meleg/DCM_MODEL/experiments/onset*/results'))
def read(name):
    with (res/name).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
onsets={(int(r['case_id']),r['dataset'],r['prediction_mode'],r['threshold']):r for r in read('onset_all_scores.csv')}
mld={(int(r['case_id']),r['dataset'],r['prediction_mode']):r for r in read('mld_neighbourhood_scores.csv')}
for r in read('onset_results.csv'):
    cid=int(r['case_id'])
    if r['mode']=='prior':continue
    line={k:r[k] for k in ('case_id','equation','subset','split','fit_seconds','max_rhat','min_bulk_ess','min_tail_ess','min_ebfmi','divergences','treedepth_hits','diagnostic_screen_pass')}
    for d in ('validation','reserve'):
        for metric in ('rmse','shape_rmse','peak_depth_mae','coverage90','interval_score90','median_residual_lag1'):
            line[d+'_'+metric]=r[d+'_conditional_'+metric]
        o=onsets[cid,d,'conditional','0.5']
        for metric in ('onset_mae_given_both_detected','both_detected','false_positives','false_negatives','predicted_left_censored','rapid_rise_mae'):
            line[d+'_'+metric]=o[metric]
        line[d+'_mld_rmse']=mld[cid,d,'conditional']['mld_neighbourhood_rmse']
    print(json.dumps(line))
for file in ('onset_planned_paired_contrasts.csv','planned_controlled_contrasts.csv'):
    print('\nCONTRASTS',file)
    for r in read(file):
        if r.get('prediction_mode','conditional')=='conditional' and r.get('threshold','0.5')=='0.5' and (int(r['case_id']) in (8,9,14,15,20,21,25,26,27)):
            print(json.dumps(r))
print('\nCAST SCHEMA')
for file in ('onset_all_cast_metrics.csv','all_validation_conditional_cast_metrics.csv'):
    with (res/file).open(encoding='utf-8-sig',newline='') as f:
        rr=csv.DictReader(f);print(file,rr.fieldnames);print(next(rr))
