"""Targeted, descriptive comparisons from every available result table."""
from pathlib import Path
import csv, json, math, statistics, collections
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'DCM_all_run_comparison_20261016'; SNAP=ROOT/'dcm_snapshot_20261016_090217'
def read(path):
    with path.open(newline='',encoding='utf-8-sig') as f:return list(csv.DictReader(f))
def n(x):
    try:return float(x)
    except:return None
rows=read(OUT/'01_all_registered_runs.csv'); ext={int(r['case_id']):r for r in rows if r['campaign']=='extensive150'}
def show(ids,fields):
    for i in ids:
        r=ext[i];print(i,' | '.join(f'{k}={r.get(k)}' for k in fields))
print('PRIORS / COMPUTE / NEW DATA')
show([2,4,6,8,10,18,20,49,55,56,57,58,59,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,88,94,95],
['prior','threads','cpus','cpu_model','data_md5','requested_core_hours','fit_hours','max_mcse_over_sd',
 'training_conditional_rmse','training_population_rmse','validation_conditional_rmse','validation_population_rmse',
 'validation_conditional_bias','validation_conditional_peak_height_mae','validation_conditional_interval_width90'])
regions=read(OUT/'all_validation_conditional_by_region.csv')
population=read(OUT/'all_validation_population_by_region.csv')
print('REGIONAL CONDITIONAL RMSE / SHAPE / BIAS / COVERAGE')
for i in [2,4,6,8,10,18,20,24,30,36,40,49,74,75,84,85,88,95]:
    rr=[r for r in regions if int(r['case_id'])==i]
    print('CASE',i,[(r['Region'],r['casts'],*[round(float(r[k]),3) for k in ['rmse','shape_rmse','bias','coverage90']]) for r in rr])
print('POPULATION PRIORS REGIONAL RMSE/BIAS')
for i in [6,74,75,76,49,77,78,79,84,85]:
    print(i,[(r['Region'],round(float(r['rmse']),3),round(float(r['bias']),3)) for r in population if int(r['case_id'])==i])
strata=read(OUT/'all_validation_conditional_strata.csv')
print('STRATA')
for i in [6,10,49,74,84,85,88,95]:
    print('CASE',i,[(r['stratum'],r['bin'],r['casts'],*[round(float(r[k]),3) for k in ['rmse','bias','coverage90','interval_width90']]) for r in strata if int(r['case_id'])==i])
print('SCALE PARAMETERS AND MAIN COEFFICIENTS')
for i in [6,49,74,75,76,77,78,79,84,85,88,94,95]:
    path=SNAP/ext[i]['run_path']/'checks'/'parameter_diagnostics.csv'
    d=read(path)
    print(i,[(r['variable'],round(float(r['mean']),4),round(float(r['sd']),4)) for r in d if r['variable'] in ['sigma','nu','b_Intercept','b_MLD_temp_MIX','b_log10_distance'] or r['variable'].startswith('sd_')])
print('MATCHED PRIORS AGAINST CUSTOM / K CHANGES')
pairs=read(OUT/'04_every_case_vs_original_k15.csv')
for r in pairs:
    if int(r['case_id']) in [2,6,8,10,18,20,49,74,84,85,95] and r['prediction_mode']=='conditional':print(r)
contrasts=read(OUT/'05_extensive_anchor_contrasts.csv')
for r in contrasts:
    if r['run_id'] in ['extensive150:074','extensive150:075','extensive150:077','extensive150:084','extensive150:085','extensive150:088']:
        print({k:v for k,v in r.items() if k in ['run_id','reference_run','same_cpu_model','same_training_data','delta_mean_cast_mse','descriptive_paired_se','fraction_casts_lower_mse','delta_validation_conditional_rmse','delta_validation_population_rmse']})
