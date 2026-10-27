"""Focused paired contrasts and uncertainty, using existing audited CSVs only."""
from pathlib import Path
import csv,json,math,collections,statistics,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'DCM_review_20261019'
EXP=ROOT/'dcm_snapshot_20261019_084423/experiments/followup150_20261018T150024Z_3898749'
def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def number(v):
    try:return float(v) if math.isfinite(float(v)) else None
    except (ValueError,TypeError):return None
def write(n,rs):
    fields=list(dict.fromkeys(k for r in rs for k in r))
    with (OUT/n).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fields);w.writeheader();w.writerows(rs)
rows=read(OUT/'all_registered_runs.csv');by={r['run_id']:r for r in rows}
ids={r['run_key']:r['run_id'] for r in rows}
metadata={r['cast_key']:r for r in read(EXP/'full_cast_metadata.csv')}
cast={}
for d in ('validation','reserve'):
    cast[d]=collections.defaultdict(dict)
    for c in read(OUT/f'all_{d}_conditional_cast_metrics.csv'):cast[d][ids[c['run_key']]][c['cast_key']]=c
weights={}
def bootstrap_weights(keys,d):
    signature=hashlib.sha256('\n'.join(keys).encode()).hexdigest()
    if signature in weights:return weights[signature]
    rng=np.random.default_rng(20261019)
    groups=collections.defaultdict(lambda:collections.defaultdict(list))
    for i,k in enumerate(keys):
        m=metadata[k];groups[m['Region']][(m['Ocean'],m['Region'],m['date'])].append(i)
    w=np.zeros((5000,len(keys)),dtype=float)
    for region,gs in groups.items():
        g=list(gs.values());counts=rng.multinomial(len(g),[1/len(g)]*len(g),size=5000)
        for j,indices in enumerate(g):w[:,indices]=counts[:,j,None]
    w/=w.sum(axis=1,keepdims=True)
    weights[signature]=(w,sum(len(gs) for gs in groups.values()))
    return weights[signature]
comparisons=[('followup150:005','followup150:001','mean_interaction_vs_original'),('followup150:003','followup150:001','mean_both_vs_original'),('followup150:005','followup150:003','mean_interaction_vs_both'),('followup150:027','followup150:005','training_composition'),('followup150:028','followup150:005','depth_basis'),('followup150:030','followup150:005','numeric_by_basis')]
for a,b in [(2,1),(4,3),(6,5),(7,1),(8,3),(9,5),(10,1),(11,3),(12,5)]:comparisons.append((f'followup150:{a:03}',f'followup150:{b:03}','seed_sampler_threads'))
for a,b in [(3,1),(5,1),(4,2),(6,2)]:comparisons.append((f'baseline150:{a:03}',f'baseline150:{b:03}','baseline_training_size'))
for a,b in [(3,2),(9,8),(18,13),(16,8),(17,8),(47,2),(49,8)]:comparisons.append((f'targeted150:{a:03}',f'targeted150:{b:03}','earlier_controlled'))
results=[];sensitivity=[]
for a,b,question in comparisons:
    for d in ('validation','reserve'):
        ca,cb=cast[d].get(a,{}),cast[d].get(b,{})
        if not ca or ca.keys()!=cb.keys():continue
        ha,hb=by[a].get(d+'_md5'),by[b].get(d+'_md5')
        if not ha or ha!=hb:continue
        keys=sorted(ca);w,groups=bootstrap_weights(keys,d)
        ea=np.array([float(ca[k]['mse']) for k in keys]);eb=np.array([float(cb[k]['mse']) for k in keys]);delta=ea-eb
        bMSE=w@delta;bRMSE=np.sqrt(w@ea)-np.sqrt(w@eb)
        r={'run_a':a,'run_b':b,'question':question,'dataset':d,'matched_casts':len(keys),'nested_date_groups':groups,'bootstrap_draws':5000,'delta_mean_cast_mse':float(delta.mean()),'delta_rmse':math.sqrt(float(ea.mean()))-math.sqrt(float(eb.mean())),'bootstrap_low_mse':float(np.quantile(bMSE,.025)),'bootstrap_high_mse':float(np.quantile(bMSE,.975)),'bootstrap_low_rmse':float(np.quantile(bRMSE,.025)),'bootstrap_high_rmse':float(np.quantile(bRMSE,.975)),'fraction_casts_a_better':float((delta<0).mean()),'a_diagnostic_flags':by[a]['review_diagnostic_flags'],'b_diagnostic_flags':by[b]['review_diagnostic_flags'],'interpretation':'Exploratory region-stratified Ocean/Region/date cluster bootstrap; reused cohorts; no model-selection adjustment'}
        for metric in ('shape_rmse','peak_depth_mae','peak_height_mae','coverage90','interval_score90','mae'):
            r['delta_'+metric]=number(by[a].get(d+'_conditional_'+metric))-number(by[b].get(d+'_conditional_'+metric))
        results.append(r)
        order=np.argsort(eb)[::-1]
        for remove in (0,1,5):
            keep=np.ones(len(keys),bool);keep[order[:remove]]=False
            sensitivity.append({'run_a':a,'run_b':b,'dataset':d,'removed_reference_b_highest_mse_casts':remove,'retained_casts':int(keep.sum()),'rmse_a':float(np.sqrt(ea[keep].mean())),'rmse_b':float(np.sqrt(eb[keep].mean())),'delta_mean_mse':float(delta[keep].mean()),'fraction_casts_a_better':float((delta[keep]<0).mean()),'interpretation':'Influence sensitivity only; the same reference-defined casts omitted from both models; full-data scores remain primary'})
write('focused_date_bootstrap_contrasts.csv',results);write('paired_influence_sensitivity.csv',sensitivity)

# Relate elapsed time to allocated CPU time, without assuming identical node load.
compute=[]
for a,b in [(10,1),(11,3),(12,5),(7,1),(8,3),(9,5)]:
    ra,rb=by[f'followup150:{a:03}'],by[f'followup150:{b:03}']
    compute.append({'run_a':ra['run_id'],'run_b':rb['run_id'],'threads_a':ra['threads'],'threads_b':rb['threads'],'delta_a':ra['adapt_delta'],'delta_b':rb['adapt_delta'],'fit_hours_a':ra['fit_hours'],'fit_hours_b':rb['fit_hours'],'speedup_b_over_a':number(rb['fit_hours'])/number(ra['fit_hours']),'allocated_core_hours_a':ra['requested_core_hours'],'allocated_core_hours_b':rb['requested_core_hours'],'core_hour_ratio_a_over_b':number(ra['requested_core_hours'])/number(rb['requested_core_hours']),'same_cpu_model':ra['cpu_model']==rb['cpu_model'],'host_a':ra['hostname'],'host_b':rb['hostname'],'a_divergences':ra['divergences'],'b_divergences':rb['divergences'],'a_validation_rmse':ra['validation_conditional_rmse'],'b_validation_rmse':rb['validation_conditional_rmse']})
write('compute_contrasts.csv',compute)

# Check predictor support and nested-date sharing, without inspecting outcomes.
splits=read(EXP/'split_membership.csv');groups=collections.defaultdict(list)
for m in splits:groups[m['dataset']].append(m)
support=[]
for train in ('primary_p10_d150','primary_p10_alt_d150','primary_p20_d150'):
    tr=groups[train];dates={(m['Ocean'],m['Region'],m['date']) for m in tr}
    for evaluation in ('primary_validation_d150','primary_reserve_d150'):
        for region in sorted({m['Region'] for m in groups[evaluation]}):
            te=[m for m in groups[evaluation] if m['Region']==region];tt=[m for m in tr if m['Region']==region]
            z={'training_dataset':train,'evaluation_dataset':evaluation,'Region':region,'training_casts':len(tt),'evaluation_casts':len(te),'unseen_date_casts':sum((m['Ocean'],m['Region'],m['date']) not in dates for m in te)}
            for predictor in ('log10_distance','MLD_temp_MIX'):
                v=[float(m[predictor]) for m in tt];e=[float(m[predictor]) for m in te]
                z[predictor+'_training_min']=min(v);z[predictor+'_training_max']=max(v);z[predictor+'_evaluation_min']=min(e);z[predictor+'_evaluation_max']=max(e);z[predictor+'_outside_training_range']=sum(x<min(v) or x>max(v) for x in e)
            support.append(z)
write('predictor_support.csv',support)
print('FOCUSED CONTRASTS',len(results))
for r in results:
    if r['question'] in ('mean_interaction_vs_original','mean_interaction_vs_both','training_composition','baseline_training_size'):print(json.dumps(r))
print('COMPUTE',json.dumps(compute))
