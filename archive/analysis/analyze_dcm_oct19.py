"""Audit existing local result tables and compare all DCM runs; no model execution.

All numerical analyses read the uploaded snapshot. No network access, model
fitting, chain parsing, or archive extraction occurs in this script.
"""
from pathlib import Path
import csv,collections,math,statistics,json,hashlib
import numpy as np

ROOT=Path(__file__).resolve().parent
SNAP=ROOT/'dcm_snapshot_20261019_084423'
EXP=SNAP/'experiments/followup150_20261018T150024Z_3898749'
RESULTS=EXP/'results'
OUT=ROOT/'DCM_review_20261019'
OUT.mkdir(exist_ok=True)

def read(p,sep=','):
    if not p.exists():return []
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter=sep))
def num(v):
    try:
        z=float(v);return z if math.isfinite(z) else None
    except (ValueError,TypeError):return None
def truth(v):return str(v).lower()=='true'
def write(n,rs):
    keys=list(dict.fromkeys(k for r in rs for k in r))
    with (OUT/n).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,keys);w.writeheader();w.writerows(rs)
def local(p):return SNAP/str(p).split('/ibex/user/meleg/DCM_MODEL/')[-1]

rows=read(RESULTS/'all_results.csv')
assert len(rows)==270 and len({r['run_key'] for r in rows})==270
configs={r['run_path']:r for r in read(ROOT/'DCM_review_20261018/config_metadata.csv')}
for r in rows:
    path=str(r.get('run_dir','')).split('/ibex/user/meleg/DCM_MODEL/')[-1]
    for k,v in configs.get(path,{}).items():
        if r.get(k) in (None,'','NA') and v not in ('','NA'):r[k]=v
    r['short_campaign']=r['campaign'].split('_')[0]
    cid=next((str(r[k]) for k in ('case_id','case') if r.get(k) not in (None,'','NA')),Path(path).name)
    r['run_id']=r['short_campaign']+':'+(cid.zfill(3) if cid.isdigit() else cid)
    r['prior_only']=r.get('mode')=='prior'
    r['fit_hours']=num(r.get('fit_seconds'))/3600 if num(r.get('fit_seconds')) is not None else None
    r['configured_samples_per_chain']=num(r.get('samples')) or ((num(r.get('iter')) or 0)-(num(r.get('warmup')) or 0))
    vals=[num(r.get(k)) for k in ('max_rhat','min_bulk_ess','min_tail_ess','divergences','min_ebfmi')]
    screen=(vals[0]<1.01 and vals[1]>=400 and vals[2]>=400 and vals[3]==0 and vals[4]>=.3) if all(v is not None for v in vals) else None
    r['review_convergence_screen']=screen
    r['review_diagnostic_flags']='; '.join(k for k,bad in [('Rhat>=1.01',vals[0] is not None and vals[0]>=1.01),('bulkESS<400',vals[1] is not None and vals[1]<400),('tailESS<400',vals[2] is not None and vals[2]<400),('divergences',vals[3] is not None and vals[3]>0),('E-BFMI<0.3',vals[4] is not None and vals[4]<.3)] if bad)
    r['zero_tree_and_convergence']=screen is True and num(r.get('treedepth_hits'))==0
by={r['run_key']:r for r in rows}
assert len({r['run_id'] for r in rows})==len(rows)
short={r['run_id']:r for r in rows}
cast={};regions={};strata={}
audit=[]
for dataset in ('training','validation','reserve'):
    for mode in ('conditional','population'):
        prefix=f'{dataset}_{mode}'
        cs=read(RESULTS/f'all_{prefix}_cast_metrics.csv');cast[prefix]=collections.defaultdict(dict)
        for c in cs:
            k=c['run_key'];assert k in by
            assert c['cast_key'] not in cast[prefix][k]
            cast[prefix][k][c['cast_key']]=c
        for k,casts in cast[prefix].items():
            r=by[k];value=math.sqrt(statistics.mean(float(c['mse']) for c in casts.values()))
            reported=num(r.get(prefix+'_rmse'));assert reported is not None
            diff=abs(value-reported);assert diff<1e-9,(k,prefix,diff)
            audit.append({'run_key':k,'dataset':dataset,'mode':mode,'casts':len(casts),'reported_rmse':reported,'recomputed_rmse':value,'absolute_discrepancy':diff})
            r[prefix+'_cast_signature']=hashlib.sha256('\n'.join(sorted(casts)).encode()).hexdigest()
            errors=sorted([float(c['mse']) for c in casts.values()],reverse=True)
            r[prefix+'_top1_mse_share']=errors[0]/sum(errors)
            r[prefix+'_top5_mse_share']=sum(errors[:5])/sum(errors)
            r[prefix+'_median_cast_rmse']=statistics.median(math.sqrt(float(c['mse'])) for c in casts.values())
            r[prefix+'_trim_top1_rmse']=math.sqrt(statistics.mean(errors[1:]))
            r[prefix+'_trim_top5_rmse']=math.sqrt(statistics.mean(errors[5:]))
        regions[prefix]=read(RESULTS/f'all_{prefix}_by_region.csv')
        strata[prefix]=read(RESULTS/f'all_{prefix}_strata.csv')
        for k in by:
            rg=[x for x in regions[prefix] if x['run_key']==k]
            if rg:by[k][prefix+'_macro_region_rmse']=math.sqrt(statistics.mean(float(x['rmse'])**2 for x in rg))
        write('all_'+prefix+'_cast_metrics.csv',cs)
        write('all_'+prefix+'_by_region.csv',regions[prefix])
        write('all_'+prefix+'_strata.csv',strata[prefix])
assert len(audit)==860,(len(audit),'RMSE audit count')
write('rmse_audit.csv',audit)

# Keep validation cohorts separate using both saved data hashes and cast identities.
for r in rows:
    h=r.get('validation_md5');sig=r.get('validation_conditional_cast_signature')
    r['validation_cohort']='not_scored' if not sig else 'primary_md5_verified' if h=='950184ac377851361fd4af891249ca22' else 'date_blocked_md5_verified' if h=='7a1cb728145f25ada97e9cfdcad8ea98' else 'other_hash_or_membership_only'

splits=read(EXP/'split_membership.csv');sets=collections.defaultdict(set)
for s in splits:sets[s['dataset']].add(s['cast_key'])
split_audit=[]
for n,s in sets.items():split_audit.append({'check':'dataset size','dataset':n,'casts':len(s),'passed':True})
reserve_names=[n for n in sets if 'reserve' in n]
assert reserve_names
reserve=sets[reserve_names[0]]
assert all(sets[n]==reserve for n in reserve_names),reserve_names
split_audit.append({'check':'primary and date-blocked reserve identities','dataset':'; '.join(reserve_names),'casts':len(reserve),'passed':True})
for n,s in sets.items():
    if n in reserve_names:continue
    overlap=len(reserve&s);assert overlap==0,(n,overlap)
    split_audit.append({'check':'reserve overlap','dataset':n,'overlap_casts':overlap,'passed':True})
write('split_audit.csv',split_audit)

# Annotate influential casts with original metadata, retaining all observed data.
metadata={r['cast_key']:r for r in read(EXP/'full_cast_metadata.csv')}
influence=[];contributions=[]
for dataset in ('validation','reserve'):
    prefix=dataset+'_conditional'
    for r in rows:
        cs=cast[prefix].get(r['run_key'],{})
        if not cs:continue
        total=sum(float(c['mse']) for c in cs.values());n=len(cs)
        for rank,c in enumerate(sorted(cs.values(),key=lambda x:float(x['mse']),reverse=True),1):
            influence.append({**metadata.get(c['cast_key'],{}),**c,'run_id':r['run_id'],'dataset':dataset,'mse_rank':rank,'fraction_total_mse':float(c['mse'])/total,'magnitude_component_mse':float(c['mse'])-float(c['shape_mse'])})
        for region in sorted({c['Region'] for c in cs.values()}):
            cr=[c for c in cs.values() if c['Region']==region];mse=sum(float(c['mse']) for c in cr)
            contributions.append({'run_id':r['run_id'],'dataset':dataset,'Region':region,'casts':len(cr),'regional_rmse':math.sqrt(mse/len(cr)),'fraction_total_mse':mse/total,'contribution_overall_mse':mse/n,'mean_cast_bias':statistics.mean(float(c['bias']) for c in cr),'shape_rmse':math.sqrt(statistics.mean(float(c['shape_mse']) for c in cr))})
write('cast_influence.csv',influence);write('regional_error_contributions.csv',contributions)

# Recompute all saved paired contrasts from aligned individual casts.
pairs=read(RESULTS/'all_matched_cast_comparisons.csv')
for p in pairs:
    ca=cast[p['dataset']+'_conditional'][p['run_a']];cb=cast[p['dataset']+'_conditional'][p['run_b']]
    assert ca.keys()==cb.keys()
    d=statistics.mean(float(ca[k]['mse'])-float(cb[k]['mse']) for k in ca)
    assert abs(d-float(p['delta_mean_cast_mse']))<1e-9
write('all_matched_cast_comparisons.csv',pairs)

# Planned contrasts include the campaign's region-stratified date bootstrap.
for name in ('planned_controlled_contrasts.csv','planned_regional_contrasts.csv','all_predictive_batch_sensitivity.csv','all_prior_predictive_quantiles.csv'):
    write(name,read(RESULTS/name))

coverage=[]
for campaign in sorted({r['short_campaign'] for r in rows}):
    rs=[r for r in rows if r['short_campaign']==campaign];post=[r for r in rs if truth(r.get('fit_completed')) and not r['prior_only']]
    coverage.append({'campaign':campaign,'registered':len(rs),'posterior_complete':len(post),'prior_complete':sum(truth(r.get('fit_completed')) and r['prior_only'] for r in rs),'primary_scored':sum(r['validation_cohort']=='primary_md5_verified' for r in post),'date_scored':sum(r['validation_cohort']=='date_blocked_md5_verified' for r in post),'reserve_scored':sum(num(r.get('reserve_conditional_rmse')) is not None for r in post),'convergence_pass':sum(r['review_convergence_screen'] is True for r in post),'zero_tree_and_convergence':sum(r['zero_tree_and_convergence'] for r in post)})
write('coverage.csv',coverage)
write('all_registered_runs.csv',rows)
write('completed_posterior.csv',[r for r in rows if truth(r.get('fit_completed')) and not r['prior_only']])
write('prior_runs.csv',[r for r in rows if r['prior_only']])
write('followup_all_cases.csv',[r for r in rows if r['short_campaign']=='followup150'])
meta={'runs':len(rows),'completed':sum(truth(r.get('fit_completed')) for r in rows),'posterior_complete':sum(truth(r.get('fit_completed')) and not r['prior_only'] for r in rows),'prior_complete':sum(truth(r.get('fit_completed')) and r['prior_only'] for r in rows),'rmse_aggregates':len(audit),'maximum_rmse_discrepancy':max(x['absolute_discrepancy'] for x in audit),'paired_contrasts_verified':len(pairs)}
(OUT/'audit.json').write_text(json.dumps(meta,indent=2))
print(json.dumps(meta));print(json.dumps(coverage))
for run_id in ('followup150:001','followup150:003','followup150:005','followup150:027'):
    r=short[run_id];print('CANDIDATE',json.dumps({k:r.get(k) for k in ('run_id','validation_conditional_rmse','reserve_conditional_rmse','validation_conditional_top5_mse_share','reserve_conditional_top5_mse_share','validation_conditional_trim_top5_rmse','reserve_conditional_trim_top5_rmse')}))
    print('RESERVE REGIONS',json.dumps([x for x in contributions if x['run_id']==run_id and x['dataset']=='reserve']))
    print('TOP RESERVE',json.dumps([{k:x.get(k) for k in ('cast_key','Region','fraction_total_mse','mse','bias','shape_mse','observed_peak_depth','peak_height_abs_error')} for x in influence if x['run_id']==run_id and x['dataset']=='reserve' and x['mse_rank']<=5]))
