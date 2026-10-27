"""Audit all saved DCM results and calculate matched, exploratory contrasts.

No fitted model is loaded or rerun. All new calculations use the uploaded CSVs.
The older review supplies missing configuration metadata, never updated scores.
"""
from pathlib import Path
import csv, collections, json, math, statistics, hashlib
import numpy as np

ROOT=Path(__file__).resolve().parent
SNAP=ROOT/'dcm_snapshot_20261021_124243'
RES=next(SNAP.glob('ibex/user/meleg/DCM_MODEL/experiments/onset*/results'))
OUT=ROOT/'DCM_review_20261021'
OUT.mkdir(exist_ok=True)

def read(path, sep=','):
    with Path(path).open(encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f,delimiter=sep))

def num(x):
    try:
        a=float(x)
        return a if math.isfinite(a) else None
    except (ValueError,TypeError):return None

def write(name, rows):
    rows=list(rows)
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with (OUT/name).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)

def mean(values):
    a=[num(x) for x in values];a=[x for x in a if x is not None]
    return statistics.mean(a) if a else None

rows=read(RES/'all_results.csv')
assert len(rows)==304 and len({r['run_key'] for r in rows})==304
old={r['run_key']:r for r in read(ROOT/'DCM_review_20261019/all_registered_runs.csv')}
config_keys=('formula','configured_k','max_treedepth','cmdstan_version','k','allocated_cpus','input_md5')
for r in rows:
    r['short_campaign']=r['campaign'].split('_')[0]
    cid=r.get('case_id')
    if cid in (None,'NA',''):cid=r.get('case')
    if cid in (None,'NA',''):cid=r['run_dir'].split('/')[-1]
    r['run_id']=r['short_campaign']+':'+(str(cid).zfill(3) if str(cid).isdigit() else str(cid))
    for k in config_keys:
        if r.get(k) in (None,'NA','') and old.get(r['run_key'],{}).get(k) not in (None,'NA',''):
            r[k]=old[r['run_key']][k]
    r['evaluation_only_reuse']=r.get('execution')=='reuse'
    r['unique_fit_identity']=r['source_run'] if r['evaluation_only_reuse'] else r['run_dir']
    r['fit_hours']=num(r.get('fit_seconds'))/3600 if num(r.get('fit_seconds')) is not None else None
    vals=[num(r.get(k)) for k in ('max_rhat','min_bulk_ess','min_tail_ess','divergences','min_ebfmi')]
    r['review_convergence_screen']=(vals[0]<1.01 and vals[1]>=400 and vals[2]>=400 and vals[3]==0 and vals[4]>=.3) if all(v is not None for v in vals) else None
    r['review_zero_tree_screen']=r['review_convergence_screen'] is True and num(r.get('treedepth_hits'))==0
    r['review_diagnostic_flags']='; '.join(k for k,bad in [('Rhat>=1.01', vals[0] is not None and vals[0]>=1.01),('bulkESS<400',vals[1] is not None and vals[1]<400),('tailESS<400',vals[2] is not None and vals[2]<400),('divergences',vals[3] is not None and vals[3]>0),('E-BFMI<0.3',vals[4] is not None and vals[4]<.3),('tree_hits',num(r.get('treedepth_hits')) is not None and num(r.get('treedepth_hits'))>0)] if bad)
by={r['run_key']:r for r in rows}
byid={r['run_id']:r for r in rows}
cast={}; audit=[]
for d in ('training','validation','reserve'):
    for mode in ('conditional','population'):
        prefix=d+'_'+mode
        cr=read(RES/('all_'+prefix+'_cast_metrics.csv'))
        groups=collections.defaultdict(dict)
        for c in cr:
            assert c['cast_key'] not in groups[c['run_key']]
            groups[c['run_key']][c['cast_key']]=c
        cast[prefix]=groups
        for key,cs in groups.items():
            r=by[key]
            expected=math.sqrt(statistics.mean(float(c['mse']) for c in cs.values()))
            reported=num(r.get(prefix+'_rmse'));assert reported is not None
            delta=abs(expected-reported);assert delta<1e-9,(key,prefix,delta)
            audit.append({'run_key':key,'run_id':r['run_id'],'target':prefix,'casts':len(cs),'reported_rmse':reported,'recomputed_rmse':expected,'absolute_difference':delta})
            signature=hashlib.sha256('\n'.join(sorted(cs)).encode()).hexdigest()
            r[prefix+'_cohort_signature']=signature
            errors=sorted([float(c['mse']) for c in cs.values()],reverse=True)
            r[prefix+'_top1_mse_share']=errors[0]/sum(errors)
            r[prefix+'_top5_mse_share']=sum(errors[:5])/sum(errors)
            r[prefix+'_median_cast_rmse']=statistics.median(math.sqrt(e) for e in errors)
            regional=collections.defaultdict(list)
            for c in cs.values():regional[c['Region']].append(float(c['mse']))
            r[prefix+'_macro_region_rmse']=math.sqrt(statistics.mean(statistics.mean(x) for x in regional.values()))
        write('all_'+prefix+'_cast_metrics.csv',cr)
        for suffix in ('by_region','strata'):
            write('all_'+prefix+'_'+suffix+'.csv',read(RES/('all_'+prefix+'_'+suffix+'.csv')))
assert len(audit)==1136
write('rmse_audit.csv',audit)

scores=read(RES/'onset_all_scores.csv')
mld=read(RES/'mld_neighbourhood_scores.csv')
onset_cast=read(RES/'onset_all_cast_metrics.csv')
onset={(int(r['case_id']),r['dataset'],r['prediction_mode'],r['threshold']):r for r in scores}
onset_cast_index=collections.defaultdict(dict)
for c in onset_cast:
    k=(int(c['case_id']),c['dataset'],c['prediction_mode'],c['threshold'])
    assert c['cast_key'] not in onset_cast_index[k]
    onset_cast_index[k][c['cast_key']]=c

# Independently verify each onset MAE and all missing/censored status counts.
onset_audits=[]
for key,s in onset.items():
    cs=list(onset_cast_index[key].values())
    assert len(cs)==int(s['casts'])
    for side in ('observed','predicted'):
        for status,column in [('detected','detected'),('left_censored','left_censored'),('not_reached_in_window','not_reached'),('gap_censored','gap_censored'),('insufficient_support','insufficient_support')]:
            assert sum(c[side+'_status']==status for c in cs)==int(s[side+'_'+column])
    ma=mean(c['onset_abs_error'] for c in cs)
    assert abs(ma-float(s['onset_mae_given_both_detected']))<1e-10
    onset_audits.append({'case_id':key[0],'dataset':key[1],'prediction_mode':key[2],'threshold':key[3],'casts':len(cs),'recomputed_mae':ma,'absolute_difference':abs(ma-float(s['onset_mae_given_both_detected']))})
write('onset_score_audit.csv',onset_audits)

for m in mld:
    r=byid['onset150:'+m['case_id'].zfill(3)]
    r[m['dataset']+'_'+m['prediction_mode']+'_mld_neighbourhood_rmse']=m['mld_neighbourhood_rmse']
for s in scores:
    if s['threshold']!='0.5':continue
    r=byid['onset150:'+s['case_id'].zfill(3)]
    for k,v in s.items():
        if k not in ('case_id','equation','subset','split','dataset','prediction_mode','threshold'):
            r[s['dataset']+'_'+s['prediction_mode']+'_onset50_'+k]=v
write('all_registered_runs.csv',rows)
write('onset_all_threshold_scores.csv',scores)
write('onset_all_cast_metrics.csv',onset_cast)
write('mld_neighbourhood_scores.csv',mld)
write('saved_onset_paired_contrasts.csv',read(RES/'onset_planned_paired_contrasts.csv'))
write('saved_controlled_contrasts.csv',read(RES/'planned_controlled_contrasts.csv'))

# Use the cast metadata embedded in onset scores for reproducible date clusters.
meta={c['cast_key']:c for c in onset_cast}
weight_cache={}
def bootstrap(keys):
    signature=hashlib.sha256('\n'.join(keys).encode()).hexdigest()
    if signature in weight_cache:return weight_cache[signature]
    clusters=collections.defaultdict(lambda:collections.defaultdict(list))
    for i,k in enumerate(keys):
        c=meta[k];clusters[c['Region']][(c['Ocean'],c['Region'],c['date'])].append(i)
    rng=np.random.default_rng(20261021)
    w=np.zeros((5000,len(keys)))
    for groups in clusters.values():
        group_list=list(groups.values())
        counts=rng.multinomial(len(group_list),[1/len(group_list)]*len(group_list),size=5000)
        for j,indices in enumerate(group_list):w[:,indices]=counts[:,j,None]
    w/=w.sum(axis=1,keepdims=True)
    result=(w,sum(len(g) for g in clusters.values()))
    weight_cache[signature]=result
    return result

comparison_pairs=[(a,b) for a,b in [(2,1),(4,3),(6,5),(7,1),(8,1),(9,1),(10,1),(11,1),(12,1),(13,3),(14,3),(15,3),(16,3),(17,3),(18,3),(19,5),(20,5),(21,5),(22,5),(23,5),(24,5),(9,8),(15,14),(21,20),(25,8),(26,8),(27,8),(28,11),(29,11),(30,11)]]
contrasts=[];onset_contrasts=[];regional_contrasts=[];influence=[]
for a,b in comparison_pairs:
    ra=byid[f'onset150:{a:03}'];rb=byid[f'onset150:{b:03}']
    for d in ('validation','reserve'):
        for mode in ('conditional','population'):
            prefix=d+'_'+mode
            ca,cb=cast[prefix][ra['run_key']],cast[prefix][rb['run_key']]
            assert ca.keys()==cb.keys()
            keys=sorted(ca);w,clusters=bootstrap(keys)
            ea=np.array([float(ca[k]['mse']) for k in keys]);eb=np.array([float(cb[k]['mse']) for k in keys])
            rmse_difference=np.sqrt(w@ea)-np.sqrt(w@eb)
            delta=ea-eb
            out={'candidate':ra['run_id'],'reference':rb['run_id'],'candidate_equation':ra['equation'],'reference_equation':rb['equation'],'training_subset':ra['subset'],'split':ra['split'],'dataset':d,'prediction_mode':mode,'matched_casts':len(keys),'nested_date_clusters':clusters,'bootstrap_replicates':5000,'delta_rmse':float(np.sqrt(ea.mean())-np.sqrt(eb.mean())),'rmse_difference_ci_low':float(np.quantile(rmse_difference,.025)),'rmse_difference_ci_high':float(np.quantile(rmse_difference,.975)),'fraction_casts_lower_mse':float((delta<0).mean()),'candidate_convergence_pass':ra['review_convergence_screen'],'reference_convergence_pass':rb['review_convergence_screen'],'interpretation':'Exploratory paired bootstrap of Ocean/Region/date clusters within regions; reused development holdouts, no search adjustment'}
            for metric in ('shape_rmse','peak_depth_mae','mae','coverage90','interval_score90','mld_neighbourhood_rmse'):
                na,nb=num(ra.get(prefix+'_'+metric)),num(rb.get(prefix+'_'+metric))
                out['delta_'+metric]=na-nb if na is not None and nb is not None else None
            contrasts.append(out)
            for region in sorted({c['Region'] for c in ca.values()}):
                ii=[i for i,k in enumerate(keys) if ca[k]['Region']==region]
                regional_contrasts.append({**{k:out[k] for k in ('candidate','reference','dataset','prediction_mode')},'Region':region,'casts':len(ii),'rmse_candidate':float(np.sqrt(ea[ii].mean())),'rmse_reference':float(np.sqrt(eb[ii].mean())),'delta_mean_cast_mse':float(delta[ii].mean()),'fraction_casts_lower_mse':float((delta[ii]<0).mean())})
            for n in (0,1,5):
                order=np.argsort(eb)[::-1];keep=np.ones(len(keys),dtype=bool);keep[order[:n]]=False
                influence.append({'candidate':ra['run_id'],'reference':rb['run_id'],'dataset':d,'prediction_mode':mode,'reference_highest_error_casts_removed':n,'retained_casts':int(keep.sum()),'candidate_rmse':float(np.sqrt(ea[keep].mean())),'reference_rmse':float(np.sqrt(eb[keep].mean())),'delta_rmse':float(np.sqrt(ea[keep].mean())-np.sqrt(eb[keep].mean())),'interpretation':'Influence sensitivity only; identical reference-defined removals; full-data scores remain primary'})
            for threshold in ('0.25','0.5','1'):
                oa=onset_cast_index[a,d,mode,threshold];ob=onset_cast_index[b,d,mode,threshold]
                shared=sorted(k for k in keys if num(oa[k]['onset_abs_error']) is not None and num(ob[k]['onset_abs_error']) is not None)
                ww,nc=bootstrap(shared)
                eoa=np.array([float(oa[k]['onset_abs_error']) for k in shared]);eob=np.array([float(ob[k]['onset_abs_error']) for k in shared])
                difference=eoa-eob;dist=ww@difference
                onset_contrasts.append({'candidate':ra['run_id'],'reference':rb['run_id'],'dataset':d,'prediction_mode':mode,'threshold':threshold,'all_matched_casts':len(keys),'jointly_detected_casts':len(shared),'candidate_mae_same_casts':float(eoa.mean()),'reference_mae_same_casts':float(eob.mean()),'delta_mae_same_casts':float(difference.mean()),'mae_difference_ci_low':float(np.quantile(dist,.025)),'mae_difference_ci_high':float(np.quantile(dist,.975)),'fraction_joint_casts_lower_error':float((difference<0).mean()),'interpretation':'Conditional on observed and both candidate/reference expected curves detecting onset; counts and censoring reported separately; exploratory clustered bootstrap'})
write('focused_profile_contrasts.csv',contrasts)
write('focused_onset_contrasts.csv',onset_contrasts)
write('focused_regional_contrasts.csv',regional_contrasts)
write('paired_influence_sensitivity.csv',influence)

coverage=[]
for campaign in sorted({r['short_campaign'] for r in rows}):
    rs=[r for r in rows if r['short_campaign']==campaign]
    posterior=[r for r in rs if r['mode']=='posterior' and r['status']=='complete']
    coverage.append({'campaign':campaign,'registered':len(rs),'complete_posterior_evaluations':len(posterior),'complete_prior':sum(r['mode']=='prior' and r['status']=='complete' for r in rs),'incomplete':sum(r['status']!='complete' for r in rs),'evaluation_only_reuse':sum(r['evaluation_only_reuse'] for r in rs),'convergence_pass':sum(r['review_convergence_screen'] is True for r in posterior),'zero_tree_and_convergence':sum(r['review_zero_tree_screen'] for r in posterior),'validation_scored':sum(num(r.get('validation_conditional_rmse')) is not None for r in posterior),'reserve_scored':sum(num(r.get('reserve_conditional_rmse')) is not None for r in posterior)})
write('campaign_coverage.csv',coverage)
write('onset_campaign_comparison.csv',[r for r in rows if r['short_campaign']=='onset150'])
write('prior_checks.csv',read(RES/'prior_checks.csv'))
write('all_prior_predictive_quantiles.csv',read(RES/'all_prior_predictive_quantiles.csv'))
for name in ('all_matched_cast_comparisons.csv','planned_regional_contrasts.csv','all_predictive_batch_sensitivity.csv'):
    write(name,read(RES/name))

posterior=[r for r in rows if r['mode']=='posterior' and r['status']=='complete']
unique={r['unique_fit_identity']:r for r in posterior}
summary={'registered':len(rows),'complete_registrations':sum(r['status']=='complete' for r in rows),'complete_posterior_evaluations':len(posterior),'distinct_completed_posterior_fits':len(unique),'complete_prior_fits':sum(r['status']=='complete' and r['mode']=='prior' for r in rows),'incomplete':sum(r['status']!='complete' for r in rows),'reused_control_evaluations':sum(r['evaluation_only_reuse'] for r in rows),'rmse_audits':len(audit),'maximum_rmse_discrepancy':max(r['absolute_difference'] for r in audit),'onset_score_audits':len(onset_audits),'distinct_posterior_fits_convergence_pass':sum(r['review_convergence_screen'] is True for r in unique.values()),'distinct_posterior_fits_zero_tree_pass':sum(r['review_zero_tree_screen'] for r in unique.values()),'snapshot_archive_sha256':'aa21a660aaa005ec74fd259e5459c8c0fae464210eba20195e4f699780204643'}
(OUT/'audit_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print('SUMMARY',json.dumps(summary))
print('COVERAGE',json.dumps(coverage))
print('ONSET FOCUSED COMPARISONS')
for r in onset_contrasts:
    if r['prediction_mode']=='conditional' and r['threshold']=='0.5' and r['candidate'] in ('onset150:008','onset150:009','onset150:014','onset150:015','onset150:020','onset150:021') and r['reference'] in ('onset150:001','onset150:003','onset150:005'):
        print(json.dumps(r))
print('PROFILE FOCUSED COMPARISONS')
for r in contrasts:
    if r['prediction_mode']=='conditional' and r['candidate'] in ('onset150:008','onset150:009','onset150:014','onset150:015','onset150:020','onset150:021') and r['reference'] in ('onset150:001','onset150:003','onset150:005'):
        print(json.dumps(r))
