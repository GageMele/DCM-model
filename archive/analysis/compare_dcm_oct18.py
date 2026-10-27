"""Build matched comparisons, retain all cases, and quantify regional tradeoffs."""
from review_dcm_oct18 import ROOT,SNAP,OUT,read,num,write
import collections,math,statistics,json,itertools,csv

rows=read(OUT/'all_registered_runs.csv')
configs={r['run_path']:r for r in read(OUT/'config_metadata.csv')}
for r in rows:
    for k,v in configs.get(r['run_path'],{}).items():
        if r.get(k) in (None,'','NA') and v not in ('','NA'):r[k]=v
    if not r.get('treedepth'):r['treedepth']=r.get('max_treedepth')
    r['comparison_cohort']='primary_md5_verified' if r.get('validation_md5')=='950184ac377851361fd4af891249ca22' else 'date_blocked_md5_verified' if r.get('validation_md5')=='7a1cb728145f25ada97e9cfdcad8ea98' else 'primary_membership_only' if r.get('validation_casts_review')=='108' else 'not_scored'
by={r['run_id']:r for r in rows}
casts=collections.defaultdict(dict)
for c in read(OUT/'all_validation_conditional_cast_metrics.csv'):casts[c['run_id']][c['cast_key']]=c
regional=collections.defaultdict(list)
for c in read(OUT/'all_validation_conditional_by_region.csv'):regional[c['run_id']].append(c)
for r in rows:
    cs=casts[r['run_id']];rg=regional[r['run_id']]
    if cs:
        errors=sorted([float(c['mse']) for c in cs.values()],reverse=True)
        r['top_one_cast_mse_fraction']=errors[0]/sum(errors)
        r['top_five_cast_mse_fraction']=sum(errors[:5])/sum(errors)
    if rg:
        r['macro_region_rmse']=math.sqrt(statistics.mean(float(x['rmse'])**2 for x in rg))
        r['macro_region_coverage90']=statistics.mean(float(x['coverage90']) for x in rg)
        r['worst_region_rmse']=max(float(x['rmse']) for x in rg)
        r['min_region_coverage90']=min(float(x['coverage90']) for x in rg)
write('all_registered_runs.csv',rows)
post=[r for r in rows if r['fit_completed']=='True' and r['prior_only']=='False']
write('all_completed_posterior.csv',post)

metrics=('rmse','mae','bias','shape_rmse','peak_depth_mae','peak_height_mae','coverage90','interval_width90','interval_score90','positive_auc_mae')
def paired(a,b):
    ra,rb=by[a],by[b];ca,cb=casts[a],casts[b]
    if not ca or ca.keys()!=cb.keys():return {}
    ha,hb=ra.get('validation_md5'),rb.get('validation_md5')
    strict=bool(ha and hb and ha==hb)
    delta=[float(ca[k]['mse'])-float(cb[k]['mse']) for k in ca]
    out={'run_a':a,'run_b':b,'matched_casts':len(ca),'validation_md5_matches':strict,
      'cohort_a':ra['comparison_cohort'],'cohort_b':rb['comparison_cohort'],
      'a_convergence':ra['harmonized_convergence_screen'],'b_convergence':rb['harmonized_convergence_screen'],
      'delta_mean_cast_mse':statistics.mean(delta),'descriptive_paired_se':statistics.stdev(delta)/math.sqrt(len(delta)),
      'fraction_casts_a_better':sum(x<0 for x in delta)/len(delta)}
    for m in metrics:
        x,y=num(ra.get('validation_conditional_'+m)),num(rb.get('validation_conditional_'+m))
        out['delta_'+m]=x-y if x is not None and y is not None else None
    out['delta_macro_region_rmse']=(num(ra.get('macro_region_rmse'))-num(rb.get('macro_region_rmse'))) if num(ra.get('macro_region_rmse')) is not None and num(rb.get('macro_region_rmse')) is not None else None
    return out

# Every pair with identical validation-data hash; no accuracy-based selection.
groups=collections.defaultdict(list)
for r in post:
    if r.get('validation_md5') and casts[r['run_id']]:groups[r['validation_md5']].append(r['run_id'])
pairs=[]
for ids in groups.values():
    for a,b in itertools.combinations(ids,2):
        p=paired(a,b)
        if p:pairs.append(p)
write('all_matched_pairs.csv',pairs)

controlled=[]
controls=('equation','variance','cast','family','autocor','prior','subset','split','depth_max','k_depth','k_mld','k_distance','warmup','samples','adapt_delta','treedepth','threads','seed')
for r in rows:
    if r['campaign']!='targeted150' or r['prior_only']=='True':continue
    ref=num(r.get('reference_case'))
    refid='targeted150:'+str(int(ref)).zfill(3) if ref else 'baseline150:001' if r['case_id']=='61' else None
    if not refid:continue
    b=by[refid]
    out={'run_a':r['run_id'],'run_b':refid,'block':r.get('block'),'status_a':r['status'],'status_b':b['status'],
      'changed_settings':'; '.join(k+': '+str(b.get(k))+' -> '+str(r.get(k)) for k in controls if r.get(k)!=b.get(k)),
      'same_cpu_model':r.get('cpu_model')==b.get('cpu_model'),'same_training_data_md5':r.get('data_md5')==b.get('data_md5')}
    out.update(paired(r['run_id'],refid))
    for m in ('fit_hours','requested_core_hours','min_bulk_ess','max_rhat','divergences','treedepth_hits'):
        x,y=num(r.get(m)),num(b.get(m));out['delta_'+m]=x-y if x is not None and y is not None else None
    controlled.append(out)
write('targeted_controlled_contrasts.csv',controlled)

# Compare every primary evaluated model with the historical-equation baseline.
reference=[]
for r in post:
    p=paired(r['run_id'],'baseline150:001')
    if p:reference.append(p)
write('all_vs_original_equation_baseline.csv',reference)

# Explicit hypothesis panels retain every member rather than selecting by score.
panels={
 'residual_structure':range(1,13), 'mean_equations':range(13,19),
 'depth_basis':(2,8,19,20,21,22), 'priors':(2,8,23,24,25,26,27,28),
 'sampler_additive':(2,29,30,31,32,33,34,35,36,37),
 'sampler_both':(8,38,39,40,41,42,43,44,45,46),
 'training_size':(2,8,47,48,49,50), 'date_generalization':(2,8,51,52),
 'residual_dependence':(59,60,62,63), 'historical_depth':(61,)}
panelrows=[]
for panel,ids in panels.items():
    for cid in ids:panelrows.append({'question_panel':panel,**by['targeted150:'+str(cid).zfill(3)]})
write('targeted_question_panels.csv',panelrows)

# Retain scientifically useful parameters and the worst diagnostic per fit.
parameters=[];worst=[];priors=[]
for r in rows:
    run=SNAP/r['run_path']
    p=run/'checks/parameter_diagnostics.csv'
    if not p.exists():p=run/'checks_and_plots/parameter_diagnostics.csv'
    ds=read(p)
    eligible=[x for x in ds if num(x.get('rhat')) is not None]
    if eligible:
        for x in sorted(eligible,key=lambda x:float(x['rhat']),reverse=True)[:3]:worst.append({'run_id':r['run_id'],**x})
    for x in ds:
        v=x.get('variable','')
        if v.startswith(('b_','sd_','sds_','ar[')) or v in ('sigma','nu'):
            parameters.append({'run_id':r['run_id'],**x})
    for x in read(run/'resolved_priors.csv'):priors.append({'run_id':r['run_id'],**x})
write('parameters.csv',parameters);write('worst_parameter_diagnostics.csv',worst);write('resolved_priors.csv',priors)

print('MATCHED PAIRS',len(pairs),'MAIN COHORT',len(groups.get('950184ac377851361fd4af891249ca22',[])))
for panel,ids in panels.items():
    print('PANEL',panel)
    for cid in ids:
        r=by['targeted150:'+str(cid).zfill(3)]
        print(cid,*[f'{k}={r.get(k)}' for k in ('fit_hours','requested_core_hours','macro_region_rmse','min_region_coverage90','top_five_cast_mse_fraction')])
print('KEY PAIRED DIFFERENCES')
for a,b in [(2,1),(8,7),(18,8),(13,2),(16,2),(17,8),(47,2),(49,8),(59,62),(60,63),(61,0)]:
    aid='targeted150:'+str(a).zfill(3);bid='targeted150:'+str(b).zfill(3) if b else 'baseline150:001'
    print(json.dumps(paired(aid,bid)))
print('REGIONAL COVERAGE / ERROR')
for cid in (1,2,3,8,13,18,47,49):
    print(cid,[(x['Region'],x['casts'],round(float(x['rmse']),3),round(float(x['coverage90']),3)) for x in regional['targeted150:'+str(cid).zfill(3)]])
