"""Preserve scientific diagnostics and explain influential dates from saved tables."""
from pathlib import Path
import csv,collections,json,math,statistics,re
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'DCM_review_20261019'
SNAP=ROOT/'dcm_snapshot_20261019_084423'
EXP=SNAP/'experiments/followup150_20261018T150024Z_3898749'
def read(p,sep=','):
    if not p.exists():return []
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter=sep))
def write(n,rs):
    keys=list(dict.fromkeys(k for r in rs for k in r))
    with (OUT/n).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,keys);w.writeheader();w.writerows(rs)
rows=read(OUT/'all_registered_runs.csv');by={r['run_id']:r for r in rows};ids={r['run_key']:r['run_id'] for r in rows}

# Independently reconcile the collector's registry with campaign manifests.
recorded=set()
manifest_provenance=[]
for manifest in SNAP.rglob('manifest.tsv'):
    if 'scheduler_backup_' in str(manifest):continue
    records=read(manifest,'\t')
    manifest_provenance.append({'manifest':str(manifest.relative_to(SNAP)),'source':'current snapshot','registrations':len(records)})
    for r in records:
        if r.get('run_dir'):recorded.add(r['run_dir'])
current_manifest_count=len(recorded)
OLD=ROOT/'dcm_snapshot_20261018_105731'
for manifest in OLD.rglob('manifest.tsv'):
    if 'scheduler_backup_' in str(manifest) or (SNAP/manifest.relative_to(OLD)).exists():continue
    records=read(manifest,'\t')
    manifest_provenance.append({'manifest':str(manifest.relative_to(OLD)),'source':'previous snapshot: frozen campaign registration only','registrations':len(records)})
    for r in records:
        if r.get('run_dir'):recorded.add(r['run_dir'])
listed={r['run_dir'] for r in rows}
assert len(recorded)==269 and recorded<=listed,(len(recorded),len(recorded-listed))
assert len(listed-recorded)==1
write('registry_reconciliation.csv',[{'current_manifest_registered':current_manifest_count,'earlier_frozen_manifest_registered':len(recorded)-current_manifest_count,'additional_pilot':len(listed-recorded),'total_unique_runs':len(listed),'passed':True}])
write('manifest_provenance.csv',manifest_provenance)

# Retain all scientific coefficient diagnostics, with the worst parameters too.
parameters=[];worst=[];priors=[];phases=[];times=[]
for r in rows:
    p=SNAP/r['run_dir'].split('/ibex/user/meleg/DCM_MODEL/')[-1]
    provenance='current snapshot'
    ps=read(p/'checks/parameter_diagnostics.csv') or read(p/'checks_and_plots/parameter_diagnostics.csv')
    if not ps:
        old=OLD/r['run_dir'].split('/ibex/user/meleg/DCM_MODEL/')[-1]
        ps=read(old/'checks/parameter_diagnostics.csv') or read(old/'checks_and_plots/parameter_diagnostics.csv')
        provenance='previous snapshot: unchanged completed fit'
    for s in ps:
        v=s.get('variable','')
        if (v.startswith(('b_','sd_','sds_')) and not v.startswith('b_Intercept_')) or v in ('sigma','nu','b_Intercept'):parameters.append({'run_id':r['run_id'],'evidence_source':provenance,**s})
    def value(s,k):
        try:return float(s.get(k,''))
        except (TypeError,ValueError):return -math.inf if k=='rhat' else math.inf
    if ps:
        for k,rev in [('rhat',True),('ess_bulk',False),('ess_tail',False)]:
            for s in sorted(ps,key=lambda s:value(s,k),reverse=rev)[:3]:worst.append({'run_id':r['run_id'],'criterion':k,'evidence_source':provenance,**s})
    priors.extend({'run_id':r['run_id'],**s} for s in read(p/'resolved_priors.csv'))
    for folder in ('checks','checks_and_plots'):
        phases.extend({'run_id':r['run_id'],**s} for s in read(p/folder/'sampler_summary_by_phase.csv'))
        times.extend({'run_id':r['run_id'],**s} for s in read(p/folder/'chain_times.csv'))
write('parameters.csv',parameters);write('worst_parameters.csv',worst);write('resolved_priors.csv',priors);write('sampler_phases.csv',phases);write('chain_times.csv',times)

splits=read(EXP/'split_membership.csv');train={n:{(r['Ocean'],r['Region'],r['date']) for r in splits if r['dataset']==n} for n in ('primary_p10_d150','primary_p10_alt_d150')}
influence=read(OUT/'cast_influence.csv')
date_rows=[]
for dataset in ('validation','reserve'):
    ca={r['cast_key']:r for r in influence if r['run_id']=='followup150:027' and r['dataset']==dataset}
    cb={r['cast_key']:r for r in influence if r['run_id']=='followup150:005' and r['dataset']==dataset}
    dg=collections.defaultdict(list)
    for k,c in cb.items():dg[(c['Ocean'],c['Region'],c['date'])].append(k)
    total=sum(float(ca[k]['mse'])-float(cb[k]['mse']) for k in cb)
    for d,keys in dg.items():
        delta=sum(float(ca[k]['mse'])-float(cb[k]['mse']) for k in keys)
        date_rows.append({'dataset':dataset,'Ocean':d[0],'Region':d[1],'date':d[2],'casts':len(keys),'seen_date_original_training':d in train['primary_p10_d150'],'seen_date_alternative_training':d in train['primary_p10_alt_d150'],'sum_mse_change_alternative_minus_original':delta,'contribution_to_overall_mse_change':delta/len(cb),'fraction_net_mse_change':delta/total,'total_casts':len(cb)})
write('training_composition_date_contributions.csv',date_rows)

# Classification agreement is not specificity if every held-out profile has a DCM.
confusion=[]
for dataset in ('validation','reserve'):
    cs=read(OUT/f'all_{dataset}_conditional_cast_metrics.csv');groups=collections.defaultdict(list)
    for c in cs:groups[ids[c['run_key']]].append(c)
    for rid,rs in groups.items():
        z={'run_id':rid,'dataset':dataset,'casts':len(rs)}
        for o,p,label in [(True,True,'true_positive'),(True,False,'false_negative'),(False,True,'false_positive'),(False,False,'true_negative')]:z[label]=sum((c['observed_dcm_in_range']=='TRUE')==o and (c['predicted_dcm_in_range']=='TRUE')==p for c in rs)
        z['observed_negative_casts']=z['false_positive']+z['true_negative'];confusion.append(z)
write('dcm_confusion.csv',confusion)

# Pointwise reference profiles are population means at recorded regional settings.
profiles=[];examples=[]
for rid in ('followup150:001','followup150:003','followup150:005','followup150:027'):
    r=by[rid];run=SNAP/r['run_dir'].split('/ibex/user/meleg/DCM_MODEL/')[-1]
    profiles.extend({'run_id':rid,**p} for p in read(run/'checks/population_reference_profiles.csv'))
    for dataset in ('validation','reserve'):
        examples.extend({'run_id':rid,'dataset':dataset,**p} for p in read(run/f'checks/{dataset}_conditional_example_profiles.csv'))
write('selected_population_reference_profiles.csv',profiles);write('selected_observed_predicted_profiles.csv',examples)

controls=[]
for r in rows:
    if r['short_campaign'] not in ('targeted150','followup150'):continue
    ref=r.get('reference_case')
    if ref in ('NA','',None):continue
    base=by.get(r['short_campaign']+':'+str(int(float(ref))).zfill(3))
    if not base:continue
    keys=('equation','variance','cast','family','autocor','subset','split','prior','depth_max','chains','warmup','samples','adapt_delta','treedepth','threads','k_depth','k_mld','k_distance','sigma_smooth_sd','seed')
    controls.append({'run_a':r['run_id'],'run_b':base['run_id'],'block':r.get('block'),'status_a':r['status'],'status_b':base['status'],'changed_settings':'; '.join(f'{k}: {base.get(k)} -> {r.get(k)}' for k in keys if r.get(k)!=base.get(k)),'same_training_md5':r.get('data_md5')==base.get('data_md5'),'same_cpu_model':r.get('cpu_model')==base.get('cpu_model')})
write('all_planned_control_settings.csv',controls)
print(json.dumps({'parameters':len(parameters),'priors':len(priors),'date_contribution_rows':len(date_rows),'profiles':len(profiles),'example_rows':len(examples),'planned_controls':len(controls)}))
print('MOST INFLUENTIAL COMPOSITION DATES',json.dumps(sorted([r for r in date_rows if r['dataset']=='reserve'],key=lambda r:r['sum_mse_change_alternative_minus_original'])[:4]))
