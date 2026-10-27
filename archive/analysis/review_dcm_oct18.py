"""Audit the October 18 DCM snapshot without fitting or modifying models."""
from pathlib import Path
import csv, json, tarfile, shutil, collections, math, statistics, hashlib

ROOT=Path(__file__).resolve().parent
ARCHIVE=ROOT/'dcm_oct18_input/dcm_targeted_comparison_20261018T105731Z_3066645.tar.gz'
SNAP=ROOT/'dcm_snapshot_20261018_105731'
OUT=ROOT/'DCM_review_20261018'
OUT.mkdir(exist_ok=True)

def read(path,sep=','):
    if not path.exists() or not path.stat().st_size:return []
    with path.open(newline='',encoding='utf-8-sig') as f:return list(csv.DictReader(f,delimiter=sep))

def num(v):
    try:
        x=float(v);return x if math.isfinite(x) else None
    except (TypeError,ValueError):return None

def write(name,rows):
    cols=list(dict.fromkeys(k for r in rows for k in r))
    with (OUT/name).open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows(rows)

def extract():
    SNAP.mkdir(exist_ok=True)
    with tarfile.open(ARCHIVE,'r:gz') as t:
        members=t.getmembers()
        for m in members:
            p=(SNAP/m.name).resolve()
            if not p.is_relative_to(SNAP.resolve()) or not m.isfile():raise ValueError(m.name)
        for m in members:
            p=SNAP/m.name;p.parent.mkdir(parents=True,exist_ok=True)
            with t.extractfile(m) as src,p.open('wb') as dst:shutil.copyfileobj(src,dst)
    print('Extracted',len(members),'files')

def inventory():
    for p in sorted(SNAP.rglob('manifest.tsv')):
        if 'scheduler_backup_' in str(p):continue
        print(p.relative_to(SNAP),len(read(p,'\t')))
    for n in ('all_results.csv','case_status.csv','test_summary.csv','diagnostic_summary.csv','validation_conditional_cast_metrics.csv'):
        ps=list(SNAP.rglob(n));print(n,len(ps))
        for p in ps[:2]:
            rows=read(p);print(str(p.relative_to(SNAP)),len(rows),list(rows[0]) if rows else [])

def audit():
    # Resolve each fitted run once, even when a later campaign rescored it.
    registry={}
    for manifest in sorted(SNAP.rglob('manifest.tsv')):
        if 'scheduler_backup_' in str(manifest):continue
        base=manifest.parent
        campaign=base.name.split('_')[0]
        cases={r['case_id']:r for r in read(base/'cases.tsv','\t')}
        for record in read(manifest,'\t'):
            path=record['run_dir'].split('/ibex/user/meleg/DCM_MODEL/')[-1]
            cid=record.get('case_id',record.get('case',''))
            registry[path]={**cases.get(cid,{}),**record,'campaign':campaign,'case_id':cid,'run_id':campaign+':'+str(cid).zfill(3),'run_path':path}
    # Include unmanifested original pilot/threading fits without duplication.
    for config in SNAP.rglob('run_config.rds'):
        path=config.parent.relative_to(SNAP).as_posix()
        if '/legacy_evaluations/' in path or 'scheduler_backup_' in path:continue
        if path not in registry:
            parts=Path(path).parts
            camp=parts[1].split('_')[0] if parts[0] in ('experiments','benchmarks') else 'pilot'
            registry[path]={'campaign':camp,'case_id':parts[-1],'run_id':camp+':'+parts[-1],'run_path':path}
    targeted=next(p.parent for p in SNAP.rglob('existing_evaluations.tsv') if 'scheduler_backup_' not in str(p))
    evaluations={r['source_run'].split('/ibex/user/meleg/DCM_MODEL/')[-1]:r for r in read(targeted/'existing_evaluations.tsv','\t')}
    rows=[];tables=collections.defaultdict(list);verified=0;maxdiff=0
    for path,reg in registry.items():
        run=SNAP/path;r=reg.copy()
        for n in ('diagnostic_summary.csv','test_summary.csv'):
            s=read(run/n)
            if s:r.update(s[0])
        r['fit_completed']=(run/'FIT_COMPLETE.txt').exists()
        r['checks_completed']=(run/'CHECKS_COMPLETE.txt').exists()
        r['status']='complete' if r['checks_completed'] else 'awaiting_checks' if r['fit_completed'] else 'fit_failed' if (run/'FIT_FAILED.txt').exists() else 'incomplete'
        r['prior_only']=r.get('mode')=='prior'
        evaluation=evaluations.get(path)
        eval_run=None
        if evaluation:
            ep=evaluation['run_dir'].split('/ibex/user/meleg/DCM_MODEL/')[-1];eval_run=SNAP/ep
            r['common_evaluation_path']=ep
            r['common_evaluation_complete']=(eval_run/'CHECKS_COMPLETE.txt').exists()
            r['common_evaluation_failed']=(eval_run/'CHECKS_FAILED.txt').exists()
            es=read(eval_run/'test_summary.csv')
            if es:
                for k,v in es[0].items():
                    if k.startswith(('validation_','training_','prediction_','raw_training_','evaluation_type')):r[k]=v
        r['fit_hours']=(num(r.get('fit_seconds')) or 0)/3600 if num(r.get('fit_seconds')) is not None else None
        vals=[num(r.get(k)) for k in ('max_rhat','min_bulk_ess','min_tail_ess','divergences','min_ebfmi')]
        r['harmonized_convergence_screen']=(vals[0]<=1.01 and vals[1]>=400 and vals[2]>=400 and vals[3]==0 and vals[4]>=.3) if all(x is not None for x in vals) else None
        r['diagnostic_reason']='; '.join(k for k,fail in [('Rhat>1.01',vals[0] is not None and vals[0]>1.01),('bulkESS<400',vals[1] is not None and vals[1]<400),('tailESS<400',vals[2] is not None and vals[2]<400),('divergences',vals[3] is not None and vals[3]>0),('E-BFMI<0.3',vals[4] is not None and vals[4]<.3)] if fail)
        draws=num(r.get('retained_draws'));hits=num(r.get('treedepth_hits'))
        r['tree_hit_fraction_review']=hits/draws if draws and hits is not None else None
        r['evaluation_scope']='prior_predictive' if r['prior_only'] else 'common_holdout' if num(r.get('validation_conditional_rmse')) is not None else 'legacy_in_sample_only' if num(r.get('training_rmse')) is not None else 'no_predictive_scores'
        # An evaluation may be active while its summary exists; retain that fact.
        metric_run=eval_run if eval_run and (eval_run/'test_summary.csv').exists() else run
        for folder in ('checks','checks_and_plots'):
            for p in (metric_run/folder).glob('*.csv'):
                name=p.name
                if any(x in name for x in ('by_region','strata','cast_metrics','prior_predictive_quantiles','example_profiles','population_reference_profiles','population_expected_peak_summary','population_ocean_profiles_equal_region_weight','prediction_scores')):
                    for x in read(p):tables[name].append({'run_id':r['run_id'],**x})
            for p in (run/folder).glob('*.csv'):
                if p.name in ('sampler_summary_by_phase.csv','chain_times.csv'):
                    for x in read(p):tables[p.name].append({'run_id':r['run_id'],**x})
        # Recompute cast-weighted RMSE from every available saved cast table.
        for dataset in ('training','validation'):
            for mode in ('conditional','population'):
                cs=read(metric_run/'checks'/f'{dataset}_{mode}_cast_metrics.csv')
                reported=num(r.get(f'{dataset}_{mode}_rmse'))
                if cs and reported is not None:
                    assert len({x['cast_key'] for x in cs})==len(cs),(r['run_id'],dataset)
                    value=math.sqrt(statistics.mean(float(x['mse']) for x in cs))
                    diff=abs(value-reported);assert diff<1e-9,(r['run_id'],value,reported)
                    verified+=1;maxdiff=max(maxdiff,diff)
                    if dataset=='validation' and mode=='conditional':
                        r['validation_cast_signature_review']=hashlib.sha256('\n'.join(sorted(x['cast_key'] for x in cs)).encode()).hexdigest()
                        r['validation_casts_review']=len(cs)
        rows.append(r)
    rows.sort(key=lambda r:(r['campaign'],num(r['case_id']) or 0,str(r['case_id'])))
    assert len({r['run_id'] for r in rows})==len(rows)
    write('all_registered_runs.csv',rows)
    write('all_completed_posterior.csv',[r for r in rows if r['fit_completed'] and not r['prior_only']])
    write('all_prior_runs.csv',[r for r in rows if r['prior_only']])
    for n,rs in tables.items():write('all_'+n,rs)
    coverage=[]
    for c in sorted({r['campaign'] for r in rows}):
        rs=[r for r in rows if r['campaign']==c];ps=[r for r in rs if r['fit_completed'] and not r['prior_only']]
        coverage.append(dict(campaign=c,registered=len(rs),fits_complete=sum(r['fit_completed'] for r in rs),posterior_complete=len(ps),prior_complete=sum(r['fit_completed'] and r['prior_only'] for r in rs),common_holdout=sum(r['evaluation_scope']=='common_holdout' for r in ps),convergence_pass=sum(r['harmonized_convergence_screen'] is True for r in ps),zero_tree_and_convergence=sum(r['harmonized_convergence_screen'] is True and num(r.get('treedepth_hits'))==0 for r in ps)))
    write('campaign_coverage.csv',coverage)
    print('COVERAGE',json.dumps(coverage))
    print('VERIFIED RMSE',verified,'max error',maxdiff)
    print('COHORTS',collections.Counter((r.get('validation_md5'),r.get('validation_cast_signature_review'),r.get('validation_casts_review')) for r in rows if r['evaluation_scope']=='common_holdout'))
    print('EVALUATION REGISTRY',len(evaluations),'complete',sum(r.get('common_evaluation_complete',False) for r in rows),'failed',sum(r.get('common_evaluation_failed',False) for r in rows))
    print('TARGETED')
    fields=('case_id','equation','variance','cast','prior','subset','fit_hours','max_rhat','min_bulk_ess','min_tail_ess','divergences','treedepth_hits','harmonized_convergence_screen','validation_conditional_rmse','validation_conditional_shape_rmse','validation_conditional_peak_depth_mae','validation_conditional_coverage90','validation_conditional_interval_score90')
    for r in rows:
        if r['campaign']=='targeted150':print(' | '.join(str(r.get(k,'')) for k in fields))
    (OUT/'audit_verification.json').write_text(json.dumps({'aggregate_rmse_verified':verified,'maximum_discrepancy':maxdiff,'run_count':len(rows)},indent=2))

if __name__=='__main__':
    import sys
    if '--extract' in sys.argv:extract()
    if '--audit' in sys.argv:audit()
    else:inventory()
