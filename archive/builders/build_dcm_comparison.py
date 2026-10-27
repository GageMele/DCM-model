"""Compare every recorded DCM run without error-based selection.

Uses saved summaries and cast metrics only. No fitting or new prediction.
"""
from pathlib import Path
import csv, json, math, statistics, collections, hashlib, html

ROOT = Path(__file__).resolve().parent
SNAP = ROOT / 'dcm_snapshot_20261016_090217'
OUT = ROOT / 'DCM_all_run_comparison_20261016'
OUT.mkdir(exist_ok=True)
CONFIGS = json.loads((ROOT / 'snapshot_configs.json').read_text(encoding='utf-8'))

def read_csv(path, sep=','):
    if not path.exists(): return []
    with path.open(newline='', encoding='utf-8-sig') as handle:
        return list(csv.DictReader(handle, delimiter=sep))

def num(value):
    try:
        v = float(value)
        return v if math.isfinite(v) else None
    except (ValueError, TypeError): return None

def yes(value):
    return value is True or str(value).upper() == 'TRUE'

def write_csv(name, rows, first=()):
    cols = list(first)
    for row in rows:
        for key in row:
            if key not in cols: cols.append(key)
    with (OUT / name).open('w', newline='', encoding='utf-8-sig') as handle:
        w = csv.DictWriter(handle, fieldnames=cols, extrasaction='ignore')
        w.writeheader(); w.writerows(rows)

def short_campaign(path):
    s = path.parts
    return s[1].split('_')[0] if s[0] in ('experiments', 'benchmarks') else 'pilot'

def local_path(value):
    marker = '/ibex/user/meleg/DCM_MODEL/'
    return value.split(marker, 1)[-1] if marker in value else value

def coerce(rows):
    for row in rows:
        for key, val in list(row.items()):
            if val in ('NA', 'NaN', 'Inf', '-Inf', ''): row[key] = None
            elif str(val).upper() in ('TRUE', 'FALSE'): row[key] = yes(val)
            elif isinstance(val, str) and num(val) is not None:
                v = num(val); row[key] = int(v) if v.is_integer() else v
    return rows

def uid(campaign, case):
    return campaign + ':' + (f'{int(case):03d}' if num(case) is not None else str(case))

registry = {}
for manifest in sorted(SNAP.rglob('manifest.tsv')):
    campaign = short_campaign(manifest.relative_to(SNAP))
    cases = {str(x['case_id']): x for x in read_csv(manifest.parent / 'cases.tsv', '\t')}
    for record in read_csv(manifest, '\t'):
        relative = local_path(record['run_dir'])
        case = record.get('case_id', record.get('case'))
        row = dict(cases.get(str(case), {})); row.update(record)
        row.update(campaign=campaign, case_id=case, run_path=relative, run_id=uid(campaign, case))
        registry[relative] = row
for config_path, config in CONFIGS.items():
    relative = config_path.split(SNAP.name + '/', 1)[-1]
    if relative not in registry:
        p = Path(relative); campaign = short_campaign(p)
        registry[relative] = dict(campaign=campaign, case_id=p.name, run_path=relative,
                                  run_id=uid(campaign, p.name))
    row = registry[relative]
    for key, value in config.items():
        if key == 'case' and isinstance(value, list) and value:
            for k, v in value[0].items(): row[k] = v
        elif key not in ('case', 'mode') and not isinstance(value, (dict, list)):
            row[key] = value
        elif key == 'k': row['configured_k'] = str(value)
    row['config_mode'] = config.get('mode')

all_rows = []
for relative, row in registry.items():
    run = SNAP / relative
    for filename in ('diagnostic_summary.csv', 'test_summary.csv'):
        rows = read_csv(run / filename)
        if rows: row.update(rows[0])
    fit = (run / 'FIT_COMPLETE.txt').exists()
    checks = (run / 'CHECKS_COMPLETE.txt').exists()
    row['fit_completed'] = fit; row['checks_completed'] = checks
    row['status'] = 'complete' if fit and checks else 'fit_complete_checks_unfinished' if fit else 'not_complete_in_snapshot'
    row['prior_only'] = row.get('mode') == 'prior' or row.get('config_mode') == 'prior'
    row['evaluation_scope'] = ('prior_predictive' if row['prior_only'] else
        'common_108_cast_validation' if 'validation_conditional_rmse' in row else
        'bounded_training_check' if 'training_rmse' in row else 'no_predictive_scores')
    row['fit_hours'] = num(row.get('fit_seconds')) / 3600 if num(row.get('fit_seconds')) is not None else None
    row['configured_retained_per_chain'] = (num(row['iter']) - num(row['warmup'])) if num(row.get('iter')) is not None and num(row.get('warmup')) is not None else None
    # Recompute the same convergence screen across historical and new summaries.
    required = ('max_rhat', 'min_bulk_ess', 'min_tail_ess', 'divergences', 'min_ebfmi')
    vals = {k: num(row.get(k)) for k in required}
    row['harmonized_convergence_screen'] = (vals['max_rhat'] <= 1.01 and vals['min_bulk_ess'] >= 400
        and vals['min_tail_ess'] >= 400 and vals['divergences'] == 0 and vals['min_ebfmi'] >= .3
        ) if all(v is not None for v in vals.values()) else None
    hits = num(row.get('treedepth_hits')); draws = num(row.get('retained_draws'))
    row['harmonized_tree_hit_fraction'] = hits / draws if hits is not None and draws else None
    row['harmonized_screen_and_zero_tree_hits'] = (row['harmonized_convergence_screen'] and hits == 0
        ) if row['harmonized_convergence_screen'] is not None and hits is not None else None
    prior_file = run / 'resolved_priors.csv'
    row['resolved_prior_file_md5_review'] = hashlib.md5(prior_file.read_bytes()).hexdigest() if prior_file.exists() else None
    all_rows.append(row)
all_rows = coerce(all_rows)
all_rows.sort(key=lambda x: (x['campaign'], num(x['case_id']) or 0, str(x['case_id'])))
BY_ID = {x['run_id']: x for x in all_rows}
FIRST = ('run_id', 'campaign', 'case_id', 'label', 'status', 'evaluation_scope')
completed = [r for r in all_rows if r['fit_completed']]
post = [r for r in completed if not r['prior_only']]
ext = [r for r in post if r['campaign'] == 'extensive150']
write_csv('01_all_registered_runs.csv', all_rows, FIRST)
write_csv('02_all_completed_posterior_runs.csv', post, FIRST)
write_csv('03_prior_predictive_runs.csv', [r for r in all_rows if r['prior_only']], FIRST)

# Keep every row of every available evaluation/diagnostic subgroup table.
long_tables = collections.defaultdict(list)
for row in all_rows:
    run = SNAP / row['run_path']
    for folder in ('checks', 'checks_and_plots'):
        for path in (run / folder).glob('*.csv'):
            name = path.name
            if any(token in name for token in ('by_region', 'strata', 'cast_metrics', 'prediction_scores',
                'prior_predictive_quantiles', 'sampler_summary_by_phase', 'chain_times')):
                for values in read_csv(path):
                    values.update(run_id=row['run_id'], campaign=row['campaign'], case_id=row['case_id'])
                    long_tables[name].append(values)
for name, rows in long_tables.items(): write_csv('all_' + name, rows, ('run_id', 'campaign', 'case_id'))

counts = []
for campaign in sorted({r['campaign'] for r in all_rows}):
    rs = [r for r in all_rows if r['campaign'] == campaign]
    ps = [r for r in rs if r['fit_completed'] and not r['prior_only']]
    counts.append(dict(campaign=campaign, registered=len(rs), complete=sum(r['fit_completed'] for r in rs),
        complete_posterior=len(ps), prior_complete=sum(r['fit_completed'] and r['prior_only'] for r in rs),
        convergence_screen_pass=sum(r['harmonized_convergence_screen'] is True for r in ps),
        zero_tree_screen_pass=sum(r['harmonized_screen_and_zero_tree_hits'] is True for r in ps),
        posterior_with_validation=sum(r['evaluation_scope'] == 'common_108_cast_validation' for r in ps)))
write_csv('00_campaign_coverage.csv', counts)

# Cast-paired comparisons for all available extensive posterior cases vs the original k=15.
def cast_data(case, dataset='validation', mode='conditional'):
    row = BY_ID.get(uid('extensive150', case))
    if not row: return {}
    rows = coerce(read_csv(SNAP / row['run_path'] / 'checks' / f'{dataset}_{mode}_cast_metrics.csv'))
    return {r['cast_key']: r for r in rows}

pairwise = []
def paired(case, ref, dataset='validation', mode='conditional'):
    a, b = cast_data(case, dataset, mode), cast_data(ref, dataset, mode)
    if not a or a.keys() != b.keys(): return None
    delta = [a[k]['mse'] - b[k]['mse'] for k in a]
    return dict(case_id=case, reference_case=ref, dataset=dataset, prediction_mode=mode,
        matched_casts=len(delta), delta_mean_cast_mse=statistics.mean(delta),
        descriptive_paired_se=statistics.stdev(delta)/math.sqrt(len(delta)),
        delta_rmse=math.sqrt(statistics.mean(r['mse'] for r in a.values()))-math.sqrt(statistics.mean(r['mse'] for r in b.values())),
        fraction_casts_lower_mse=sum(d < 0 for d in delta)/len(delta),
        case_convergence=BY_ID[uid('extensive150',case)]['harmonized_convergence_screen'],
        reference_convergence=BY_ID[uid('extensive150',ref)]['harmonized_convergence_screen'])
for row in ext:
    for mode in ('conditional','population'):
        result = paired(row['case_id'], 4, mode=mode)
        if result: pairwise.append(result)
write_csv('04_every_case_vs_original_k15.csv', pairwise)

# Explicit within-campaign contrasts; retain every completed member, including diagnostic failures.
contrasts = []
fields = ('equation','subset','prior','mode','family','chains','iter','warmup','adapt_delta',
          'treedepth','threads','k_depth','k_mld','k_distance','seed')
metrics = ('fit_seconds','requested_core_hours','min_bulk_ess','max_rhat','divergences','treedepth_hits',
    'validation_conditional_rmse','validation_population_rmse','validation_conditional_shape_rmse',
    'validation_conditional_peak_depth_mae','validation_conditional_coverage90',
    'validation_conditional_interval_width90','training_conditional_rmse','training_population_rmse')
for row in ext:
    if row.get('equation') not in ('additive','both_scaled'): continue
    ref = BY_ID[uid('extensive150',6 if row['equation']=='additive' else 49)]
    if row['case_id'] == ref['case_id']: continue
    rec = dict(run_id=row['run_id'], reference_run=ref['run_id'], block=row['block'],
        changed_settings='; '.join(f'{k}: {ref.get(k)} -> {row.get(k)}' for k in fields if row.get(k)!=ref.get(k)),
        same_cpu_model=row.get('cpu_model')==ref.get('cpu_model'), same_training_data=row.get('data_md5')==ref.get('data_md5'),
        case_convergence=row['harmonized_convergence_screen'], reference_convergence=ref['harmonized_convergence_screen'])
    for metric in metrics:
        a,b = num(row.get(metric)),num(ref.get(metric))
        rec['delta_'+metric] = a-b if a is not None and b is not None else None
    p = paired(row['case_id'],ref['case_id'])
    if p: rec.update({k:v for k,v in p.items() if k in ('delta_mean_cast_mse','descriptive_paired_se','fraction_casts_lower_mse')})
    contrasts.append(rec)
write_csv('05_extensive_anchor_contrasts.csv', contrasts)

def show(rows, fields):
    for r in rows:
        print(' | '.join(f'{k}={r.get(k)}' for k in fields))

print('COVERAGE'); show(counts, counts[0].keys())
print('ALL EXTENSIVE POSTERIOR RESULTS')
show(ext, ('case_id','equation','block','subset','prior','k_depth','fit_hours','max_rhat','min_bulk_ess',
    'divergences','treedepth_hits','harmonized_convergence_screen','validation_conditional_rmse',
    'validation_population_rmse','validation_conditional_shape_rmse','validation_conditional_peak_depth_mae',
    'validation_conditional_coverage90','validation_conditional_interval_width90','validation_conditional_median_residual_lag1'))
print('LEGACY COMPLETED POSTERIOR RESULTS')
show([r for r in post if r['campaign']!='extensive150'], ('run_id','label','equation','subset','prior','fit_hours',
    'max_rhat','min_bulk_ess','divergences','treedepth_hits','harmonized_convergence_screen','training_rmse','training_ppc_95_coverage'))
print('TOTAL',len(all_rows),'REGISTERED;',len(completed),'COMPLETE;',len(post),'POSTERIOR;',len(ext),'VALIDATION')
print('SAMPLE TABLE SCHEMAS')
for name, rows in long_tables.items(): print(name,len(rows),list(rows[0]) if rows else [])
