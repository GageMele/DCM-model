"""Print schemas and small previews of existing CSV result tables."""
from pathlib import Path
import csv, collections, json
root = Path(__file__).resolve().parent / 'dcm_snapshot_20261021_124243'
res = next(root.glob('ibex/user/meleg/DCM_MODEL/experiments/onset*/results'))
for name in ('all_results.csv','onset_results.csv','onset_all_scores.csv','mld_neighbourhood_scores.csv','onset_planned_paired_contrasts.csv','planned_controlled_contrasts.csv','all_prior_predictive_quantiles.csv','score_audit.csv'):
    with (res/name).open(encoding='utf-8-sig',newline='') as f: rows=list(csv.DictReader(f))
    print('\nTABLE', name, 'ROWS',len(rows),'FIELDS', list(rows[0]) if rows else [])
    if name=='all_results.csv':
        print('STATUSES',collections.Counter((r['campaign'].split('_')[0],r['status']) for r in rows))
        print('MODES',collections.Counter(r.get('mode') for r in rows))
    elif name=='onset_all_scores.csv':
        print('GROUPS',collections.Counter((r.get('dataset'),r.get('mode'),r.get('method_id'),r.get('metric')) for r in rows))
        for r in rows[:3]: print(json.dumps(r))
    else:
        for r in rows[:2]: print(json.dumps(r))
for p in root.glob('experiments/onset*/*'):
    if p.is_file(): print('CAMPAIGN FILE',p.name,p.stat().st_size)
