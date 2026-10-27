"""Independently verify every saved pair on the same individual casts."""
from pathlib import Path
import csv,collections,math,json,statistics
root=Path(__file__).resolve().parent/'DCM_review_20261021'
def read(name):
    with (root/name).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
tables={}
for d in ('validation','reserve'):
    groups=collections.defaultdict(dict)
    for r in read('all_'+d+'_conditional_cast_metrics.csv'):groups[r['run_key']][r['cast_key']]=r
    tables[d]=groups
max_error=0.;count=0
for r in read('all_matched_cast_comparisons.csv'):
    a,b=tables[r['dataset']][r['run_a']],tables[r['dataset']][r['run_b']]
    assert a.keys()==b.keys() and len(a)==int(r['matched_casts'])
    keys=list(a)
    ea=[float(a[k]['mse']) for k in keys];eb=[float(b[k]['mse']) for k in keys]
    vals={'delta_mean_cast_mse':statistics.mean(x-y for x,y in zip(ea,eb)),
          'delta_rmse':math.sqrt(statistics.mean(ea))-math.sqrt(statistics.mean(eb)),
          'fraction_casts_a_better':sum(x<y for x,y in zip(ea,eb))/len(keys),
          'delta_shape_rmse':math.sqrt(statistics.mean(float(a[k]['shape_mse']) for k in keys))-math.sqrt(statistics.mean(float(b[k]['shape_mse']) for k in keys))}
    for label,field in [('delta_mae','mae'),('delta_peak_depth_mae','peak_depth_abs_error'),('delta_interval_score90','interval_score90')]:
        vals[label]=statistics.mean(float(a[k][field])-float(b[k][field]) for k in keys)
    for k,v in vals.items():
        diff=abs(v-float(r[k]));max_error=max(diff,max_error)
        assert diff<1e-8,(r,k,v)
    count+=1
result={'saved_pairwise_comparisons_verified':count,'maximum_pairwise_discrepancy':max_error}
(root/'pairwise_audit.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result))
