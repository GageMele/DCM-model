"""Summarize every model family and verify saved pairwise score differences."""
from pathlib import Path
import csv,collections,statistics,math,json
ROOT=Path(__file__).resolve().parent; OUT=ROOT/'DCM_review_20261021'
def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def num(x):
    try:
        v=float(x);return v if math.isfinite(v) else None
    except (ValueError,TypeError):return None
def write(name,rs):
    rs=list(rs); fields=list(dict.fromkeys(k for r in rs for k in r))
    with (OUT/name).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fields);w.writeheader();w.writerows(rs)
rows=read(OUT/'all_registered_runs.csv'); by={r['run_key']:r for r in rows}
post=[r for r in rows if r['mode']=='posterior' and r['status']=='complete']
groups=collections.defaultdict(list)
for r in post:groups[(r['equation'],r.get('variance','NA'),r.get('family','NA'))].append(r)
families=[]
for (eq,variance,family),rs in sorted(groups.items()):
    errors=[num(r.get('validation_conditional_rmse')) for r in rs]
    errors=[v for v in errors if v is not None]
    families.append({'equation':eq,'variance':variance,'family':family,'completed_evaluations':len(rs),'distinct_completed_fits':len({r['unique_fit_identity'] for r in rs}),'convergence_screen_pass':sum(r['review_convergence_screen']=='True' for r in rs),'zero_tree_screen_pass':sum(r['review_zero_tree_screen']=='True' for r in rs),'observed_validation_rmse_min':min(errors) if errors else None,'observed_validation_rmse_max':max(errors) if errors else None,'interpretation':'Descriptive range across settings/subsets; different cohorts are not ranked against each other'})
write('model_family_summary.csv',families)

# The newest export already contains all registered pairwise comparisons.
# Compare the stored aggregate differences with the source model summaries.
pair=read(OUT/'all_matched_cast_comparisons.csv')
print('PAIR FIELDS',list(pair[0]))
print('PAIR EXAMPLE',json.dumps(pair[0]))

print('TOP PRIMARY SCORES (descriptive only, eligible training >=10%, unique fits)')
rank=[]; seen=set()
primary=next(r['validation_conditional_cohort_signature'] for r in post if r['run_id']=='onset150:001')
for r in sorted(post,key=lambda r:num(r.get('validation_conditional_rmse')) or float('inf')):
    if r.get('validation_conditional_cohort_signature')!=primary or r.get('subset') not in ('p10','p10_alt','p20','p50'):continue
    if r['unique_fit_identity'] in seen:continue
    seen.add(r['unique_fit_identity']);rank.append(r)
for r in rank[:14]:print(' | '.join(str(r.get(k)) for k in ('run_id','equation','subset','variance','validation_conditional_rmse','reserve_conditional_rmse','review_diagnostic_flags')))
print('HISTORICAL SMALL-SUBSET COUNTS',collections.Counter(r['subset'] for r in rows))

# Retain exact primary and reserve error contributions, including the extreme
# casts: these are influence diagnostics, not reasons to omit observations.
cast=read(OUT/'all_reserve_conditional_cast_metrics.csv')
selected={r['run_key'] for r in rows if r['run_id'] in ('onset150:001','onset150:008','onset150:009','onset150:014','onset150:015')}
influence=[]
for key in selected:
    cs=[c for c in cast if c['run_key']==key];total=sum(float(c['mse']) for c in cs)
    for i,c in enumerate(sorted(cs,key=lambda c:float(c['mse']),reverse=True),1):
        influence.append({'run_id':by[key]['run_id'],'error_rank':i,**c,'share_of_total_cast_mse':float(c['mse'])/total})
write('reserve_cast_influence.csv',influence)
print('INFLUENCE FIRST',json.dumps([r for r in influence if r['error_rank']==1]))

# Keep the original, full-data 250 m fit as a historical comparator. Its
# in-training example errors do not enter held-out model rankings.
historic=read(ROOT/'original_model_review/historical_summary.csv')
for r in historic:r.update({'run_id':'historical_original_full_250m','k_depth':15,'adapt_delta':.97,'max_treedepth':14,'training_membership':'Not reconstructed for an independent held-out comparison','comparison_role':'Historical diagnostic/profile reference, separate from campaign registry'})
write('historical_original_model.csv',historic)
