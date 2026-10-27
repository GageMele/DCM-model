"""Audit the uploaded DCM snapshot and produce reproducible comparison tables.

Usage: python analyze_shape_review.py EXTRACTED_EXPERIMENT OUTPUT_DIRECTORY
Requires numpy, pandas and matplotlib. Does not fit or modify any model.
"""
from pathlib import Path
import sys, json, hashlib, zipfile
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root = Path(sys.argv[1])
out = Path(sys.argv[2]); out.mkdir(parents=True, exist_ok=True)
s = pd.read_csv(root/'results/shape_results.csv')
a = pd.read_csv(root/'results/all_results.csv')
cs = pd.read_csv(root/'results/all_shape_cast_metrics.csv')
ss = pd.read_csv(root/'results/all_shape_shape_cast_metrics.csv')
complete = s[s.status.eq('complete') & s['mode'].eq('posterior')].copy()
assert a.run_key.is_unique and s.case_id.is_unique
assert len(complete) == 67
# Recalculate all available equal-cast aggregate metrics, without changing data.
errors=[]
for (cid, dataset, mode), d in cs.groupby(['case_id','dataset','prediction_mode']):
    row=s.set_index('case_id').loc[cid]
    assert d.cast_key.is_unique
    metrics={'rmse':np.sqrt(d.mse.mean()), 'shape_rmse':np.sqrt(d.shape_mse.mean()),
             'mae':d.mae.mean(), 'peak_depth_mae':d.peak_depth_abs_error.mean(),
             'coverage90':d.coverage90.mean(), 'interval_score90':d.interval_score90.mean()}
    for name,value in metrics.items():
        old=row[f'{dataset}_{mode}_{name}']
        if np.isfinite(old):errors.append(abs(value-old))
assert max(errors)<1e-10

# Retain every registration; completed posterior runs are a separate view.
a.to_csv(out/'all_586_registrations.csv', index=False)
s.to_csv(out/'all_282_new_cases.csv',index=False)
fields=['case_id','label','block','equation','subset','split','selection_seed','prior','k_depth','k_profile',
        'adapt_delta','treedepth','warmup','samples','fit_seconds','requested_core_hours',
        'max_rhat','min_bulk_ess','min_tail_ess','min_ebfmi','divergences','treedepth_hits',
        'convergence_screen_pass','diagnostic_screen_pass']
for ds in ['validation','reserve']:
    for mode in ['conditional','population']:
        fields += [f'{ds}_{mode}_{name}' for name in ['rmse','mae','shape_rmse','peak_depth_mae',
          'rise_slope_mae','rise_slope_pairs','rise_width_mae','rise_width_pairs','onset50_onset_mae_given_both_detected',
          'onset50_both_detected','coverage90','interval_score90','median_residual_lag1']]
complete[fields].to_csv(out/'all_67_completed_posterior_runs.csv',index=False)
s.query('mode=="prior"').to_csv(out/'12_prior_predictive_checks.csv',index=False)

# Cache the per-cast +0.5 onset results. Censored/missing detections stay missing.
onset={}
for cid in complete.case_id:
    for ds in ['validation','reserve']:
        for mode in ['conditional','population']:
            p=root/f'runs/case_{cid:03d}/checks/{ds}_{mode}_onset_cast_metrics.csv'
            x=pd.read_csv(p);x=x[x.threshold.eq(.5)].set_index('cast_key')
            assert x.index.is_unique
            onset[cid,ds,mode]=x

def pair(cid, reference, dataset, mode):
    """Compare two models only when their evaluation casts are identical."""
    x=cs.query('case_id==@cid and dataset==@dataset and prediction_mode==@mode').set_index('cast_key')
    y=cs.query('case_id==@reference and dataset==@dataset and prediction_mode==@mode').set_index('cast_key')
    if len(x)==0 or set(x.index)!=set(y.index):return None
    y=y.loc[x.index]
    result={'case_id':cid,'reference_case':reference,'dataset':dataset,'prediction_mode':mode,
            'casts':len(x),'delta_rmse':np.sqrt(x.mse.mean())-np.sqrt(y.mse.mean()),
            'fraction_casts_lower_mse':np.mean(x.mse<y.mse),
            'delta_shape_rmse':np.sqrt(x.shape_mse.mean())-np.sqrt(y.shape_mse.mean()),
            'delta_peak_depth_mae':(x.peak_depth_abs_error-y.peak_depth_abs_error).mean()}
    for key,column in [('shape','rise_slope_abs_error'),('shape','rise_width_abs_error'),('onset','onset_abs_error')]:
        if key=='shape':
            xx=ss.query('case_id==@cid and dataset==@dataset and prediction_mode==@mode').set_index('cast_key')
            yy=ss.query('case_id==@reference and dataset==@dataset and prediction_mode==@mode').set_index('cast_key')
        else:xx=onset[cid,dataset,mode];yy=onset[reference,dataset,mode]
        assert set(xx.index)==set(yy.index)
        yy=yy.loc[xx.index];good=np.isfinite(xx[column])&np.isfinite(yy[column])
        result[column+'_common_casts']=int(good.sum())
        result['delta_'+column+'_common_casts']=(xx.loc[good,column]-yy.loc[good,column]).mean()
    return result

# Every eligible completed new case is compared with its planned equation and
# settings references. These are descriptive paired comparisons, not rankings.
results=[]
for _,row in complete.iterrows():
    refs=set(int(row[n]) for n in ['equation_reference_case','settings_reference_case'] if pd.notna(row[n]))
    for reference in refs:
        if reference==row.case_id or reference not in set(complete.case_id):continue
        for ds in ['validation','reserve']:
            for mode in ['conditional','population']:
                result=pair(int(row.case_id),reference,ds,mode)
                if result:results.append(result)
pairs=pd.DataFrame(results);pairs.to_csv(out/'paired_comparisons_common_casts.csv',index=False)
# Check saved planned contrasts independently, including comparable-cast counts.
planned=pd.read_csv(root/'results/planned_paired_contrasts.csv')
merged=pairs.merge(planned.rename(columns={'mode':'prediction_mode'}),on=['case_id','reference_case','dataset','prediction_mode'],suffixes=('_audit','_saved'))
assert np.max(np.abs(merged.delta_rmse_audit-merged.delta_rmse_saved))<1e-10
assert np.max(np.abs(merged.delta_shape_rmse_audit-merged.delta_shape_rmse_saved))<1e-10

# Four selections of the same training fraction; sampler seed remains fixed.
repeat_ids=[19,22,25,28,31,34]+list(range(197,203))+list(range(205,211))+list(range(213,219))
rep=s[s.case_id.isin(repeat_ids)].copy();rep.to_csv(out/'four_20percent_cast_selections.csv',index=False)
labels={'original_z':'Baseline depth','interaction_z':'Shared MLD × distance','mld_vary_shared':'Shared MLD–depth',
        'mld_halfshift':'Half shift','mld_shift':'Full shift','mld_vary_region':'Regional MLD–depth'}
order=list(labels)
colors=['#475569','#9ca3af','#8b5cf6','#0891b2','#2563eb','#059669']
fig,axs=plt.subplots(1,2,figsize=(12,5),layout='constrained')
for ax,ds in zip(axs,['validation','reserve']):
    for seed,grp in rep.groupby('selection_seed'):
        grp=grp.set_index('equation').reindex(order)
        ax.plot(np.arange(6),grp[f'{ds}_conditional_rmse'],color='#d1d5db',lw=1,zorder=1)
    for j,(eq,color) in enumerate(zip(order,colors)):
        vals=rep.loc[rep.equation.eq(eq),f'{ds}_conditional_rmse']
        ax.scatter(np.full(len(vals),j),vals,color=color,s=40,zorder=2)
    ax.set_xticks(range(6),[labels[e].replace(' ','\n',1) for e in order],fontsize=9)
    ax.set_ylabel('Equal-cast RMSE (chl_norm units)');ax.set_title(ds.capitalize()+' cohort (108 casts)')
    ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.18)
fig.suptitle('20% training: results depend on which casts are selected',fontsize=14)
fig.savefig(out/'cast_selection_comparison.png',dpi=160);plt.close(fig)

# Shape illustration uses the same primary 20% training casts, fixed regional
# distance, and zero date/ENSO effects. Dashed lines flag unsupported combinations.
fig,axs=plt.subplots(2,3,figsize=(12,8),sharey=True,layout='constrained')
levels=['Shallow MLD (10%)','Median MLD','Deep MLD (90%)']
palette=['#2563eb','#059669','#e11d48']
for col,(cid,title) in enumerate([(19,'Baseline depth equation'),(25,'Half shift'),(129,'Regional MLD–depth; RS k=15')]):
    data=pd.read_csv(root/f'runs/case_{cid:03d}/checks/mld_shape_slices.csv')
    for i,region in enumerate(['Central RS','South RS']):
        ax=axs[i,col]
        for level,color in zip(levels,palette):
            x=data[data.Region.eq(region)&data.mld_level.eq(level)].sort_values('depth')
            ok=x.joint_support.astype(bool)
            ax.plot(x['mean'].where(ok),x.depth,color=color,lw=1.8,label=f'{level.split(" MLD")[0]}: {x.MLD_temp_MIX.iloc[0]:.0f} m')
            ax.plot(x['mean'].where(~ok),x.depth,color=color,lw=1,ls='--')
            ax.fill_betweenx(x.depth.to_numpy(),x.lower90.to_numpy(),x.upper90.to_numpy(),where=ok.to_numpy(),color=color,alpha=.09)
        ax.set_ylim(150,10);ax.set_title(title if i==0 else '',fontsize=11)
        ax.set_xlabel('Expected chl_norm');ax.set_ylabel(region+'\nDepth (m)' if col==0 else 'Depth (m)')
        ax.grid(alpha=.15);ax.legend(fontsize=8,loc='lower right');ax.spines[['top','right']].set_visible(False)
fig.suptitle('MLD can change the regional profile shape, not only its depth',fontsize=14)
fig.savefig(out/'mld_shape_comparison.png',dpi=160);plt.close(fig)

# A broad bootstrap diagnostic describes uncertainty across casts. It is not
# an independent final-test interval and does not account for shared dates.
bootstrap=[]
for cid,ref in [(7,1),(8,1),(25,19),(34,19),(129,34)]:
    for ds in ['validation','reserve']:
        x=cs.query('case_id==@cid and dataset==@ds and prediction_mode=="conditional"').set_index('cast_key')
        y=cs.query('case_id==@ref and dataset==@ds and prediction_mode=="conditional"').set_index('cast_key').loc[x.index]
        rng=np.random.default_rng(123);idx=rng.integers(0,len(x),(5000,len(x)))
        boot=np.sqrt(x.mse.to_numpy()[idx].mean(1))-np.sqrt(y.mse.to_numpy()[idx].mean(1))
        low,high=np.quantile(boot,[.025,.975])
        bootstrap.append({'case_id':cid,'reference_case':ref,'dataset':ds,'delta_rmse':np.sqrt(x.mse.mean())-np.sqrt(y.mse.mean()),'cast_bootstrap_q025':low,'cast_bootstrap_q975':high,'note':'Paired cast bootstrap; shared dates and model selection not accounted for'})
pd.DataFrame(bootstrap).to_csv(out/'exploratory_cast_bootstrap.csv',index=False)

# Quantify underprediction of sharp rises using exactly the same eligible casts.
sharp=[]
for ds in ['validation','reserve']:
    part=ss.query('case_id in [19,25,34,129] and dataset==@ds and prediction_mode=="conditional"')
    valid_sets=[set(x.loc[x.rise_slope_abs_error.notna(),'cast_key'])for _,x in part.groupby('case_id')]
    common=set.intersection(*valid_sets)
    for cid,x in part[part.cast_key.isin(common)].groupby('case_id'):
        sharp.append({'dataset':ds,'case_id':cid,'common_casts':len(x),'observed_median_max_slope':x.observed_max_rise_slope.median(),'predicted_median_max_slope':x.predicted_max_rise_slope.median(),'median_slope_ratio':(x.predicted_max_rise_slope/x.observed_max_rise_slope).median()})
pd.DataFrame(sharp).to_csv(out/'sharpness_common_casts.csv',index=False)
# Validate that the same primary evaluation casts were used in all four seeds.
for ds in ['validation','reserve']:
    sets=[set(x.cast_key)for _,x in cs.query('case_id in @repeat_ids and dataset==@ds and prediction_mode=="conditional"').groupby('case_id')]
    assert all(x==sets[0]for x in sets)

audit={'registrations':len(a),'completed_registrations':int(a.status.eq('complete').sum()),
       'new_cases':len(s),'completed_new_posterior':len(complete),'completed_new_prior':12,
       'unique_registration_keys':bool(a.run_key.is_unique),'aggregate_values_recomputed':len(errors),
       'maximum_aggregate_discrepancy':max(errors),'paired_comparisons_recomputed':len(pairs),
       'matched_saved_comparison_rows':len(merged),'same_primary_validation_and_reserve_casts_across_four_selections':True,
       'scope':'Archive snapshot only; no refits; no R/Stan rerun; older rows reviewed from aggregate summaries'}
(out/'audit.json').write_text(json.dumps(audit,indent=2))
print(json.dumps(audit,indent=2));print(pd.DataFrame(sharp).round(4).to_string(index=False))
