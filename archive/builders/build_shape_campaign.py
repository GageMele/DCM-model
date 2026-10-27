"""Build the frozen design and self-contained Ibex MLD-shape test package.

This is local packaging only: it neither connects to Ibex nor submits jobs.
Training-selection seeds and MCMC seeds are distinct experimental factors.
"""
from pathlib import Path
import csv, collections, hashlib, json, shutil, zipfile, base64

ROOT=Path(__file__).resolve().parent
PKG=ROOT/'DCM_Ibex_shape_tests'
PKG.mkdir(exist_ok=True)

CORE=['original_z','interaction_z','mld_halfshift','mld_shift',
      'mld_vary_shared','mld_vary_region']
TENSORS=['tensor_region','tensor_shift_region','tensor_halfshift_region']
cases=[]
def add(eq,fraction=20,k=7,km=4,block='size_basis',selection='frozen',
        selection_seed=0,split='primary',prior='custom_z',mode='posterior',
        delta=.97,warmup=1000,samples=2000,mcmc_seed=1234,k_profile='uniform',
        threads=None):
    if fraction==100:selection='all';selection_seed=0
    threads=threads or (8 if fraction>=75 and mode=='posterior' else 4)
    tensor=eq in TENSORS
    memory=(256 if tensor else 192) if fraction>=75 else ((160 if tensor else 96) if fraction>=40 else (96 if tensor else 64))
    if mode=='prior':memory=64
    ident=len(cases)+1
    dataset=f'{split}_p{fraction:02d}_{selection}' + (f'_{selection_seed}' if selection_seed else '')
    role='training_reconstruction' if fraction==100 else 'heldout'
    cases.append(dict(case_id=ident,label=f's{ident:03d}_{block}_{eq}_p{fraction:02d}_k{k}_{k_profile}',
      block=block,equation=eq,subset=f'p{fraction:02d}',fraction=fraction/100,
      split=split,dataset_id=dataset,selection=selection,selection_seed=selection_seed,
      prior=prior,mode=mode,variance='region',family='gaussian',cast='none',autocor='none',
      depth_max=150,k_depth=k,k_mld=km,k_distance=5,k_profile=k_profile,
      chains=4,cores=4,threads=threads,cpus=4*threads,memory_gb=memory,
      iter=warmup+samples,warmup=warmup,samples=samples,adapt_delta=delta,treedepth=14,
      seed=mcmc_seed,refresh=100,evaluation_role=role,fit_days=7))

# Nested training sizes retain the existing 10/20/50% cast sets as anchors.
for fraction in (10,20,40,50,75):
    for eq in CORE:
        for k in (7,11,15):add(eq,fraction,k)
for eq in CORE:add(eq,100,7,block='full_data')
for eq in ('original_z','mld_halfshift','mld_vary_region'):
    add(eq,100,15,block='full_data_high_k')
for fraction in (20,50):
    for eq in CORE:
        for k in (9,13):add(eq,fraction,k,block='intermediate_k')

# Give Central and South RS k=15 while holding the other regional bases at 7.
for fraction in (20,50,75):
    for eq in CORE:add(eq,fraction,7,block='regional_k',k_profile='rs15_other7')

# These surfaces let MLD alter rise steepness, width, peak and other shape.
for fraction in (10,20,40,50,75):
    for eq in TENSORS:
        for k in (7,11,15):add(eq,fraction,k,block='mld_surface')
for eq in ('tensor_region','tensor_halfshift_region'):
    add(eq,100,7,block='full_data_surface')
for eq in ('tensor_region','tensor_halfshift_region'):
    for k in (7,15):
        for km in (5,7):add(eq,20,k,km=km,block='mld_basis')

# Repeat model comparisons on three NEW 20% cast selections, holding MCMC seed.
for selection_seed in (20261023,20261024,20261025):
    for eq in CORE+['tensor_region','tensor_halfshift_region']:
        add(eq,20,7,block='cast_selection',selection='repeat',selection_seed=selection_seed)

# Test an explicit sharpening/flattening hypothesis, alongside flexible models.
for fraction in (20,50):
    for eq in ('warp_half_pos025','warp_half_pos050','warp_half_neg025'):
        for k in (7,15):add(eq,fraction,k,block='mld_width')

# Change one sampling control at a time relative to the corresponding anchor.
for eq in ('mld_halfshift','mld_shift','mld_vary_region','tensor_halfshift_region'):
    add(eq,block='delta95',delta=.95)
    add(eq,block='delta99',delta=.99)
    add(eq,block='warmup2000',warmup=2000)
    add(eq,block='samples4000',samples=4000)
for eq in ('mld_vary_region','tensor_halfshift_region'):
    for prior in ('b_tight','scales_halfnormal'):add(eq,block='prior_sensitivity',prior=prior)
for eq in ('original_z','mld_halfshift','mld_shift','mld_vary_region','tensor_halfshift_region','warp_half_pos025'):
    for seed in (1234,20261022):
        add(eq,block='prior_predictive',mode='prior',warmup=500,samples=500,mcmc_seed=seed)

# These validation casts are excluded by entire Ocean/Region/date groups.
for selection_seed in (20261022,20261023):
    for eq in CORE+TENSORS:
        add(eq,block='date_blocked',split='date_blocked',selection='repeat',selection_seed=selection_seed)

assert len(cases)==282
assert len({r['label'] for r in cases})==len(cases)
signatures=[tuple((k,r[k]) for k in r if k not in ('case_id','label','block')) for r in cases]
assert len(set(signatures))==len(cases),'Duplicate experimental settings'
assert all(r['fraction']>=.1 and r['iter']==r['warmup']+r['samples'] and r['cpus']==r['chains']*r['threads'] for r in cases)
for r in cases:
    # Launch useful 20% comparisons early; leave very large fits behind them.
    r['priority']=(0 if r['mode']=='prior' else 10 if r['fraction']==.2 and r['k_depth']==7 else 20 if r['fraction']<=.2 else 30 if r['fraction']<=.5 else 40)
    r['equation_reference_case']='NA';r['settings_reference_case']='NA'
    if r['mode']=='posterior':
        same=lambda q:all(q[k]==r[k] for k in ('dataset_id','mode','seed'))
        anchors=[q for q in cases if same(q) and q['equation']=='original_z' and
                 q['k_depth']==r['k_depth'] and q['k_profile']==r['k_profile'] and
                 q['prior']=='custom_z' and q['block'] in ('size_basis','full_data','full_data_high_k','intermediate_k','regional_k','cast_selection','date_blocked')]
        if anchors:r['equation_reference_case']=anchors[0]['case_id']
        base_eq='mld_halfshift' if r['equation'].startswith('warp_') else r['equation']
        bases=[q for q in cases if same(q) and q['equation']==base_eq and q['k_depth']==7 and
               q['k_profile']=='uniform' and q['k_mld']==4 and q['prior']=='custom_z' and
               q['adapt_delta']==.97 and q['warmup']==1000 and q['samples']==2000]
        if bases:r['settings_reference_case']=bases[0]['case_id']
with (PKG/'cases.tsv').open('w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,list(cases[0]),delimiter='\t');w.writeheader();w.writerows(cases)
summary={'cases':len(cases),'posterior_fits':sum(r['mode']=='posterior' for r in cases),
 'prior_predictive_fits':sum(r['mode']=='prior' for r in cases),
 'full_data_posterior_fits':sum(r['mode']=='posterior' and r['fraction']==1 for r in cases),
 'fractions':dict(collections.Counter(r['subset'] for r in cases)),
 'depth_basis_dimensions':sorted({r['k_depth'] for r in cases}),
 'blocks':dict(collections.Counter(r['block'] for r in cases))}
(PKG/'DESIGN.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')

# Reuse read-only sampler/diagnostic utilities, not any old scheduler or scorer.
for name in ('common.R','live_diagnostics.awk'):
    shutil.copyfile(ROOT/'DCM_Ibex_onset_tests'/name,PKG/name)
shutil.copyfile(ROOT/'DCM_Ibex_onset_support_repair_r1/onset_metrics.R',PKG/'onset_metrics.R')
legacy=(ROOT/'DCM_Ibex_onset_tests/predictive.R').read_text(encoding='utf-8')
helpers=legacy[legacy.index('trap_area <-'):legacy.index('plot_validation_examples<-')]
(PKG/'profile_helpers.R').write_text('# Existing cast-weighted profile metrics, reused without changing observations.\n'+helpers,encoding='utf-8')
checker=(ROOT/'DCM_Ibex_onset_tests/03_check_case.R').read_text(encoding='utf-8')
checker=checker[:checker.index('  trace_pars <-')]
checker=checker.replace('"campaign.R"','"models.R"')
checker=checker.replace('run <- file.path(experiment, "runs", sprintf("case_%02d", case$case_id))','run <- run_directory(experiment, case)')
checker=checker.replace('rhat <= 1.01','rhat < 1.01')
checker+='''  source(file.path(experiment,"scoring.R"))
  result<-score_predictions(fit,experiment,run,case,result)
  write.csv(result,file.path(run,"test_summary.csv"),row.names=FALSE)
  if(case$mode=="posterior")tryCatch(plot_shape_sensitivity(fit,experiment,run,case),
    error=function(e)writeLines(conditionMessage(e),file.path(run,"FIGURES_FAILED.txt")))
  writeLines(c("Complete; all diagnostic warnings remain visible.",
    "Full-data reconstruction is never labelled held-out validation.",
    "Conditional and zero-group-effect population predictions have distinct interpretations.",
    "Onset/width errors retain their comparable-case counts and censoring."),
    file.path(run,"CHECKS_COMPLETE.txt"))
},error=function(e){writeLines(conditionMessage(e),file.path(run,"CHECKS_FAILED.txt"));stop(e)})
'''
(PKG/'03_check_case.R').write_text(checker,encoding='utf-8')
print(json.dumps(summary,indent=2))
