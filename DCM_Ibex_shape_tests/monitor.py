#!/usr/bin/env python3
"""Read all campaign manifests and saved diagnostics; submit nothing."""
from pathlib import Path
import csv,collections,json,subprocess,sys
root=Path(sys.argv[1]);diagnostics='--diagnostics' in sys.argv[2:]
def read(path,delimiter=','):
    try:
        with path.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter=delimiter))
    except (OSError,csv.Error):return []
p=subprocess.run(['squeue','--me','--noheader','--format=%i|%T|%R'],stdout=subprocess.PIPE,universal_newlines=True)
live={x[0]:(x[1],x[2]) for x in (line.split('|',2) for line in p.stdout.splitlines()) if len(x)==3}
for parent in ('experiments','benchmarks'):
    for manifest in sorted((root/parent).glob('*/manifest.tsv')):
        items=read(manifest,'\t');counts=collections.Counter()
        if not items:continue
        print('\n'+manifest.parent.name)
        for row in items:
            run=Path(row.get('run_dir',''))
            status='complete' if (run/'CHECKS_COMPLETE.txt').exists() else 'checks_failed' if (run/'CHECKS_FAILED.txt').exists() else 'fit_failed' if (run/'FIT_FAILED.txt').exists() else 'awaiting_checks' if (run/'FIT_COMPLETE.txt').exists() else 'incomplete'
            counts[status]+=1
            if diagnostics:
                records=read(run/'test_summary.csv') or read(run/'diagnostic_summary.csv');d=records[0] if records else {}
                job=row.get('check_job') if status=='awaiting_checks' else row.get('fit_job','')
                state,reason=live.get(job,('',''))
                print('case',row.get('case_id',row.get('case','?')),status,state,
                    'Rhat='+d.get('max_rhat','NA'),'div='+d.get('divergences','NA'),
                    'tree='+d.get('treedepth_hits','NA'),'valRMSE='+d.get('validation_conditional_rmse','NA'),
                    'role='+d.get('evaluation_role','legacy'),reason)
        print(' '.join(str(n)+' '+s for s,n in sorted(counts.items())))
pointer=root/'logs/latest_shape150_experiment.txt'
if pointer.exists():
    exp=Path(pointer.read_text().strip());path=exp/'dispatch_status.json'
    if path.exists():print('\nNEW CAMPAIGN DISPATCHER\n'+path.read_text())
    for p in (exp/'scheduler').glob('*.rejected.json'):print('SUBMISSION REJECTED',p.name,json.loads(p.read_text())['reason'])
    if (exp/'DISPATCHER_RESTART_NEEDED.txt').exists():print((exp/'DISPATCHER_RESTART_NEEDED.txt').read_text())
