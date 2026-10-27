#!/usr/bin/env python3
"""Resumable, bounded Slurm dispatch for the frozen MLD-shape campaign.

Only this file submits jobs. Preparation precedes fitting; checks run after
saved fits exist. Accepted model fits are never automatically resubmitted.
"""
import argparse, collections, contextlib, csv, datetime, hashlib, json, os
from pathlib import Path
import re, shutil, subprocess, sys, time

TERMINAL={'COMPLETED','FAILED','CANCELLED','TIMEOUT','OUT_OF_MEMORY','NODE_FAIL',
          'PREEMPTED','BOOT_FAIL','DEADLINE','REVOKED'}
class QueueLimit(Exception):pass
class UncertainSubmission(Exception):pass
class DefiniteRejection(Exception):pass

def atomic_text(path,text):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_name(path.name+'.tmp.'+str(os.getpid()))
    temp.write_text(text,encoding='utf-8');os.replace(str(temp),str(path))
def atomic_json(path,value):atomic_text(path,json.dumps(value,indent=2)+'\n')
def read_json(path,default=None):
    return json.loads(Path(path).read_text(encoding='utf-8')) if Path(path).exists() else default
def read_cases(experiment):
    with (Path(experiment)/'cases.tsv').open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def runner(args):
    p=subprocess.run(args,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
    return p.returncode,p.stdout,p.stderr
def read_command(args,run=runner):
    code,out,err=run(args)
    if code:raise RuntimeError('Read-only Slurm query failed: '+' '.join(args[:2])+' '+err.strip())
    return out

@contextlib.contextmanager
def lock(path):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('a+b') as f:
        if os.name=='nt':
            import msvcrt
            f.seek(0,2)
            if f.tell()==0:f.write(b'0');f.flush()
            f.seek(0);msvcrt.locking(f.fileno(),msvcrt.LK_NBLCK,1)
            try:yield
            finally:f.seek(0);msvcrt.locking(f.fileno(),msvcrt.LK_UNLCK,1)
        else:
            import fcntl
            fcntl.flock(f.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
            try:yield
            finally:fcntl.flock(f.fileno(),fcntl.LOCK_UN)

def live_jobs(run=runner):
    out=read_command(['squeue','--me','--noheader','--format=%i|%j|%T'],run)
    result={}
    for line in out.splitlines():
        fields=line.strip().split('|')
        if len(fields)>=3 and fields[0].isdigit():result[fields[0]]={'name':fields[1],'state':fields[2]}
    return result

class Ledger:
    def __init__(self,experiment,run=runner):
        self.exp=Path(experiment);self.run=run
        self.directory=self.exp/'scheduler';self.directory.mkdir(exist_ok=True)
        self.meta=read_json(self.exp/'campaign.json')
        self.token=self.meta['token']
    def record(self,stage):return self.directory/(stage+'.json')
    def name(self,stage):return 'dcm-s150-'+self.token+'-'+stage
    def lookup(self,name):
        found={j for j,v in live_jobs(self.run).items() if v['name']==name}
        start=self.meta['created_utc'][:10]+'T00:00:00'
        out=read_command(['sacct','-X','--noheader','--parsable2','--starttime='+start,
            '--name='+name,'--format=JobIDRaw,JobName%100'],self.run)
        for line in out.splitlines():
            f=line.strip().split('|')
            if len(f)>=2 and f[0].isdigit() and f[1]==name:found.add(f[0])
        return sorted(found)
    def ensure(self,stage,args):
        path=self.record(stage);intent=self.directory/(stage+'.intent.json')
        saved=read_json(path)
        if saved:return saved['job_id']
        if intent.exists():
            found=self.lookup(self.name(stage))
            if len(found)==1:
                atomic_json(path,{'job_id':found[0],'stage':stage,'name':self.name(stage),'recovered_utc':utc()})
                return found[0]
            raise UncertainSubmission('Uncertain earlier submission for '+stage+
                '; accepted jobs will not be duplicated. Inspect '+str(intent))
        request=['sbatch','--parsable','--job-name='+self.name(stage)]+list(args)
        atomic_json(intent,{'stage':stage,'name':self.name(stage),'requested_utc':utc(),'args':request})
        code,out,err=self.run(request)
        accepted=re.fullmatch(r'\s*([0-9]+)(?:;[^\s]+)?\s*',out)
        if code==0 and accepted:
            job=accepted.group(1)
            atomic_json(path,{'job_id':job,'stage':stage,'name':self.name(stage),'accepted_utc':utc()})
            return job
        message=(out+'\n'+err).strip()
        if 'sbatch: error:' in message and not accepted:
            # Explicit rejection means no job was accepted; retrying a queue
            # capacity rejection cannot duplicate a model fit.
            intent.unlink()
            quota=re.search(r'MaxSubmit|maximum number of jobs|too many jobs|QOSMaxSubmit|AssocMaxSubmit',message,re.I)
            if quota:raise QueueLimit(message)
            atomic_json(self.directory/(stage+'.rejected.json'),{'utc':utc(),'reason':message,'args':request})
            raise DefiniteRejection(message)
        raise UncertainSubmission('sbatch outcome uncertain for '+stage+': '+message)
    def records(self):
        return {p.stem:read_json(p) for p in self.directory.glob('*.json')
                if not p.name.endswith(('.intent.json','.rejected.json'))}
    def accounting(self,records):
        ids=sorted({r['job_id'] for r in records.values() if 'job_id' in r})
        if not ids:return {}
        out=read_command(['sacct','-X','--noheader','--parsable2',
            '--starttime='+self.meta['created_utc'][:10]+'T00:00:00',
            '--jobs='+','.join(ids),'--format=JobIDRaw,State'],self.run)
        result={}
        for line in out.splitlines():
            f=line.strip().split('|')
            if len(f)>=2 and f[0].isdigit():result[f[0]]=f[1].split()[0].rstrip('+')
        return result

def partition_time(partition,requested,run=runner):
    out=read_command(['scontrol','show','partition',partition,'-o'],run)
    match=re.search(r'\bMaxTime=(\S+)',out)
    if not match or match.group(1) in ('UNLIMITED','INFINITE'):return requested
    def seconds(value):
        day,clock=value.split('-',1) if '-' in value else ('0',value)
        fields=list(map(int,clock.split(':')))
        if len(fields)==2:fields=[0]+fields
        return int(day)*86400+fields[0]*3600+fields[1]*60+fields[2]
    return match.group(1) if seconds(match.group(1))<seconds(requested) else requested

def job_args(ledger,stage,case=None,dependency=None):
    meta=ledger.meta;root=meta['root'];exp=str(ledger.exp)
    args=['--partition='+meta['partition'],'--chdir='+root,'--nodes=1','--ntasks=1',
          '--export=ALL,DCM_ROOT='+root+',EXPERIMENT_DIR='+exp+',STAGE='+stage,
          '--output='+str(ledger.exp/'logs/%x-%j.out'),'--error='+str(ledger.exp/'logs/%x-%j.err')]
    if case is not None:
        args[4]+=',CASE_ID='+str(case['case_id'])
    if stage=='prepare':cpu,mem,wall=1,96,meta['prepare_time']
    elif stage=='controller':cpu,mem,wall=1,4,meta['controller_time']
    elif stage=='check':cpu,mem,wall=4,64,meta['check_time']
    else:cpu,mem,wall=int(case['cpus']),int(case['memory_gb']),meta['prior_time'] if case['mode']=='prior' else meta['fit_time']
    args+=['--cpus-per-task='+str(cpu),'--mem='+str(mem)+'G','--time='+wall,'--hint=nomultithread']
    if meta['constraint'] and stage=='fit':args+=['--constraint='+meta['constraint']]
    if dependency:args+=['--dependency=afterany:'+dependency,'--kill-on-invalid-dep=yes']
    args+=[str(ledger.exp/'worker.slurm')]
    return args

def case_status(exp,case,records,states,preflight):
    cid=int(case['case_id']);run=exp/'runs'/('case_%03d'%cid)
    if preflight.get(cid) is False:return 'preflight_failed'
    if (run/'CHECKS_COMPLETE.txt').exists():return 'complete'
    if (run/'CHECKS_FAILED.txt').exists():return 'checks_failed'
    if (run/'FIT_FAILED.txt').exists():return 'fit_failed'
    if (run/'FIT_COMPLETE.txt').exists():
        rec=records.get('check_%03d'%cid)
        if rec and states.get(rec['job_id']) in TERMINAL:return 'checks_failed'
        return 'checking' if rec else 'awaiting_checks'
    rec=records.get('fit_%03d'%cid)
    if rec and states.get(rec['job_id']) in TERMINAL:return 'fit_failed'
    return 'fitting' if rec else 'not_submitted'

def dispatch_once(ledger,max_active,max_cpus):
    cases=read_cases(ledger.exp);records=ledger.records();live=live_jobs(ledger.run)
    states=ledger.accounting(records);states.update({j:v['state'] for j,v in live.items()})
    with (ledger.exp/'preflight_cases.csv').open(encoding='utf-8',newline='') as f:
        preflight={int(r['case_id']):r['valid'].upper()=='TRUE' for r in csv.DictReader(f)}
    resources={};known={}
    for c in cases:
        for stage in ('fit','check'):
            tag=stage+'_%03d'%int(c['case_id']);r=records.get(tag)
            if r:resources[r['job_id']]=int(c['cpus']) if stage=='fit' else 4;known[r['job_id']]=tag
    # Reserve slots for accepted jobs that have not yet appeared in squeue or
    # accounting; absence from a query is not evidence that a job finished.
    active=set()
    for job,tag in known.items():
        stage,cid=tag.split('_');run=ledger.exp/'runs'/('case_%03d'%int(cid))
        marked=any((run/name).exists() for name in
            (('FIT_COMPLETE.txt','FIT_FAILED.txt') if stage=='fit' else ('CHECKS_COMPLETE.txt','CHECKS_FAILED.txt')))
        if job in live or (states.get(job) not in TERMINAL and not marked):active.add(job)
    cpus=sum(resources[j] for j in active)
    pending=[];status={}
    for c in cases:
        cid=int(c['case_id']);s=case_status(ledger.exp,c,records,states,preflight);status[cid]=s
        if s in ('fit_failed','checks_failed'):
            marker='FIT_FAILED.txt' if s=='fit_failed' else 'CHECKS_FAILED.txt'
            run=ledger.exp/'runs'/('case_%03d'%cid)
            if not (run/marker).exists():
                atomic_text(run/marker,'Accepted Slurm stage ended without a completion marker; inspect its log. No automatic refit.\n')
        if s=='awaiting_checks':pending.append((0,int(c['priority']),cid,'check',c))
        elif s=='not_submitted':pending.append((1,int(c['priority']),cid,'fit',c))
    submitted=[]
    for _,_,cid,stage,c in sorted(pending):
        if len(active)>=max_active:break
        cpu=int(c['cpus']) if stage=='fit' else 4
        if cpus+cpu>max_cpus:continue
        tag=stage+'_%03d'%cid
        if (ledger.directory/(tag+'.rejected.json')).exists():continue
        try:job=ledger.ensure(tag,job_args(ledger,stage,c))
        except QueueLimit as e:
            atomic_text(ledger.exp/'QUEUE_CAPACITY_WAIT.txt',utc()+' '+str(e)+'\n');break
        except DefiniteRejection as e:
            print('SUBMISSION REJECTED',tag,str(e),flush=True);continue
        active.add(job);cpus+=cpu;submitted.append((cid,stage,job))
        print('QUEUED',stage,cid,c['label'],job,flush=True)
    records=ledger.records()
    # All cases stay visible, including unsubmitted, failed and full-data cases.
    fields=['case_id','label','fit_job','check_job','cpus','run_dir']
    path=ledger.exp/'manifest.tsv';temp=path.with_name(path.name+'.tmp')
    with temp.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader()
        for c in cases:
            cid=int(c['case_id']);w.writerow(dict(case_id=cid,label=c['label'],
                fit_job=records.get('fit_%03d'%cid,{}).get('job_id','NA'),
                check_job=records.get('check_%03d'%cid,{}).get('job_id','NA'),cpus=c['cpus'],
                run_dir=str(ledger.exp/'runs'/('case_%03d'%cid))))
    os.replace(str(temp),str(path))
    counts=collections.Counter(status.values())
    atomic_json(ledger.exp/'dispatch_status.json',{'utc':utc(),'statuses':dict(counts),
        'active_model_jobs':len(active),'allocated_model_cpus':cpus,'submitted_this_pass':submitted})
    finished=all(s in ('complete','fit_failed','checks_failed','preflight_failed') for s in status.values())
    return finished

def launch(root,package,source,new=False):
    root=Path(root).resolve();package=Path(package).resolve();source=Path(source).resolve()
    for p in ('ibex_environment.sh','full_input_250m.rds','scaling.rds','region_levels.rds'):
        if not (source/p).is_file():raise RuntimeError('Missing prepared follow-up source file: '+str(source/p))
    pointer=root/'logs/latest_shape150_experiment.txt';pointer.parent.mkdir(exist_ok=True)
    signature=hashlib.sha256((package/'SHA256SUMS.txt').read_bytes()).hexdigest()
    with lock(root/'logs/shape150_launch.lock'):
        exp=Path(pointer.read_text().strip()) if pointer.exists() and not new else None
        if exp is not None:
            if read_json(exp/'campaign.json')['package_sha256']!=signature:
                raise RuntimeError('Latest campaign uses another frozen package. Set SHAPE_NEW=1 only to request a new campaign.')
        else:
            stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
            exp=root/'experiments'/('shape150_'+stamp+'_'+str(os.getpid()))
            exp.mkdir(parents=True);(exp/'logs').mkdir();(exp/'scheduler').mkdir()
            for p in package.iterdir():
                if p.is_file():shutil.copyfile(str(p),str(exp/p.name))
            shutil.copyfile(str(source/'ibex_environment.sh'),str(exp/'ibex_environment.sh'))
            atomic_text(exp/'source_experiment.txt',str(source)+'\n')
            partition=os.getenv('SHAPE_PARTITION','batch')
            token=hashlib.sha256(str(exp).encode()).hexdigest()[:10]
            meta=dict(root=str(root),created_utc=utc(),token=token,package_sha256=signature,
                partition=partition,constraint=os.getenv('SHAPE_CONSTRAINT','intel&cascadelake'),
                fit_time=partition_time(partition,os.getenv('SHAPE_FIT_TIME','7-00:00:00')),
                prior_time=partition_time(partition,'04:00:00'),prepare_time=partition_time(partition,'08:00:00'),
                check_time=partition_time(partition,'2-00:00:00'),controller_time=partition_time(partition,'7-00:00:00'))
            atomic_json(exp/'campaign.json',meta);atomic_text(pointer,str(exp)+'\n')
        ledger=Ledger(exp)
        limits=read_json(exp/'dispatch_limits.json',{'max_active':48,'max_cpus':1024})
        if 'SHAPE_MAX_ACTIVE' in os.environ:limits['max_active']=int(os.environ['SHAPE_MAX_ACTIVE'])
        if 'SHAPE_MAX_CPUS' in os.environ:limits['max_cpus']=int(os.environ['SHAPE_MAX_CPUS'])
        if min(limits.values())<1:raise RuntimeError('Dispatch limits must be positive')
        atomic_json(exp/'dispatch_limits.json',limits)
        if (exp/'CAMPAIGN_FINISHED.txt').exists():
            print('Campaign already finished; accepted fits were preserved:',exp);return exp
        if (exp/'PREPARATION_FAILED.txt').exists():raise RuntimeError('Existing preparation failed; inspect '+str(exp/'PREPARATION_FAILED.txt'))
        prep=ledger.ensure('preparation',job_args(ledger,'prepare'))
        live=live_jobs();records=ledger.records();states=ledger.accounting(records)
        controllers=[(s,r) for s,r in records.items() if s.startswith('controller_')]
        if any(r['job_id'] in live for _,r in controllers):
            print('Existing dispatcher remains queued/running; no duplicate submitted.')
        else:
            generation=max([int(s.split('_')[1]) for s,_ in controllers]+[0])+1
            stage='controller_%03d'%generation
            # Use dependencies only while the preparation job is visible. A
            # completed marker is authoritative even after accounting expiry.
            dependency=prep if not (exp/'PREPARATION_COMPLETE.txt').exists() else None
            controller=ledger.ensure(stage,job_args(ledger,'controller',dependency=dependency))
            print('Preparation job:',prep,'Dispatcher job:',controller)
        print('Experiment:',exp)
        print('All 282 cases are frozen. The dispatcher submits at most',limits['max_active'],'model/check jobs concurrently.')
        print('Monitor: bash '+str(exp/'monitor_all_dcm_v8.sh')+' --diagnostics')
    return exp

def collect(exp):
    p=subprocess.run(['Rscript','--vanilla',str(Path(exp)/'04_collect_results.R'),str(exp)])
    if p.returncode:print('Collector returned',p.returncode,'; completed fits are preserved.',flush=True)

def control(experiment,once=False):
    ledger=Ledger(experiment)
    if not (ledger.exp/'PREPARATION_COMPLETE.txt').exists():
        raise RuntimeError('Preparation did not complete. No model fits were submitted by this dispatcher.')
    with lock(ledger.exp/'scheduler/dispatcher.lock'):
        cycles=0;started=time.monotonic()
        while True:
            limits=read_json(ledger.exp/'dispatch_limits.json',{'max_active':48,'max_cpus':1024})
            finished=dispatch_once(ledger,int(limits['max_active']),int(limits['max_cpus']))
            if cycles%5==0 or finished:collect(ledger.exp)
            if finished:
                atomic_text(ledger.exp/'CAMPAIGN_FINISHED.txt',utc()+' All accepted stages terminal; failures remain visible.\n');return
            if once:return
            # Leave time to hand over before this lightweight controller ends.
            limit=os.getenv('SLURM_JOB_END_TIME','')
            expiring=limit.isdigit() and int(limit)-time.time()<300
            if expiring or time.monotonic()-started>6*86400:
                generation=max([int(s.split('_')[1]) for s in ledger.records() if s.startswith('controller_')]+[0])+1
                stage='controller_%03d'%generation
                try:job=ledger.ensure(stage,job_args(ledger,'controller',dependency=os.getenv('SLURM_JOB_ID')))
                except QueueLimit:
                    # Requeue only the lightweight controller, never a model.
                    # This does not need a new submitted-job slot.
                    current=os.getenv('SLURM_JOB_ID','')
                    code,out,err=ledger.run(['scontrol','requeue',current]) if current.isdigit() else (1,'','Missing controller job ID')
                    if code:atomic_text(ledger.exp/'DISPATCHER_RESTART_NEEDED.txt',utc()+' No capacity for successor and controller requeue was unavailable; rerun queue_shape_dcm.sh. '+err+'\n')
                    return
                print('Successor dispatcher:',job,flush=True);return
            cycles+=1;time.sleep(45)

def status(experiment):
    exp=Path(experiment);print('Experiment:',exp)
    s=read_json(exp/'dispatch_status.json',{})
    print(json.dumps(s,indent=2))
    rejected=list((exp/'scheduler').glob('*.rejected.json'))
    if rejected:
        print('Rejected submissions (settings were not silently changed):')
        for p in rejected:print(p.name,read_json(p)['reason'])

def main():
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='command',required=True)
    launchp=sub.add_parser('launch');launchp.add_argument('--root',required=True);launchp.add_argument('--package',required=True);launchp.add_argument('--source',required=True);launchp.add_argument('--new',action='store_true')
    for name in ('control','status'):
        s=sub.add_parser(name);s.add_argument('experiment')
        if name=='control':s.add_argument('--once',action='store_true')
    a=p.parse_args()
    if a.command=='launch':launch(a.root,a.package,a.source,a.new)
    elif a.command=='control':control(a.experiment,a.once)
    else:status(a.experiment)
if __name__=='__main__':main()
