"""Offline tests with fake Slurm; no real submissions or network access."""
import csv,importlib.util,json,os,tempfile,unittest
from pathlib import Path
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('shape_scheduler',HERE/'scheduler.py')
s=importlib.util.module_from_spec(spec);spec.loader.exec_module(s)
class Fake:
    def __init__(self):self.calls=[];self.next=100;self.live={};self.account={};self.reject=None;self.uncertain=False
    def __call__(self,args):
        self.calls.append(list(args))
        if args[0]=='squeue':return 0,''.join(f'{j}|{n}|RUNNING\n' for j,n in self.live.items()),''
        if args[0]=='sacct':
            if any(x.startswith('--name=') for x in args):
                name=next(x.split('=',1)[1] for x in args if x.startswith('--name='))
                return 0,''.join(f'{j}|{n}\n' for j,n in self.live.items() if n==name),''
            return 0,''.join(f'{j}|{v}\n' for j,v in self.account.items()),''
        if args[0]=='scontrol':return 0,'PartitionName=batch MaxTime=7-00:00:00',''
        if args[0]=='sbatch':
            if self.reject:return 1,'','sbatch: error: '+self.reject
            if self.uncertain:return 1,'','Connection lost before acknowledgement'
            self.next+=1;job=str(self.next);self.live[job]=next(x.split('=',1)[1] for x in args if x.startswith('--job-name='))
            return 0,job+';ibex\n',''
        raise AssertionError(args)
class SchedulerTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='shape_scheduler_test_',dir=str(HERE.parent))
        self.exp=Path(self.temp.name).resolve()
        assert os.path.commonpath([str(self.exp),str(HERE.parent.resolve())])==str(HERE.parent.resolve())
        self.fake=Fake();s.atomic_json(self.exp/'campaign.json',dict(token='test',root=str(HERE.parent),
          created_utc='2026-10-22T00:00:00+00:00',partition='batch',constraint='intel&cascadelake',
          prepare_time='08:00:00',controller_time='7-00:00:00',check_time='2-00:00:00',fit_time='7-00:00:00',prior_time='04:00:00'))
        self.ledger=s.Ledger(self.exp,self.fake)
        with (HERE/'cases.tsv').open(encoding='utf-8') as f:self.cases=list(csv.DictReader(f,delimiter='\t'))[:3]
        with (self.exp/'cases.tsv').open('w',encoding='utf-8',newline='') as f:
            w=csv.DictWriter(f,list(self.cases[0]),delimiter='\t');w.writeheader();w.writerows(self.cases)
        (self.exp/'preflight_cases.csv').write_text('case_id,valid,reason\n1,TRUE,\n2,TRUE,\n3,TRUE,\n')
    def tearDown(self):
        # The absolute temporary target was checked inside the workspace above.
        self.temp.cleanup()
    def submitted(self):return [c for c in self.fake.calls if c[0]=='sbatch']
    def test_idempotent_accepted_fit(self):
        a=self.ledger.ensure('fit_001',[]);self.assertEqual(a,self.ledger.ensure('fit_001',[]));self.assertEqual(len(self.submitted()),1)
    def test_lost_ack_recovery(self):
        self.fake.uncertain=True
        with self.assertRaises(s.UncertainSubmission):self.ledger.ensure('fit_001',[])
        self.fake.live['907']=self.ledger.name('fit_001');self.fake.uncertain=False
        self.assertEqual(self.ledger.ensure('fit_001',[]),'907');self.assertEqual(len(self.submitted()),1)
    def test_unresolved_ack_not_retried(self):
        self.fake.uncertain=True
        with self.assertRaises(s.UncertainSubmission):self.ledger.ensure('fit_001',[])
        self.fake.uncertain=False
        with self.assertRaises(s.UncertainSubmission):self.ledger.ensure('fit_001',[])
        self.assertEqual(len(self.submitted()),1)
    def test_explicit_quota_rejection(self):
        self.fake.reject='QOSMaxSubmitJobPerUserLimit'
        with self.assertRaises(s.QueueLimit):self.ledger.ensure('fit_001',[])
        self.assertFalse((self.exp/'scheduler/fit_001.intent.json').exists())
        self.fake.reject=None;self.assertEqual(self.ledger.ensure('fit_001',[]),'101')
    def test_environment_and_resources(self):
        args=s.job_args(self.ledger,'fit',self.cases[0])
        self.assertIn('CASE_ID=1',next(x for x in args if x.startswith('--export=')))
        self.assertNotIn('CASE_ID',next(x for x in args if x.startswith('--output=')))
        self.assertIn('--cpus-per-task=16',args)
    def test_slot_limit_and_repeat(self):
        s.dispatch_once(self.ledger,2,512);self.assertEqual(len(self.submitted()),2)
        s.dispatch_once(self.ledger,2,512);self.assertEqual(len(self.submitted()),2)
    def test_cpu_limit(self):
        s.dispatch_once(self.ledger,24,16);self.assertEqual(len(self.submitted()),1)
    def test_missing_squeue_row_reserves_slot(self):
        s.dispatch_once(self.ledger,1,512);self.fake.live={}
        s.dispatch_once(self.ledger,1,512);self.assertEqual(len(self.submitted()),1)
    def test_timeout_never_refits(self):
        job=self.ledger.ensure('fit_001',[]);self.fake.live={};self.fake.account[job]='TIMEOUT'
        s.dispatch_once(self.ledger,1,512)
        self.assertTrue((self.exp/'runs/case_001/FIT_FAILED.txt').exists())
        self.assertEqual(sum('--job-name='+self.ledger.name('fit_001') in x for x in self.submitted()),1)
    def test_checks_prioritized(self):
        job=self.ledger.ensure('fit_001',[]);self.fake.live={};self.fake.account[job]='COMPLETED'
        s.atomic_text(self.exp/'runs/case_001/FIT_COMPLETE.txt','complete');s.dispatch_once(self.ledger,1,512)
        self.assertIn('--job-name='+self.ledger.name('check_001'),self.submitted()[-1])
    def test_partition_time(self):self.assertEqual(s.partition_time('batch','14-00:00:00',self.fake),'7-00:00:00')
class DesignTests(unittest.TestCase):
    def test_design(self):
        with (HERE/'cases.tsv').open(encoding='utf-8') as f:c=list(csv.DictReader(f,delimiter='\t'))
        self.assertEqual(len(c),282)
        self.assertEqual({float(r['fraction']) for r in c},{.1,.2,.4,.5,.75,1.})
        self.assertEqual({int(r['k_depth']) for r in c},{7,9,11,13,15})
        full=[r for r in c if float(r['fraction'])==1]
        self.assertEqual(len(full),11);self.assertTrue(all(r['evaluation_role']=='training_reconstruction' for r in full))
        repeat=[r for r in c if r['block']=='cast_selection']
        self.assertEqual({r['selection_seed'] for r in repeat},{'20261023','20261024','20261025'})
        self.assertEqual({r['seed'] for r in repeat},{'1234'})
        self.assertTrue(any(r['k_profile']=='rs15_other7' for r in c))
if __name__=='__main__':unittest.main(verbosity=2)
