"""Build a self-contained installer; validate its bytes and frozen design offline.

This script never starts R, Stan or network clients. Ibex runs the R checks.
"""
from pathlib import Path,PurePosixPath
import base64,csv,hashlib,io,json,tarfile,textwrap,collections,math

ROOT=Path(__file__).resolve().parent
PKG=ROOT/'DCM_Ibex_followup_tests'
def read(path,sep=','):
    with path.open(newline='',encoding='utf-8-sig') as f:return list(csv.DictReader(f,delimiter=sep))
cases=read(PKG/'cases.tsv','\t')
assert len(cases)==34 and len({c['label'] for c in cases})==34
assert sum(c['mode']=='posterior' for c in cases)==30
assert [int(c['case_id']) for c in cases]==list(range(1,35))
assert all(int(c['cpus'])==int(c['chains'])*int(c['threads']) for c in cases)
assert all(c['subset'] in ('p10','p10_alt','p20') for c in cases)
assert all(int(c['iter'])==int(c['warmup'])+int(c['samples']) for c in cases)
assert all(int(c['treedepth'])==14 for c in cases)
assert {int(c['seed']) for c in cases[:6]}=={1234,20261018}
assert collections.Counter(c['equation'] for c in cases[:6])==dict(original_z=2,both_z=2,interaction_z=2)
assert collections.Counter(c['block'] for c in cases)==dict(stability=6,delta=3,threads=3,training_size=3,sigma_control=2,sigma_depth=2,sigma_prior=2,date_blocked=3,training_composition=3,depth_basis=1,effect_basis_low=1,effect_basis_high=1,prior_sigma=4)
for a,b in [(17,16),(20,19)]:
    x,y=cases[a-1],cases[b-1]
    assert all(x[k]==y[k] for k in ('equation','samples','warmup','adapt_delta','seed','treedepth','threads','sigma_smooth_sd'))
    assert x['variance']=='region_depth' and y['variance']=='region'
for a,b in [(18,17),(21,20)]:
    x,y=cases[a-1],cases[b-1]
    assert all(x[k]==y[k] for k in ('equation','variance','samples','warmup','adapt_delta','seed','treedepth','threads'))
    assert float(x['sigma_smooth_sd'])==.25 and float(y['sigma_smooth_sd'])==.5

# Inspect the completed data-only audit rather than rerunning the blocked R runtime.
marker=ROOT/'validation_followup150/latest_data_check.txt'
AUDIT=Path(marker.read_text(encoding='utf-8').strip())
membership=read(AUDIT/'split_membership.csv')
sets=collections.defaultdict(set);dates=collections.defaultdict(set)
for x in membership:
    sets[x['dataset']].add(x['cast_key'])
    dates[x['dataset']].add((x['Ocean'],x['Region'],x['date']))
expected={'primary_p10_d150':108,'primary_p10_alt_d150':108,'primary_p20_d150':216,'primary_p50_d150':540,
 'primary_validation_d150':108,'date_blocked_p10_d150':108,'date_blocked_validation_d150':117,
 'primary_reserve_d150':108,'date_blocked_reserve_d150':108}
assert {k:len(v) for k,v in sets.items()}==expected
assert not sets['primary_p10_d150']&sets['primary_p10_alt_d150']
assert sets['primary_p10_d150']|sets['primary_p10_alt_d150']==sets['primary_p20_d150']
excluded=set().union(*(sets[k] for k in ['primary_p50_d150','date_blocked_p10_d150','primary_validation_d150','date_blocked_validation_d150']))
assert not sets['primary_reserve_d150']&excluded
assert sets['primary_reserve_d150']==sets['date_blocked_reserve_d150']
assert not dates['date_blocked_p10_d150']&dates['date_blocked_validation_d150']
for r in read(AUDIT/'source_checksums.csv'):
    assert hashlib.md5((AUDIT/r['file']).read_bytes()).hexdigest()==r['md5']
    assert hashlib.md5(Path(r['source']).read_bytes()).hexdigest()==r['md5']

# Do not duplicate an existing posterior case with the same statistical settings.
fields=('equation','subset','split','prior','mode','variance','cast','family','autocor','depth_max','chains','warmup','samples','adapt_delta','treedepth','threads','k_depth','k_mld','k_distance','seed','sigma_smooth_sd')
def canonical(r):
    z=dict(r)
    defaults={'split':'primary','variance':'constant','cast':'none','family':'gaussian','autocor':'none','depth_max':'150','sigma_smooth_sd':'.5'}
    for k,v in defaults.items():
        if z.get(k) in ('','NA',None):z[k]=v
    if not z.get('samples') and z.get('iter') and z.get('warmup'):z['samples']=str(int(z['iter'])-int(z['warmup']))
    result=[]
    for k in fields:
        value=z.get(k,'')
        try:value=float(value)
        except (ValueError,TypeError):pass
        result.append(value)
    return tuple(result)
old=[]
for p in (ROOT/'dcm_snapshot_20261018_105731/experiments').glob('*/cases.tsv'):old+=read(p,'\t')
old_signatures={canonical(r) for r in old}
assert all(canonical(c) not in old_signatures for c in cases if c['mode']=='posterior')

# Check delimiters/strings without interpreting or executing R.
def check_r_delimiters(path):
    text=path.read_text(encoding='utf-8-sig');stack=[];quote=None;escaped=False;comment=False
    for ch in text:
        if comment:
            if ch=='\n':comment=False
            continue
        if quote:
            if escaped:escaped=False
            elif ch=='\\':escaped=True
            elif ch==quote:quote=None
            continue
        if ch=='#':comment=True
        elif ch in ('"',"'",'`'):quote=ch
        elif ch in '({[':stack.append(ch)
        elif ch in ')}]':assert stack and stack.pop()=={')':'(',']':'[','}':'{'}[ch],path
    assert not quote and not stack,path
for p in PKG.glob('*.R'):check_r_delimiters(p)

validation=f'''# Validation status

The installer contains34 cases:30 posterior and4 prior-only. Python checks passed for case IDs, unique labels, resource arithmetic, block counts, sampler iterations, matched variance controls, prior controls, and lack of exact posterior duplicates among{len(old)} earlier registered cases.

Existing completed data-only audit files were independently inspected: original108 and alternate108 training casts form the216-cast pool; reserve108 is outside primary540, date-blocked training108, and both previous holdouts108/117; date groups remain disjoint for date-blocked training/validation. All copied source files retain identical MD5 checksums. The completed data-only R audit compared raw response/predictor rows before writing its completion marker; it did not fit models.

R delimiters/strings and installer payload/checksums were checked without executing R. These are static checks, not a complete R parser or a substitute for running the included tests.

Automatic approval review blocked local R scheduling/reporting validation after detecting an attempted HTTPS connection to an external Visual Studio telemetry endpoint. It was not retried or bypassed. Windows also prevented local Bash startup with an OS named-object permission error. Local R scheduler tests and full Bash syntax validation are therefore not claimed as passed.

On Ibex, the installer runs bash -n and all base-R scheduling/reporting recovery tests before submitting jobs. The preparation job then constructs every actual formula/prior/data combination before releasing fits through afterok dependencies. Failure prevents fitting. These runtime gates are required, not optional.

Local data audit datasets: {json.dumps(expected,sort_keys=True)}

Aggregate requested capacity if every fitting job ran simultaneously: {sum(int(c['cpus']) for c in cases)} CPUs and{sum(int(c['memory_gb']) for c in cases)}GB. Slurm controls actual concurrency. No runtime/total-cost forecast is implied.
'''
(PKG/'VALIDATION.md').write_text(validation,encoding='utf-8',newline='\n')

for p in PKG.iterdir():
    if p.is_file() and p.suffix in ('.R','.slurm','.sh','.awk','.md','.tsv','.csv'):
        p.write_text(p.read_text(encoding='utf-8-sig'),encoding='utf-8',newline='\n')
files=[p for p in sorted(PKG.iterdir()) if p.is_file() and p.name!='SHA256SUMS.txt']
(PKG/'SHA256SUMS.txt').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.name+'\n' for p in files),encoding='ascii',newline='\n')
files.append(PKG/'SHA256SUMS.txt')
buffer=io.BytesIO()
with tarfile.open(fileobj=buffer,mode='w:gz',format=tarfile.USTAR_FORMAT) as t:
    for p in files:
        info=tarfile.TarInfo(PKG.name+'/'+p.name);data=p.read_bytes();info.size=len(data);info.mode=0o644;info.mtime=1792310400
        t.addfile(info,io.BytesIO(data))
payload=buffer.getvalue();digest=hashlib.sha256(payload).hexdigest()
header='''#!/bin/bash
# Install and submit34 focused DCM tests; repeating resumes existing registrations.
# Use --extract-only to install for inspection without queueing.
set -euo pipefail
bash -n "$0"
mode="${1:-submit}"
[[ "$mode" == submit || "$mode" == --extract-only ]] || { echo 'Usage: bash queue_followup_dcm.sh [--extract-only]' >&2; exit 2; }
root="${DCM_ROOT:-/ibex/user/meleg/DCM_MODEL}"
[[ -d "$root" ]] || { echo "Project directory missing: $root" >&2; exit 1; }
cd -- "$root"
root="$PWD"
package="$root/DCM_Ibex_followup_tests"
mkdir -p logs
for command in base64 tar sha256sum; do command -v "$command" >/dev/null || { echo "Missing $command" >&2; exit 1; }; done
if [[ "$mode" == submit ]]; then command -v sbatch >/dev/null || { echo 'Run this on an Ibex login node.' >&2; exit 1; }; fi
payload=$(mktemp "$root/logs/followup_payload.XXXXXXXX.tar.gz")
trap 'rm -f -- "$payload"' EXIT
base64 --decode > "$payload" <<'__DCM_FOLLOWUP_PAYLOAD__'
'''
footer='''
__DCM_FOLLOWUP_PAYLOAD__
printf '%s  %s\\n' '__DIGEST__' "$payload" | sha256sum --check --status
if [[ -L "$package" ]]; then echo 'Package path is a symbolic link; refusing extraction.' >&2; exit 1; fi
if [[ ! -e "$package" ]]; then
    staging=$(mktemp -d "$root/.dcm_followup_install.XXXXXXXX")
    tar --extract --gzip --file "$payload" --directory "$staging" --no-same-owner
    mv -- "$staging/DCM_Ibex_followup_tests" "$package"
    rmdir -- "$staging"
fi
[[ -d "$package" && -f "$package/SHA256SUMS.txt" ]] || { echo 'Existing package is incomplete; preserve it for inspection before installing.' >&2; exit 1; }
expected=$(tar -xOzf "$payload" DCM_Ibex_followup_tests/SHA256SUMS.txt | sha256sum | cut -d' ' -f1)
actual=$(sha256sum "$package/SHA256SUMS.txt" | cut -d' ' -f1)
[[ "$expected" == "$actual" ]] || { echo 'Existing package version differs; no files overwritten or jobs submitted.' >&2; exit 1; }
(cd -- "$package"; sha256sum --check --quiet SHA256SUMS.txt; for file in *.sh *.slurm; do bash -n "$file"; done)
printf '\\nInstalled and verified: %s\\nPlan: %s/PLAN.md\\n' "$package" "$package"
if [[ "$mode" == --extract-only ]]; then echo 'No jobs submitted. To queue: bash DCM_Ibex_followup_tests/submit_tests.sh'; exit 0; fi
bash "$package/submit_tests.sh"
'''.replace('__DIGEST__',digest)
installer=header+'\n'.join(textwrap.wrap(base64.b64encode(payload).decode('ascii'),76))+footer
destination=ROOT/'queue_followup_dcm.sh';destination.write_text(installer,encoding='ascii',newline='\n')
saved=destination.read_text(encoding='ascii');start=saved.index("<<'__DCM_FOLLOWUP_PAYLOAD__'\n")+len("<<'__DCM_FOLLOWUP_PAYLOAD__'\n");end=saved.index('\n__DCM_FOLLOWUP_PAYLOAD__\n',start)
decoded=base64.b64decode(saved[start:end]);assert hashlib.sha256(decoded).hexdigest()==digest
with tarfile.open(fileobj=io.BytesIO(decoded),mode='r:gz') as t:
    members=t.getmembers();assert len(members)==len(files)
    for member in members:
        p=PurePosixPath(member.name);assert not p.is_absolute() and '..' not in p.parts and p.parts[0]==PKG.name and member.isfile()
        assert t.extractfile(member).read()==(ROOT/member.name).read_bytes()
assert b'\r\n' not in destination.read_bytes()
print('PASS:34 cases; all controls/resources; frozen split independence and source hashes; zero exact old posterior duplicates; balanced R delimiters; safe archive paths; exact payload/source checksums and Unix newlines.')
print('INSTALLER',destination.name,destination.stat().st_size,'bytes; FILES',len(files),'PAYLOAD_SHA256',digest)
