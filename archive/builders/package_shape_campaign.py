"""Create one LF-terminated, self-extracting launcher for upload to Ibex."""
from pathlib import Path
import base64,hashlib,io,json,textwrap,zipfile,ast
ROOT=Path(__file__).resolve().parent;PKG=ROOT/'DCM_Ibex_shape_tests'
files=[p for p in PKG.iterdir() if p.is_file() and p.name!='SHA256SUMS.txt']
for p in files:
    if p.suffix in ('.py','.R','.sh','.slurm','.awk','.tsv','.json','.md','.txt'):
        p.write_bytes(p.read_bytes().replace(b'\r\n',b'\n'))
    if p.suffix=='.py':ast.parse(p.read_text(encoding='utf-8'),filename=p.name)
checksums=''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.name+'\n' for p in sorted(files))
(PKG/'SHA256SUMS.txt').write_bytes(checksums.encode())
buffer=io.BytesIO()
with zipfile.ZipFile(buffer,'w',zipfile.ZIP_DEFLATED) as z:
    for p in sorted(files+[PKG/'SHA256SUMS.txt']):z.writestr(p.name,p.read_bytes())
payload=buffer.getvalue();digest=hashlib.sha256(payload).hexdigest()
encoded='\n'.join(textwrap.wrap(base64.b64encode(payload).decode(),76))
launcher='''#!/bin/bash
# Self-contained MLD-shape campaign. Run on Ibex from DCM_MODEL.
set -euo pipefail
export LC_ALL=C PYTHONUTF8=1
export DCM_ROOT="${DCM_ROOT:-/ibex/user/meleg/DCM_MODEL}"
mkdir -p "$DCM_ROOT/logs"
cd -- "$DCM_ROOT"
command -v python3 >/dev/null || { echo 'Python 3 is required on Ibex.' >&2; exit 1; }
package=$(python3 - "$DCM_ROOT" <<'PY_INSTALL'
from pathlib import Path,PurePosixPath
import base64,hashlib,io,os,zipfile
import sys
root=Path(sys.argv[1]).resolve()
payload=base64.b64decode("""__PAYLOAD__""")
digest=hashlib.sha256(payload).hexdigest()
if digest!="__DIGEST__":raise SystemExit('Embedded package checksum failed')
package=root/'DCM_Ibex_shape_tests'/'releases'/digest[:16]
package.mkdir(parents=True,exist_ok=True)
with zipfile.ZipFile(io.BytesIO(payload)) as z:
    for member in z.infolist():
        name=PurePosixPath(member.filename)
        if name.is_absolute() or '..' in name.parts or len(name.parts)!=1:
            raise SystemExit('Unsafe embedded package member')
        path=package/name.name
        if os.path.commonpath([str(path.resolve()),str(package.resolve())])!=str(package.resolve()):
            raise SystemExit('Embedded path escapes package directory')
        data=z.read(member)
        if path.exists() and path.read_bytes()!=data:
            raise SystemExit('Existing frozen release differs: '+str(path))
        if not path.exists():path.write_bytes(data)
for line in (package/'SHA256SUMS.txt').read_text().splitlines():
    expected,name=line.split('  ',1)
    if hashlib.sha256((package/name).read_bytes()).hexdigest()!=expected:
        raise SystemExit('Package file checksum failed: '+name)
stable=root/'DCM_Ibex_shape_tests'
for name in ('monitor_all_dcm_v8.sh','collect_results.sh','export_results.sh'):
    text='#!/bin/bash\\nset -euo pipefail\\nroot="${DCM_ROOT:-/ibex/user/meleg/DCM_MODEL}"\\nexperiment=$(cat "$root/logs/latest_shape150_experiment.txt")\\nexec bash "$experiment/'+name+'" "$@"\\n'
    path=stable/name
    if path.exists() and path.read_text()!=text:
        raise SystemExit('Existing stable wrapper differs; no file overwritten: '+str(path))
    if not path.exists():path.write_text(text)
print(package)
PY_INSTALL
)
source_experiment="${SHAPE_SOURCE:-$(cat "$DCM_ROOT/logs/latest_followup150_experiment.txt")}"
source "$source_experiment/ibex_environment.sh"
# The tests use fake Slurm responses; they submit no real model jobs.
python3 "$package/test_scheduler.py"
arguments=(launch --root "$DCM_ROOT" --package "$package" --source "$source_experiment")
if [[ "${SHAPE_NEW:-0}" == 1 ]]; then arguments+=(--new); fi
python3 "$package/scheduler.py" "${arguments[@]}"
'''.replace('__PAYLOAD__',encoded).replace('__DIGEST__',digest)
path=ROOT/'queue_shape_dcm.sh';path.write_bytes(launcher.encode('utf-8'))
compile(launcher.split("<<'PY_INSTALL'\n",1)[1].split('\nPY_INSTALL',1)[0],'embedded_installer','exec')
with zipfile.ZipFile(io.BytesIO(payload)) as z:assert z.testzip() is None
print(json.dumps({'launcher':str(path),'bytes':path.stat().st_size,'package_sha256':digest,'files':len(files)+1}))
