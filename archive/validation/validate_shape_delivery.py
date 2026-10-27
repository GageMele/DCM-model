"""Validate the delivered installer and frozen design without running R/Slurm."""
from pathlib import Path
import ast,base64,contextlib,csv,hashlib,io,json,os,re,sys,tempfile,zipfile
root=Path(__file__).resolve().parent
text=(root/'queue_shape_dcm.sh').read_text(encoding='utf-8')
installer=text.split("<<'PY_INSTALL'\n",1)[1].split('\nPY_INSTALL',1)[0]
ast.parse(installer)
payload=base64.b64decode(re.search(r'payload=base64.b64decode\("""(.*?)"""\)',installer,re.S).group(1))
digest=hashlib.sha256(payload).hexdigest()
assert 'if digest!="'+digest+'"' in installer
with zipfile.ZipFile(io.BytesIO(payload)) as z:
    assert z.testzip() is None
    manifest=z.read('SHA256SUMS.txt').decode()
    for line in manifest.splitlines():
        h,name=line.split('  ',1);assert hashlib.sha256(z.read(name)).hexdigest()==h
    cases=list(csv.DictReader(io.StringIO(z.read('cases.tsv').decode()),delimiter='\t'))
    assert len(cases)==282
    assert sum(float(c['fraction'])==1 and c['mode']=='posterior' for c in cases)==11
    assert all(c['evaluation_role']=='training_reconstruction' for c in cases if float(c['fraction'])==1)
    for name in z.namelist():
        if name.endswith('.py'):ast.parse(z.read(name).decode(),filename=name)
        if name.endswith(('.R','.sh','.slurm','.py')):assert b'\r\n' not in z.read(name)
    assert "args[4]+=',CASE_ID='" in z.read('scheduler.py').decode()
    members=len(z.namelist())
with tempfile.TemporaryDirectory(prefix='shape_installer_check_',dir=str(root)) as tmp:
    target=Path(tmp).resolve()
    assert os.path.commonpath([str(target),str(root.resolve())])==str(root.resolve())
    before=sys.argv[:]
    try:
        sys.argv=['embedded_installer',str(target)]
        for repeat in range(2):
            with contextlib.redirect_stdout(io.StringIO()):exec(compile(installer,'embedded_installer','exec'),{})
        stable=target/'DCM_Ibex_shape_tests'
        for name in ('monitor_all_dcm_v8.sh','collect_results.sh','export_results.sh'):
            wrapper=(stable/name).read_text()
            assert wrapper.startswith('#!/bin/bash\nset -euo pipefail\n')
            assert 'latest_shape150_experiment.txt' in wrapper
    finally:sys.argv=before
result={'installer_bytes':len(text.encode()),'package_members':members,'cases':len(cases),
 'full_data_cases':11,'zip_and_file_hashes_valid':True,'python_sources_parse':True,
 'installer_first_and_repeat_execution':'PASS in an isolated temporary workspace',
 'R_and_Stan_validation':'Scheduled on Ibex before fitting; not executed locally'}
(root/'DCM_Ibex_shape_tests/DELIVERY_VALIDATION.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
