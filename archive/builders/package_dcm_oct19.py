"""Bundle audited data and annotated R plotting code, excluding raw models."""
from pathlib import Path
import zipfile,json
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'DCM_review_20261019'
dest=OUT/'DCM_comparison_data_and_R_figures_20261019.zip'
files=sorted([p for p in OUT.iterdir() if p.is_file() and (p.suffix in ('.csv','.svg','.R','.md','.json')) and p.name!='audit_console.txt'])
with zipfile.ZipFile(dest,'w',zipfile.ZIP_DEFLATED) as z:
    for p in files:z.write(p,'DCM_review_20261019/'+p.name)
with zipfile.ZipFile(dest) as z:
    assert z.testzip() is None
    assert 'DCM_review_20261019/reproduce_figures.R' in z.namelist()
    assert 'DCM_review_20261019/all_registered_runs.csv' in z.namelist()
print(json.dumps({'file':str(dest),'files':len(files),'bytes':dest.stat().st_size}))
