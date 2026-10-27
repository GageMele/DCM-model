"""Bundle only local reviewed tables, figures and annotated reproducibility code."""
from pathlib import Path
import zipfile,json
root=Path(__file__).resolve().parent; out=root/'DCM_review_20261021'
archive=out/'DCM_comparison_tables_and_code_20261021.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
    for path in sorted(out.iterdir()):
        if path.is_file() and path.suffix in ('.csv','.json','.svg','.R','.md') and path.name!='pairwise_audit.json':
            z.write(path,'DCM_review_20261021/'+path.name)
    for name in ('analyze_dcm_oct21.py','enrich_dcm_oct21.py','deliver_dcm_oct21.py','inspect_onset_archive_20261021.py','make_dcm_workbook.py'):
        z.write(root/name,'audit_code/'+name)
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None
    assert 'DCM_review_20261021/all_registered_runs.csv' in z.namelist()
print(json.dumps({'archive':str(archive),'bytes':archive.stat().st_size}))
