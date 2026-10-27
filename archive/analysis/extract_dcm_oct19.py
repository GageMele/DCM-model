"""Extract the user-provided snapshot after checking every archive path."""
from pathlib import Path
import tarfile,shutil,collections,json
ROOT=Path(__file__).resolve().parent
ARCHIVE=ROOT/'dcm_oct19_input/dcm_followup_comparison_20261019T084423Z_2726681.tar.gz'
DEST=ROOT/'dcm_snapshot_20261019_084423'
DEST.mkdir(exist_ok=True)
with tarfile.open(ARCHIVE,'r:gz') as t:
    members=t.getmembers()
    for m in members:
        p=(DEST/m.name).resolve()
        if not m.isfile() or not p.is_relative_to(DEST.resolve()): raise ValueError(m.name)
    for m in members:
        p=DEST/m.name;p.parent.mkdir(parents=True,exist_ok=True)
        with t.extractfile(m) as src,p.open('wb') as dst:shutil.copyfileobj(src,dst)
print(json.dumps({'files':len(members),'bytes':sum(m.size for m in members),'destination':str(DEST)}))
for p in sorted(DEST.rglob('all_results.csv')):print(p.relative_to(DEST))
for p in sorted(DEST.rglob('split_membership.csv')):print(p.relative_to(DEST))
