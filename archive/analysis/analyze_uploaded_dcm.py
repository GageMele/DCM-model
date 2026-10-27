"""Inventory the user-supplied DCM snapshot; never modify original model runs."""
from pathlib import Path
import csv, json, tarfile, collections

ROOT = Path(__file__).resolve().parent
ARCHIVE = ROOT / 'all_results_input' / 'dcm_all_results_20261016T090217Z_81828.tar.gz'
DEST = ROOT / 'dcm_snapshot_20261016_090217'

def extract():
    DEST.mkdir(exist_ok=True)
    with tarfile.open(ARCHIVE, 'r:gz') as archive:
        members = archive.getmembers()
        for member in members:
            target = (DEST / member.name).resolve()
            if not target.is_relative_to(DEST.resolve()) or not member.isfile():
                raise ValueError('Unsafe archive member: ' + member.name)
        for member in members:
            target = DEST / member.name
            target.parent.mkdir(parents=True, exist_ok=True)
            with archive.extractfile(member) as source, target.open('wb') as output:
                while chunk := source.read(1024 * 1024):
                    output.write(chunk)
    print('Extracted files:', len(members))
    print('Campaign files:', collections.Counter('/'.join(m.name.split('/')[:2]) for m in members))
    manifests = sorted(DEST.rglob('manifest.tsv'))
    for path in manifests:
        with path.open(newline='', encoding='utf-8-sig') as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        print(path.relative_to(DEST), len(rows), list(rows[0]) if rows else [])
    for filename in ['FIT_COMPLETE.txt', 'test_summary.csv', 'diagnostic_summary.csv', 'CHECKS_COMPLETE.txt']:
        print(filename, len(list(DEST.rglob(filename))))

if __name__ == '__main__':
    extract()
