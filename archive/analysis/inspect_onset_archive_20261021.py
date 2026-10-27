"""Inspect and safely unpack the uploaded result archive; never execute its code."""
from pathlib import Path, PurePosixPath
import collections, json, tarfile, hashlib

root = Path(__file__).resolve().parent
archive = root / 'dcm_oct21_input/dcm_onset_compact_20261021T124243Z_2673664.tar.gz'
target = root / 'dcm_snapshot_20261021_124243'
target.mkdir(exist_ok=True)
inventory = []
with tarfile.open(archive, 'r:gz') as source:
    for member in source:
        if member.isdir():
            continue
        name = PurePosixPath(member.name)
        if not member.isfile() or name.is_absolute() or '..' in name.parts:
            raise ValueError('Unsafe archive member: ' + member.name)
        destination = target.joinpath(*name.parts)
        if not destination.resolve().is_relative_to(target.resolve()):
            raise ValueError('Path escapes the result folder')
        destination.parent.mkdir(parents=True, exist_ok=True)
        content = source.extractfile(member).read()
        destination.write_bytes(content)
        inventory.append({'name': member.name, 'bytes': len(content)})
(target / 'archive_inventory.json').write_text(json.dumps(inventory, indent=2), encoding='utf-8')
print('Archive SHA256:', hashlib.sha256(archive.read_bytes()).hexdigest())
print('Files:', len(inventory), 'Unpacked bytes:', sum(x['bytes'] for x in inventory))
print('Directory groups:', collections.Counter('/'.join(x['name'].split('/')[:2]) for x in inventory))
for row in inventory:
    if '/results/' in row['name'] or row['name'].endswith(('README.md', 'all_results.csv', 'cases.csv')):
        print(row['bytes'], row['name'])
