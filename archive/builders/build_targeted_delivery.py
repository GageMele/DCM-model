"""Build and verify the single-file Ibex installer, without running it."""
from pathlib import Path, PurePosixPath
import base64
import hashlib
import io
import shutil
import tarfile
import textwrap

root = Path(__file__).resolve().parent
package = root / 'DCM_Ibex_targeted_tests'
review = package / 'original_review'
review.mkdir(exist_ok=True)
for source in sorted((root / 'original_model_review').glob('*')):
    if source.is_file() and source.suffix in {'.csv', '.png', '.rds'}:
        shutil.copyfile(source, review / source.name)
shutil.copyfile(root / 'validation_targeted150/actual/split_counts.csv',
                package / 'verified_local_split_counts.csv')

# Text entry points use Unix newlines even when built on Windows.
for path in package.iterdir():
    if path.is_file() and path.suffix in {'.R', '.sh', '.slurm', '.awk', '.md', '.tsv'}:
        text = path.read_text(encoding='utf-8-sig')
        path.write_text(text, encoding='utf-8', newline='\n')
files = [p for p in sorted(package.rglob('*')) if p.is_file() and p.name != 'SHA256SUMS.txt']
assert files and all(not p.is_symlink() for p in files)
checksums = ''.join(f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(package).as_posix()}\n' for p in files)
(package / 'SHA256SUMS.txt').write_text(checksums, encoding='ascii', newline='\n')
files.append(package / 'SHA256SUMS.txt')

buffer = io.BytesIO()
with tarfile.open(fileobj=buffer, mode='w:gz', format=tarfile.USTAR_FORMAT) as archive:
    for path in files:
        name = f'{package.name}/{path.relative_to(package).as_posix()}'
        contents = path.read_bytes()
        info = tarfile.TarInfo(name)
        info.size = len(contents)
        info.mode = 0o644
        info.mtime = 1792108800
        archive.addfile(info, io.BytesIO(contents))
payload = buffer.getvalue()
digest = hashlib.sha256(payload).hexdigest()
encoded = '\n'.join(textwrap.wrap(base64.b64encode(payload).decode('ascii'), 76))
header = r'''#!/bin/bash
# Self-contained DCM follow-up: 57 posterior + 6 prior-predictive fits.
# Default: install and submit. --extract-only installs for review without submission.
set -euo pipefail
bash -n "$0"
mode="${1:-submit}"
if [[ "$mode" != submit && "$mode" != --extract-only ]]; then
    echo "Usage: bash queue_targeted_dcm.sh [--extract-only]" >&2
    exit 2
fi
root="${DCM_ROOT:-/ibex/user/meleg/DCM_MODEL}"
[[ -d "$root" ]] || { echo "Project directory missing: $root" >&2; exit 1; }
cd -- "$root"
root="$PWD"
package="$root/DCM_Ibex_targeted_tests"
if [[ -e "$package" || -L "$package" ]]; then
    echo "Package already exists: $package" >&2
    echo "To monitor: bash DCM_Ibex_targeted_tests/monitor_all_dcm_v5.sh --diagnostics" >&2
    echo "To add old-fit predictions only: bash DCM_Ibex_targeted_tests/rescore_completed.sh" >&2
    echo "To intentionally start a fresh 63-fit campaign: bash DCM_Ibex_targeted_tests/submit_tests.sh" >&2
    exit 1
fi
for command in base64 tar sha256sum; do
    command -v "$command" >/dev/null || { echo "Missing command: $command" >&2; exit 1; }
done
if [[ "$mode" == submit ]]; then
    command -v sbatch >/dev/null || { echo "Run this installer on an Ibex login node." >&2; exit 1; }
    if [[ -z "${SOURCE_EXPERIMENT:-}" && ! -f logs/latest_extensive150_experiment.txt ]]; then
        echo "Missing frozen extensive150 source pointer." >&2
        exit 1
    fi
fi
mkdir -p logs
payload=$(mktemp "$root/logs/targeted_payload.XXXXXXXX.tar.gz")
trap 'rm -f -- "$payload"' EXIT
base64 --decode > "$payload" <<'__DCM_TARGETED_PAYLOAD__'
'''
footer = r'''
__DCM_TARGETED_PAYLOAD__
printf '%s  %s\n' '__PAYLOAD_SHA256__' "$payload" | sha256sum --check --status
tar --extract --gzip --file "$payload" --directory "$root" --no-same-owner
(
    cd -- "$package"
    sha256sum --check --quiet SHA256SUMS.txt
    for file in *.sh *.slurm; do bash -n "$file"; done
)
printf '\nInstalled: %s\n' "$package"
printf 'Plan: %s/PLAN.md\nOriginal fit review: %s/ORIGINAL_REVIEW.md\n' "$package" "$package"
if [[ "$mode" == --extract-only ]]; then
    echo "Extraction complete. No jobs submitted."
    echo "To submit: bash DCM_Ibex_targeted_tests/submit_tests.sh"
    exit 0
fi
bash "$package/submit_tests.sh"
'''.replace('__PAYLOAD_SHA256__', digest)
installer = header + encoded + footer
destination = root / 'queue_targeted_dcm.sh'
destination.write_text(installer, encoding='ascii', newline='\n')

# Independently decode the deliverable and compare every archived byte/hash.
saved = destination.read_text(encoding='ascii')
start = saved.index("<<'__DCM_TARGETED_PAYLOAD__'\n") + len("<<'__DCM_TARGETED_PAYLOAD__'\n")
end = saved.index('\n__DCM_TARGETED_PAYLOAD__\n', start)
decoded = base64.b64decode(saved[start:end], validate=False)
assert hashlib.sha256(decoded).hexdigest() == digest
with tarfile.open(fileobj=io.BytesIO(decoded), mode='r:gz') as archive:
    members = archive.getmembers()
    assert len(members) == len(files)
    for member in members:
        name = PurePosixPath(member.name)
        assert not name.is_absolute() and '..' not in name.parts
        assert name.parts[0] == package.name and member.isfile()
        assert archive.extractfile(member).read() == (root / member.name).read_bytes()
assert b'\r\n' not in destination.read_bytes()
print(f'PASS: {len(files)} archive files, exact payload checksum, safe archive paths, source bytes and Unix newlines.')
print(f'Installer: {destination.name}; {destination.stat().st_size:,} bytes; payload SHA256 {digest}')
