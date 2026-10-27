"""Package tested scheduler fixes in one self-contained, non-refitting repair."""
from pathlib import Path, PurePosixPath
import base64
import hashlib
import io
import tarfile
import textwrap

root = Path(__file__).resolve().parent
package = root / 'DCM_Ibex_targeted_tests'
names = ['submit.R', 'submit_existing.R', 'test_submission.R', 'test_existing.R',
         'test_dependency_recovery.R', 'submit_tests.sh', 'repair_existing_campaign.R',
         'DEPENDENCY_FIX_20261016.md']
contents = {n: (package / n).read_bytes() for n in names}
contents['PATCH_SHA256SUMS.txt'] = ''.join(
    f'{hashlib.sha256(contents[n]).hexdigest()}  {n}\n' for n in names).encode('ascii')
buffer = io.BytesIO()
with tarfile.open(fileobj=buffer, mode='w:gz', format=tarfile.USTAR_FORMAT) as archive:
    for name, data in contents.items():
        info = tarfile.TarInfo(name)
        info.size = len(data)
        info.mode = 0o644
        info.mtime = 1792144800
        archive.addfile(info, io.BytesIO(data))
payload = buffer.getvalue()
digest = hashlib.sha256(payload).hexdigest()
encoded = '\n'.join(textwrap.wrap(base64.b64encode(payload).decode('ascii'), 76))
header = r'''#!/bin/bash
# Repair interrupted evaluation scheduling. This script submits no model fits.
# Default target is the existing experiment from the supplied terminal log.
set -euo pipefail
export LC_ALL=C
bash -n "$0"
root="${DCM_ROOT:-/ibex/user/meleg/DCM_MODEL}"
cd -- "$root"
root=$(pwd -P)
experiment="${EXPERIMENT_DIR:-$root/experiments/targeted150_20261016T104550Z_296126}"
experiment=$(cd -- "$experiment" && pwd -P)
case "$experiment" in
  "$root"/experiments/targeted150_*) ;;
  *) echo "Unexpected experiment location: $experiment" >&2; exit 1 ;;
esac
for file in PREPARATION_COMPLETE.txt manifest.tsv cases.tsv existing_evaluations.tsv ibex_environment.sh; do
    [[ -f "$experiment/$file" ]] || { echo "Required prepared-campaign file missing: $file" >&2; exit 1; }
done
for command in sbatch squeue base64 tar sha256sum; do
    command -v "$command" >/dev/null || { echo "Missing command: $command" >&2; exit 1; }
done
lock="$experiment/.dependency_repair.lock"
mkdir -- "$lock" || { echo "Another repair may be running: $lock" >&2; exit 1; }
trap 'rmdir -- "$lock" 2>/dev/null || true' EXIT
staging=$(mktemp -d "$root/logs/dcm_dependency_repair.XXXXXXXX")
base64 --decode > "$staging/patch.tar.gz" <<'__DCM_DEPENDENCY_PATCH__'
'''
footer = r'''
__DCM_DEPENDENCY_PATCH__
printf '%s  %s\n' '__PAYLOAD_SHA256__' "$staging/patch.tar.gz" | sha256sum --check --status
tar --extract --gzip --file "$staging/patch.tar.gz" --directory "$staging" --no-same-owner
(
    cd -- "$staging"
    sha256sum --check --quiet PATCH_SHA256SUMS.txt
    bash -n submit_tests.sh
)
# The scientific specification is read only for the scheduler's fixture tests.
cp -- "$experiment/campaign.R" "$staging/campaign.R"
source "$experiment/ibex_environment.sh"
Rscript --vanilla - "$staging" "$experiment" <<'RS_PREFLIGHT'
a<-commandArgs(trailingOnly=TRUE)
for(f in list.files(a[1],pattern="\\.R$",full.names=TRUE))parse(f)
m<-read.delim(file.path(a[2],"manifest.tsv"),stringsAsFactors=FALSE)
stopifnot(nrow(m)==63L,setequal(m$case_id,1:63),
  length(unique(c(m$fit_job,m$check_job)))==126L,
  all(grepl("^[0-9]+$",c(m$fit_job,m$check_job))))
cat("PASS: existing 63 fit/check pairs verified; no replacement fitting jobs needed.\n")
RS_PREFLIGHT
for test in test_submission.R test_existing.R test_dependency_recovery.R; do
    Rscript --vanilla "$staging/$test" "$staging"
done
stamp=$(date -u +%Y%m%dT%H%M%SZ)_$$
targets=("$experiment")
if [[ -d "$root/DCM_Ibex_targeted_tests" ]]; then
    [[ ! -L "$root/DCM_Ibex_targeted_tests" ]] || { echo "Installed package must be a directory, not a symlink." >&2; exit 1; }
    targets+=("$root/DCM_Ibex_targeted_tests")
fi
for target in "${targets[@]}"; do
    backup="$target/scheduler_backup_$stamp"
    mkdir -- "$backup"
    for file in submit.R submit_existing.R test_submission.R test_existing.R test_dependency_recovery.R submit_tests.sh repair_existing_campaign.R DEPENDENCY_FIX_20261016.md; do
        if [[ -f "$target/$file" ]]; then cp -p -- "$target/$file" "$backup/$file"; fi
        cp -- "$staging/$file" "$target/$file"
    done
    printf 'Scheduler backup: %s\n' "$backup"
done
export DCM_ROOT="$root" EXPERIMENT_DIR="$experiment"
Rscript --vanilla "$experiment/repair_existing_campaign.R" "$experiment" "$root"
printf '\nRepair finished. No model fits were submitted or cancelled.\n'
printf 'Monitor: bash DCM_Ibex_targeted_tests/monitor_all_dcm_v5.sh --diagnostics\n'
'''.replace('__PAYLOAD_SHA256__', digest)
text = header + encoded + footer
destination = root / 'repair_targeted_dcm_dependencies.sh'
destination.write_text(text, encoding='ascii', newline='\n')

# Verify the actual delivered payload and safe, ordinary file paths.
saved = destination.read_text(encoding='ascii')
start = saved.index("<<'__DCM_DEPENDENCY_PATCH__'\n") + len("<<'__DCM_DEPENDENCY_PATCH__'\n")
end = saved.index('\n__DCM_DEPENDENCY_PATCH__\n', start)
decoded = base64.b64decode(saved[start:end], validate=False)
assert hashlib.sha256(decoded).hexdigest() == digest
with tarfile.open(fileobj=io.BytesIO(decoded), mode='r:gz') as archive:
    assert len(archive.getmembers()) == len(contents)
    for member in archive.getmembers():
        path = PurePosixPath(member.name)
        assert member.isfile() and not path.is_absolute() and '..' not in path.parts
        assert archive.extractfile(member).read() == contents[member.name]
assert saved.count('__DCM_DEPENDENCY_PATCH__') == 2
assert b'\r\n' not in destination.read_bytes()
print(f'PASS: {len(contents)} patch files; exact payload SHA256; safe paths; Unix newlines.')
print(f'{destination.name}: {destination.stat().st_size:,} bytes; payload {digest}')
