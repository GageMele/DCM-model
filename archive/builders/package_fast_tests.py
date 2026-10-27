"""Package only the new exploratory scripts, with Linux-compatible text files."""
from pathlib import Path
import hashlib
import zipfile

root = Path.cwd()
package = root / "DCM_Ibex_fast_tests"
files = sorted(path for path in package.iterdir() if path.is_file())
assert all(path.suffix in {".R", ".sh", ".slurm", ".tsv", ".md", ".txt"} for path in files)
for path in files:
    # Preserve the generated code, while removing Windows newlines/BOM if present.
    value = path.read_text(encoding="utf-8-sig").replace("\r\n", "\n")
    path.write_bytes(value.encode("utf-8"))
archive = root / "DCM_Ibex_fast_tests_150m.zip"
with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as output:
    for path in files:
        output.write(path, arcname=f"DCM_Ibex_fast_tests/{path.name}")
with zipfile.ZipFile(archive) as check:
    assert check.testzip() is None
    assert len(check.namelist()) == len(files)
    assert all(name.startswith("DCM_Ibex_fast_tests/") for name in check.namelist())
print(f"Created {archive.name}: {len(files)} files; {archive.stat().st_size} bytes")
print(f"SHA256: {hashlib.sha256(archive.read_bytes()).hexdigest()}")
