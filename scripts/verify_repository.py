#!/usr/bin/env python3
"""Offline integrity checks. Does not access Ibex, fit models, or submit jobs."""
import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path
import re
from urllib.parse import unquote

def read_csv(path, delimiter=","):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))

def verify(root):
    errors = []
    checked_links = 0
    docs = [root / "README.md", root / "CONTRIBUTING.md", root / "CITATION.md"]
    docs += sorted((root / "docs").glob("*.md"))
    for path in docs:
        text = path.read_text(encoding="utf-8")
        # Ignore code examples when checking Markdown destinations.
        text = re.sub(r"```.*?```", "", text, flags=re.S)
        for link in re.findall(r"!?\[[^\]]*\]\(([^)]+)\)", text):
            if re.match(r"(https?://|mailto:|#)", link):
                continue
            target = unquote(link.split("#", 1)[0])
            if not target:
                continue
            checked_links += 1
            if not (path.parent / target).exists():
                errors.append(f"Missing link: {path.relative_to(root)} -> {link}")
    registry = read_csv(root / "results/20261025/all_586_registrations.csv")
    if len(registry) != 586:
        errors.append("Expected 586 registry rows")
    if len({r["run_key"] for r in registry}) != len(registry):
        errors.append("Duplicate registry keys")
    current = read_csv(root / "results/20261025/all_282_new_cases.csv")
    if len(current) != 282:
        errors.append("Expected 282 shape rows")
    counts = collections.Counter((r.get("mode"), r.get("status")) for r in current)
    if counts[("posterior", "complete")] != 67 or counts[("prior", "complete")] != 12:
        errors.append("Unexpected current completed mode counts")
    all_cases = read_csv(root / "provenance/all_cases.tsv", "\t")
    if len(all_cases) != 581:
        errors.append("Expected 581 frozen campaign cases")
    if len({(r["package"], r["case_id"]) for r in all_cases}) != len(all_cases):
        errors.append("Duplicate package/case identifiers")
    shape_cases = read_csv(root / "DCM_Ibex_shape_tests/cases.tsv", "\t")
    if sum(r["mode"] == "prior" for r in shape_cases) != 12:
        errors.append("Expected 12 shape prior cases")
    if sum(float(r["fraction"]) == 1 for r in shape_cases) != 11:
        errors.append("Expected 11 full-data shape cases")
    if any(float(r["fraction"]) < .1 for r in shape_cases if r["mode"] == "posterior"):
        errors.append("A new shape posterior trains on less than 10%")
    inventory = read_csv(root / "provenance/source_inventory.csv")
    for row in inventory:
        path = root / row["repository_path"]
        if not path.is_file():
            errors.append("Missing copied source: " + row["repository_path"])
        elif hashlib.sha256(path.read_bytes()).hexdigest() != row["sha256"]:
            errors.append("Changed copied source: " + row["repository_path"])
    checksums = root / "provenance/repository_checksums.json"
    if checksums.exists():
        for relative, expected in json.loads(checksums.read_text(encoding="utf-8")).items():
            path = root / relative
            if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != expected:
                errors.append("Repository hash mismatch: " + relative)
    live = (root / "evidence/shape_20261025/campaign/04_collect_results.R").read_text(encoding="utf-8")
    if 'id_text<-as.character(id)' not in live or 'else id_text' not in live:
        errors.append("Live collector identifier repair missing")
    return {
        "status": "PASS" if not errors else "FAIL",
        "checked_local_links": checked_links,
        "source_files_verified": len(inventory),
        "registered_rows": len(registry),
        "frozen_campaign_cases": len(all_cases),
        "current_shape_cases": len(current),
        "completed_shape_posterior": counts[("posterior", "complete")],
        "completed_shape_prior": counts[("prior", "complete")],
        "errors": errors,
        "scope": "Offline documentation, source integrity, configuration/count checks only; no new R/Stan/HPC execution or full statistical re-audit."
    }

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    report = verify(args.root.resolve())
    result = json.dumps(report, indent=2) + "\n"
    if args.report:
        args.report.write_text(result, encoding="utf-8")
    print(result)
    raise SystemExit(0 if report["status"] == "PASS" else 1)

if __name__ == "__main__":
    main()
