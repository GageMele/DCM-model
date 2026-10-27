#!/usr/bin/env python3
"""Capture the authenticated user's complete live Slurm queue, without changing jobs."""
import argparse
import csv
import datetime as dt
import getpass
import json
from pathlib import Path
import subprocess
import tarfile
import time

FIELDS = ["job_id", "job_name", "state", "partition", "user", "elapsed",
          "time_limit", "start_time", "cpus", "minimum_memory", "nodes",
          "node_list", "reason"]
FORMAT = "%i|%j|%T|%P|%u|%M|%l|%S|%C|%m|%D|%N|%R"

def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()

def capture(user):
    command = ["squeue", "--array", "--user", user, "--states=all",
               "--noheader", "--format", FORMAT]
    started = now()
    result = subprocess.run(command, text=True, capture_output=True,
                            check=True, timeout=60)
    ended = now()
    rows = []
    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        values = line.strip().split("|")
        if len(values) != len(FIELDS):
            raise RuntimeError("Unexpected delimiter/column count; raw queue must be inspected.")
        rows.append(dict(zip(FIELDS, values)))
    return started, ended, command, result.stdout, rows

def manifest_links(root):
    links = {}
    for directory in ("experiments", "benchmarks", "runs"):
        parent = root / directory
        if not parent.exists():
            continue
        for name in ("manifest.tsv", "existing_evaluations.tsv"):
            for path in sorted(parent.glob("*/" + name)):
                with path.open(encoding="utf-8-sig", newline="") as handle:
                    for row in csv.DictReader(handle, delimiter="\t"):
                        for key, value in row.items():
                            if not value or not (key.endswith("_job") or key in ("job", "job_id")):
                                continue
                            if not value.isdigit():
                                continue
                            record = {
                                "campaign": path.parent.name, "stage_field": key,
                                "case_id": row.get("case_id", ""),
                                "label": row.get("label", row.get("source_key", "")),
                                "manifest": str(path),
                                "run_dir": row.get("run_dir", ""),
                            }
                            links.setdefault(value, []).append(record)
    return links

def write_csv(path, columns, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--user", default=getpass.getuser())
    parser.add_argument("--project-root", type=Path, default=Path("/ibex/user/meleg/DCM_MODEL"))
    parser.add_argument("--output-root", type=Path)
    args = parser.parse_args()
    # A failed squeue call raises before any empty-success snapshot is written.
    start, end, command, raw, rows = capture(args.user)
    links = manifest_links(args.project_root)
    output_root = args.output_root or args.project_root / "logs"
    output_root.mkdir(parents=True, exist_ok=True)
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = output_root / ("ibex_jobs_" + stamp + "_" + str(time.time_ns()))
    out.mkdir()
    (out / "squeue_raw.txt").write_text(raw, encoding="utf-8")
    write_csv(out / "all_queued_jobs.csv", FIELDS, rows)
    running = [r for r in rows if r["state"] == "RUNNING"]
    write_csv(out / "running_jobs.csv", FIELDS, running)
    mapped = []
    for row in rows:
        # Array task IDs remain distinct; manifests normally use scalar IDs.
        for match in links.get(row["job_id"], []):
            mapped.append({**row, **match})
    write_csv(out / "dcm_manifest_matches.csv", FIELDS + [
        "campaign", "stage_field", "case_id", "label", "manifest", "run_dir"], mapped)
    metadata = {
        "query_started_utc": start, "query_finished_utc": end, "user": args.user,
        "command": command, "queued_rows": len(rows), "running_rows": len(running),
        "project_root": str(args.project_root),
        "scope": "All visible current queue states for this user; array tasks expanded. "
                 "Not cluster-wide jobs for other users. Snapshot may change after capture.",
        "mutations": "Only local snapshot files written; no Slurm submissions, cancellations or updates.",
    }
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    lines = ["# Live Ibex job snapshot", "", "Captured UTC: " + end,
             "", "All queued jobs: " + str(len(rows)) + "; RUNNING: " + str(len(running)),
             "", "| Job ID | Name | State | CPUs | Elapsed | Nodes |",
             "|---|---|---|---:|---|---|"]
    for row in rows:
        lines.append("| " + " | ".join(row[k].replace("|", "/") for k in
                     ("job_id", "job_name", "state", "cpus", "elapsed", "node_list")) + " |")
    (out / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    archive = out.with_suffix(".tar.gz")
    with tarfile.open(archive, "w:gz") as handle:
        handle.add(out, arcname=out.name)
    print(json.dumps(metadata, indent=2))
    print("Snapshot:", out)
    print("Upload this archive:", archive)

if __name__ == "__main__":
    main()
