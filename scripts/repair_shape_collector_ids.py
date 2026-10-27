#!/usr/bin/env python3
"""Inspect or apply the documented identifier-only shape collector repair.

Usage:
  python repair_shape_collector_ids.py /absolute/04_collect_results.R
  python repair_shape_collector_ids.py /absolute/04_collect_results.R --apply

Default is inspection. --apply edits only the exact known old expression,
creates a unique backup, and preserves line endings. It never submits jobs.
"""
import argparse
import datetime
from pathlib import Path

OLD = '    key<-paste0(basename(campaign),":",sprintf("%03d",as.integer(id)))'
NEW = (
    '    id_text<-as.character(id)\n'
    '    id_key<-if(grepl("^[0-9]+$",id_text))sprintf("%03d",as.integer(id_text))else id_text\n'
    '    key<-paste0(basename(campaign),":",id_key)'
)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("collector", type=Path)
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    path = args.collector.resolve(strict=True)
    raw = path.read_bytes()
    text = raw.decode("utf-8")
    normalized = text.replace("\r\n", "\n")
    if NEW in normalized:
        print("Already repaired:", path)
        return
    if normalized.count(OLD) != 1:
        raise SystemExit("No unique exact legacy pattern; inspect manually. No file changed.")
    print("Legacy identifier conversion found:", path)
    if not args.apply:
        print("Inspection only. Use --apply to back up and apply the documented correction.")
        return
    ending = "\r\n" if b"\r\n" in raw else "\n"
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    backup = path.with_name(path.name + ".before_legacy_id_fix." + stamp)
    with backup.open("xb") as handle:
        handle.write(raw)
    updated = text.replace(OLD, NEW.replace("\n", ending), 1).encode("utf-8")
    # Atomic publication in the same directory, retaining the original backup.
    temporary = path.with_name(path.name + ".repair." + stamp)
    with temporary.open("xb") as handle:
        handle.write(updated)
    temporary.replace(path)
    print("Repaired:", path)
    print("Backup:", backup)

if __name__ == "__main__":
    main()
