#!/bin/bash
# Diagnose why the shape150 controller stopped. Read-only: submits, cancels and edits nothing.
set -uo pipefail
root=/ibex/user/meleg/DCM_MODEL
cd "$root" || exit 1
experiment=$(cat logs/latest_shape150_experiment.txt)
echo "=== experiment: $experiment"

echo; echo "=== 1. CONTROLLER JOB ACCOUNTING"
sacct -j 54250665 --format=JobID,JobName%40,State,ExitCode,Start,End,Elapsed,Timelimit,MaxRSS,ReqMem -P
echo "--- any other controller/check/prep jobs since 22 Oct:"
sacct -u "$USER" -X -n -P --starttime=2026-10-22 --format=JobID,JobName,State,Elapsed,ExitCode \
  | grep -iE "controller|prepare|collect" | tail -n 20

echo; echo "=== 2. CONTROLLER LOG TAILS"
for f in $(ls "$experiment"/logs/*54250665* "$root"/logs/*54250665* 2>/dev/null); do
  echo "--- $f"; tail -n 40 "$f"
done

echo; echo "=== 3. SCHEDULER DIRECTORY (newest 15 files) AND LOCKS"
ls -lt --time-style=+%Y-%m-%dT%H:%M "$experiment/scheduler" 2>/dev/null | head -n 16
find "$experiment" -maxdepth 2 -name "*.lock*" -printf "%TY-%Tm-%TdT%TH:%TM %p\n" 2>/dev/null

echo; echo "=== 4. WHAT ARE THE PENDING JOBS? (name pattern, reason, count)"
squeue -u "$USER" -h -t PENDING -o "%j|%r" | sed -E 's/[0-9]+/N/g' | sort | uniq -c | sort -rn | head -n 15

echo; echo "=== 5. CASE ID LISTS"
python3 - "$experiment" <<'PY'
import csv, glob, os, sys
exp = sys.argv[1]
with open(os.path.join(exp, "cases.tsv")) as fh:
    ids = [r["case_id"] for r in csv.DictReader(fh, delimiter="\t")]
groups = {"fit_done_unchecked": [], "not_finished": [], "no_run_dir": [], "check_failed": [], "fit_failed": []}
for c in ids:
    d = os.path.join(exp, "runs", "case_%03d" % int(c))
    has = lambda f: os.path.exists(os.path.join(d, f))
    if not os.path.isdir(d): groups["no_run_dir"].append(c)
    elif has("CHECKS_COMPLETE.txt"): continue
    elif has("CHECKS_FAILED.txt"): groups["check_failed"].append(c)
    elif has("FIT_FAILED.txt"): groups["fit_failed"].append(c)
    elif has("FIT_COMPLETE.txt"): groups["fit_done_unchecked"].append(c)
    else: groups["not_finished"].append(c)
for k, v in groups.items():
    print("%s (%d): %s" % (k, len(v), " ".join(v)))
PY
echo; echo "=== END"
