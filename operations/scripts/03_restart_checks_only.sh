#!/bin/bash
# Restart the shape150 controller in "checks only" mode.
# 1) Marks every not-yet-submitted FIT as deferred (scheduler/fit_NNN.rejected.json),
#    which the frozen scheduler already honours: it skips those stages. Reversible:
#    deleting a marker re-enables that fit. No running job is touched.
# 2) Submits one new lightweight controller (1 CPU, 4G). It will score the 28
#    finished-but-unchecked fits, keep tracking the 19 running fits, and refresh summaries.
set -euo pipefail
root=/ibex/user/meleg/DCM_MODEL
cd "$root"
experiment=$(cat logs/latest_shape150_experiment.txt)
echo "experiment: $experiment"

python3 - "$experiment" <<'PY'
import sys, json, datetime
from pathlib import Path
exp = Path(sys.argv[1]); sys.path.insert(0, str(exp))
import scheduler as s

ledger = s.Ledger(exp)
records = ledger.records()
live = s.live_jobs()

# Safety: refuse if any controller is already queued or running.
running = [(k, r["job_id"]) for k, r in records.items() if k.startswith("controller_") and r["job_id"] in live]
if running:
    sys.exit("A controller is already queued/running: %s. Nothing changed." % running)

# 1) Defer every fit that was never submitted.
now = datetime.datetime.now(datetime.timezone.utc).isoformat()
deferred = []
for c in s.read_cases(exp):
    cid = int(c["case_id"]); tag = "fit_%03d" % cid
    run = exp / "runs" / ("case_%03d" % cid)
    if tag in records or run.exists():
        continue
    deferred.append(cid)
if len(deferred) != 141:
    sys.exit("Expected 141 unsubmitted fits, found %d; nothing changed. Paste this output to Claude." % len(deferred))
for cid in deferred:
    tag = "fit_%03d" % cid
    marker = ledger.directory / (tag + ".rejected.json")
    if not marker.exists():
        s.atomic_json(marker, {"utc": now, "reason": "DEFERRED by user 2026-10-27 (trimmed plan); delete this file to re-enable", "args": []})
print("deferred fits (%d): %s" % (len(deferred), " ".join(map(str, deferred))))

# 2) Submit the next controller generation, exactly as the launcher would.
generation = max([int(k.split("_")[1]) for k in records if k.startswith("controller_")] + [0]) + 1
stage = "controller_%03d" % generation
job = ledger.ensure(stage, s.job_args(ledger, "controller"))
print("submitted", stage, "job", job)
PY

echo; echo "--- queue check (controller should appear):"
sleep 5
squeue -u "$USER" -o "%.10i %.45j %.9T %.10M" | grep -iE "JOBID|controller|check" || true
