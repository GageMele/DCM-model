#!/bin/bash
# DCM status snapshot, read-only: no jobs submitted, no fits loaded.
# Prints a compact summary to paste back, and writes a full CSV of every
# completed shape150 case to logs/ for later download.
set -uo pipefail
root=/ibex/user/meleg/DCM_MODEL
cd "$root" || exit 1
experiment=$(cat logs/latest_shape150_experiment.txt)
stamp=$(date -u +%Y%m%dT%H%M%SZ)
out="$root/logs/status_snapshot_$stamp.csv"

echo "=== SNAPSHOT $stamp  experiment: $(basename "$experiment")"

echo; echo "=== 1. QUEUE (all my jobs, by state)"
squeue -u "$USER" -h -o "%T" | sort | uniq -c
echo "--- running DCM fits (case, elapsed, limit)"
squeue -u "$USER" -h -t RUNNING -o "%j %M %l" | grep -i dcm | sed -E 's/.*(fit|check)_?([0-9]+)/\1_\2/' | sort

echo; echo "=== 2. SLURM ACCOUNTING since 2026-10-22 (DCM jobs)"
sacct -u "$USER" -X -n -P --starttime=2026-10-22 --format=JobID,JobName,State,Elapsed \
  | grep -i dcm > /tmp/dcm_sacct_$$.txt
cut -d'|' -f3 /tmp/dcm_sacct_$$.txt | awk '{print $1}' | sort | uniq -c
echo "--- not COMPLETED/RUNNING/PENDING:"
awk -F'|' '$3!~/COMPLETED|RUNNING|PENDING/' /tmp/dcm_sacct_$$.txt | tail -n 40
rm -f /tmp/dcm_sacct_$$.txt

echo; echo "=== 3. SHAPE150 CASES + NEW RESULTS SINCE 25 OCT"
python3 - "$experiment" "$out" <<'PY'
import csv, glob, os, sys, datetime
exp, out = sys.argv[1], sys.argv[2]
cut = datetime.datetime(2026, 10, 25, 11, 15).timestamp()
counts, rows, fails = {}, [], []
for d in sorted(glob.glob(os.path.join(exp, "runs", "case_*"))):
    case = os.path.basename(d)[5:]
    has = lambda f: os.path.exists(os.path.join(d, f))
    if has("CHECKS_COMPLETE.txt"): s = "checked"
    elif has("CHECKS_FAILED.txt"): s = "check_failed"
    elif has("FIT_FAILED.txt"): s = "fit_failed"
    elif has("FIT_COMPLETE.txt"): s = "fit_done_unchecked"
    else: s = "not_finished"
    counts[s] = counts.get(s, 0) + 1
    if "failed" in s: fails.append(case + ":" + s)
    p = os.path.join(d, "test_summary.csv")
    if os.path.exists(p):
        with open(p) as fh:
            for r in csv.DictReader(fh):
                r["status"] = s
                r["new_since_oct25"] = os.path.getmtime(p) > cut
                rows.append(r)
print("status counts:", counts)
print("failed:", " ".join(fails) or "none")
if rows:
    keys = []
    for r in rows:
        for k in r:
            if k not in keys: keys.append(k)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
    print("full table:", out, "(%d rows)" % len(rows))
def f(r, k, fmt="%.3f"):
    try: return fmt % float(r.get(k, ""))
    except ValueError: return "-"
cols = ["case", "equation", "subset", "split", "seed", "prior", "k", "kmld", "delta",
        "div", "th", "rhat", "essB", "hours", "valRMSE", "resRMSE", "valPop", "resPop", "valCov90", "valOnset"]
print("\nnew completed cases (since 25 Oct 11:15 UTC):")
print("\t".join(cols))
for r in sorted(rows, key=lambda r: int(r["case_id"])):
    if not r["new_since_oct25"]: continue
    print("\t".join([r["case_id"], r.get("equation", ""), r.get("subset", ""), r.get("split", ""),
        r.get("selection_seed", ""), r.get("prior", ""), r.get("k_depth", "") + "/" + r.get("k_profile", ""),
        r.get("k_mld", ""), r.get("adapt_delta", ""), f(r, "divergences", "%.0f"), f(r, "treedepth_hits", "%.0f"),
        f(r, "max_rhat", "%.4f"), f(r, "min_bulk_ess", "%.0f"), f({"h": str(float(r["fit_seconds"]) / 3600) if r.get("fit_seconds", "").replace(".", "", 1).isdigit() else ""}, "h", "%.1f"),
        f(r, "validation_conditional_rmse"), f(r, "reserve_conditional_rmse"),
        f(r, "validation_population_rmse"), f(r, "reserve_population_rmse"),
        f(r, "validation_conditional_coverage90"),
        f(r, "validation_conditional_onset50_onset_mae_given_both_detected", "%.1f")]))
PY

echo; echo "=== 4. SCHEDULER"
python3 "$experiment/scheduler.py" status "$experiment" 2>&1 | tail -n 25
echo; echo "=== END"
