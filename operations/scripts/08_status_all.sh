#!/bin/bash
# Read-only status of the round3 and bins1m campaigns (plus a shape150 count).
# Prints per-case status and scores, controller state and the tail of any failed log.
root=/ibex/user/meleg/DCM_MODEL
cd "$root" || exit 1
for name in round3 bins1m; do
  exp=$(cat logs/latest_${name}_experiment.txt)
  tok=$(python3 -c "import json;print(json.load(open('$exp/campaign.json'))['token'])")
  echo; echo "================ $name  $(basename "$exp")  token $tok"
  echo "--- jobs (sacct, since 9 Nov)"
  sacct -u "$USER" -X -n -P --starttime=2026-11-09 --format=JobID,JobName,State,Elapsed | grep "$tok" \
    | sed "s/dcm-s150-$tok-//" | awk -F'|' '{printf "%s %s %s %s\n",$2,$1,$3,$4}' | sort
  python3 - "$exp" <<'PY'
import csv, os, sys, glob
exp = sys.argv[1]
with open(os.path.join(exp, "cases.tsv")) as fh: cases = list(csv.DictReader(fh, delimiter="\t"))
def f(r, k, fmt="%.3f"):
    try: return fmt % float(r.get(k, ""))
    except (TypeError, ValueError): return "-"
print("\n--- cases")
print("\t".join(["case", "label", "status", "h", "rhat", "div", "valRMSE", "resRMSE", "valCov90", "valLPD", "valELPD", "dbRPS/ dcmMAE", "onset"]))
for c in cases:
    d = os.path.join(exp, "runs", "case_%03d" % int(c["case_id"]))
    has = lambda n: os.path.exists(os.path.join(d, n))
    s = ("checked" if has("CHECKS_COMPLETE.txt") else "check_failed" if has("CHECKS_FAILED.txt") else
         "fit_failed" if has("FIT_FAILED.txt") else "fit_done" if has("FIT_COMPLETE.txt") else
         "fitting" if os.path.isdir(d) else "not_started")
    r = {}
    for n in ("test_summary.csv", "diagnostic_summary.csv"):
        if has(n):
            with open(os.path.join(d, n)) as fh: r = next(csv.DictReader(fh), {})
            break
    try: h = "%.1f" % (float(r["fit_seconds"]) / 3600)
    except (KeyError, ValueError): h = "-"
    extra = f(r, "validation_conditional_dcm_depth_mae_m", "%.1fm") if c["equation"] == "dcm_depth" else f(r, "validation_conditional_rps")
    print("\t".join([c["case_id"], c["label"][5:45], s, h, f(r, "max_rhat", "%.4f"), f(r, "divergences", "%.0f"),
        f(r, "validation_conditional_rmse"), f(r, "reserve_conditional_rmse"), f(r, "validation_conditional_coverage90"),
        f(r, "validation_conditional_lpd_per_row"), f(r, "validation_conditional_elpd_date_joint_sum", "%.0f"), extra,
        f(r, "validation_conditional_onset50_onset_mae_given_both_detected", "%.1f")]))
    for n in ("FIT_FAILED.txt", "CHECKS_FAILED.txt", "FIGURES_FAILED.txt"):
        if has(n): print("    " + n + ": " + open(os.path.join(d, n)).read().strip().replace("\n", " | ")[:400])
PY
  echo "--- controller log tail"
  ls -t "$exp"/logs/*controller*.err 2>/dev/null | head -1 | xargs -r tail -n 8
  ls -t "$exp"/logs/*controller*.out 2>/dev/null | head -1 | xargs -r tail -n 5
  [ -f "$exp/DISPATCHER_RESTART_NEEDED.txt" ] && { echo "!! DISPATCHER_RESTART_NEEDED:"; cat "$exp/DISPATCHER_RESTART_NEEDED.txt"; }
done
echo; echo "================ shape150"
s=$(cat logs/latest_shape150_experiment.txt)
for st in CHECKS_COMPLETE CHECKS_FAILED FIT_FAILED FIT_COMPLETE; do
  printf "%s: %s\n" $st $(ls "$s"/runs/case_*/$st.txt 2>/dev/null | wc -l); done
echo "run dirs: $(ls -d "$s"/runs/case_* | wc -l) of $(($(wc -l < "$s/cases.tsv")-1)) cases"
echo "=== END"
