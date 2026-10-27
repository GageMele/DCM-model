"""Offline import of the supplied 2026-10-27 queue snapshot; no Slurm calls."""
import collections, csv, hashlib, io, json, pathlib, tarfile
ROOT = pathlib.Path("DCM-model").resolve()
ARCHIVE = pathlib.Path("github_handover_inputs/jobs_20261027/ibex_jobs_20261027T104740Z_1793098060392524450.tar.gz").resolve()
BASE = pathlib.Path("operations/jobs/20261027T104740Z")
EXPECTED = "d49aae3377f68af41d2fa664872ca1c4061a570a6c7d2ec46aa05cbdc010ddcc"
ARCHIVE_ROOT = ARCHIVE.name.removesuffix(".tar.gz")
CHANGED, checks = [], []
def check(name, condition):
    checks.append({"check": name, "passed": bool(condition)})
    if not condition: raise RuntimeError(name)
def sha(data): return hashlib.sha256(data).hexdigest()
def write(rel, data):
    if isinstance(data, str): data = data.encode("utf-8")
    path = ROOT / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)
    rel = pathlib.Path(rel).as_posix()
    if rel not in CHANGED: CHANGED.append(rel)
def write_json(rel, obj): write(rel, json.dumps(obj, indent=2, ensure_ascii=False) + "\n")
def write_csv(rel, fields, rows):
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
    writer.writeheader(); writer.writerows(rows)
    write(rel, stream.getvalue())
def read_csv_bytes(data, sep=","):
    return list(csv.DictReader(io.StringIO(data.decode("utf-8-sig")), delimiter=sep))
def markdown(headers, rows):
    escape = lambda x: str(x).replace("|", "\\|").replace("\n", " ")
    return "| " + " | ".join(headers) + " |\n|" + "|".join("---" for _ in headers) + "|\n" + "".join("| " + " | ".join(map(escape, row)) + " |\n" for row in rows)
check("Original archive SHA256", sha(ARCHIVE.read_bytes()) == EXPECTED)
expected_files = {"README.md","all_queued_jobs.csv","dcm_manifest_matches.csv","metadata.json","running_jobs.csv","squeue_raw.txt"}
raw = {}
with tarfile.open(ARCHIVE) as handle:
    for member in handle:
        path = pathlib.PurePosixPath(member.name)
        if member.isdir():
            check("Safe archive directory: " + member.name, path.parts == (ARCHIVE_ROOT,))
            continue
        check("Safe regular archive member: " + member.name,
              member.isfile() and len(path.parts) == 2 and path.parts[0] == ARCHIVE_ROOT and path.name in expected_files and path.name not in raw)
        raw[path.name] = handle.extractfile(member).read()
check("Exactly six expected archive files", set(raw) == expected_files)
queue = read_csv_bytes(raw["all_queued_jobs.csv"])
running = read_csv_bytes(raw["running_jobs.csv"])
matches = read_csv_bytes(raw["dcm_manifest_matches.csv"])
meta = json.loads(raw["metadata.json"])
fields = list(queue[0])
queue_by_id = {r["job_id"]: r for r in queue}
match_by_id = collections.defaultdict(list)
for row in matches: match_by_id[row["job_id"]].append(row)
check("Metadata equals earlier reported terminal JSON", meta == json.loads((ROOT/BASE/"reported_capture_metadata.json").read_bytes()))
check("335 unique whole-account job IDs", len(queue) == len(queue_by_id) == meta["queued_rows"] == 335)
check("20 RUNNING rows exactly match whole-account RUNNING subset", running == [r for r in queue if r["state"] == "RUNNING"] and len(running) == meta["running_rows"] == 20)
check("All rows belong to captured account", all(r["user"] == meta["user"] for r in queue))
raw_values = [line.strip().split("|") for line in raw["squeue_raw.txt"].decode().splitlines() if line.strip()]
check("All raw squeue rows have 13 columns", len(fields)==13 and all(len(row)==13 for row in raw_values))
check("Raw squeue exactly equals complete CSV, including order", [dict(zip(fields,row)) for row in raw_values] == queue)
check("Manifest matches reproduce exact queue fields", all(m["job_id"] in queue_by_id and all(m[f] == queue_by_id[m["job_id"]][f] for f in fields) for m in matches))
check("19 unique DCM matches, one association each", len(matches) == len(match_by_id) == 19 and all(len(v) == 1 for v in match_by_id.values()))
check("Only observed dcm-prefixed jobs are the 19 manifest matches", {r["job_id"] for r in queue if r["job_name"].startswith("dcm-")} == set(match_by_id))
cases_rel = "evidence/shape_20261025/campaign/cases.tsv"
case_bytes = (ROOT/cases_rel).read_bytes()
case_rows = read_csv_bytes(case_bytes, "\t")
cases = {r["case_id"]:r for r in case_rows}
check("Archived case matrix matches preserved delivered matrix", case_bytes == (ROOT/"DCM_Ibex_shape_tests/cases.tsv").read_bytes())
campaign = "shape150_20261022T102851Z_3426481"
check("All DCM matches belong to shape fit jobs and exact frozen case labels",
      all(m["campaign"] == campaign and m["stage_field"] == "fit_job" and m["state"] == "RUNNING" and m["case_id"] in cases and m["label"] == cases[m["case_id"]]["label"] and m["job_name"] == "dcm-s150-77caaba437-fit_" + m["case_id"].zfill(3) for m in matches))
check("Recorded fit CPU, memory and walltime agree with frozen plan", all(
    m["cpus"] == cases[m["case_id"]]["cpus"] and m["minimum_memory"] == cases[m["case_id"]]["memory_gb"]+"G" and m["time_limit"] == cases[m["case_id"]]["fit_days"]+"-00:00:00" for m in matches))
state_counts = collections.Counter(r["state"] for r in queue)
reason_counts = collections.Counter(r["reason"] for r in queue if r["state"]=="PENDING")
check("Exact recorded state totals", dict(state_counts) == {"PENDING":315,"RUNNING":20})
check("Exact recorded pending reasons", dict(reason_counts) == {"(Dependency)":310,"(Priority)":1,"(DependencyNeverSatisfied)":4})
for name, data in sorted(raw.items()): write(BASE/"raw"/name, data)
classified = []
for row in queue:
    links = match_by_id.get(row["job_id"], [])
    classified.append(dict(row,
        dcm_classification="manifest_matched" if links else "not_matched",
        dcm_evidence="exact_job_id_in_captured_dcm_manifest_matches" if links else "no_captured_manifest_match",
        campaign=";".join(m["campaign"] for m in links), stage_field=";".join(m["stage_field"] for m in links),
        case_id=";".join(m["case_id"] for m in links), label=";".join(m["label"] for m in links)))
write_csv(BASE/"classified_queue.csv", list(classified[0]), classified)
enriched = [dict(m, **{"plan_"+k:v for k,v in cases[m["case_id"]].items()}) for m in sorted(matches,key=lambda x:int(x["case_id"]))]
write_csv(BASE/"dcm_running_with_settings.csv", list(enriched[0]), enriched)
pending = [r for r in queue if r["state"] == "PENDING"]
write_csv(BASE/"pending_jobs.csv", fields, pending)
unmatched = [r for r in running if r["job_id"] not in match_by_id]
write_csv(BASE/"running_without_dcm_match.csv", fields, unmatched)
def name_group(name):
    for prefix in ("dcm-","C_D1_","VR_","V_"):
        if name.startswith(prefix): return prefix
    return name
groups = collections.Counter((name_group(r["job_name"]),r["state"]) for r in queue)
summary = {
 "query_started_utc":meta["query_started_utc"],"query_finished_utc":meta["query_finished_utc"],
 "account":meta["user"],"queue_rows":len(queue),"state_counts":dict(state_counts),
 "dcm_manifest_matched_rows":len(matches),"dcm_manifest_matched_unique_jobs":len(match_by_id),
 "dcm_running_fits":len(matches),"dcm_pending_matches":0,"running_without_dcm_match":len(unmatched),
 "visible_dcm_controllers":sum(r["job_name"].startswith("dcm-") and "controller" in r["job_name"] for r in queue),
 "visible_dcm_checks":sum(r["job_name"].startswith("dcm-") and "-check_" in r["job_name"] for r in queue),
 "dcm_running_cpus":sum(int(r["cpus"]) for r in matches),"whole_account_running_cpus":sum(int(r["cpus"]) for r in running),
 "pending_reason_counts":dict(reason_counts),
 "name_prefix_state_counts":[{"name_prefix":g,"state":s,"count":n} for (g,s),n in sorted(groups.items())],
 "running_dcm_equations":dict(collections.Counter(cases[m["case_id"]]["equation"] for m in matches)),
 "running_dcm_training_fractions":dict(collections.Counter(cases[m["case_id"]]["fraction"] for m in matches)),
 "limitations":[
    "Dated queue observation, not a live feed or complete Slurm accounting history.",
    "DCM classification uses captured exact manifest associations; prefixes do not establish other project ownership.",
    "plan_* settings come from the frozen case matrix; no new run_config, diagnostics, predictions, ledger or completion markers were supplied.",
    "No DCM controller is visible; this does not establish its terminal state, cause, or whether a restart is needed.",
    "Absence of other jobs does not establish completed or failed models.",
    "start_time is as emitted by Slurm without an offset; query start/end explicitly record UTC.",
    "CPU and minimum-memory fields are scheduler allocations/requests, not observed utilization.",
    "The October 25 scientific result snapshot remains unchanged."
 ]}
write_json(BASE/"summary.json", summary)
write_json(BASE/"provenance.json", {
 "source":"User-provided Ibex queue capture archive",
 "archive_filename":ARCHIVE.name,"archive_bytes":ARCHIVE.stat().st_size,"archive_sha256":EXPECTED,
 "original_archive_root":ARCHIVE_ROOT,"reported_remote_path":meta["project_root"]+"/logs/"+ARCHIVE.name,
 "capture_script_repository_path":"scripts/capture_ibex_jobs.py",
 "capture_script_repository_sha256":sha((ROOT/"scripts/capture_ibex_jobs.py").read_bytes()),
 "raw_files":{name:{"bytes":len(data),"sha256":sha(data)} for name,data in sorted(raw.items())},
 "case_matrix_repository_path":cases_rel,"case_matrix_sha256":sha(case_bytes),
 "case_join":"Captured campaign name + string case_id + exact label; plan_* fields from October 25 frozen matrix.",
 "raw_preservation":"All six archive file contents retained without byte changes; source tar.gz hash recorded above.",
 "derived_files":["classified_queue.csv","dcm_running_with_settings.csv","pending_jobs.csv","running_without_dcm_match.csv","summary.json","validation.json"],
 "rebuild_source":"build_inventory.py",
 "mutations":"Local dated evidence only; no SSH, Slurm submission, cancellation, update, model execution or GitHub mutation."
})
write_json(BASE/"validation.json", {"status":"PASS","checks":checks,"verified_raw_files":6,
 "queue_rows":335,"running_rows":20,"dcm_fits_joined":19,"source_archive_sha256":EXPECTED})
write(BASE/"build_inventory.py", pathlib.Path(__file__).read_bytes())
print(json.dumps({"written":CHANGED,"summary":summary},indent=2))
