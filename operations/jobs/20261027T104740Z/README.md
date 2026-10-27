# Queue capture reported on 27 October 2026

The user ran scripts/capture_ibex_jobs.py on Ibex and provided its successful terminal output. The query ended at 10:47:40 UTC (13:47:40 Saudi time), reporting **335 total queue entries and 20 RUNNING jobs** for user meleg, with array tasks expanded.

[reported_capture_metadata.json](reported_capture_metadata.json) transcribes the JSON shown in that terminal output. It has not been independently read from the generated archive.

**Per-job archive not yet received.** The reported remote path is:

/ibex/user/meleg/DCM_MODEL/logs/ibex_jobs_20261027T104740Z_1793098060392524450.tar.gz

Expected contents are all_queued_jobs.csv, running_jobs.csv, dcm_manifest_matches.csv, squeue_raw.txt, metadata.json, and README.md. Do not fabricate these files from aggregate counts or infer that every running job is a DCM fit. The 315 non-RUNNING entries may include multiple states.

Keep the [October 25 inventory](../20261025/shape_all_282_job_registrations.csv) as a historical record. Once the archive arrives, preserve its bytes/checksum, compare metadata to this reported summary, and publish the full account inventory with DCM matches and explicit timestamps. No fits or Slurm settings were changed by this documentation update.
