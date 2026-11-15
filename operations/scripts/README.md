# Ibex helper scripts used from 27 Oct to 15 Nov 2026

These scripts were written in the project thread and run by Gage on the Ibex login node. The usual pattern was: download the script, `scp` it to `/ibex/user/meleg/DCM_MODEL/`, run it there with output teed to `logs/`, and paste the output back. The full dated log is [RUN_LOG_20261027_20261115.md](../RUN_LOG_20261027_20261115.md).

| Script | Effect | Used |
|---|---|---|
| [01_status_snapshot.sh](01_status_snapshot.sh) | Read-only. Prints shape150 controller/queue state, per-case markers and new results | 27 Oct |
| [02_controller_diagnosis.sh](02_controller_diagnosis.sh) | Read-only. Shows why the shape150 controller died (sacct and log tails) and who owns the pending jobs | 27 Oct |
| [03_restart_checks_only.sh](03_restart_checks_only.sh) | Defers the 141 never-submitted shape150 fits with reversible `scheduler/fit_NNN.rejected.json` markers, then submits a new controller (`controller_002`, job 54668377) so that 28 finished fits get checked | 27 Oct |
| [08_status_all.sh](08_status_all.sh) | Read-only. Prints the round3 and bins1m job states (sacct), a per-case score table, failure text, controller log tails and shape150 marker counts | 12–15 Nov |

The deferral markers were removed later on 27 Oct (16:18 UTC), which returned all 141 fits to the shape150 queue.

## Typical commands

```powershell
# Windows PowerShell: copy a script up
scp $HOME\Downloads\08_status_all.sh meleg@ilogin.ibex.kaust.edu.sa:/ibex/user/meleg/DCM_MODEL/
```

```bash
# Ibex login node
cd /ibex/user/meleg/DCM_MODEL && bash 08_status_all.sh 2>&1 | tee logs/status_all.txt
squeue -u $USER -o "%.18i %.45j %.9T %.10M" | grep dcm-
```

The campaign installers are kept with their packages:
- [DCM_Ibex_next150](../../DCM_Ibex_next150/): scripts 04 and 05;
- [DCM_Ibex_round3](../../DCM_Ibex_round3/): script 06;
- [DCM_Ibex_bins1m](../../DCM_Ibex_bins1m/): script 07;
- [DCM_Ibex_plots_top5](../../DCM_Ibex_plots_top5/): script 09.
