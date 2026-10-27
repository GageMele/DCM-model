# DCM MLD-onset campaign

Use the single-file installer queue_onset_dcm.sh. It installs into the new
DCM_Ibex_onset_tests directory and queues 28 new fits plus six existing-fit
evaluations. Prior campaigns and pointer files remain untouched. Read PLAN.md
for the scientific design and VALIDATION.md for checks and limits.

Run from the Ibex project directory:

    bash queue_onset_dcm.sh

To install for inspection without queueing:

    bash queue_onset_dcm.sh --extract-only

Repeating submission resumes the frozen campaign. No automatic refits or job
cancellations occur. Check results/onset_all_scores.csv,
results/onset_planned_paired_contrasts.csv and results/mld_neighbourhood_scores.csv
alongside all-run predictive/diagnostic tables. All completed or failed cases remain
visible. Lower conditional onset error alone is insufficient for model selection.
