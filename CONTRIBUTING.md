# Maintaining the project handover

Keep scientific history and operational recovery usable for a new researcher.

1. Read docs/START_HERE.md and the latest dated review before editing.
2. Preserve raw fitting observations, whole-cast membership and frozen experiment identity. Record proposed scientific changes separately.
3. Use campaign-qualified case IDs. Keep all completed, failed, prior-only, reused and incomplete entries.
4. Add new dated results rather than overwriting older snapshots. Update README's evidence date and docs/NEXT_STEPS.md.
5. Explain the question a new case tests, its direct reference case, and every changed equation/prior/data/sampler/resource setting.
6. Keep model fitting, checking, reporting and publication separate. A reporting repair should not start replacement fitting jobs.
7. Verify local links and table/count integrity with `python scripts/verify_repository.py`. Full R/model checks require the matching environment and actual data.
8. Preserve original package hashes. If source changes, version the changed package explicitly and distinguish it from archived executed code.
9. Do not commit raw research data, fitted objects, chain parameter CSVs, credentials or software installations. Compact results/metadata and documented diagnostic excerpts are intentional exceptions to the raw-output boundary.
10. Do not silently assign a licence or claim final model validity. Discuss data rights and publication choices with the owner.

Historical scripts in archive/ are preserved for interpretation. Their absolute paths, older assumptions and previously reported validation limits should be made explicit if adapted.

The repository is a continuation record, not an instruction to automatically execute every queue launcher.
