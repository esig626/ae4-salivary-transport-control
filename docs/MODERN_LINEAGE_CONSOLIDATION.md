# Modern lineage consolidation

This integration branch was created from `main` at `66f25d7eef0e66aa896276d56c4394d9099ba108` and selectively imports the scientifically active clean-workspace lineage through Task 38 from commit `d4bce14c5caed4e05587e237f002814f0f27e226`.

The purpose is to preserve the complete historical repository on `main` while restoring the modern production source, focused tests, compact reference checkpoints, and Task 30-38 reports/results without importing the clean workspace's mass deletions of archive, manuscript, figures, earlier analyses, or historical task records.

Imported material:

- `src/modern_full_model/` files from the Task 38 clean lineage, overlaid file-by-file so unrelated source files already on `main` are preserved.
- focused tests for Tasks 28, 30, 31 and 36.
- compact `reference/` checkpoints and provenance files required by the repaired model.
- Task 28-31 and 37-38 prompts used to reproduce the active work.
- analysis/report directories for Tasks 30, 31, 36, 37 and 38.
- result directories for Tasks 30, 31, 37 and 38.
- the small execution helpers used by the repaired lineage.

Not imported:

- any deletion from the clean workspace;
- changes to historical `archive/`, manuscript, figures, old analysis/results, or legacy model files;
- Task 39 scientific changes, because Task 39 had not yet executed at consolidation time.

Task 39 should be applied later as a normal incremental change on top of the consolidated production lineage after its source-fixed validation is complete.
