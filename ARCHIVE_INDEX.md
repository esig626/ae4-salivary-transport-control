# Previous phase preservation

Repository: `esig626/ae4-salivary-transport-control`.

The previous work is archived by reference, not by compressing or uploading
its outputs again. No history was rewritten and no existing branch was deleted.

| Snapshot branch | Commit | What it preserves |
| --- | --- | --- |
| `archive/pre-phase2-main-20260913` | `b9d4e769747e636ad79d44a23c64c729f31b506d` | Published main through the routed WT work, Task 25 text results and Task 26 diagnostic record |
| `archive/pre-phase2-task27-20260913` | `b88ca91740459b10efa55905c57bfd1dd9cbb91f` | Last remotely visible Task 27 branch, which contains the inverse task prompt |

All other historical task branches remain where they were. In particular,
closed pooled model studies are not silently imported into the active model.
`main` and its historical `archive/` tree remain unchanged.

IMPORTANT: these snapshots preserve only published Git objects. The last
visible Task 27 remote contained no numerical search results. They do not
back up unpushed files or running processes in a separate Codex workspace.
Do not delete that workspace or claim that its search cache is archived.

The new branch inherits the main commit as its parent but constructs a reduced
working tree. Old results are available through the snapshot branches and Git
history, not scanned by the new task. This is a clean working tree, not a
history purge, repository replacement or an automatic merge into main.

The source of every imported file is recorded in
`reference/import_manifest.json`; imported scientific files are initially
reused byte for byte. New instructions and setup checks are labelled separately.
