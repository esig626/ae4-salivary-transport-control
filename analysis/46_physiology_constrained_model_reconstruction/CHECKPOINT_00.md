# Checkpoint 00 — branch bootstrap and recovery record

- Milestone: 00, branch bootstrap and recovery record.
- Branch: `analysis/task-46-physiology-constrained-model-reconstruction`.
- Branch head before milestone: `850706dd6b64d68793d7bca89adc65679c73d4b0`.

## Completed

Read AGENTS.md first, then all three mandatory Task 46 specifications and the
user's stricter checkpoint directive. Verified the remote branch list had no
Task 46 branch/checkpoint to resume. Created an isolated Task 46 branch from the
current main commit. Task 43 and Task 44 are fixed read-only inputs at
`547f113d9123ab1976774639a69483faf40cef34` and
`e5fa9840bf96be3147c6118daee4942468c78f8b` respectively. Their remote heads
match those pins and both commit objects exist locally. Inspected their source
trees and completion records in the requested order. No numerical work or
model changes have been performed.

## Verified

The clone started clean and the Task 46 branch is at its stated base.
All three required prompt blobs and AGENTS.md are pinned. Remote reference
identities and base tree are recorded. Independent read-only protocol review
confirmed the user directive supersedes the manuscript active-task assignment,
requires separate 00–06 checkpoints, and permits this branch base. The local
publication guard checks branch identity, source commit availability, prompt
hashes and that every change is confined to the Task 46 directory.

## Products

- `PROTOCOL.md`: scope, scientific rules, evidence discipline, milestone order,
  publication and recovery procedure.
- `source_pins.json`: base, companion sources, specification blob identities and
  initial remote references.
- `verify_checkpoint.py`: read-only preservation/publication guard.
- `CURRENT_STATUS.md` and this immutable `CHECKPOINT_00.md`.

## Scientific conclusions so far

None newly established. Task 44's failure-baseline findings are inputs to be
reproduced, not independently verified Task 46 results. No model is accepted.

## Unresolved issues

Task 44 reproduction, compensation diagnosis, evidence split and every
scientific reconstruction/validation stage remain pending. Evidence screening
versus unseen final holdouts will be declared in milestone 03. Runtime currently
has Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0 and Matplotlib 3.10.8; SymPy is not
available in the active runtime. Resolve dependencies only if required by the
pinned replay. LaTeX and PDF-rendering executables are available.

## Exact next milestone

Milestone 01: reproduce the relevant Task 44 baseline without scientific
changes, including WT, 5% and null; stationary/600 s deficits; local and
finite-time NKCC compensation; pH/Na/Cl/volume; eigenvalues and slow modes.
Do not begin causal diagnosis or model search until checkpoint 01 is published.

## Resume commands

From the repository root, first fetch and reconcile only by fast-forward:

```bash
git fetch origin
git switch analysis/task-46-physiology-constrained-model-reconstruction
git merge --ff-only origin/analysis/task-46-physiology-constrained-model-reconstruction
python analysis/46_physiology_constrained_model_reconstruction/verify_checkpoint.py --published
git show HEAD:analysis/46_physiology_constrained_model_reconstruction/CURRENT_STATUS.md
```

Read pinned Task 43/44 source reports and run instructions without changing them:

```bash
git show 547f113d9123ab1976774639a69483faf40cef34:analysis/43_parameter_provenance/report.tex
git show e5fa9840bf96be3147c6118daee4942468c78f8b:analysis/44_full_system_mathematics/README.md
git show e5fa9840bf96be3147c6118daee4942468c78f8b:analysis/44_full_system_mathematics/run_analysis.py
```

Then implement a Task 46 replay wrapper that imports unchanged pinned source,
redirects only output destinations, and writes `output/checkpoint_01_*`.
No replay script or numerical result is claimed to exist at checkpoint 00.

## Publication cleanliness and identity

This checkpoint is published only if the commit leaves a clean working tree
and the remote branch equals that commit, independently checked with
`verify_checkpoint.py --published` immediately after push. Pre-commit changes
are precisely this milestone's Task 46 payload; protected files are unchanged.
The containing Git commit is the checkpoint SHA (avoids a self-referential SHA).
Publication failure prevents progression; do not treat local existence as
remote publication. No force push or merge is authorised or used.
