# Current status: Task 46G local solver verification complete

* Branch: `analysis/task-46-physiology-constrained-model-reconstruction`.
* Preceding branch SHA: `3598bd866c33d41d4a6914bda56d17c73900e82a`.
* Current immutable checkpoint: `CHECKPOINT_46G_LOCAL_CORE.md`.
* The containing commit identifies this checkpoint after publication.

## Completed and verified

Read Task 46G and the local core documentation. Installed its declared
dependencies and passed the supplied smoke test. The independent analytical
control passed two FBA objectives and reference/fast/prepared FVA checks across
three regions. All 90 comparisons with analytical endpoints and all 60
comparisons between reference and fast paths passed, with zero absolute
endpoint error. Thirty explicitly specified complete states independently
verify endpoint attainability.

All solver imports resolved within the required local vendor directory. Core
file hashes remained unchanged. No CarbonScope repository was accessed or
modified. No AE4 scientific CBM optimisation, new kinetics or phenotype fitting
was performed. The full Task 46 reconstruction is not complete.

Read `LOCAL_CBM_VERIFICATION.md` and
`output/checkpoint_46g_local_core_verification.json` for the exact control,
versions, tolerances, imports, hashes and results.

## Preserved scientific state

Checkpoint 01 remains the authoritative Task 44 failure replay, published at
`83e69048a329135b95ef00f007610610ff8afc6a`. Its scripts and outputs are
unchanged. Its report and records retain the known thermodynamic, protocol
and mathematical limitations. No AE4 conclusion follows from this control.

Main, manuscript, Tasks 1–45, production source, previous outputs and the
local core are unchanged. The Task 46 branch remains unmerged.

## Next scientific milestone

Freeze the construction/holdout evidence split and implement Task 46E CBM
milestone 1: one curated network, transport ledger, exact S matrix, charge and
conservation audit, provenance of bounds, rank, nullspace dimension and shared
WT/KO capacities. Publish it before phenotype optimisation. Do not rerun the
old dynamic model merely to begin this stage. Keep all AE4 biology above the
local generic solver layer, without an external CarbonScope dependency.

Subsequent CBM milestones cover WT feasibility/FVA, knockout feasibility,
the justified NKCC constraint, targeted coupling, sparse structural repair and
the decision report. Each requires its own remote checkpoint. New kinetics
remain gated on the published CBM decision.

## Resume and verify publication

```bash
git fetch origin
git pull --ff-only
python analysis/46_physiology_constrained_model_reconstruction/verify_checkpoint.py --published
```

Publication is complete only when the remote branch matches this checkpoint's
commit and the working tree is clean. `PROTOCOL.md` records the connected
GitHub publication route if shell write authentication is unavailable.
