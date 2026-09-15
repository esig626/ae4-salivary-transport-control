# Task 46G checkpoint: local CBM core verified

* Branch: `analysis/task-46-physiology-constrained-model-reconstruction`.
* Parent: `3598bd866c33d41d4a6914bda56d17c73900e82a`.
* User authority: Task 46G at main commit
  `1148397ca7add328b1e3da71864c56e5382ae324`.
* Local core provenance: CarbonScope source commit
  `a11e25f176a2cc70de22c3f36a1c0e623248393a`.

Read the local README, provenance and third party notices, installed its
requirements, and passed its smoke test. An independent analytical control
verified two FBA objectives, three FVA regions, 30 feasible endpoint witnesses,
90 comparisons against analytical endpoints and 60 comparisons between
reference and fast paths. There were zero failures and zero absolute endpoint
error. Every core import resolved inside the local vendor directory, whose
source files remained unchanged.

## Products

* `verify_local_cbm_core.py`.
* `LOCAL_CBM_VERIFICATION.md` with the derivation and exact ranges.
* `output/checkpoint_46g_local_core_verification.json` with full evidence.
* Updated execution protocol and current status.

No solver extension was required. The work used no CarbonScope repository
access. New tracked work is confined to Task 46. Main, manuscript, Tasks 1–45,
production source, the existing core and previous outputs are unchanged.
The branch has not been merged.

The first preservation guard could not find pinned Task 43/44 commit objects
in the new AE4 clone. The exact AE4 objects were fetched and the unchanged
guard passed. No prior numerical model or scientific replay was executed.

## Interpretation and next step

Task 46G's local solver prerequisite passes. No AE4 scientific CBM calculation,
kinetic candidate or phenotype fitting was performed. The full Task 46
reconstruction remains open. Next comes the evidence split and Task 46E
network construction checkpoint. The later addenda supersede the historical
instruction to extend baseline forensics before starting the CBM screen.
Keep all AE4 network, paired WT/KO, coupling and sparse repair code above the
local generic solver layer. Publish each CBM milestone before continuing.

## Publication and resumption

The containing commit identifies this immutable checkpoint. Publication is
complete only when the remote Task 46 branch matches that commit and the local
checkout is clean. Never force push or write to CarbonScope.

```bash
git fetch origin
git pull --ff-only
python analysis/46_physiology_constrained_model_reconstruction/verify_checkpoint.py --published
PYTHONDONTWRITEBYTECODE=1 python analysis/46_physiology_constrained_model_reconstruction/verify_local_cbm_core.py --output /tmp/ae4_task46g_verification.json
```

Use a new output filename if that scratch verification record already exists.
