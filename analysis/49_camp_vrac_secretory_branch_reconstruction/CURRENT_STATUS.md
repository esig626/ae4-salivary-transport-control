# Task 49 current status

## Canonical branch

`analysis/task-49-camp-vrac-secretory-branch-reconstruction`

## Scientific parent

Merged Tasks 47-48:

`5ca70c61e8ccc49f084bac1b91ea05c612d197e8`

## Current state

Task 49 has been staged but no Task 49 scientific inference, parameter calibration or production trajectory has been run.

The task is restricted to one previously untested, independently source-backed architectural reconstruction:

- retain the existing CCh/Ca -> TMEM16A branch;
- add the documented beta/cAMP-associated, volume-regulated apical chloride conductance;
- add the independently documented beta/cAMP input to NKCC1 regulation without replacing the Palk/Benjamin core law;
- do not alter AE4 stoichiometry, NHE1, NBC, pump, K channels, Ca amplitude or Task 48 genotype resting-state machinery.

The binding no-repeat ledger is `EXCLUSION_LEDGER.md`.

## Publication discipline

Milestones 49A through 49F must each be committed and pushed immediately. Every milestone must update this file and write a machine-readable checkpoint receipt under `output/checkpoints/`. The orchestrator must verify the remote SHA before dependent work begins.

Only the orchestrator may modify the canonical branch, integrate worker outputs, run authoritative inference/production trajectories, or accept scientific changes.

## Next step

Launch Codex from the current remote Task 49 head. Read `CODEX_START_HERE.md` and execute `prompts/49_camp_vrac_secretory_branch_reconstruction.md` using up to five parallel workers for milestone 49A.
