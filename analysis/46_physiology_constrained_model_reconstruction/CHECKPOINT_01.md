# Checkpoint 01 — Task 44 failure baseline replay

- Milestone: 01, reproduce the unchanged relevant Task 44 baseline.
- Branch: `analysis/task-46-physiology-constrained-model-reconstruction`.
- Branch head before milestone: `6b8fc26df59da943e2656cf3237ee603311c00a7`.
- Previous remotely verified checkpoint: 00 at that same SHA.

## Completed

Reproduced WT, 5% AE4 and AE4-null 600 s production trajectories, the C0
constant-input equilibrium continuation, eigenvalues/slow modes and local
compensation derivative. Recomputed cumulative secretion and the two distinct
60–600 s compensation measures. Retained pH, Na, K, Cl, carbon, alkalinity,
cell/lumen volumes/compositions, electrical closure and AE4 regulation.
No scientific parameter or equation changed.

Read the Task 43 and Task 44 pinned reports, PDFs and numerical source context.
All original sources remain in clean detached checkouts. The replay calls
unchanged source functions and writes exclusively to Task 46 paths.

## Verified

The primary comparison record has 1,794 numerical comparisons, all bit-identical
to the corresponding pinned references. The three production cases each pass
934 full inherited state/conservation checks; all eight C0 roots pass full
physiology gates and have negative-real-part spectra. Four independently
solved nearby roots verify the local WT IFT derivative. A separate independent
verifier checks saved dense states, volume extrema, reconstructed quadrature,
root gates and source preservation; its exact counts and result are recorded
in `output/checkpoint_01_independent_verification.json`: 67,067 numerical
comparisons (67,060 exact, seven quadrature-roundoff differences), 3,886 logical
checks and 1,806 fresh saved-state evaluations, with zero failures.

One preliminary WT attempt failed only during JSON serialisation of a NumPy
boolean. Its CSV outputs and failure record are preserved. The fixed run used
four integrations; including the preliminary WT attempt, this milestone used
five integrations, twelve stationary solves, no optimisation, no parameter
search and no scientific retry/rescue.

## Files and results

- `replay_baseline.py`, `verify_baseline.py`.
- `BASELINE_REPRODUCTION.md`.
- `output/checkpoint_01_replay/`: trajectories, states, summaries, roots, IFT,
  compensation and primary comparison record.
- `output/checkpoint_01_independent_verification.json`.
- `output/checkpoint_01_execution.json`.
- `output/checkpoint_01_attempt01_failure.json` and
  `output/checkpoint_01_baseline/`: preserved implementation failure outputs.
- Updated `PROTOCOL.md` describes the verified GitHub-integration publication
  route after shell Git lacked write authentication.

## Scientific conclusions so far

The failure baseline is reproduced: cumulative deficits are 3.417577608% (5%)
and 3.857486297% (null); the stationary null deficit is 0.9461727943%.
The local independent-verifier NKCC chloride compensation is 92.90107418%;
finite-time replacement is 88.50703510%; relative integrated NKCC increase is
23.16344951%. WT/5%/null slowest decay times are
345.81857188/442.08904428/459.91603104 s. These are separate local and finite-time
numerical results. No new model has been proposed or accepted.

## Unresolved issues

Causal diagnosis remains pending. The inherited thermodynamic opposition of
C0's substituted source convention, the shared WT resting initial condition
for mutants, isolated-assay protocol mismatch and lack of global/finite-range
robustness proofs remain explicit. No evidence split or model search has yet
been performed. The baseline is not a successful physiological reconstruction.

## Exact next milestone

Milestone 02: causal diagnosis of excessive NKCC compensation using the pinned
Task 43 provenance and Task 44 derivatives, state/network interactions and
chloride storage. Rank causes in Markdown and machine-readable form. Do not
change equations or begin evidence/model search before checkpoint 02 publication.

## Exact resume commands

From the Task 46 repository root:

```bash
git fetch origin
git switch analysis/task-46-physiology-constrained-model-reconstruction
git merge --ff-only origin/analysis/task-46-physiology-constrained-model-reconstruction
python analysis/46_physiology_constrained_model_reconstruction/verify_checkpoint.py --published
git show HEAD:analysis/46_physiology_constrained_model_reconstruction/CURRENT_STATUS.md
git show 547f113d9123ab1976774639a69483faf40cef34:analysis/43_parameter_provenance/output/task44_sensitivity_crosswalk.csv
git show e5fa9840bf96be3147c6118daee4942468c78f8b:analysis/44_full_system_mathematics/output/expression_sensitivities.csv
git show e5fa9840bf96be3147c6118daee4942468c78f8b:src/modern_full_model/nkcc1_palk2010.py
```

Read `BASELINE_REPRODUCTION.md`, then implement the new Task 46 diagnosis
script/report; no diagnosis script is claimed to exist at this checkpoint.
Do not rerun completed baseline calculations merely to resume milestone 02.
For an explicitly requested independent replay, use a new scratch destination:

```bash
PYTHONDONTWRITEBYTECODE=1 python analysis/46_physiology_constrained_model_reconstruction/replay_baseline.py --source /absolute/path/to/clean-pinned-source44 --output /absolute/path/to/new-empty-replay-output
```

The source checkout must be exactly
`e5fa9840bf96be3147c6118daee4942468c78f8b`; the output directory must not exist.

## Publication cleanliness and identity

Checkpoint 00 was independently confirmed clean and remotely published at
`6b8fc26df59da943e2656cf3237ee603311c00a7` before this work began. This checkpoint
is published only if its independently reviewed tree is committed, the remote
Task 46 ref contains that commit, and the local checkout is fast-forwarded and
clean. The post-publication guard is mandatory before milestone 02.
The containing commit identifies checkpoint 01. Earlier checkpoint files and
outputs are unchanged. No force push, history rewrite or branch merge is used.
