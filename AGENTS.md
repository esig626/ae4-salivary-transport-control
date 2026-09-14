# Repository state after Task 40

`main` is the canonical scientific lineage.

The production model is consolidated through Task 40. There is currently no active numbered task branch and no instruction to continue work on a completed task branch.

Before starting new scientific work:

1. Read `README.md` for the current scientific interpretation.
2. Read `docs/MODERN_LINEAGE_CONSOLIDATION.md` for the modern lineage repair and consolidation history.
3. Read the most recent task report under `analysis/40_ae4_equal_cation_routing/` and its machine readable outputs under `results/40_ae4_equal_cation_routing/`.
4. Create a new branch from current `main` for any new numbered scientific task.

## Canonical current result

Task 40 tests fixed 50:50 Na/K routing of the inherited AE4 cation component while preserving the total AE4 cycle and every other inherited scientific mechanism.

One WT resting solve and all three 600 s trajectories passed the declared numerical, physiological and conservation gates without retries or tuning.

At 600 s, AE4 loss gives final cumulative secretion deficits of about 3.42% at AE4 = 0.05 and 3.86% at AE4 = 0.00. NKCC1 compensation is about +20.82% and +23.16%, respectively. The current model therefore remains far from the held out experimental secretion reduction of about 35%.

This is a negative model result, not permission to tune to the held out phenotype.

## Scientific constraints

* Treat the approximately 35% experimental AE4 loss secretion reduction as held out phenotype information, not a calibration target.
* Do not reinterpret a failed or weak genotype phenotype as permission for post hoc parameter fitting.
* Preserve exact provenance for any source fixed transporter law, parameter or experimental quantity.
* Keep genotype specific compensation separate from constitutive WT model calibration unless a future task explicitly introduces and sources such regulation.
* Preserve conservation, charge accounting, resting closure and the declared production protocol unless a future task explicitly changes one of them.
* Distinguish numerical failure, structural failure and biological mismatch in all reports.

## Repository constraints

* `archive/` is historical and provenance material. Do not rewrite or delete it.
* Do not delete earlier numbered analyses or results merely because later tasks supersede their interpretation.
* New scientific changes should be made on a new branch from current `main`, with a task specific prompt, tests, compact results and a final report.
* Keep temporary computation products out of the repository unless they are required for reproducibility or audit.
* Use the connected GitHub integration for remote writes. Do not add PATs, SSH credentials or local authentication material to the repository.
* Do not merge a new scientific branch until its declared checks and result audit are complete.

## Current main lineage

The latest scientific merge is Task 40, PR #32. Repository housekeeping after that merge updates only top level project guidance and does not alter the Task 40 scientific outputs.
