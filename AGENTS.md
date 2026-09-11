# AGENTS.md

These instructions apply to all work in this repository.

## Current scientific objective

The active task is Task 19. Execute:

`prompts/19_ae4_loss_nak_pump_coupling.md`

Work only on branch:

`codex/task-19-ae4-loss-nak-pump-coupling`

Tracking issue: #17.

Repository inheritance is merged Task 18 on `main`, merge commit:

`52a56cc8e8788cdb5cae60dcbaa75257de2878dd`.

## Scientific question

Task 18 showed that changing WT chloride allocation does not repair the wrong-direction AE4-loss phenotype. Task 19 tests one narrow alternative mechanism:

> Does reducing effective Na/K-ATPase support together with AE4 near-loss make AE4 loss secretion-limiting?

This is a diagnostic sufficiency test, not evidence that AE4 experimentally regulates the pump.

## Primary model scope

Use the inherited baseline chloride-allocation model only.

Do not use Task 18's 0.10 or 0.30 rebalanced WT parameterizations in the primary Task 19 factorial panel. Those conditions require large capacity shifts and would confound this mechanism diagnostic.

Retain all ten inherited roots, both AE4 routing families, and calcium 0.10, 0.25 and 0.50 uM.

## Fixed factorial panel

Use exactly:

AE4 expression:

- 1.0
- 0.05

Na/K-pump capacity scale:

- 1.0
- 0.50
- 0.10

No intermediate pump scales, threshold searches or phenotype-guided refinement.

Reuse exact hash-valid scale-1.0 results wherever possible.

## Pump intervention

Scale the complete Na/K-ATPase capacity while preserving:

- 3Na:2K stoichiometry;
- pump law and direction;
- apical/basolateral pump partition;
- every non-pump parameter.

If total capacity plus fraction is used, scale total capacity and keep the fraction fixed. If separate effective capacities are used, multiply both by the same declared scale.

Do not change pump affinity constants or topology.

## Hard freeze outside AE4 expression and pump scale

Do not modify:

- AE4 transport law, routing or thermodynamics;
- chloride allocation;
- NKCC1;
- AE2;
- NHE1;
- CaCC;
- K channels;
- paracellular pathways;
- acid-base chemistry;
- water/lumen transport;
- bath;
- geometry;
- regulation;
- calcium inputs;
- solver or conservation tolerances.

Do not refit any capacity after changing AE4 expression or pump scale.

## Mandatory controls

Pump impairment alone can reduce secretion. Therefore for pump scales 0.50 and 0.10, evaluate both AE4 expression 1.0 and AE4 expression 0.05.

Do not call the mechanism successful merely because `Q(AE4=0.05,pump<1)` is below normal WT flow.

Report:

- matched AE4 effect at each pump scale;
- pump-only effect in AE4-intact cells;
- pump effect in AE4-low cells;
- AE4 × pump interaction;
- coupled-hypothesis comparison against normal WT.

Use the exact definitions in the Task 19 prompt.

## Numerical/root discipline

Use connected resting-state continuation from the inherited root.

For each pump scale below 1, first continue the AE4-intact state to the target pump capacity, then continue AE4 expression to 0.05 at fixed pump scale.

Internal continuation points are numerical only, not extra scientific pump conditions.

Do not rescue failed roots by changing any other parameter or jumping to a disconnected root.

## Interpretation discipline

The experimental approximately 35% AE4-loss secretion reduction is context only after the complete factorial results are frozen.

Do not choose a pump scale based on closeness to experiment.

Do not fit a continuous AE4-to-pump coupling law in Task 19.

A real mechanism signal requires distinguishing an AE4 × pump interaction from the trivial effect of pump inhibition itself.

## General discipline

Do not edit `archive/`.

Do not draft manuscript text.

Do not perform unrelated calibration, model reduction, GSPT or identifiability work.

Commit and push completed Task 19 work to `codex/task-19-ae4-loss-nak-pump-coupling`.

After Task 19 is complete and the required tests pass, open a pull request to `main` and merge the completed task using a merge commit. If a genuine scientific/code failure or merge conflict is exposed, report it rather than forcing the merge.
