# AGENTS.md

These instructions apply to all work in this repository.

## Current scientific objective

The active task is Task 16. The executable prompt is:

`prompts/16_wt_constrained_ae4_chloride_allocation.md`

Work on branch:

`codex/task-16-wt-chloride-allocation`

Tracking issue: #10.

## Task 15 status and explicit override

Task 15 completed its scientific audit locally but its large audit package did not successfully advance the remote Task 15 branch before the Codex session froze. Do not wait for that remote push and do not rerun Task 15.

The Task 15 scientific result is accepted for purposes of Task 16:

- no clerical or algebraic AE4 sign error was found;
- no whole-cell AE4 source-mapping error was found;
- the 0.05 AE4 expression factor is applied exactly once to the whole-cell AE4 source;
- AE4 net loads chloride at every audited frozen trajectory point;
- the Na and K AE4 branches oppose each other throughout the audited trajectories;
- WT branch cancellation is 96.32% to 99.22%, median 98.24%;
- 5% AE4 branch cancellation is 60.55% to 88.48%, median 75.87%;
- the large mixed-cation slip/cancellation is a property of the current model architecture and is not established native physiology;
- Task 15 did not establish that slip causes the wrong secretion phenotype;
- all 187 frozen scientific inputs were confirmed unchanged;
- the scientific model and `archive/` were unchanged.

These conclusions are diagnostic context only. They do not alter any scientific model parameter or equation.

**This section explicitly supersedes the `Mandatory Task 15 prerequisite` section in the Task 16 prompt.** Do not require a final Task 15 remote SHA, do not require Task 15 artifacts to be present on the Task 16 branch, and do not stop because the Task 15 branch itself did not finish pushing.

Base Task 16 scientific calculations on the unchanged completed Task 14B scientific state, commit:

`4d4403f279fc36aec940c9e4759486d02f07e3f6`

plus the accepted Task 15 diagnostic conclusions above.

Where Task 16 needs chloride partition, AE4 branch cancellation, or flux values numerically, evaluate those quantities directly from the unchanged baseline model or hash-valid saved Task 14/14B trajectories. Do not rerun the entire Task 15 audit.

## Scientific state motivating Task 16

Task 14B found that reducing AE4 expression to 5% increased total 0 to 600 s secretion in all 30 root/calcium cases by 8.23% to 21.66%, despite the observed AE4-loss phenotype having the opposite direction.

Task 15 then established that this is not explained by a simple sign or source-scaling mistake. It also found very large opposing Na/K AE4 branch currents with small net productive Cl loading.

Task 16 tests a conservative hypothesis before changing the AE4 microscopic architecture:

> Does the current model assign AE4 too little of the WT basolateral chloride-loading burden, and can WT-consistent capacity reallocation from NKCC1 toward AE4 restore a biologically sensible AE4 role?

## Phenotype firewall

The AE4-null phenotype is already known, so Task 16 is not blind. Nevertheless, candidate construction and acceptance must remain independent of that phenotype.

Do not use the approximately 35% experimental AE4-null secretion reduction, the Task 14B wrong-direction result, or any genotype secretion result to:

- choose AE4/NKCC1 capacity scales;
- choose a target chloride share;
- accept or reject a WT candidate;
- choose a root;
- choose a calcium value;
- choose a stopping rule;
- rank WT candidates.

First construct and freeze the complete WT-admissible candidate set. Only then evaluate 5% AE4 and AE2 loss.

## Allowed scientific changes

Task 16 permits exactly two primary allocation axes:

1. a multiplicative scale on complete AE4 carrier amount / active AE4 capacity;
2. a multiplicative scale on complete NKCC1 capacity.

These must be derived by the predeclared reference-state chloride-allocation rule in the Task 16 prompt.

Do not scale the AE4 chloride source alone.

An AE4 capacity change must scale the entire coupled AE4 source vector and preserve Na/K, Cl, HCO3, TIC, alkalinity, electroneutrality, local detailed balance and reversal.

An NKCC1 capacity change must scale the complete Na:K:2Cl cycle.

## Hard freeze outside the allocation axis

Keep fixed:

- all AE4 stoichiometry and source signs;
- the shared-carrier AE4 graph;
- AE4 Na/K routing within each retained routing family;
- AE4 rate-ratio/detailed-balance construction;
- beta/cAMP/PKA regulation and `R1_G125_P21_PROSE_TREFERENCE`;
- AE2 and NHE1 parameters;
- pump/K/CaCC topology and parameters;
- paracellular pathways;
- acid-base chemistry and buffers;
- water, lumen and outflow parameters;
- bath and geometry;
- calcium values `0.10`, `0.25`, `0.50 uM`;
- CCh + IPR protocol;
- solver and conservation tolerances.

Do not introduce genotype-specific compensation.

Do not alter the AE4 mixed-cation mechanism or suppress slip in Task 16. If simple capacity allocation fails, that becomes the reason for a later mechanism task.

## WT allocation panel

Use exactly three scientific allocation conditions:

- inherited baseline;
- 0.10;
- 0.30.

Do not evaluate 0.025, 0.05, 0.20, 0.40, midpoint refinements, threshold searches, or any additional phenotype-guided share.

These are sensitivity coordinates, not measured physiological fractions. Baseline anchors the inherited model, 0.10 tests a moderate AE4 loading share, and 0.30 tests a substantially larger share while NKCC1 remains dominant at the reference allocation.

At each original WT root, derive AE4 and NKCC1 scales from the inherited positive Cl-loading sources so that the AE4+NKCC1 Cl-loading sum is preserved at the reference state before re-equilibration.

Do not hide or fold a negative AE2 counterflux into the positive-loading share denominator.

Numerical continuation may use internal steps to reach the fixed 0.10 or 0.30 endpoint, but those internal points are not scientific candidates and must not trigger WT dynamics or genotype simulations.

Reuse hash-valid inherited baseline WT and genotype trajectories rather than recomputing them. Only the two modified allocation endpoints require new candidate calculations.

## WT-only admissibility

Continue each WT root only to the two modified allocation endpoints without genotype information.

Use only inherited WT physiological and numerical gates to decide admissibility. Preserve positivity, WT Cl/pH/Na/K constraints, physical domains, current closure, charge/carbon/water consistency and root-rank/residual requirements.

Then run WT dynamics for every WT-rest-admissible modified candidate at all three calcium values and apply only inherited WT/numerical gates.

The unresolved absolute one-SMG scale remains a nonblocking observation diagnostic.

No AE4-null time-course shape is a WT acceptance criterion.

## Freeze before genotype evaluation

Before any 5% AE4 or AE2 simulation under the reallocated candidates:

- save every WT candidate and decision;
- freeze every WT-admissible candidate;
- hash all parameter/root/trajectory payloads;
- commit and push the frozen candidate checkpoint.

No candidate may be dropped after seeing its genotype result.

## Post-freeze genotype evaluation

For each frozen modified candidate at target share 0.10 or 0.30:

- continue AE4 expression from 1.0 to 0.05 using the candidate parameters;
- do not reuse the old Task 14B 5% state after capacities change;
- do not attempt exact zero;
- run production Radau at all retained calcium values;
- evaluate matched AE2 loss as specified in the Task 16 prompt;
- report paired WT-relative total secretion and AE4-versus-AE2 contrast.

Reuse the inherited baseline genotype result when its hashes and scientific inputs match exactly.

Only after these results are frozen may the experimental approximately 35% reduction be shown as context. Do not treat it as an exact target or acceptance band.

Detailed minute-by-minute phenotype shape remains diagnostic only.

## Interpretation

The central question is whether increasing AE4's WT productive chloride share, while preserving WT physiology and transporter thermodynamics, changes the direction of the 5% AE4 secretion effect.

Do not call one convenient root/share/calcium combination a successful repair.

Report robustness across the complete WT-admissible ensemble, both routing families and all retained calcium values.

If a sign change occurs between baseline, 0.10 and 0.30, report only the sparse result. Do not refine the share axis in Task 16.

If increasing AE4 capacity mainly magnifies the Task 15 Na/K slip and cannot produce a WT-admissible phenotype correction, state that clearly. That would motivate revisiting the mixed-cation architecture rather than further capacity scaling.

## General discipline

Do not edit `archive/`.

Do not merge to `main`.

Do not draft manuscript text.

Do not perform model reduction, GSPT, or a new unrelated identifiability study.

Commit and push completed Task 16 work only to `codex/task-16-wt-chloride-allocation`.
