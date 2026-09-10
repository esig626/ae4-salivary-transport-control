# AGENTS.md

These instructions apply to all work in this repository.

## Current scientific objective

The active prepared task is Task 16. The executable prompt is:

`prompts/16_wt_constrained_ae4_chloride_allocation.md`

Work on branch:

`codex/task-16-wt-chloride-allocation`

Tracking issue: #10.

## Mandatory prerequisite

Task 16 was prepared while Task 15 was still running.

Do not begin Task 16 scientific calculations until the completed final commit from

`codex/task-15-ae4-sign-slip-audit`

has been fetched, verified, and incorporated into the Task 16 branch.

Record the final Task 15 SHA and classification in the Task 16 manifest.

If Task 15 finds a clerical/algebraic AE4 sign error, a whole-cell AE4 source-mapping or expression-scaling error, or net reverse AE4 transport in the frozen physiological states, stop Task 16. The chloride-allocation hypothesis is then not the correct next step.

## Scientific state motivating Task 16

Task 14B found that reducing AE4 expression to 5% increased total 0 to 600 s secretion in all 30 root/calcium cases by 8.23% to 21.66%, despite the observed AE4-loss phenotype having the opposite direction.

Preliminary Task 15 work has found no simple sign/source-scaling error and has identified large opposing Na/K AE4 branch currents with small net productive Cl loading. Task 15 must finish before this preliminary result is treated as final.

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

Use the predeclared positive AE4+NKCC1 loading-share panel from the Task 16 prompt:

- inherited baseline;
- 0.025;
- 0.05;
- 0.10;
- 0.20;
- 0.30;
- 0.40.

These are sensitivity coordinates, not measured physiological fractions.

At each original WT root, derive AE4 and NKCC1 scales from the inherited positive Cl-loading sources so that the AE4+NKCC1 Cl-loading sum is preserved at the reference state before re-equilibration.

Do not hide or fold a negative AE2 counterflux into the positive-loading share denominator.

## WT-only admissibility

Continue each WT root through the allocation panel without genotype information.

Use only inherited WT physiological and numerical gates to decide admissibility. Preserve positivity, WT Cl/pH/Na/K constraints, physical domains, current closure, charge/carbon/water consistency and root-rank/residual requirements.

Then run WT dynamics for every WT-rest-admissible candidate at all three calcium values and apply only inherited WT/numerical gates.

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

For each frozen candidate:

- continue AE4 expression from 1.0 to 0.05 using the candidate parameters;
- do not reuse the old Task 14B 5% state after capacities change;
- do not attempt exact zero;
- run production Radau at all retained calcium values;
- evaluate matched AE2 loss as specified in the Task 16 prompt;
- report paired WT-relative total secretion and AE4-versus-AE2 contrast.

Only after these results are frozen may the experimental approximately 35% reduction be shown as context. Do not treat it as an exact target or acceptance band.

Detailed minute-by-minute phenotype shape remains diagnostic only.

## Interpretation

The central question is whether increasing AE4's WT productive chloride share, while preserving WT physiology and transporter thermodynamics, changes the direction of the 5% AE4 secretion effect.

Do not call one convenient root/share/calcium combination a successful repair.

Report robustness across the complete WT-admissible ensemble, both routing families and all retained calcium values.

If increasing AE4 capacity mainly magnifies the Task 15 Na/K slip and cannot produce a WT-admissible phenotype correction, state that clearly. That would motivate revisiting the mixed-cation architecture rather than further capacity scaling.

## General discipline

Do not edit `archive/`.

Do not merge to `main`.

Do not draft manuscript text.

Do not perform model reduction, GSPT, or a new unrelated identifiability study.

Commit and push completed Task 16 work only to `codex/task-16-wt-chloride-allocation`.
