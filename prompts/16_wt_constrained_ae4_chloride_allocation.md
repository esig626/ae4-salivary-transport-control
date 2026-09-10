# Task 16: WT-constrained AE4 chloride-allocation analysis

## Mission

Task 14B showed that reducing AE4 expression to 5% increased total secretion in every retained root/calcium case. Task 15 was then designed to determine whether this wrong-direction phenotype came from a sign/source-mapping defect or from the current mixed-cation AE4 architecture.

Task 16 asks a narrower calibration-structure question:

> Does the modern model assign AE4 too small a productive share of WT basolateral chloride loading, and can that chloride-loading burden be redistributed toward AE4 while preserving WT physiology and all transporter stoichiometry/thermodynamics?

This is **not** a fit to the known AE4-loss phenotype.

The approximately 35% experimental AE4-null secretion reduction is known. It must not enter candidate generation, parameter selection, admissibility, root selection, calcium selection, stopping rules, or ranking before the WT-admissible candidate set is frozen.

Work on branch:

`codex/task-16-wt-chloride-allocation`

Tracking issue: #10.

Do not edit `archive/`. Do not merge to `main`.

## Mandatory Task 15 prerequisite

This prompt was prepared before Task 15 finished.

Before any Task 16 scientific calculation:

1. fetch the completed remote branch `codex/task-15-ae4-sign-slip-audit`;
2. identify its final pushed commit;
3. verify that the required Task 15 final report and result artifacts exist;
4. incorporate that final Task 15 commit into the current Task 16 branch before proceeding;
5. record the Task 15 final SHA and classification in the Task 16 manifest.

If Task 15 concludes either

- `CLERICAL OR ALGEBRAIC AE4 SIGN ERROR FOUND`, or
- `WHOLE-CELL AE4 SOURCE MAPPING OR EXPRESSION-SCALING ERROR FOUND`,

stop Task 16. A chloride-allocation study is not the correct next step until the exact defect is repaired in a separate controlled task.

If Task 15 finds the AE4 net transport direction is reversed in the frozen physiological states, stop and report that Task 16's premise is invalid.

Proceed only if Task 15 establishes that the implemented signs/source mapping are correct and the present model is scientifically meaningful enough to test a WT chloride-allocation hypothesis.

## Scientific principle

Do **not** increase the intracellular chloride source of AE4 by itself.

AE4 is an electroneutral coupled transporter. Any increase in AE4 capacity must scale the complete AE4 cycle and therefore preserve its linked Na/K, Cl, HCO3, TIC and alkalinity sources, its local detailed balance, and its reversal surfaces.

Likewise, any NKCC1 adjustment must scale the complete Na:K:2Cl cycle.

The primary Task 16 intervention is therefore a **capacity reallocation** between AE4 and NKCC1, not a hand-edited chloride source.

The existing AE4 Na/K routing, stoichiometry and microscopic graph remain fixed in this task. Task 16 is intended to answer whether simple WT chloride-load allocation is sufficient before changing the AE4 architecture itself.

## Hard freeze apart from the declared allocation axis

Keep fixed exactly as inherited from the final Task 15 scientific state:

- all ten retained root families/seeds;
- AE4 reaction stoichiometry and source matrix;
- AE4 common-Cl/shared-carrier graph;
- AE4 Na/K mixed-bath routing for each retained routing family;
- AE4 forward/reverse rate ratios and detailed-balance construction;
- AE4 beta/cAMP/PKA regulation and `R1_G125_P21_PROSE_TREFERENCE`;
- AE2 law and capacity;
- NHE1 law and capacity;
- CaCC and calcium gate;
- Na/K pump topology and parameters;
- K-channel topology and parameters;
- paracellular transport;
- acid-base chemistry and finite buffers;
- water/hydraulic model and lumen outflow;
- bath composition and geometry;
- calcium values `0.10`, `0.25`, `0.50 uM`;
- CCh + IPR protocol and 600 s duration;
- solver and conservation tolerances.

The **only scientific parameter changes licensed in the primary allocation panel** are:

1. a multiplicative scale on AE4 carrier amount / total active AE4 capacity, applied to the complete AE4 transporter;
2. a multiplicative scale on NKCC1 capacity, applied to the complete NKCC1 cycle.

Do not change AE4 Na/K routing, common-Cl attempt-rate ratios, stoichiometry, cooperative-gate structure, regulation, CaCC, water transport, or any other parameter in this task.

Do not introduce genotype-specific compensation.

## Stage A: quantify the inherited WT chloride-loading partition

Using the final Task 15 state, evaluate every original WT root at rest and report the component chloride sources.

At minimum record:

- AE4 net intracellular Cl source `A0`;
- NKCC1 intracellular Cl source `N0 = 2 * J_NKCC1`;
- AE2 intracellular Cl source `E0` separately, with its sign;
- total basolateral Cl source;
- AE4 productive fraction/cancellation fraction from Task 15;
- intracellular Na, K, Cl, HCO3, pH and volume.

Do not combine a negative AE2 counterflux into a positive-loading share denominator.

Define the primary positive-loading pool as

`L0 = A0 + N0`

only when both `A0 > 0` and `N0 > 0`.

Define inherited AE4 positive-loading share

`f_AE4_0 = A0 / (A0 + N0)`.

If either AE4 or NKCC1 is not positive-loading in a root, flag the root and do not force the allocation formula onto it.

## Stage B: predeclare the allocation panel

Before solving any modified WT equilibrium, write and commit a machine-readable allocation contract.

Use the same target-share panel for every eligible root:

- inherited baseline share;
- `0.025`;
- `0.05`;
- `0.10`;
- `0.20`;
- `0.30`;
- `0.40`.

These values are **sensitivity coordinates**, not measured physiological fractions and not phenotype-fit targets.

Do not add or remove a share value after looking at any 5% AE4 phenotype result.

For each eligible root and each non-baseline target share `f`, derive scales at the original WT reference state by keeping the AE4+NKCC1 positive Cl-loading sum fixed:

`A_target = f * L0`

`N_target = (1 - f) * L0`

`lambda_AE4 = A_target / A0`

`mu_NKCC1 = N_target / N0`.

Implement `lambda_AE4` by scaling the complete AE4 carrier amount/capacity.

Implement `mu_NKCC1` by scaling the complete NKCC1 capacity.

This construction is a reference-state partition transformation. It does not assert that the realized share after re-equilibration remains exactly `f`.

Require both capacity scales to be finite and nonnegative. Record their magnitudes transparently. Do not impose an arbitrary hidden upper bound simply because a large scale looks uncomfortable; let WT physiological constraints determine admissibility and discuss implausibly large scales separately.

## Stage C: WT-only continuation and admissibility

For each original root, continue the WT resting equilibrium from the inherited allocation toward increasing target AE4 shares in the predeclared order.

No AE4-loss phenotype may be used in continuation, root selection or acceptance.

Use the inherited root-distance/continuation discipline where applicable. Preserve alternate roots for audit if encountered. Do not jump to a disconnected root merely because it later gives a preferred genotype result.

For every candidate WT equilibrium require the inherited WT resting gates, including at minimum:

- positivity;
- intracellular Cl gate;
- intracellular pH gate;
- intracellular Na and K gates;
- physical volume/domain gates;
- current closure;
- charge consistency;
- carbon/alkalinity accounting;
- full resting numerical rank / root residual requirements.

Use the existing source-supported WT constraints already encoded in the repository. Do not invent new numerical WT acceptance bands for Task 16.

For every WT-rest-admissible candidate, report the **realized** resting chloride partition after re-equilibration:

`f_AE4_realized = A / (A + N)`

when `A > 0` and `N > 0`, together with AE2's signed contribution.

Also report actual AE4 Na/K branch cancellation, gross branch turnover and net productive Cl loading, because increasing total AE4 capacity while leaving the slip geometry unchanged may create large cation traffic.

## Stage D: WT dynamic validation only

Before any reduced-AE4 simulation, run the WT dynamic protocol for every WT-rest-admissible candidate at all three fixed calcium values `0.10`, `0.25`, `0.50 uM`.

Use production Radau and inherited numerical/conservation standards. Use BDF confirmation according to a predeclared representative rule if the existing framework requires it.

Candidate WT dynamic admissibility must be based only on WT evidence and inherited numerical/physiological gates.

The inherited absolute whole-gland 9 to 10 uL/min mapping remains a **nonblocking observation-scale diagnostic**.

Do not require the WT model to reproduce any minute-by-minute AE4-null time signature. That is irrelevant to WT candidate admission.

Record, but do not silently veto on, changes in per-cell WT secretion relative to the original model unless an existing source-supported WT gate is violated.

## Stage E: freeze the WT-admissible candidate set

This is the critical firewall.

Before running any new 5% AE4 or AE2 genotype simulation under the reallocated candidates:

1. write the full candidate table;
2. identify every WT-admissible candidate using only Stages A-D;
3. retain every admissible root/share/calcium combination according to the predeclared rules;
4. write a frozen-candidate manifest with hashes of all parameter payloads, roots, WT trajectories and selection decisions;
5. commit and push this checkpoint.

The candidate set must not depend on the known direction or magnitude of the AE4-loss phenotype.

No candidate may be dropped after genotype evaluation because its result is inconvenient.

## Stage F: post-freeze genotype evaluation

Only after the WT-admissible checkpoint is pushed may the reallocated models be evaluated for genotype discrimination.

For every frozen WT-admissible model family:

### AE4 near-loss

Continue AE4 expression from `1.0` to `0.05` under that candidate's reallocated parameters using the same target-free numerical continuation principle as before.

Do not use the old Task 14B 5% resting state after changing capacities.

Do not attempt exact zero in Task 16.

For every valid 5% state and each calcium value, run production Radau and compute

`R_AE4_5pct = total_0_600_AE4_5pct / total_0_600_WT`

`D_AE4_5pct = 1 - R_AE4_5pct`.

### AE2 loss

Also evaluate exact AE2 deletion or the inherited validated AE2 deletion procedure for the same frozen candidate set, without changing any other parameter.

Report

`R_AE2 = total_0_600_AE2 / total_0_600_WT`

and the AE4-versus-AE2 contrast.

No genotype-specific rescaling is allowed.

## Stage G: interpretation after the frozen test

After all frozen genotype results exist, report how the 5% AE4 effect varies with the **realized WT AE4 chloride-loading share**.

At minimum show:

- target share;
- realized resting share;
- AE4 and NKCC1 capacity scales;
- WT resting physiology;
- WT total secretion;
- WT AE4 cancellation/productive fraction;
- 5% AE4 secretion ratio/reduction;
- AE2 deletion secretion ratio/reduction;
- routing family;
- calcium value;
- numerical status.

Analyze whether increasing WT AE4 productive chloride share changes the sign of the 5% AE4 secretion phenotype.

If a sign change occurs, report where it occurs in the predeclared panel. Do not call that share "the correct value" merely because it gives the desired phenotype.

Only after all results are frozen may the experimental approximately 35% reduction be shown as contextual magnitude. It is not an acceptance band, optimization target or ranking criterion.

Detailed minute-by-minute phenotype shape remains diagnostic only, not an acceptance criterion.

## What this task can and cannot establish

Task 16 can establish whether **simple capacity allocation** between the existing AE4 mechanism and NKCC1 can reconcile WT physiology with the direction of the AE4-loss secretion effect.

It cannot establish that a particular AE4 share is measured in vivo unless supported independently by experiment.

It cannot justify changing the AE4 Na/K routing or eliminating cation slip. Those are separate mechanism questions.

If no WT-admissible allocation produces reduced secretion after 5% AE4, the next task should examine the AE4 mixed-cation architecture itself rather than continuing to scale capacity.

If only extreme allocation factors survive or if increasing AE4 share creates pathological Na/K slip, state that plainly.

## Required classifications

Choose exactly one primary classification:

1. `WT-CONSTRAINED REALLOCATION FEASIBLE; AE4 NEAR-LOSS ROBUSTLY REDUCES SECRETION`
2. `WT-CONSTRAINED REALLOCATION FEASIBLE; AE4 NEAR-LOSS STILL INCREASES SECRETION`
3. `WT-CONSTRAINED REALLOCATION CANNOT INCREASE AE4 SHARE WITHOUT BREAKING WT PHYSIOLOGY`
4. `REALLOCATION RESULT DEPENDS STRONGLY ON SHARE OR ROUTING FAMILY`
5. `TASK 16 INVALIDATED BY TASK 15 SIGN OR SOURCE-MAPPING DEFECT`
6. `TASK 16 NUMERICALLY INCONCLUSIVE`

Do not classify as success merely because one root/share/calcium pair gives reduced secretion. Robustness across the frozen WT-admissible ensemble matters.

## Required artifacts

Create at minimum:

- `results/16_wt_chloride_allocation/allocation_contract.json`
- `results/16_wt_chloride_allocation/inherited_partition.csv`
- `results/16_wt_chloride_allocation/wt_rest_candidates.csv`
- `results/16_wt_chloride_allocation/wt_dynamic_candidates.csv`
- `results/16_wt_chloride_allocation/frozen_candidate_manifest.json`
- `results/16_wt_chloride_allocation/ae4_5pct_results.csv`
- `results/16_wt_chloride_allocation/ae2_results.csv`
- `results/16_wt_chloride_allocation/allocation_vs_genotype.csv`
- `results/16_wt_chloride_allocation/test_results.txt`
- `analysis/16_wt_chloride_allocation/method.md`
- `analysis/16_wt_chloride_allocation/wt_admissibility.md`
- `analysis/16_wt_chloride_allocation/genotype_interpretation.md`
- `analysis/16_wt_chloride_allocation/final_answer.md`

The final answer must state plainly:

- the inherited AE4 share of positive AE4+NKCC1 WT Cl loading;
- how much AE4/NKCC1 capacity reallocation is required to realize larger shares;
- which target shares remain WT-admissible;
- whether those allocations worsen or reduce the Task 15 Na/K cancellation problem;
- whether 5% AE4 then decreases or increases secretion;
- whether AE2 remains comparatively neutral;
- whether the result is robust across roots, routing families and calcium values;
- whether simple chloride-load allocation is sufficient or whether the AE4 mixed-cation mechanism itself must be revisited.

## Final actions

Run focused tests and inherited tests sufficient to establish conservation, source scaling, candidate-selection independence from genotype outcome, and complete ensemble accounting.

Commit and push the completed Task 16 work to `codex/task-16-wt-chloride-allocation`.

Do not edit `archive/`.

Do not merge to `main`.

Do not draft manuscript text.
