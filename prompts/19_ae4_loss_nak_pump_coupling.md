# Task 19: AE4-loss / Na-K pump coupling diagnostic

## Mission

Task 18 showed that reallocating more WT chloride loading to AE4 does not repair the wrong-direction AE4-loss secretion phenotype. The next mechanistic hypothesis is narrower:

> Does AE4 loss become secretion-limiting if loss of AE4 function also reduces effective Na/K-ATPase support?

This is a diagnostic mechanism test, not a fit to the experimental ~35% secretion reduction.

Work on branch:

`codex/task-19-ae4-loss-nak-pump-coupling`

Tracking issue: #17.

Repository inheritance is merged Task 18 on `main`, merge commit:

`52a56cc8e8788cdb5cae60dcbaa75257de2878dd`.

Use the inherited baseline chloride-allocation model only for this task. Do not use the Task 18 0.10 or 0.30 rebalanced WT parameterizations in the primary experiment. The point is to isolate the AE4 × pump interaction without importing extreme rebalancing changes.

Do not edit `archive/`.

## Scientific question

The current model keeps Na/K-pump capacity fully available after AE4 is reduced to 5%, and AE4 near-loss therefore increases secretion in the inherited baseline model.

Test whether reducing effective Na/K-pump capacity together with AE4 loss changes that direction.

A reduction in pump capacity is itself expected to affect secretion. Therefore the experiment must include matched WT-AE4 pump-reduction controls. A lower secretion value in the combined arm is not evidence for an AE4-specific mechanism unless the matched factorial comparisons are also reported.

## Fixed sparse factorial panel

Use exactly:

AE4 expression:

- `1.0`
- `0.05`

Na/K-pump capacity scale:

- `1.0`
- `0.50`
- `0.10`

Calcium:

- `0.10 uM`
- `0.25 uM`
- `0.50 uM`

Retain all ten inherited roots and both AE4 routing families.

Do not add intermediate pump scales, threshold searches, or phenotype-guided refinements.

Where exact hash-valid scale-1.0 WT and AE4-0.05 resting states/trajectories already exist from the completed baseline work, reuse them rather than recomputing them.

The only new scientific pump scales are 0.50 and 0.10.

## Pump scaling rule

For a pump scale `s`, multiply the complete effective Na/K-ATPase capacity by `s` while preserving:

- 3 Na : 2 K stoichiometry;
- pump direction and thermodynamic law;
- the inherited apical/basolateral pump partition;
- all other transporter parameters.

If the production model stores total capacity plus an apical fraction, scale total capacity by `s` and keep the fraction fixed.

If it stores separate effective apical and basolateral capacities, multiply both by exactly `s`.

Do not change pump affinity constants, reversal logic, membrane topology, or stoichiometry.

## Hard freeze outside AE4 expression and pump scale

Keep fixed:

- inherited baseline chloride allocation;
- AE4 carrier architecture and Na/K routing;
- AE4 thermodynamic ratios;
- NKCC1;
- AE2;
- NHE1;
- CaCC;
- K channels;
- paracellular pathways;
- acid-base chemistry and buffers;
- water and lumen transport;
- bath composition;
- geometry;
- beta/cAMP/PKA regulation;
- calcium protocol;
- solver and conservation tolerances.

Do not re-fit any capacity after changing AE4 expression or pump scale.

Do not introduce genotype-specific compensation other than the explicitly declared pump scale.

## Stage A: freeze the experiment contract

Before running new trajectories, create and commit a machine-readable contract containing:

- exact branch and baseline commit;
- ten root IDs;
- both AE4 routing families;
- AE4 expressions `{1.0, 0.05}`;
- pump scales `{1.0, 0.50, 0.10}`;
- calcium values `{0.10, 0.25, 0.50}`;
- precise rule for pump scaling;
- continuation/root-selection rule;
- numerical gates;
- definitions of all primary comparison metrics;
- declaration that no pump scale may be added after results are seen.

Do not use the experimental phenotype to define the contract.

## Stage B: resting-state continuation

For pump scales 0.50 and 0.10, obtain connected resting states for both AE4-expression arms.

Use deterministic continuation from the corresponding inherited root. Preserve root identity and report failures rather than jumping to a disconnected root because its phenotype is convenient.

A recommended numerical path is:

1. from `(AE4=1.0, pump=1.0)` continue pump scale to the target with AE4 fixed at 1.0;
2. from that matched pump-reduced WT state continue AE4 expression from 1.0 to 0.05 at the same fixed pump scale.

Internal continuation steps are numerical only and are not additional scientific pump conditions.

Do not change any other parameter to rescue a failed root.

For every resting state record at minimum:

- Na_i;
- K_i;
- Cl_i;
- pH_i;
- cell volume;
- membrane voltages;
- total/apical/basolateral pump fluxes;
- AE4 Na and K branch fluxes and net Cl source;
- NKCC1, NHE1 and AE2 fluxes;
- full resting residual/conservation diagnostics.

## Stage C: dynamic simulations

For every valid factorial state run the 600 s CCh + IPR production protocol at all three calcium values using production Radau.

Use inherited numerical and conservation gates.

Reuse exact scale-1 trajectories when valid and hash-matched.

Do not use absolute whole-gland flow scale or exact minute-by-minute AE4-null shape as acceptance gates.

## Stage D: mandatory factorial comparisons

For every root and calcium, define total 0-600 s secretion `Q(e,s)` where `e` is AE4 expression and `s` is pump scale.

Report all of the following.

### 1. Matched AE4 effect at each pump scale

`R_AE4(s) = Q(0.05,s) / Q(1.0,s)`

This asks whether AE4 near-loss is harmful when pump capacity is held at the same level in the matched control.

### 2. Pump-only effect in AE4-intact cells

`R_pump_WT(s) = Q(1.0,s) / Q(1.0,1.0)`

This quantifies how much secretion falls simply because the pump has been impaired.

### 3. Pump effect in AE4-low cells

`R_pump_AE4(s) = Q(0.05,s) / Q(0.05,1.0)`

### 4. Interaction ratio

For `s < 1`, define

`I(s) = R_AE4(s) / R_AE4(1.0)`.

`I(s) < 1` means pump reduction makes AE4 loss more deleterious relative to its matched WT control.

### 5. Coupled-hypothesis phenotype

Also report

`R_coupled(s) = Q(0.05,s) / Q(1.0,1.0)`.

This is the phenotype predicted by a hypothetical rule in which pump support is reduced specifically when AE4 is nearly absent.

Do not call `R_coupled(s) < 1` an AE4-specific success by itself. The pump-only and interaction results must be shown beside it.

## Stage E: interpretation discipline

The experimental approximately 35% AE4-loss secretion reduction may be shown only after all simulations and comparisons are frozen.

It is context only.

Do not select 0.50 or 0.10 because one is closer to 0.65.

Do not fit a continuous AE4-to-pump coupling law in Task 19.

Do not claim that AE4 experimentally regulates Na/K-ATPase merely because this intervention improves the phenotype. Task 19 tests model sufficiency of such a coupling, not its biological truth.

A useful mechanism signal would require more than severe pump inhibition simply reducing flow. Report whether the matched AE4 ratio and interaction change materially.

## Computational budget

This task is intentionally bounded.

The full scientific panel contains:

- 10 roots;
- 2 AE4 expressions;
- 3 pump scales;
- 3 calcium values.

Scale-1.0 outputs should be reused wherever hashes match, so only pump scales 0.50 and 0.10 require new scientific trajectories.

Do not expand the panel.

## Required classifications

Choose exactly one primary classification:

1. `PUMP IMPAIRMENT MAKES AE4 NEAR-LOSS ROBUSTLY HARMFUL WITH A NONTRIVIAL INTERACTION`
2. `PUMP IMPAIRMENT LOWERS SECRETION BUT DOES NOT CREATE AN AE4-SPECIFIC INTERACTION`
3. `AE4 NEAR-LOSS REMAINS WRONG-DIRECTION AT ALL TESTED PUMP SCALES`
4. `RESULT DEPENDS STRONGLY ON ROUTING FAMILY, ROOT OR CALCIUM`
5. `SEVERE PUMP IMPAIRMENT BREAKS PHYSIOLOGY OR NUMERICS BEFORE THE HYPOTHESIS CAN BE TESTED`
6. `TASK 19 NUMERICALLY INCONCLUSIVE`

Do not classify the task as mechanistic success solely because `Q(0.05,0.10)` is smaller than normal WT secretion.

## Required artifacts

Create at minimum:

- `results/19_ae4_loss_nak_pump_coupling/contract.json`
- `results/19_ae4_loss_nak_pump_coupling/resting_states.csv`
- `results/19_ae4_loss_nak_pump_coupling/dynamic_results.csv`
- `results/19_ae4_loss_nak_pump_coupling/factorial_comparisons.csv`
- `results/19_ae4_loss_nak_pump_coupling/state_flux_diagnostics.csv`
- `results/19_ae4_loss_nak_pump_coupling/test_results.txt`
- `analysis/19_ae4_loss_nak_pump_coupling/method.md`
- `analysis/19_ae4_loss_nak_pump_coupling/mechanism_interpretation.md`
- `analysis/19_ae4_loss_nak_pump_coupling/final_answer.md`

The final answer must state plainly:

- whether 50% or 10% pump capacity permits valid resting/dynamic solutions;
- whether matched AE4 loss becomes secretion-reducing at either pump level;
- how much secretion is reduced by pump impairment alone;
- whether there is a genuine AE4 × pump interaction;
- whether the coupled AE4-low/pump-low state reproduces the experimental direction;
- whether results are robust across roots, routing families and calcium;
- what happens to intracellular Na/K/Cl and the relevant transporter fluxes;
- whether the result motivates a biologically explicit AE4-pump coupling model or rejects this hypothesis.

## Tests and final actions

Run focused tests proving:

- exact pump scaling on both membranes;
- preserved 3Na:2K stoichiometry;
- no changes to non-pump parameters;
- correct reuse of scale-1 results;
- complete factorial accounting;
- conservation and numerical gates;
- no phenotype-guided pump-scale selection.

Commit and push the completed Task 19 branch.

After the task is complete and all focused/inherited tests required by the task pass, open a pull request from `codex/task-19-ae4-loss-nak-pump-coupling` to `main` and merge the completed task into `main` using a merge commit.

If a genuine test/CI failure or merge conflict exposes a scientific or code defect, do not hide it or force the merge. Report the defect instead.

Do not draft manuscript text.
