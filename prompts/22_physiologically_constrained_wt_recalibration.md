# Task 22: physiologically constrained WT recalibration

## Goal

Continue directly from Task 21 and determine whether the `POOLED_CATION_112_NO_SLIP` AE4 mechanism can support a complete physiologically admissible WT salivary-acinar model once the WT calibration domain is properly constrained.

Task 21 is not a mechanism failure. It found resting witnesses with measured WT intracellular Cl and pH, positive pooled AE4 affinity and positive chloride loading, but the optimizer escaped to nonphysiological luminal concentrations, enormous lumen volumes and million-fold capacity shifts. Task 22 must close those loopholes before any further genotype test.

Do not change the AE4 mechanism in this task.

## Repository inheritance

Start from Task 21 final commit:

`db734f52b8d71d670fcf793beab922967212bdd6`

Work only on:

`codex/task-22-physiological-wt-domain`

Tracking issue: #23.

Task 21 PR #22 remains unmerged. Task 22 intentionally branches from the Task 21 final state because it reuses the pooled AE4 implementation and WT calibration infrastructure.

Do not edit `archive/`.

## Mechanism freeze

Keep `POOLED_CATION_112_NO_SLIP` unchanged.

Do not:

- reintroduce independent Na/K branches or slip;
- change the 1:1:2 Cl:cation:HCO3 working stoichiometry;
- add transporters, currents, buffers or compartments;
- change bath chemistry, acid-base chemistry, source signs, geometry equations, water equations or regulation law;
- use genotype outcomes to choose WT bounds or parameters.

## Stage 0: freeze the full WT physiological domain before optimisation

Before running or inspecting any new Task 22 optimisation result, create and push:

- `results/22_physiological_wt_domain/precalibration_contract.json`
- `analysis/22_physiological_wt_domain/wt_physiology_contract.md`

The contract must contain every numerical bound below and the exact scientific/provenance rationale.

### Intracellular WT constraints

Retain the Task 21 screens unchanged:

- `[Cl-]i`: 47.10 to 53.10 mM;
- `pHi`: 6.77 to 7.05;
- `[Na+]i`: 2 to 60 mM;
- `[K+]i`: 60 to 200 mM;
- cell volume: 0.3 to 5.0 pL.

The Cl and pH bands are the predeclared ±2 reported SEM screens from the native WT measurements. The Na/K/volume intervals remain broad modelling/provenance screens, not measured uncertainty intervals.

### Luminal WT constraints

Primary acinar saliva is near isotonic and plasma-like. The Task 21 molar luminal concentrations are therefore inadmissible.

Use these broad hard screens:

- luminal Na: 100 to 200 mM;
- luminal K: 1 to 30 mM;
- luminal Cl: 80 to 180 mM;
- luminal pH: 6.8 to 8.0;
- luminal TIC: 1 to 80 mM;
- lumen volume: 0.02 to 0.50 pL.

Treat these as deliberately broad physiological/model envelopes, not exact mouse-submandibular measurements. The repository reference lumen volume of 0.10 pL is provenance only.

At minimum cite the salivary-physiology basis that primary acinar secretion is approximately isotonic/plasma-like and the published primary-fluid composition ranges already used by the project literature lineage. Do not tighten the ranges after inspecting optimisation results.

### Osmotic consistency

Use the production model's own osmolarity calculation.

At resting WT require:

- lumen osmolarity between 0.80 and 1.20 times bath osmolarity;
- intracellular osmolarity between 0.80 and 1.20 times bath osmolarity.

Do not replace this by a hand-written concentration sum if the production model already defines osmolarity including buffers/impermeant species.

### Capacity/conductance domain

Task 21's unbounded positive capacities produced million-fold excursions. Prevent that explicitly.

For every adjustable positive capacity or conductance with a positive inherited reference value, enforce:

`0.01 <= candidate/reference <= 100`

for the complete WT optimisation.

Fractions remain in `[0,1]`.

For a coordinate whose inherited reference is exactly zero:

- do not silently turn on a previously absent pathway;
- keep it at zero unless the pathway is already structurally present and Task 21 explicitly treated it as adjustable;
- if a nonzero absolute domain is scientifically necessary, declare that domain in the Stage-0 contract before optimisation and justify it from existing model provenance, not from phenotype outcome.

Do not widen any Task 22 capacity domain after seeing failures.

### Anti-cancellation screen

Record the absolute magnitude of every major resting transporter/source flux.

Define the resting positive basolateral chloride-loading pool using the same signed-ledger semantics as Tasks 18/20.

Require that no individual absolute transporter/source flux magnitude exceeds 100 times that positive chloride-loading pool.

This is a deliberately broad anti-pathology/numerical-conditioning screen. It is not a measured biological ratio and must not be fitted.

## WT admissibility

A Task 22 WT candidate is admissible only if all of the following pass simultaneously:

1. all intracellular physiology screens;
2. all luminal physiology screens;
3. cell and lumen osmotic-consistency screens;
4. complete production resting RHS within inherited tolerance;
5. membrane-current closure within inherited tolerance;
6. state-charge closure within inherited tolerance;
7. complete Na/K/Cl/TIC/alkalinity/water/lumen conservation bookkeeping;
8. all states finite and physically valid;
9. all capacities/permeabilities finite, nonnegative and inside the Stage-0 domains;
10. all fractions in `[0,1]`;
11. pooled AE4 affinity positive for the productive chloride-loading direction;
12. pooled AE4 net chloride flux strictly positive;
13. pooled no-slip Na/K sign invariant;
14. anti-cancellation screen;
15. valid 600 s production WT dynamics at calcium 0.10, 0.25 and 0.50 uM;
16. every dynamic trajectory passes the existing conservation, current and state-charge gates.

Do not relax solver/conservation tolerances to obtain feasibility.

Higher numerical precision, better variable scaling, analytic Jacobians, continuation, trust-region changes or other purely numerical improvements are allowed if they do not change the scientific equations, bounds or acceptance gates. Document them.

## Search strategy

Use:

- all ten inherited WT roots;
- every Task 21 resting witness;
- additional WT-only numerical starts if needed for robustness.

Task 21 witnesses are seeds only. Their pathological lumen values are not allowed to survive the Task 22 hard screens.

A local optimizer failure is not a structural impossibility proof.

Preserve failed attempts and distinguish:

- numerical nonconvergence;
- physiological-bound conflict;
- capacity-bound conflict;
- thermodynamic conflict;
- conservation/dynamic-gate failure;
- proven structural contradiction, if one is actually established.

## WT objective

Use WT evidence only.

Among admissible WT candidates rank lexicographically:

1. minimise normalized squared error of measured WT Cl and pH from 50.10 mM and 6.91;
2. minimise the largest absolute log-fold change of adjustable positive capacities from inherited references;
3. minimise the sum of squared log-fold capacity changes.

Do not penalise state coordinates for departing from Task 21 or legacy states beyond the declared physiological screens.

Retain every distinct admissible WT solution passing the same WT-only rule. Do not select using later genotype behaviour.

## WT checkpoint and phenotype firewall

This ordering is mandatory.

Before running or inspecting **any** Task 22 genotype calculation or historical genotype regression:

1. finish the complete Task 22 WT search;
2. write all admissible states, capacities, flux ledgers and WT trajectories;
3. verify all Task 22 gates;
4. commit the complete WT evidence set;
5. push it to `codex/task-22-physiological-wt-domain`;
6. remotely verify the pushed checkpoint SHA and manifest.

Historical genotype regression tests must not be invoked before this checkpoint.

If the admissible WT ensemble is empty, stop. Do not run AE4-low or AE2-loss.

## Held-out genotype evaluation

Only after a non-empty WT ensemble is frozen and remotely verified:

### AE4 near-loss

Set AE4 expression to 0.05 using existing genotype semantics.

Freeze every other calibrated parameter.

No compensation, no refit, no bound changes.

Find the connected genotype resting state and run 600 s at calcium 0.10, 0.25 and 0.50 uM.

Report matched `AE4_0.05 / WT` secretion ratios and resting-state shifts.

### AE2 specificity control

Evaluate AE2 loss with the same frozen WT parameterisation and no refit.

Run the same calcium panel and report matched `AE2_loss / WT` ratios.

The approximately 35% experimental AE4-loss secretion reduction remains held-out context until all Task 22 genotype outputs are frozen.

## Required outputs

At minimum create:

- `analysis/22_physiological_wt_domain/wt_physiology_contract.md`
- `analysis/22_physiological_wt_domain/methods.md`
- `analysis/22_physiological_wt_domain/final_answer.md`
- `results/22_physiological_wt_domain/precalibration_contract.json`
- complete WT search/decision table including all failed attempts;
- WT calibration table;
- luminal/intracellular physiology table;
- osmolarity diagnostics;
- capacity-fold table;
- resting flux and anti-cancellation diagnostics;
- frozen per-solution WT payloads;
- remotely verified WT checkpoint receipt;
- WT dynamic results;
- AE4 5% results only if WT succeeds;
- AE2-loss results only if WT succeeds;
- matched comparison table only if genotype evaluation is authorised;
- focused and broad validation transcripts.

Add focused tests for every new Task 22 screen and for checkpoint ordering.

## Completion

Commit and push Task 22 only to:

`codex/task-22-physiological-wt-domain`

Open a PR to `main` when complete. Do not merge automatically. Report the scientific result for review first.
