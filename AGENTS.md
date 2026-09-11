# AGENTS.md

These instructions apply to all work in this repository.

## Current scientific objective

The active task is Task 21. Execute:

`prompts/21_pooled_ae4_wt_physiology_recalibration.md`

Work only on branch:

`codex/task-21-pooled-ae4-wt-physiology-recalibration`

Tracking issue: #21.

Repository inheritance is merged Task 20 on `main`, merge commit:

`3f41579ef81f1e73cc37dc1ce4ab3380f98c4659`.

## Scientific premise

Task 20 established that the exact legacy WT states are thermodynamically incompatible with positive chloride loading under the pooled-cation 1:1:2 no-slip AE4 law. Those exact states were generated under the legacy independent Na/K branch mechanism and are no longer hard biological targets.

Task 21 therefore recalibrates the pooled/no-slip model against WT physiology rather than preserving the legacy WT state byte-for-byte.

## AE4 mechanism

Use the Task 20 `POOLED_CATION_112_NO_SLIP` implementation unchanged.

Do not:

- reintroduce independent Na/K branches;
- permit an opposing Na/K slip cycle;
- change the 1:1:2 working stoichiometry in this task;
- add a fitted cation-selectivity mechanism;
- add a new transporter or pathway.

The pooled mechanism must retain one scalar net cycle with Na and K following the same net cycle direction.

## WT data and physiology

Primary measured WT targets:

- resting intracellular chloride `50.10 ± 1.50 mM`;
- resting pH `6.91 ± 0.07`.

Use ±2 reported SEM screening bands:

- chloride `[47.10, 53.10] mM`;
- pH `[6.77, 7.05]`.

Before optimisation, predeclare broad physiological envelopes for intracellular Na, intracellular K and cell volume from the repository provenance/source ledger and existing physiological references. Record the source and rationale. Do not derive these bounds from genotype behaviour.

The legacy WT Na/K/volume values are starting points or provenance references only, not exact targets.

## WT calibration

Allow the full WT conserved state to re-equilibrate.

Allow existing uncertain capacities/conductances to recalibrate as needed, including AE4, NKCC1, NHE1, pumps, K conductance, CaCC, CO2 exchange and paracellular capacities.

Historical capacity values and old calibration boxes are reference values, not hard feasibility walls.

Every accepted WT solution must satisfy:

- full production steady-state closure;
- membrane-current closure;
- all conservation identities;
- chloride and pH screening bands;
- predeclared Na/K/volume physiological envelopes;
- finite physical states and capacities;
- positive pooled AE4 affinity for productive chloride loading;
- strictly positive pooled AE4 chloride flux;
- pooled no-slip invariant;
- valid production dynamics at calcium 0.10, 0.25 and 0.50 uM.

Do not impose the Task 18 10% or 30% chloride-share conditions.

Do not match the pooled AE4 flux to the legacy AE4 flux.

## WT-only objective and root discipline

Use WT information only.

Rank admissible WT solutions by:

1. fit to measured WT chloride and pH using their reported SEM scales;
2. minimum largest absolute log-fold capacity departure from inherited references;
3. minimum summed squared log-fold capacity departure.

Do not penalise state movement away from the legacy WT state beyond the declared physiological constraints.

Use all ten inherited WT roots as independent numerical seeds. Deduplicate identical solutions and retain every distinct admissible WT calibration. Do not select solutions using later genotype results.

## Phenotype firewall

Before any AE4-loss or AE2-loss calculation:

1. freeze the physiology contract;
2. complete the entire WT calibration ensemble;
3. run and validate all WT dynamics;
4. write all WT states, parameters and flux ledgers;
5. commit and push the complete WT set;
6. record and remotely verify the WT checkpoint SHA.

Only after this checkpoint may genotype outputs be evaluated.

## Held-out genotype tests

For every frozen WT calibration:

- reduce AE4 expression to 0.05 with all other calibrated parameters fixed;
- evaluate matched AE2 loss using existing project semantics;
- allow each genotype state to re-equilibrate only through the production equations;
- do not refit or compensate any capacity;
- run calcium 0.10, 0.25 and 0.50 uM;
- report matched secretion ratios and resting-state shifts.

The approximately 35% experimental AE4-loss secretion deficit is held-out context only after all genotype results are frozen. Never fit to it.

## General discipline

Do not edit `archive/`.

Do not draft manuscript text.

Do not perform unrelated model reduction, GSPT or identifiability work.

Commit and push completed Task 21 work only to `codex/task-21-pooled-ae4-wt-physiology-recalibration`.

After completion, open a PR to `main` but do not merge automatically. Report the Task 21 scientific result before merge.
