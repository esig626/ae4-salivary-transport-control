# AGENTS.md

These instructions apply to all work in this repository.

## Current scientific objective

The active task is Task 22. Execute:

`prompts/22_physiologically_constrained_wt_recalibration.md`

Work only on branch:

`codex/task-22-physiological-wt-domain`

Tracking issue: #23.

Repository inheritance is Task 21 final commit:

`db734f52b8d71d670fcf793beab922967212bdd6`

Task 21 PR #22 remains unmerged. Task 22 intentionally continues from that final state.

## Scientific premise

Task 21 did not disprove pooled/no-slip AE4. It found resting states with measured WT intracellular Cl and pH, positive pooled AE4 affinity and positive chloride loading, but the WT optimisation was underconstrained and escaped to absurd luminal concentrations, oversized lumen volumes and million-fold capacity changes.

Task 22 closes those calibration loopholes and retries WT before any genotype calculation.

## AE4 mechanism freeze

Use `POOLED_CATION_112_NO_SLIP` unchanged.

Do not reintroduce independent Na/K branches or slip. Do not alter the 1:1:2 working stoichiometry. Do not add pathways, transporters, currents, buffers or compartments.

## Predeclare the complete physiology domain

Before any new optimisation result is inspected, create, commit and push the Task 22 physiology contract.

Retain intracellular screens:

- Cl 47.10 to 53.10 mM;
- pH 6.77 to 7.05;
- Na 2 to 60 mM;
- K 60 to 200 mM;
- cell volume 0.3 to 5.0 pL.

Add luminal screens:

- Na 100 to 200 mM;
- K 1 to 30 mM;
- Cl 80 to 180 mM;
- pH 6.8 to 8.0;
- TIC 1 to 80 mM;
- lumen volume 0.02 to 0.50 pL.

Require model-computed cell and lumen osmolarity each to lie within 0.80 to 1.20 times bath osmolarity.

These broad luminal ranges are physiology/model screens, not exact mouse-SMG measurements. Primary acinar saliva is approximately isotonic/plasma-like; document the source basis.

## Capacity domain

For every adjustable positive capacity/conductance with a positive inherited reference, require candidate/reference in [0.01, 100].

Fractions remain in [0,1].

Do not silently activate a zero-reference pathway. Any justified nonzero absolute domain for an already-existing zero-reference coordinate must be declared before optimisation.

Do not widen capacity bounds after seeing failures.

## Anti-pathology screen

Use the same signed-ledger definition of the positive basolateral chloride-loading pool as the prior tasks.

No individual absolute resting transporter/source flux may exceed 100 times that positive loading pool.

This is a broad anti-cancellation/numerical-conditioning screen, not a fitted biological constant.

## WT admissibility

Every accepted WT solution must simultaneously satisfy:

- intracellular and luminal physiology screens;
- cell and lumen osmotic consistency;
- full resting production RHS closure;
- membrane-current closure;
- state-charge closure;
- all conservation/bookkeeping identities;
- finite physical states and bounded capacities;
- positive pooled AE4 affinity and positive chloride loading;
- pooled no-slip Na/K signs;
- anti-cancellation screen;
- valid 600 s production dynamics at calcium 0.10, 0.25 and 0.50 uM;
- all existing dynamic conservation/current/state-charge gates.

Do not relax scientific or numerical acceptance tolerances to obtain feasibility.

Purely numerical improvements such as scaling, continuation, higher internal precision or analytic Jacobians are allowed if equations, bounds and gates remain unchanged.

## Search and objective

Use all ten inherited WT roots and all Task 21 resting witnesses as numerical seeds. Additional WT-only starts are allowed for robustness.

Rank admissible solutions using WT information only:

1. fit to measured WT Cl and pH;
2. minimum largest absolute log-fold capacity departure;
3. minimum summed squared log-fold departure.

Retain every distinct admissible WT solution. Do not use genotype behaviour for selection.

## Phenotype firewall

Before any Task 22 genotype calculation OR historical genotype regression test:

1. complete the WT search;
2. validate all admissible WT dynamics;
3. write all states/parameters/flux ledgers;
4. commit and push the complete WT evidence set;
5. remotely verify the WT checkpoint SHA and manifest.

If the WT ensemble is empty, stop before genotype evaluation.

## Held-out genotype tests

Only after a non-empty frozen WT checkpoint:

- reduce AE4 expression to 0.05 with every other calibrated parameter frozen;
- evaluate matched AE2 loss with no refit;
- find connected genotype resting states through the production equations only;
- run calcium 0.10, 0.25 and 0.50 uM;
- report matched secretion ratios and resting-state changes.

The approximately 35% experimental AE4-loss deficit is held-out context only after all genotype results are frozen.

## Failure classification

If no WT succeeds, preserve all attempts and distinguish numerical nonconvergence, physiological-bound conflict, capacity-bound conflict, thermodynamic conflict, dynamic/conservation failure and genuine structural impossibility. Do not widen Task 22 bounds post hoc.

## General discipline

Do not edit `archive/`.

Do not draft manuscript text.

Commit and push completed Task 22 work only to `codex/task-22-physiological-wt-domain`.

Open a PR to `main` when complete but do not merge automatically. Report the Task 22 scientific result first.
