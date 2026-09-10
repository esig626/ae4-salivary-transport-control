# AGENTS.md

These instructions apply to all work in this repository.

## Current scientific objective

The active task is Task 17. Execute:

`prompts/17_fixed_wt_flux_repartition.md`

Work only on branch:

`codex/task-17-fixed-wt-flux-repartition`

Task 17 supersedes Task 16's moving-state allocation experiment.

## Scientific baseline

Use completed Task 14B scientific state:

`4d4403f279fc36aec940c9e4759486d02f07e3f6`

Task 15 established that AE4 signs and source mapping are correct, WT AE4 net loads chloride, and the model has large opposing Na/K branch cancellation. Task 16 changed AE4/NKCC1 capacities while freezing all other parameters and allowed WT state variables to move. Its failure does not answer the Task 17 question.

## Fixed-WT rule

For every inherited accepted WT root, the complete WT resting state is immutable.

Do not optimize, continue, or alter any WT state coordinate.

The Task 17 intervention changes only the decomposition of transport fluxes supporting that same state.

Use exactly three conditions:

- inherited baseline;
- AE4 positive chloride-loading share 0.10;
- AE4 positive chloride-loading share 0.30.

No other share values or refinements are allowed.

## Chloride repartition

Increase AE4's positive chloride-loading contribution by scaling its complete coupled transporter, never the chloride source alone.

Reduce all other positive basolateral chloride-loading pathways proportionally so the total positive basolateral chloride-loading pool at the unchanged WT state remains exactly the inherited total.

Negative chloride counterfluxes are not part of the positive-loading share denominator.

Do not alter apical CaCC simply to absorb the repartition.

## Rebalance all other equations at the same state

Because AE4 and the other chloride loaders have different stoichiometries, the chloride repartition will perturb Na, K, carbon and alkalinity balances.

Restore those balances by solving only for existing, biologically legitimate WT nuisance/capacity parameters. The net source of every conserved species at the unchanged WT state must remain exactly the inherited steady-state source, i.e. zero within inherited tolerances.

Before choosing balancing axes, compute the fixed-state source/sensitivity matrix and its rank. Use the smallest full-rank set of existing parameters that were previously unmeasured or treated as WT calibration/nuisance degrees of freedom. Do not select parameters using genotype outcomes.

Do not change AE4 stoichiometry, AE4 Na/K routing, AE4 thermodynamic ratios, calcium values, regulation, bath composition, geometry, or genotype-specific parameters.

## Phenotype firewall

The known AE4-loss phenotype must not be used to construct or select fixed-WT parameterizations.

First solve, verify, freeze, commit and push every fixed-WT feasible baseline/0.10/0.30 parameterization.

Only after that checkpoint may AE4 expression be reduced to 0.05 and AE2 loss evaluated.

## Numerical discipline

A fixed-WT candidate is valid only if the exact inherited WT state remains a steady state under the modified flux decomposition, with inherited charge, current, carbon, water and numerical tolerances satisfied.

Do not reject a target solely because one local solver fails. Distinguish structural rank deficiency, bound incompatibility and numerical failure.

Do not rerun expensive baseline simulations when hash-valid inherited outputs can be reused.

## General discipline

Do not edit `archive/`.

Do not merge to `main`.

Do not draft manuscript text.

Do not perform unrelated model reduction, GSPT or identifiability work.

Commit and push completed Task 17 work only to `codex/task-17-fixed-wt-flux-repartition`.
