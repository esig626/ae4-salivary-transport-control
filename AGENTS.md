# AGENTS.md

These instructions apply to all work in this repository.

## Current scientific objective

The active task is Task 20. Execute:

`prompts/20_ae4_pooled_cation_no_slip.md`

Work only on branch:

`codex/task-20-ae4-pooled-cation-no-slip`

Tracking issue: #19.

Repository inheritance is completed and merged Task 19 on `main`, merge commit:

`e98789b276607807f7c9c88ff53605008a19cb66`.

## Why Task 20 exists

The current physiological AE4 model permits independent Na and K branch fluxes. In the frozen WT ensemble those branches run in opposite directions and cancel by roughly 97–99%, leaving very small productive chloride transport.

Tasks 18 and 19 showed that neither chloride-load redistribution nor Na/K-pump impairment repairs the wrong-direction AE4-loss phenotype.

The 2016 AE4 experiments establish cation transport capability, electroneutrality and reversibility, but do not establish simultaneous independent opposing Na/K cycles in physiological mixed solution. Their physiological working model treats Na/K as a nonselective pooled cation pathway.

Task 20 therefore tests a less committal pooled-cation architecture.

## Mechanism rule

Do not delete the legacy `SR2_SHARED_112_QSS` model.

Add a separate pooled-cation AE4 mechanism with:

- one scalar net cycle flux;
- pooled Na+K thermodynamic activity;
- the current 1:1:2 Cl:cation:HCO3 stoichiometry retained only as a working-model assumption;
- electroneutrality and reversible local detailed balance;
- donor-side Na/K partitioning of the one net cation flux;
- no independent Na/K slip loop.

At every state the Na and K AE4 net sources must have the same sign or one must be zero.

The primary WT pooled-cation selectivity is equal Na/K weighting. Do not introduce a fitted cation-selectivity parameter in the whole-cell phenotype test.

## Evidence discipline

Do not claim that 1:1:2 stoichiometry, equal microscopic Na/K use, or donor-side partitioning are experimentally identified facts. They are the declared coarse-grained mechanism being tested.

The direct evidence to preserve is qualitative:

- Na transport support;
- K transport support;
- Cl/HCO3 exchange;
- electroneutrality;
- reversibility;
- broad monovalent-cation permissiveness.

Do not map the 2016 Hill coefficients to transported-ion count.

## Whole-cell scope

Retain all ten inherited WT conserved states as fixed physiological targets.

Use exactly two pooled-model WT loading conditions:

1. `P0`: pooled AE4 net chloride loading matched to the inherited legacy AE4 net chloride flux for that root;
2. `P10`: pooled AE4 carries 10% of the inherited positive basolateral chloride-loading pool.

Do not run a new 30% condition.

For each condition use the Task 18 fixed-WT inverse-rebalancing method. Old conductances and capacity boxes are reference values, not hard laws. Keep the exact WT conserved state fixed and solve existing uncertain capacities needed for full steady-state closure.

Use the Task 18 minimal-deviation objective and no genotype information.

## Phenotype firewall

Before AE4-loss or AE2-loss evaluation:

- finish all pooled-model WT inverse solutions;
- verify exact fixed-state preservation and full production closure;
- verify the pooled no-slip invariant;
- commit and push a frozen WT checkpoint.

Only after that checkpoint may genotype outcomes be accessed.

## Genotype evaluation

For each feasible pooled-model WT parameterization:

- run WT dynamics at calcium 0.10, 0.25 and 0.50 uM;
- reduce AE4 expression from 1.0 to 0.05 with all other parameters frozen;
- evaluate matched AE2 loss;
- do not refit either genotype;
- do not attempt exact AE4 zero.

The approximately 35% experimental AE4-loss deficit is context only after results are frozen.

Report resting Cl, pH, Na, K and major compensating fluxes as secondary diagnostics.

## Comparator discipline

Reuse hash-valid legacy-model results rather than recomputing them:

- legacy inherited allocation;
- legacy 10% AE4-loading results from Task 18.

The key comparison is whether removing independent Na/K slip changes the AE4-loss phenotype direction at matched productive chloride loading.

## General discipline

Do not edit `archive/`.

Do not draft manuscript text.

Do not perform unrelated calibration, model reduction, GSPT or identifiability work.

Commit and push completed Task 20 work only to `codex/task-20-ae4-pooled-cation-no-slip`.

After completion, open a PR to `main` but do not merge automatically. Report the Task 20 scientific result before merge.
