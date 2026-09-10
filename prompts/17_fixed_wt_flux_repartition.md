# Task 17: fixed-WT flux repartition with increased AE4 chloride share

## Mission

Task 16 changed AE4 and NKCC1 capacities and then allowed the WT state to move. That was not the intended physiological question.

Task 17 asks the correct question:

> Can the already accepted WT steady state be held fixed exactly while AE4 is assigned a larger share of the existing basolateral chloride-loading flux, with the remaining transport fluxes rebalanced so that every WT steady-state equation remains unchanged?

This is a fixed-state flux repartition problem. It is not a new WT state search and it is not a fit to the known AE4-loss phenotype.

Work on branch:

`codex/task-17-fixed-wt-flux-repartition`

Base scientific state: completed Task 14B commit

`4d4403f279fc36aec940c9e4759486d02f07e3f6`

Task 15 established the diagnostic facts that the implemented signs/source mapping are correct, WT AE4 net loads chloride, and the current mixed-cation model has large opposing Na/K branch cancellation. Task 16 showed only that a two-capacity perturbation with all other parameters frozen drives the state away from WT. Do not treat Task 16 as evidence that a larger AE4 chloride share is impossible.

Do not edit `archive/`. Do not merge to `main`.

## Core principle: keep the WT state fixed

For each of the ten accepted WT roots, let `x_WT` denote the exact inherited resting state.

In Task 17:

- `x_WT` is immutable;
- intracellular Na, K, Cl, TIC, alkalinity, pH and volume must remain exactly the inherited values;
- luminal Na, K, Cl, TIC, alkalinity, pH and volume must remain exactly the inherited values;
- no state coordinate is optimized or continued;
- the resting RHS at `x_WT` must remain zero to the inherited numerical tolerances after flux repartition;
- conservation, charge and current closure must remain valid.

The task changes only the decomposition of fluxes that support this same WT state.

## Three scientific conditions only

Use exactly:

- inherited baseline;
- AE4 positive chloride-loading share `0.10`;
- AE4 positive chloride-loading share `0.30`.

Do not add intermediate shares, threshold searches or refinements.

## Stage A: baseline flux ledger at the fixed WT state

For every WT root, evaluate the complete resting source contribution of every transporter and membrane pathway at the exact inherited state.

At minimum report source vectors for:

- AE4 Na branch;
- AE4 K branch;
- total AE4;
- NKCC1;
- AE2;
- NHE1;
- Na/K pump;
- apical and basolateral K channels;
- apical CaCC chloride flux;
- paracellular Na/K/Cl/HCO3;
- basolateral and apical CO2 exchange;
- water and lumen outflow.

Record the full independent resting RHS and verify it is the inherited zero.

Identify the set `P` of positive basolateral chloride-loading pathways other than AE4 at that exact state. Negative chloride counterfluxes are not members of `P`.

Define

`A0 = positive AE4 intracellular Cl source`

`O0 = sum of positive intracellular Cl sources from P`

`Lplus = A0 + O0`.

The target AE4 share is defined relative to this unchanged positive-loading pool.

## Stage B: impose the chloride repartition at the same state

For target share `f` in `{0.10, 0.30}` require at the exact same `x_WT`:

`A_target = f * Lplus`

and scale the complete AE4 cycle so that its intracellular chloride source equals `A_target`.

Do not scale chloride alone. The full AE4 Na/K/Cl/HCO3/TIC/alkalinity source vector must scale consistently.

All other positive basolateral chloride-loading pathways in `P` must be reduced by one common proportional factor `rho` so that

`A_target + rho * O0 = Lplus`.

Thus the total positive basolateral chloride-loading flux is unchanged and the relative proportions among the non-AE4 positive loaders are preserved.

Do not use AE2 or any negative chloride counterflux to define the positive-loading share.

Do not change apical CaCC merely to absorb the imposed repartition. The resting chloride source balance at the fixed WT state must be recovered by the prescribed loading repartition plus the permitted balancing fluxes below.

## Stage C: rebalance every other WT equation, not the state

The coupled AE4/NKCC-style repartition will generally perturb Na, K, carbon and alkalinity source components even though total positive chloride loading is fixed.

The goal is to cancel those source perturbations by adjusting existing transport capacities while leaving `x_WT` unchanged.

This is an equality-constrained parameter/flux solve at fixed state.

Before solving the target shares, construct the fixed-state sensitivity/source matrix of candidate balancing mechanisms and determine its rank.

Balancing mechanisms must satisfy all of the following:

1. they already exist in the model;
2. their capacities are unmeasured or previously treated as WT calibration/nuisance degrees of freedom;
3. they are not selected using the AE4-loss phenotype;
4. they do not change AE4 stoichiometry, AE4 Na/K routing or AE4 branch thermodynamics;
5. they do not create a new chloride-loading share outside the Stage B proportional rule.

Prefer the smallest full-rank set. Candidate balancing axes may include, where mathematically required and already present in the model:

- NHE1 capacity for Na/alkalinity balance;
- Na/K pump capacity for Na/K balance;
- K-channel conductance/capacity for K/current balance;
- neutral CO2 exchange capacity for TIC balance;
- other pre-existing WT nuisance capacities only if rank analysis proves they are required.

Do not use a common membrane scale if it necessarily changes CaCC chloride transport and thereby defeats the fixed chloride-allocation rule. If a Cl-carrying pathway must change for exact fixed-state balance, document that conflict rather than silently allowing it.

For each root and target share solve only for balancing parameter multipliers. Never optimize state coordinates.

Require the complete independent resting RHS at the inherited `x_WT` to match the baseline resting RHS within inherited tolerances. Since the baseline is a steady state, this means the modified RHS must also be zero.

Require also:

- charge conservation;
- carbon/alkalinity accounting;
- current closure;
- finite nonnegative capacities;
- no source-sign reversals that violate the retained biological direction constraints;
- all modified parameter values within previously declared admissible or source-supported bounds where such bounds already exist.

Do not invent new bounds simply to reject an inconvenient result.

## Stage D: exact interpretation of feasibility

A target share is fixed-WT feasible only if the exact same inherited state remains a valid steady state after the flux repartition and balancing solve.

Report for every root:

- target share;
- AE4 capacity multiplier;
- proportional reduction factor `rho` on other positive Cl loaders;
- every balancing capacity multiplier;
- maximum absolute/scaled RHS residual at the unchanged WT state;
- all conservation/current residuals;
- whether the fixed state is exactly retained;
- resulting AE4 Na/K cancellation and productive fraction at the unchanged state.

Do not call a target infeasible merely because a local optimizer fails. Use rank diagnostics, deterministic multistarts or an equivalent robust algebraic solve sufficient to distinguish structural inconsistency from numerical failure.

If the fixed-state equations are underdetermined, preserve the full feasible family or choose the minimum-change solution by a predeclared norm that depends only on distance from the baseline parameter vector. Do not use genotype outcomes to choose among solutions.

## Stage E: freeze fixed-WT parameterizations before genotype evaluation

Before any new AE4-loss simulation:

1. save every baseline, 0.10 and 0.30 fixed-state result;
2. identify all fixed-WT feasible parameterizations using only Stages A-D;
3. freeze and hash the complete parameter payloads;
4. commit and push the fixed-WT checkpoint.

The known AE4-loss phenotype must not be used to select among feasible WT parameterizations.

## Stage F: post-freeze genotype test

Only after the fixed-WT checkpoint is pushed may genotype simulations be performed.

For each fixed-WT feasible modified parameterization:

- keep the frozen modified parameters;
- continue AE4 expression from 1.0 to 0.05 to obtain the corresponding near-loss resting state;
- do not reuse the old Task 14B 5% state after parameters change;
- do not attempt exact zero;
- run the same 600 s CCh + IPR protocol at calcium 0.10, 0.25 and 0.50 uM;
- compute total 0-600 s secretion ratio relative to that parameterization's WT trajectory;
- evaluate matched AE2 loss as the specificity control.

Reuse baseline trajectories where hashes and parameters are unchanged.

Only after all genotype results are frozen may the experimental approximately 35% AE4-null reduction be shown as context. It is not an optimization target or acceptance band.

## Critical scientific distinction

Task 17 is not allowed to claim that AE4 can transport chloride independently of its coupled ions.

Increasing AE4 chloride loading necessarily scales its complete coupled source vector. What must remain unchanged is the **net WT balance of every other ion**, achieved by compensating changes in the rest of the network while keeping the WT state fixed.

This distinction is mandatory.

## Required classifications

Choose exactly one primary classification:

1. `FIXED-WT FLUX REPARTITION FEASIBLE; HIGHER AE4 SHARE PRESERVES THE SAME WT STATE`
2. `FIXED-WT FLUX REPARTITION FEASIBLE ONLY FOR A SUBSET OF ROOTS OR SHARES`
3. `FIXED-WT FLUX REPARTITION STRUCTURALLY INCOMPATIBLE WITH THE EXISTING NETWORK`
4. `FIXED-WT FLUX REPARTITION NUMERICALLY INCONCLUSIVE`

If fixed-WT repartition is feasible, also report the post-freeze AE4 5% secretion direction and AE2 specificity without turning those outcomes into the WT feasibility criterion.

## Required artifacts

Create at minimum:

- `results/17_fixed_wt_flux_repartition/contract.json`
- `results/17_fixed_wt_flux_repartition/baseline_flux_ledger.csv`
- `results/17_fixed_wt_flux_repartition/fixed_state_source_matrix.csv`
- `results/17_fixed_wt_flux_repartition/fixed_state_rank.json`
- `results/17_fixed_wt_flux_repartition/fixed_wt_solutions.csv`
- `results/17_fixed_wt_flux_repartition/frozen_fixed_wt_manifest.json`
- `results/17_fixed_wt_flux_repartition/ae4_5pct_results.csv`
- `results/17_fixed_wt_flux_repartition/ae2_results.csv`
- `results/17_fixed_wt_flux_repartition/test_results.txt`
- `analysis/17_fixed_wt_flux_repartition/method.md`
- `analysis/17_fixed_wt_flux_repartition/final_answer.md`

## Final actions

Run focused tests sufficient to prove that no WT state coordinate moved, the complete AE4 source vector scaled consistently, the prescribed chloride-loading repartition was respected, and every accepted fixed-WT solution satisfies the full resting balances.

Commit and push completed Task 17 work only to `codex/task-17-fixed-wt-flux-repartition`.

Do not edit `archive/`.

Do not merge to `main`.

Do not draft manuscript text.
