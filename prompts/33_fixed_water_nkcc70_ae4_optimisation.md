# Task 33: fixed-water WT optimisation with NKCC1 at 70% and frozen mechanistic NHE1

Repository: `esig626/ae4-salivary-transport-control`

## Purpose

Stop re-diagnosing already established defects. Build one WT parameterisation in which the existing routed AE4 mechanism has a substantial productive chloride-loading role, while preserving the mechanistic NHE1 repair and the accepted WT water throughput. Then freeze the WT parameterisation and test only AE4 expression 0.05 and 0.0.

This is a bounded WT optimisation task. The 5% and 0% genotype outcomes are strictly held out until the WT solution is frozen.

## Fixed model components

Start from Task 31 commit `fb4bc7d9a569051f5e41cf84837046cba6952d4f` and its accepted R09 WT state.

Do not change:

- the Cha et al. 2009 eight-state + Mod2 NHE1 equations;
- the Task 31 shared NHE1 carrier amount `2.339370005697548e-05 fmol`;
- AE4 routing topology, donor-side Na/K partition, or 1:1:2 working stoichiometry;
- AE2, NKCC1, pump, K-channel, CaCC, CO2, paracellular, acid-base, or water functional forms;
- stimulus/regulatory equations;
- bath composition or geometry.

Do not introduce a new transporter, pH clamp, hidden acid/base source, genotype-specific parameter, or signalling rescue.

## WT reference and fixed water target

Use R09 only for optimisation.

The accepted Task 31 R09 WT state is the baseline physiological reference. Do **not** force intracellular chloride to 50.1 mM. Its approximately 60 mM chloride is not a failure criterion.

Hold the WT resting ductal/luminal outflow target fixed to the accepted Task 31 value:

`q_out_WT = 0.0010757073853493032 pL/s`.

Keep all water-law parameters fixed. The optimiser must meet this water-throughput target through ionic/osmotic balance, not by changing water permeability or the outflow law. Use a tight numerical equality tolerance consistent with the stationary solver, and report the achieved relative error.

Also preserve a physiologically ordinary WT state. Use the Task 31 accepted R09 WT values as reference values, not hard equality targets, except pH should remain in the experimental WT interval 6.84--6.98, Na/K must remain positive and ordinary, volume < 3 pL, and all charge/current/carbon/water closure gates must pass.

## NKCC1 70% WT chloride-loading constraint

At the final WT resting state, require NKCC1 to supply exactly 70% of the **positive basolateral chloride-loading pool**:

`share_NKCC1 = J_NKCC1_Cl_positive / J_positive_total = 0.70`.

Use `J_NKCC1_Cl = 2 * J_NKCC1_cycle`.

Define the positive pool as the sum of positive chloride-loading contributions actually present in the model:

- `max(J_NKCC1_Cl, 0)`;
- `max(J_AE4_Cl, 0)`;
- `max(J_AE2_Cl, 0)`;
- any other implemented positive basolateral chloride loader, if one exists.

Signed counterfluxes must be reported separately and must not be hidden inside the denominator.

The 70% constraint is an explicit WT modelling constraint motivated by the experimental chloride-uptake partition. Do not use the AE4-null secretion phenotype to set it.

Do **not** separately impose an arbitrary AE4 percentage. AE4's final positive loading share is an output of the constrained WT optimisation. However, reject a solution in which AE4 is non-positive or merely negligible (<5% of the positive loading pool), because such a solution does not satisfy the purpose of this task.

## Optimisation

This task explicitly authorises one bounded low-dimensional WT optimisation. It does **not** authorise a Cartesian grid, exhaustive combinations, random/global/evolutionary/Bayesian search, Shapley analysis, or automatic expansion of the parameter set.

Before running the optimiser, write `analysis/33_fixed_water_nkcc70_ae4_optimisation/optimisation_contract.md` listing the exact free parameters, bounds, transformations and objective weights.

Use the smallest defensible free set needed to satisfy full WT closure plus the water and 70% NKCC1 constraints. Start from these existing uncertain magnitude parameters only:

1. NKCC1 capacity;
2. routed AE4 carrier amount/activity magnitude;
3. AE2 capacity;
4. total Na/K-pump capacity;
5. total K conductance magnitude;
6. apical CaCC conductance magnitude.

Keep all membrane allocation fractions fixed. Do not free NHE1 density. Do not free water, CO2, paracellular, stimulus, regulatory or AE4 cation-partition parameters.

If direct algebra removes one or more free parameters, do so. Prefer equality elimination over numerical optimisation. If a parameter is unnecessary for closure, leave it frozen rather than optimising it.

The WT objective is lexicographic / strongly constrained:

1. full stationary amount, current, charge, carbon and water closure;
2. fixed WT `q_out` target;
3. NKCC1 positive chloride-loading share = 0.70;
4. positive, non-negligible AE4 loading;
5. physiologically ordinary WT state, with pH 6.84--6.98 and volume < 3 pL;
6. among feasible solutions, minimise weighted log-deviation of free parameters from their Task 31 values and state deviation from the accepted Task 31 WT state.

Do not target knockout Cl, knockout pH, knockout secretion, or a desired AE4-null/WT ratio.

### Hard numerical budget

- one local constrained optimisation run, with at most one deterministic continuation/restart from its own last feasible point;
- at most 40 stationary solves total, including optimisation evaluations and validation;
- at most 60,000 stationary residual evaluations total;
- at most 2 stimulated WT integrations before freeze, only if needed to verify the fixed WT water target under the standard central protocol;
- after freeze, at most 4 genotype integrations total: 5% and 0% only, plus solver-method confirmations only if essential and still within this count;
- one worker, one BLAS thread;
- 45 minutes numerical execution maximum;
- checkpoint the best feasible WT point and counters continuously.

Stop at the first exhausted limit. Do not enlarge the free set after failure.

## Freeze before genotype release

When one WT parameterisation passes all gates:

1. write all free parameter values and complete WT state/flux ledger;
2. record the achieved `q_out` and NKCC1 70% constraint residuals;
3. record AE4 and AE2 signed chloride contributions and AE4 positive-pool share;
4. hash the frozen parameter payload;
5. publish/checkpoint the frozen WT text artifacts through the connected GitHub integration.

Only after that freeze may AE4 expression be changed.

## Genotype tests: exactly 5% and 0%

With every non-AE4 parameter frozen, test exactly:

- AE4 expression `0.05`;
- AE4 expression `0.00`.

WT expression `1.00` is the already frozen denominator/reference, not an additional optimisation arm.

For each of 5% and 0%:

1. solve its own REST state from the frozen WT model;
2. report admissibility, pH, Na, K, Cl, TIC/HCO3, cell/lumen volumes, voltages, osmolarities and all major transporter/water fluxes;
3. if and only if REST is admissible, run the standard central stimulation at Ca = 0.25 uM to 600 s;
4. report cumulative ductal secretion and ratio to frozen WT;
5. decompose replacement of lost AE4 chloride loading into NKCC1, AE2, altered apical delivery, paracellular return, ionic storage and water/osmotic effects.

Do not test 10%, 20%, 50%, AE2 loss, other Ca levels, other backgrounds, or any additional genotype/dose condition.

Do not refit after seeing either 5% or 0% result.

## Interpretation

The primary scientific question is now direct:

> With WT water throughput fixed, mechanistic NHE1 frozen, and NKCC1 constrained to 70% of positive WT chloride loading, can the existing AE4 architecture carry a meaningful share of productive chloride loading and does reducing AE4 to 5% or 0% reduce secretion without genotype fitting?

A successful WT optimisation does not itself establish the biological truth of the 70% partition. A failed genotype phenotype does not license further tuning in this task.

Report the 5% and 0% results exactly as obtained, including a neutral or wrong-direction effect.

## Outputs

Write compact UTF-8 outputs under:

- `analysis/33_fixed_water_nkcc70_ae4_optimisation/`
- `results/33_fixed_water_nkcc70_ae4_optimisation/`

At minimum include:

- `optimisation_contract.md`;
- `budget.json`;
- `frozen_wt_parameters.json`;
- `wt_state.json`;
- `wt_flux_ledger.csv`;
- `genotype_rest_states.csv`;
- `genotype_flux_ledger.csv`;
- `secretion_comparison.csv` if dynamics are legally run;
- `verification.json`;
- `final_answer.md`.

No binary trajectory arrays or archives in the publication deliverable.

## GitHub publication

Use the connected GitHub integration for all remote writes. Do not depend on shell Git authentication, `gh auth`, PATs or SSH credentials. Shell Git is only for local status/diff inspection. Publication failure must never cause scientific recomputation.

Do not merge to `main`.
