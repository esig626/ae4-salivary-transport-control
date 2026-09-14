# Task 40 - equal Na/K AE4 routing and frozen phenotype validation

## Objective

Test exactly one structural change to the completed Task 39 model:

replace donor-concentration-weighted AE4 cation routing with a fixed 50:50 Na/K split while preserving the inherited net AE4 cycle law, total cation export, Cl/HCO3 stoichiometry, capacity, regulation, and every non-AE4-routing mechanism.

Read `analysis/40_ae4_equal_cation_routing/source_contract.md` first and treat it as binding.

This is not a parameter search and not permission to keep changing the model until the phenotype appears.

## Parent model

Start from merged Task 39 main commit:

`afd101448763439f369f2682c467f069ea18a442`

Task 39 is the frozen comparator. Its reported AE4-loss NKCC1 compensation is approximately:

- `+29.65%` at AE4=0.05
- `+33.10%` at AE4=0.00.

Task 39 cumulative secretion deficits are approximately:

- `0.482%` at AE4=0.05
- `0.627%` at AE4=0.00.

These are comparison values only.

## Phase A - implement exactly one routing change

Implement a named selectable equal-cation routed AE4 evaluator. Prefer the smallest possible focused module or selector around the existing routed AE4 evaluator.

Let `J4` be the inherited net AE4 chloride source, positive inward.

Keep `J4` exactly unchanged.

Set:

`S_Na_AE4 = -0.5 * J4`

`S_K_AE4 = -0.5 * J4`

`S_Cl_AE4 = +J4`

`S_TIC_AE4 = -2 * J4`

`S_TA_AE4 = -2 * J4`.

Negative `J4` reverses all sources automatically.

Do not inspect donor Na/K concentrations when choosing the split.

Do not change:

- Palk/Benjamin NKCC1 law or its `alpha_eff`;
- NKCC1 activity multiplier;
- NBC;
- NHE1;
- AE4 total net cycle law, capacity, expression scaling or beta/cAMP/PKA regulation;
- AE2;
- Na/K pump;
- K channels or CaCC;
- acid/base and CO2 closure;
- paracellular transport;
- water;
- bath/geometry;
- stimulus.

Keep the old donor-weighted routed AE4 implementation selectable as a comparator but inactive in Task 40 production runs.

## Phase B - focused algebra, source and conservation tests

Before any stationary solve, verify:

1. `J4` from the equal-routing evaluator is bitwise or numerically identical to the inherited routed AE4 `J4` for the same state and regulation input.
2. Na and K cation branches are each exactly 0.5 of the one-cation cycle.
3. `S_Na + S_K = -J4`.
4. Cl, TIC and TA sources are unchanged from the inherited routed AE4 evaluator.
5. `S_Na + S_K - S_Cl - S_TA = 0` within floating-point round-off.
6. The same 50:50 source rule works for reverse `J4 < 0` without clipping or a special fitted branch.
7. AE4 expression and beta/PKA regulation scale the complete equal-routing cycle exactly as before.
8. All non-AE4-routing scientific parameters and source laws match the Task 39 parent.

Independently reproduce the algebraic preflight at the accepted Task 31 R09 state:

- donor-weighted Na fraction `0.09062450161183984`
- `J4 = 0.004937067441546469 fmol/s`
- expected equal-routing change in the old reference-state Na source `-0.0020211144444590447 fmol/s`
- expected K-source change `+0.0020211144444590447 fmol/s`
- zero direct change in Cl/TIC/TA AE4 sources
- zero charge perturbation.

If the implementation or algebra fails after one literal coding correction, stop with:

`TASK40_EQUAL_ROUTING_SOURCE_FAILED`

Do not alter a scientific parameter to fix a source test.

## Phase C - one WT resting-state solve

Exact Task 31 REST nesting is not expected because the mechanism now changes Na/K source partition at rest.

Use the accepted Task 31 R09 resting state as the sole initial guess for exactly one intended WT stationary solve with:

- Task 39 Palk NKCC1 active;
- equal 50:50 AE4 cation routing active;
- AE4 expression 1.0;
- resting stimulus/regulatory conditions;
- every other Task 39 scientific setting frozen.

No scalar calibration and no parameter fitting are permitted.

No multistart, root search family, continuation tree, alternate routing fraction, or second WT background.

One numerical retry is allowed only if the first stationary method fails for a clearly numerical reason, with the same scientific equations and parameters.

Require an admissible WT rest with the inherited physical/conservation gates. At minimum:

- finite positive tracked amounts and volumes;
- `Na_i < 40 mM`;
- `50 <= K_i <= 200 mM`;
- `30 <= Cl_i <= 80 mM`;
- `6.6 <= pH_i <= 7.3`;
- positive finite TIC/HCO3, HCO3 below 100 mM;
- cell volume `< 3 pL`;
- current, charge, carbon, water and speciation residuals within inherited production tolerances.

If no admissible WT rest is obtained, stop with:

`TASK40_EQUAL_ROUTING_WT_REST_FAILED`

Freeze the accepted Task 40 WT rest before dynamics.

## Phase D - exactly three production trajectories

From the same frozen Task 40 WT rest, run exactly:

1. WT, AE4 expression `1.0`
2. AE4 expression `0.05`
3. AE4 expression `0.00`

Use the same 600 s standard stimulation and the same production stiff solver/tolerances as Task 39.

Do not solve genotype-specific resting states.

Only `genotype.ae4_expression` changes between the three trajectories.

Apply the same physiology and conservation gates as Task 39.

If WT fails, stop with:

`TASK40_EQUAL_ROUTING_WT_FAILED`

If WT passes but 5% fails, stop with:

`TASK40_EQUAL_ROUTING_AE4_5PCT_FAILED`

If 5% passes but null fails, stop with:

`TASK40_EQUAL_ROUTING_AE4_NULL_FAILED`

If all three pass, report:

`TASK40_EQUAL_ROUTING_VALIDATION_COMPLETE`

Do not rescue a failed scientific gate.

## Required outputs

For the new WT rest, report the full core state, flux ledger, water, voltages, residuals, and the shift relative to the accepted Task 31 R09 state.

For WT, 5%, and null trajectories report at 0,60,...,600 s:

- q_out and cumulative secretion;
- Na_i, K_i, Cl_i, pH_i, HCO3_i, TIC_i, cell volume;
- NKCC1 cycle/Cl influx;
- total AE4 Cl flux and separate AE4 Na/K cation branches;
- signed and positive AE2 Cl flux;
- NBC cycle/HCO3 influx/current;
- NHE1 flux;
- total Na/K pump cycles;
- CaCC Cl export;
- paracellular Cl return;
- membrane potentials and lumen osmolarity;
- Palk activity product `X = Na_i K_i Cl_i^2`.

Integrate relevant fluxes over 60-600 s and compare each perturbation against Task 40 WT.

Explicitly answer:

1. At 600 s, is Na_i in AE4=0.05 and/or AE4=0.00 higher than WT?
2. How much does integrated positive NKCC1 chloride loading change at 5% and null relative to WT?
3. Is this compensation smaller than Task 39's `+29.65% / +33.10%`?
4. Does Cl_i remain materially lower after AE4 loss?
5. Does AE2 remain negligible as a positive loader?
6. What happens to NBC, NHE1, pump, CaCC and paracellular return?
7. Does the cumulative secretion deficit persist, grow, plateau, or collapse over 600 s?
8. What are the final cumulative secretion ratios and deficits at 5% and null?
9. As held-out context only, how does the null deficit compare with the experimental approximately 35% reduction?

Do not use any answer to trigger another mechanism or parameter change inside Task 40.

## Absolute anti-combinatorial rules

Forbidden:

- optimisation;
- routing-fraction sweep;
- testing 22%, 29%, or any split other than 50:50;
- grid/random/global/Bayesian/evolutionary search;
- multistart stationary search;
- R10;
- alternate NKCC1 law or coefficient;
- changing `alpha_eff`;
- changing NKCC1 stimulation;
- changing NBC or NHE1;
- changing AE4 total capacity or regulation;
- changing AE2, pump, channels, CO2, paracellular or water parameters;
- alternate Ca/CCh/IPR inputs;
- genotype-specific resting solves;
- adding another transporter or mechanism;
- fitting to the 35% phenotype.

## Compute budget

- one intended WT stationary solve;
- at most one purely numerical stationary retry;
- exactly three intended production integrations if the WT rest passes;
- at most one purely numerical integration retry across the entire task;
- zero optimisation calls;
- zero parameter sweeps;
- one scientific worker;
- one BLAS thread;
- maximum 15 minutes numerical execution.

## Publication

Write compact UTF-8 outputs under:

- `analysis/40_ae4_equal_cation_routing/`
- `results/40_ae4_equal_cation_routing/`

Use the connected GitHub integration for remote writes.
Do not use shell Git authentication, PATs, SSH configuration or `gh auth`.
Do not merge or modify `main`.

Stop at the first declared scientific failure or after the complete Task 40 report.
