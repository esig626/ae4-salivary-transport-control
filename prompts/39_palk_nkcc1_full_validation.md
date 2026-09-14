# Task 39 - source-fixed Palk/Benjamin NKCC1 law and complete frozen validation

## Objective

Task 38 showed that the corrected NBC/AE4 architecture is physiologically viable, but AE4 loss is almost completely rescued because the current generic NKCC1 law increases chloride loading by about 24-27%.

This task tests one specific hypothesis only:

**Does replacing the generic NKCC1 concentration-response law with the published Palk et al. / Benjamin-Johnson two-state law suppress that compensatory response while the rest of the current model is frozen?**

This is one implementation-and-validation task. It is not a search task and it is not permission to keep changing mechanisms until the phenotype appears.

Read `analysis/39_palk_nkcc1_full_validation/source_contract.md` before touching model source.

## Scientific lineage that is fixed

Start from Task 38 final commit `d4bce14c5caed4e05587e237f002814f0f27e226`.

Keep unchanged:

- accepted Task 31 R09 resting state and all 12 core coordinates;
- Cha et al. mechanistic NHE1 carrier amount `2.339370005697548e-05 fmol`;
- Task 36 stimulus-recruited electrogenic `1 Na : 2 HCO3` NBC;
- `G_NBC = 0.11570913197464398 fmol/s`;
- NBC voltage-dependent affinity and basolateral current contribution;
- NHE1 multiplier `1 + 1.3 u_sec`, therefore 1.0 at rest and 2.3 at the standard central stimulus;
- routed AE4 `1 Cl : 1 monovalent cation : 2 HCO3` stoichiometry and donor-side Na/K routing;
- existing beta/cAMP/PKA regulation of AE4;
- AE2;
- Na/K pump;
- K channels;
- CaCC;
- CO2 and acid-base closure;
- paracellular transport;
- water laws;
- bath and geometry;
- stimulus protocol;
- the existing NKCC1 stimulus-dependent N1 activity multiplier, including its resting and fully recruited values.

The only scientific model change permitted is the **NKCC1 concentration-response core law**.

## Phase A - implement the published NKCC1 law exactly

Implement a named selectable Palk/Benjamin two-state NKCC1 evaluator. Prefer a focused module such as `src/modern_full_model/nkcc1_palk2010.py` plus the smallest selector/wiring changes needed in the existing homeostasis evaluator.

The law is fixed by the source contract:

`X = Na_i * K_i * Cl_i^2`

`F_Palk = (a1 - a2_mM * X) / (a3 + a4_mM * X)`

`J_NKCC1 = alpha_eff * activity_multiplier * F_Palk`

with

- `a1 = 157.5 s^-1`
- `a2_mM = 2.0096e-5 mM^-4 s^-1`
- `a3 = 1.0306 s^-1`
- `a4_mM = 1.3852e-6 mM^-4 s^-1`.

Positive flux is one inward `1 Na : 1 K : 2 Cl` cycle.

Do not clip negative flux. The source law is reversible through its numerator zero.

Do not add voltage dependence to NKCC1. The published law is electroneutral and concentration-driven.

Do not add a new NKCC1 state. Use the published steady two-state reduction as an algebraic flux law.

Retain the current generic thermodynamic NKCC1 evaluator as a selectable comparator, but it must not be active in Task 39 production runs.

### Exact unit verification

The primary paper prints the fourth-order coefficients in its source convention, while the historical MATLAB implementation multiplies them by `1e-12` before using mM concentration variables.

Add a focused test showing that:

- source coefficients evaluated after converting mM concentrations to M; and
- the converted coefficients `2.0096e-5` and `1.3852e-6` evaluated directly on mM concentrations

produce the same shape factor to numerical precision.

Do not invent another unit convention.

### Exact one-point WT resting scale

Do not optimise `alpha_eff`.

Recompute it algebraically from the accepted Task 31 R09 resting state and resting NKCC1 cycle flux:

- Na_i `11.636125748680639 mM`
- K_i `116.76320932871538 mM`
- Cl_i `60.30496692587399 mM`
- J_NKCC1_rest `0.1281222023866353 fmol/s`

Expected audit values are recorded in the source contract:

- `X_rest = 4941065.334966328 mM^4`
- `F_rest = 7.391062769440939`
- `alpha_eff = 0.017334746894096052 fmol/s`.

Recompute these independently. If the implementation does not reproduce them within numerical precision, stop with

`PALK_NKCC1_SOURCE_REPRODUCTION_FAILED`.

Freeze/hash the derived `alpha_eff` before any stimulated trajectory. It is the only new absolute scale and is determined from WT REST only.

## Phase B - source and conservation tests

Before dynamics, verify all of the following with focused tests:

1. Direct Eq. 17 evaluation agrees with an independent literal implementation.
2. M-space and mM-space coefficient conventions agree after the exact fourth-order conversion.
3. Positive NKCC1 cycle source is exactly `(+1 Na, +1 K, +2 Cl, 0 TIC, 0 TA)`.
4. The flux is exactly zero at `X = a1/a2_mM` within floating-point tolerance.
5. Flux changes sign on opposite sides of that reversal point.
6. The inherited NKCC1 activity multiplier scales the complete Palk cycle flux linearly but does not move reversal.
7. At REST, the inherited activity multiplier remains exactly 1.0.
8. The new Palk evaluator reproduces the accepted Task 31 resting NKCC1 cycle flux exactly to production tolerance.
9. No source term violates cell charge conservation.

If these tests fail after one literal coding correction, stop with

`PALK_NKCC1_SOURCE_OR_CONSERVATION_FAILED`.

Do not change a scientific parameter to repair a failed test.

## Phase C - exact Task 31 REST nesting

Evaluate the accepted Task 31 R09 resting state directly under:

- the frozen pre-Task-39 model; and
- the Task 39 model with Palk NKCC1 active.

Do not solve another resting root.

Because the NKCC1 source signature is unchanged and `alpha_eff` is defined to match its resting cycle flux, the complete resting RHS should nest to the existing solution.

Require the same nesting tolerances used in Task 37 for:

- all amount RHS rows;
- both volume RHS rows;
- apical, basolateral and transepithelial voltages;
- water outflow;
- regulatory RHS;
- conservation diagnostics.

If REST nesting fails, stop with

`PALK_NKCC1_REST_NESTING_FAILED`.

No root solve or parameter adjustment is permitted.

## Phase D - one WT stimulated validation

Only after source tests and REST nesting pass, run exactly one intended WT trajectory:

- initial condition: accepted Task 31 R09 WT REST state;
- AE4 expression `1.0`;
- CCh `0.3 uM`;
- isoproterenol `5 uM`;
- Ca `0.25 uM`;
- beta occupancy `1.0`;
- 0-600 s;
- same production Radau settings used in Tasks 37-38.

Apply the same physiology and conservation gates as Task 37:

- finite positive tracked amounts and volumes;
- `Na_i < 40 mM`;
- `50 <= K_i <= 200 mM`;
- `30 <= Cl_i <= 80 mM`;
- `6.6 <= pH_i <= 7.3`;
- positive finite TIC/HCO3 with HCO3 below 100 mM;
- cell volume `< 3 pL`;
- all inherited conservation/current/water/speciation tolerances.

Use the same secretion activation gates as Task 37.

If WT fails a scientific gate, stop with

`PALK_NKCC1_WT_VALIDATION_FAILED`.

Do not rescue it.

### WT report

For WT report:

- cumulative 0-600 s secretion;
- mean q_out over 60-600 s;
- q_out at 600 s;
- q_out at 0,60,...,600 s;
- Na_i, K_i, Cl_i, pH_i, HCO3_i, TIC_i, cell volume at the same times;
- NKCC1 cycle and chloride-loading fluxes;
- AE4 and signed AE2 chloride fluxes;
- NBC cycle/HCO3 influx and current;
- NHE1 flux;
- total Na/K pump cycles;
- CaCC chloride export;
- paracellular chloride return;
- voltages and lumen osmolarity;
- the Palk intracellular activity product `X = Na_i K_i Cl_i^2`.

Integrate positive NKCC1, AE4 and AE2 chloride loading over 60-600 s and report their shares. Compare the NKCC1 share with the 65-75% experimental context, but do not tune anything if it lies outside that band.

Also report the ratio of Task 39 WT cumulative secretion to the frozen Task 37 WT cumulative secretion `0.986995055625111 pL`. This is a diagnostic consequence of replacing NKCC1 kinetics, not a fitting target.

## Phase E - frozen AE4 expression predictions

Only if WT passes, freeze every Task 39 scientific parameter and run exactly two additional intended trajectories, in this order:

1. AE4 expression `0.05`;
2. AE4 expression `0.00`.

For both cases:

- use exactly the same accepted Task 31 WT resting initial state;
- use exactly the same stimulus and solver as WT;
- change only `genotype.ae4_expression`;
- do not solve genotype-specific resting states;
- do not refit `alpha_eff` or anything else.

If the 5% case fails physiology/conservation, stop before the null case with

`PALK_NKCC1_AE4_5PCT_FAILED`.

If the null case fails, report

`PALK_NKCC1_AE4_NULL_FAILED`.

Otherwise finish with

`PALK_NKCC1_FULL_VALIDATION_COMPLETE`.

## Required comparison and mechanism diagnosis

For WT, 5% and null, report the same time-series and integrated flux panel.

Then compare each perturbation against Task 39 WT:

- cumulative secretion ratio and deficit;
- mean stimulated flow ratio;
- endpoint flow ratio;
- intracellular Cl difference;
- Na and K differences;
- pH difference;
- integrated positive NKCC1 chloride-loading ratio;
- integrated AE4 loading;
- signed and positive AE2 loading;
- NBC and NHE1 integrated flux ratios;
- pump-cycle ratio;
- CaCC export and paracellular return ratios;
- cumulative secretion difference versus time.

Explicitly answer these questions:

1. Does the Palk law reduce the large NKCC1 compensation seen in Task 38 (`+24.13%` at AE4=5%, `+26.89%` at null)?
2. Does AE4 loss now cause intracellular Na and/or K to rise relative to WT, as in the mechanism described by the 2018 paper, or do they still fall?
3. Does AE2 remain a negligible positive chloride loader?
4. Does intracellular chloride fall substantially after AE4 loss?
5. Does the secretion deficit grow during sustained stimulation or collapse through compensation?
6. What secretion deficit is predicted at 5% and at null?
7. How does the null result compare, as held-out validation only, with the experimental approximately 35% deficit and the original 2018 model's approximately 24% deficit?

Do not use answers to these questions to make another model change inside Task 39.

## Important interpretation rule

This Task 39 model is **not** claimed to be an exact reproduction of the entire 2018 model. It places the published Palk/Benjamin NKCC1 concentration-response core inside the current Task 38 architecture while preserving the current NKCC1 stimulation multiplier and all other current mechanisms.

Any success or failure must be attributed accordingly.

## Absolute anti-combinatorial rules

There is no parameter search in Task 39.

Forbidden:

- any optimisation call;
- any stationary/root solve;
- any grid or sweep;
- one-at-a-time parameter exploration;
- pairwise/factorial search;
- random search;
- Latin hypercube;
- multistart;
- global/evolutionary/Bayesian optimisation;
- Shapley analysis;
- alternate Palk parameter values;
- alternate NKCC1 density/prefactor values;
- fitting `alpha_eff` to stimulation, genotype or phenotype;
- modifying the inherited NKCC1 stimulation multiplier;
- trying the Palk law both with and without the stimulation multiplier;
- reverting to or mixing with the generic tanh law in production runs;
- changing NBC capacity, stoichiometry or activation;
- changing NHE1 gain or kinetics;
- changing AE4 capacity, routing, stoichiometry or regulation;
- changing AE2, pump, channels, CO2, paracellular or water parameters;
- changing Ca/CCh/IPR inputs;
- R10;
- AE2 knockout;
- any AE4 expression other than 1.0, 0.05 and 0.00;
- automatic addition of another mechanism after failure.

If the source-fixed law fails, report the failure and stop.

## Hard compute budget

- zero stationary solves;
- zero optimisation calls;
- exactly three intended production integrations: WT, 5%, null;
- maximum one numerical retry across the entire task, only for a clearly numerical solver failure and only using an already-declared production stiff solver;
- maximum four integration attempts total;
- one scientific worker;
- one BLAS thread;
- maximum 15 minutes numerical execution.

A spare retry is not permission to test another scientific condition.

## Output

Write compact UTF-8 results under

- `analysis/39_palk_nkcc1_full_validation/`
- `results/39_palk_nkcc1_full_validation/`

At minimum publish:

- source/equation verification;
- unit-conversion test evidence;
- derived `alpha_eff` audit;
- REST nesting audit;
- WT verification and time series;
- 5% verification and time series if reached;
- null verification and time series if reached;
- integrated flux comparison;
- NKCC1 compensation comparison against Task 38;
- final concise scientific interpretation;
- numerical budget.

Do not publish binary trajectory archives.

## GitHub publication

Use the connected GitHub integration for remote writes.
Do not depend on shell Git authentication, PATs, SSH setup or `gh auth`.
Shell Git is for local status/diff inspection only.
Do not merge or modify `main`.
Publication problems must not trigger scientific recomputation.

Stop when Task 39 is complete or the first scientific stop condition is reached.
