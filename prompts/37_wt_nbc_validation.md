# Task 37: WT-only validation of the corrected electrogenic NBC architecture

## Purpose

Validate the Task 36 corrected whole-cell architecture before any AE4 perturbation is allowed.

This is **not** a calibration task, **not** a search task, and **not** a genotype task.

The only scientific question is:

> With the accepted Task 31 R09 resting state, the fixed Task 36 stimulus-recruited electrogenic 1 Na : 2 HCO3 NBC surrogate, the fixed Cha et al. NHE1 model, and all other inherited mechanisms unchanged, does one standard WT CCh + IPR stimulation produce a physically admissible secreting trajectory with a plausible stimulated chloride-loading partition?

If the answer is no, stop and report the first failing balance/gate. Do not rescue the model inside this task.

---

## Fixed starting point

Use the Task 36 branch state inherited from commit:

`3ae2b2d5cb28e51f34d020040e8f2d7c70af5911`

Use the accepted Task 31 R09 WT resting state and mechanistic NHE1 carrier amount:

`2.339370005697548e-05 fmol`

The accepted Task 31 R09 resting outflow is:

`q_out_rest = 0.0010757073853493032 pL/s`

Do not solve a new resting state.

Do not use Task 33 optimised parameter vectors.

---

## Fixed Task 36 mechanism

The active NBC surrogate is:

`Na_o + 2 HCO3_o <-> Na_i + 2 HCO3_i`

Positive inward cycle source signature in conserved cell coordinates:

`(+1 Na, 0 K, 0 Cl, +2 TIC, +2 TA)`

It carries net `-1` charge-equivalent into the cell per inward cycle and therefore contributes explicitly to basolateral current closure.

The inward affinity is:

`A_NBC = log((Na_o HCO3_o^2)/(Na_i HCO3_i^2)) + V_b/(RT/F)`

The flux law is:

`J_NBC = u_sec * G_NBC * tanh(A_NBC/2)`

with:

`u_sec = clip((Ca - 0.058)/(0.25 - 0.058), 0, 1)`

The fixed structural NBC capacity is:

`G_NBC = 0.11570913197464398 fmol/s`

The same secretory coordinate applies the already-declared stimulated NHE1 multiplier:

`G_NHE = 1 + 1.3 u_sec`

so NHE1 is exactly `1.0x` at REST and `2.3x` at the standard central stimulated condition.

Do not alter any of these equations, constants, stoichiometries or stimulus mappings.

---

## Phase A: focused implementation verification

Before any integration, run only the focused Task 36 NBC tests plus any directly required Task 31 mechanistic-NHE regression tests needed to establish nesting.

Verify explicitly:

1. NBC source stoichiometry is exactly `(+1 Na, +2 TIC, +2 TA)` per inward cycle and zero direct K/Cl source.
2. NBC transported charge is exactly `-1` charge-equivalent into the cell per inward cycle.
3. The NBC current sign used in basolateral current closure is consistent with that transported charge.
4. The 1:2 NBC affinity includes the basolateral voltage term with the correct sign.
5. The reversible flux changes sign when the full electrochemical affinity changes sign.
6. At `Ca=0.058 uM`, NBC activation is exactly zero and NHE1 multiplier is exactly one.
7. At `Ca=0.25 uM`, NBC activation is exactly one and NHE1 multiplier is exactly 2.3.
8. Charge accounting, carbon accounting and TIC/TA source bookkeeping close to numerical round-off for the isolated NBC contribution.
9. At the accepted Task 31 R09 WT resting state, the Task 36 wrapper reproduces the Task 31 resting model without a new solve.

For REST nesting, compare the accepted state under the Task 31 model and under the Task 36 NBC wrapper at the same rest stimulus. Require:

- same state vector, because no solve is permitted;
- NBC cycle flux exactly zero;
- NHE1 multiplier exactly one;
- max absolute difference in amount RHS rows <= `1e-12 fmol/s`;
- max absolute difference in volume RHS rows <= `1e-14 pL/s`;
- max absolute difference in membrane potentials <= `1e-12 V`;
- absolute difference in resting `q_out` <= `1e-12 pL/s`.

If any focused test or REST nesting check fails, make only the smallest literal implementation correction needed to match the already-declared Task 36 equations. Do not change a scientific parameter, target, stoichiometry, capacity, stimulus, transporter law or model family.

After one such correction, rerun only the focused tests once.

If the focused tests or REST nesting still fail, STOP and report:

`TASK37_IMPLEMENTATION_OR_REST_NESTING_FAILED`

Do not integrate WT.

---

## Phase B: exactly one WT stimulation

Only if Phase A passes, run one intended WT trajectory.

Initial condition:

- accepted Task 31 R09 WT resting state

Genotype:

- WT only, AE4 expression `1.00`

Stimulus:

- CCh `0.3 uM`
- isoproterenol `5 uM`
- Ca `0.25 uM`
- existing standard beta/IPR input
- onset `t=0`
- duration `0--600 s`

Use the existing production stiff solver settings. Prefer the declared production Radau solver.

No stationary solve is allowed before the trajectory.

### Numerical retry rule

Allow at most one retry, and only if the first WT integration fails for a clearly numerical solver reason.

The retry may change only the numerical integration method/settings already declared as production alternatives, for example Radau to the existing production BDF configuration.

It may not change any scientific parameter or input.

---

## WT physical admissibility gates

The complete accepted trajectory must have:

- finite states and observables throughout;
- all tracked amounts and volumes positive;
- `Na_i < 40 mM`;
- `50 <= K_i <= 200 mM`;
- `30 <= Cl_i <= 80 mM`;
- `6.6 <= pH_i <= 7.3`;
- positive intracellular TIC and HCO3;
- `HCO3_i < 100 mM`;
- cell volume `< 3 pL`;
- existing current, charge, carbon and water conservation diagnostics within their declared production tolerances.

If the WT trajectory violates a physical or conservation gate, STOP and report:

`WT_NBC_TRAJECTORY_PHYSIOLOGY_FAILED`

Do not change anything and do not start another mechanism.

---

## WT secretion gate

Using

`q_out_rest = 0.0010757073853493032 pL/s`,

require:

- mean `q_out` over `60--600 s` > `1.10 * q_out_rest`;
- cumulative outflow over `0--600 s` > `1.10 * (600 * q_out_rest)`;
- `q_out(600 s) >= q_out_rest`.

If any of these fail, STOP and report:

`WT_NBC_SECRETION_GATE_FAILED`

Do not change any parameter.

---

## Required WT readout

For the successful WT trajectory, report at `t = 0, 60, 120, 180, 240, 300, 360, 420, 480, 540, 600 s`:

- `q_out`;
- intracellular Na, K, Cl;
- intracellular TIC and HCO3;
- intracellular pH;
- cell volume;
- lumen osmolarity;
- apical, basolateral and transepithelial voltages;
- NBC activation fraction;
- NBC cycle flux;
- NBC current;
- NHE1 multiplier;
- NHE1 flux;
- NKCC1 chloride influx;
- AE4 chloride influx;
- signed AE2 chloride flux;
- CaCC chloride export;
- paracellular chloride return.

Integrate over `60--600 s`:

- positive NKCC1 chloride loading;
- positive AE4 chloride loading;
- positive AE2 chloride loading;
- signed AE2 chloride flux;
- NBC cycles;
- NBC bicarbonate-equivalent influx;
- NHE1 flux;
- CaCC chloride export;
- paracellular chloride return;
- water outflow.

Define the positive basolateral chloride-loading pool as the sum of only positive NKCC1, AE4 and AE2 chloride-loading contributions.

Report each positive share.

Compare the NKCC1 share with the experimental context band `65--75%`.

This band is **validation only**. Do not tune anything to enter it.

Also report whether AE4 is a substantial positive stimulated chloride loader rather than the ~2% resting contribution previously seen.

---

## Interpretation rules

Do not fit to any knockout phenotype.

Do not infer that the model explains AE4 knockout yet. Task 37 is WT only.

Classify the outcome as one of:

- `WT_NBC_VALIDATION_PASS`
- `WT_NBC_VALIDATION_PASS_PARTITION_OUTSIDE_CONTEXT`
- one of the explicit failure labels above.

A pass requires physical admissibility, conservation closure and the secretion gate.

The chloride partition is reported separately as validation context.

---

## Absolute anti-combinatorial rules

There is **no search** in Task 37.

Forbidden:

- any optimisation;
- any stationary root solve;
- any parameter sweep;
- any NBC capacity change;
- any alternative NBC stoichiometry;
- any alternative NBC law;
- any alternative NHE1 multiplier;
- any alternative NHE1 model;
- any change to AE4 capacity, routing, stoichiometry or regulation;
- any change to NKCC1 capacity or stimulation;
- any change to AE2;
- any change to pump, K, CaCC, CO2, paracellular or water parameters;
- any change to Ca, CCh or isoproterenol;
- any alternative beta/PKA family;
- any genotype perturbation;
- AE4 `5%`;
- AE4 `0%`;
- AE2 knockout;
- R10;
- grids;
- one-at-a-time sweeps;
- pairwise/factorial searches;
- random search;
- Latin hypercube;
- multistart;
- global/evolutionary/Bayesian optimisation;
- Shapley analysis;
- automatic escalation to another mechanism after failure.

If WT fails, report the failure and stop.

---

## Hard compute budget

- zero stationary solves;
- zero optimisation calls;
- exactly one intended WT dynamic integration;
- at most one numerical WT retry for solver failure only;
- maximum two integration attempts total;
- one worker;
- one BLAS thread;
- maximum 10 minutes numerical execution;
- no genotype integrations.

STOP immediately when an applicable limit is reached.

---

## Required outputs

Publish compact UTF-8 outputs under:

`analysis/37_wt_nbc_validation/`

and

`results/37_wt_nbc_validation/`

At minimum:

- `focused_test.log`
- `rest_nesting.json`
- `wt_timeseries.csv`
- `wt_integrated_fluxes.csv`
- `wt_chloride_partition.json`
- `wt_verification.json`
- `budget.json`
- `final_answer.md`

Do not publish binary trajectories or archives.

Use the connected GitHub integration for publication.

Do not depend on shell Git authentication.

Do not merge or modify `main`.

The task ends after the WT report. Do not run AE4 5% or 0% in Task 37.
