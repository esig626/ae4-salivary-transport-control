# Task 33: fixed-water WT optimisation with NKCC1 at 70% and frozen mechanistic NHE1

Repository: `esig626/ae4-salivary-transport-control`
Branch: `codex/task-33-fixed-water-nkcc70-ae4-optimisation`

## Purpose

Stop re-diagnosing already established defects. Build one WT parameterisation in which the existing routed AE4 mechanism has a substantial productive chloride-loading role, while preserving the mechanistic NHE1 repair, the accepted WT resting water throughput, and a genuinely calcium-activated secretion response. Then freeze the WT parameterisation and test only AE4 expression 0.05 and 0.0.

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

## WT physiological guardrails and fixed resting water target

Use R09 only for optimisation.

The accepted Task 31 R09 WT state is the baseline physiological reference. Intracellular chloride is **not** forced to the 2015 central value 50.1 mM, but it is also not unconstrained.

Use this hard WT chloride guardrail:

`45 mM <= Cl_i <= 65 mM`.

This range deliberately contains both the measured WT value around 50 mM and the accepted Task 31 WT state around 60 mM. It is a physiological guardrail, not a point-fitting target. Any optimiser candidate outside this range is infeasible and must not be rescued by changing the range.

Also require:

- `6.84 <= pH_i <= 6.98`;
- finite positive concentrations;
- ordinary intracellular sodium, with `Na_i <= 30 mM`;
- cell volume `< 3 pL`;
- all existing charge, current, carbon, conservation and water closure gates.

Hold the WT **resting** ductal/luminal outflow target fixed to the accepted Task 31 value:

`q_out_rest_WT = 0.0010757073853493032 pL/s`.

Keep all water-law parameters fixed. The optimiser must meet this resting throughput through ionic/osmotic balance, not by changing water permeability or the outflow law. Report the achieved absolute and relative error.

Do not create a new water target from genotype results.

## Mandatory WT calcium-activated secretion gate

A stationary WT fit is **not sufficient**. Before the WT parameterisation can be frozen, the same candidate must show that the standard central calcium stimulus actually activates secretion.

For a WT candidate that has already passed all stationary constraints, run the existing standard central stimulation at:

- `Ca = 0.25 uM`;
- `t = 0--600 s`;
- all existing stimulus/regulatory equations unchanged.

Use the fixed resting outflow above as the no-stimulus reference. No separate no-stimulus integration is needed because the accepted stationary state defines the resting reference.

The WT candidate passes the dynamic activation gate only if all of the following hold:

1. the stimulated integration completes numerically and remains physically admissible;
2. time-mean `q_out` over `60--600 s` is at least `1.10 * q_out_rest_WT`;
3. cumulative ductal secretion over `0--600 s` is at least `1.10 * (600 s * q_out_rest_WT)`;
4. the stimulated response is not merely a single numerical spike followed by sustained secretion below rest: `q_out(600 s) >= q_out_rest_WT`;
5. concentrations and compartment volumes remain finite and positive.

The 10% activation margin is a modelling guardrail to exclude a trivially non-responsive or numerically indistinguishable secretory system. It is **not** fitted to the AE4-null phenotype and is not claimed as a measured fold-change.

Do not alter calcium amplitude, stimulus equations, water parameters, or the 10% gate after seeing results.

Only stationary-feasible WT candidates may be dynamically tested. Maximum WT dynamic gate evaluations are specified below. If no candidate passes within the predeclared budget, Task 33 stops.

## NKCC1 70% WT chloride-loading constraint

At the final WT resting state, require NKCC1 to supply 70% of the **positive basolateral chloride-loading pool**:

`share_NKCC1 = J_NKCC1_Cl_positive / J_positive_total = 0.70`.

Use `J_NKCC1_Cl = 2 * J_NKCC1_cycle`.

Define the positive pool as the sum of positive chloride-loading contributions actually present in the model:

- `max(J_NKCC1_Cl, 0)`;
- `max(J_AE4_Cl, 0)`;
- `max(J_AE2_Cl, 0)`;
- any other already-implemented positive basolateral chloride loader, if one exists.

Signed counterfluxes must be reported separately and must not be hidden inside the denominator.

Treat `0.70` as a constrained WT target with a numerical tolerance of `+/- 0.01`. Do not use the AE4-null phenotype to set or relax this constraint.

Do **not** separately impose an arbitrary AE4 percentage. AE4's final positive loading share is an output. However, reject a solution in which AE4 is non-positive or negligible (<5% of the positive loading pool), because such a solution does not satisfy the purpose of this task.

## Exactly one bounded WT optimisation

This task authorises one bounded low-dimensional local WT optimisation. It does **not** authorise a search programme.

Before any numerical optimisation, write and checkpoint:

`analysis/33_fixed_water_nkcc70_ae4_optimisation/optimisation_contract.md`

The contract must declare:

- the exact free parameters;
- fixed numerical bounds for every free parameter;
- parameter transformations;
- objective/constraint definitions and weights;
- the single optimisation algorithm;
- the single allowed deterministic restart rule;
- all numerical counters and stop conditions;
- the stationary-feasibility gate before dynamics;
- the WT calcium-activation gate exactly as written above.

Once written, **the contract is immutable for the rest of Task 33**. Do not change bounds, free parameters, objective weights, algorithm, chloride range, water target, NKCC1 target, calcium input or dynamic activation gate after seeing results.

Use the smallest defensible free set needed. The only parameters that may be considered are these existing uncertain magnitude parameters:

1. NKCC1 capacity;
2. routed AE4 carrier amount/activity magnitude;
3. AE2 capacity;
4. total Na/K-pump capacity;
5. total K conductance magnitude;
6. apical CaCC conductance magnitude.

Keep all membrane allocation fractions fixed. Do not free NHE1 density. Do not free water, CO2, paracellular, stimulus, regulatory or AE4 cation-partition parameters.

Prefer algebraic elimination. If a parameter is unnecessary, leave it frozen. Do not run alternative optimisations over different subsets.

### Parameter magnitude guardrails

Use Task 31 R09 values as the reference magnitudes. Unless a tighter source-based range is already encoded in the repository, enforce these hard multiplicative bounds:

- NKCC1 capacity: `0.25x` to `2x` Task 31;
- AE4 activity/carrier magnitude: `0.1x` to `50x` Task 31;
- AE2 capacity: `0.25x` to `4x` Task 31;
- total Na/K-pump capacity: `0.5x` to `2x` Task 31;
- total K conductance magnitude: `0.5x` to `2x` Task 31;
- apical CaCC conductance magnitude: `0.5x` to `2x` Task 31.

Do not widen these bounds during the task. A feasible solution sitting on a bound must be reported explicitly as boundary-dependent.

### WT constraints/objective

Treat the following as hard feasibility constraints before any soft objective:

1. full stationary amount/current/charge/carbon/water closure;
2. fixed WT resting `q_out` target;
3. NKCC1 positive chloride-loading share `0.70 +/- 0.01`;
4. `45 <= Cl_i <= 65 mM`;
5. `6.84 <= pH_i <= 6.98`;
6. positive, non-negligible AE4 loading;
7. `Na_i <= 30 mM`, positive concentrations and cell volume `<3 pL`;
8. after stationary feasibility is established, the mandatory WT Ca = 0.25 uM dynamic activation gate above.

Among WT candidates that pass both stationary and dynamic gates, minimise weighted log-deviation of free parameters from Task 31 and state deviation from the accepted Task 31 R09 WT state.

Do not target knockout chloride, knockout pH, knockout secretion, or a desired 5%/WT or 0%/WT secretion ratio.

## Absolute anti-combinatorial stop rules

These are hard stops, not suggestions.

- Maximum **one** local constrained optimisation run.
- Maximum **one** deterministic restart, and only from that run's own best feasible/least-infeasible point.
- Maximum **20 distinct WT parameter vectors** evaluated in total, including the initial vector and restart evaluations.
- Maximum **3 WT dynamic activation integrations**, and only for stationary-feasible candidate vectors.
- Never enumerate a Cartesian product of parameter values.
- Never perform one-at-a-time factor sweeps across all parameters.
- Never try multiple optimiser families.
- Never loop over alternative free-parameter subsets.
- Never widen parameter bounds, chloride range, water target, NKCC1 tolerance or calcium activation threshold after a failure.
- Never launch a second-stage search because a genotype result is disappointing.
- No random search, Latin hypercube, multistart cloud, global/evolutionary/Bayesian optimisation, pairwise scan, factorial design, Shapley analysis or automatic singles-to-pairs escalation.

Numerical ceilings across the whole task:

- at most 30 stationary solver calls total;
- at most 40,000 stationary residual evaluations total;
- at most 20 WT parameter-vector evaluations total;
- at most 3 pre-freeze WT dynamic integrations total;
- at most 20 minutes numerical execution for WT optimisation plus its dynamic gate checks;
- one worker and one BLAS thread.

If no WT solution has passed **both stationary and calcium-activation gates** when any one of those ceilings is reached, STOP TASK 33 and report `WT_OPTIMISATION_FAILED_WITHIN_PREDECLARED_BOUNDS`. Do not diagnose another mechanism, alter the contract or keep searching.

## Freeze before genotype release

When one WT parameterisation passes all stationary hard constraints **and** the WT calcium-activation gate:

1. write all free parameter values and the complete WT state/flux ledger;
2. record achieved resting `q_out`, chloride, pH and NKCC1-share residuals;
3. record WT dynamic mean q_out, q_out(600 s), cumulative 600 s secretion and all activation ratios relative to rest;
4. record AE4 and AE2 signed chloride contributions and AE4 positive-pool share;
5. hash the frozen parameter payload;
6. checkpoint/publish the frozen WT text artifacts through the connected GitHub integration.

Only after that freeze may AE4 expression be changed.

The already accepted pre-freeze WT dynamic gate trajectory becomes the frozen WT secretion denominator. Do not rerun it merely to create another copy.

## Genotype tests: exactly 5% and 0%

With every non-AE4 parameter frozen, test exactly:

- AE4 expression `0.05`;
- AE4 expression `0.00`.

For each genotype:

1. solve its own REST state from the frozen WT model;
2. report admissibility, pH, Na, K, Cl, TIC/HCO3, cell/lumen volumes, voltages, osmolarities and major transporter/water fluxes;
3. if and only if REST is admissible, run the same standard central stimulation at Ca = 0.25 uM to 600 s;
4. report cumulative ductal secretion and ratio to the frozen WT trajectory;
5. decompose replacement of lost AE4 chloride loading into NKCC1, AE2, altered apical delivery, paracellular return, ionic storage and water/osmotic effects.

Do not test 10%, 20%, 50%, AE2 loss, another Ca level, R10, another background or another genotype.

Do not refit after seeing either genotype.

If 5% fails to produce an admissible REST, report that result and continue only to the predeclared 0% REST test. If 0% also fails, stop. Do not create intermediate AE4 levels to bridge the solver.

Maximum post-freeze genotype dynamic integrations: at most one 5% and one 0% trajectory, for **2 additional integrations total**. No solver-method trajectory duplications.

## Interpretation

The question is direct:

> With WT resting water throughput fixed, mechanistic NHE1 frozen, NKCC1 constrained to 70% of positive WT chloride loading, WT chloride kept within 45--65 mM, and Ca = 0.25 uM required to produce a real stimulated secretion response, can the existing AE4 architecture carry a meaningful share of productive chloride loading and does reducing AE4 to 5% or 0% reduce secretion without genotype fitting?

Report the result exactly as obtained, including neutral, wrong-direction, boundary-dependent or non-closing outcomes. A disappointing outcome does not license another search inside this task.

## Outputs

Write compact UTF-8 outputs under:

- `analysis/33_fixed_water_nkcc70_ae4_optimisation/`
- `results/33_fixed_water_nkcc70_ae4_optimisation/`

At minimum include:

- `optimisation_contract.md`;
- `budget.json`;
- `frozen_wt_parameters.json` if WT succeeds;
- `wt_state.json`;
- `wt_flux_ledger.csv`;
- `wt_dynamic_activation.json`;
- `genotype_rest_states.csv` if WT succeeds;
- `genotype_flux_ledger.csv` if WT succeeds;
- `secretion_comparison.csv` if genotype dynamics run;
- `verification.json`;
- `final_answer.md`.

No binary trajectory arrays, archives or giant search tables.

## GitHub publication

Use the connected GitHub integration for all remote writes. Do not depend on shell Git authentication, `gh auth`, PATs or SSH credentials. Shell Git is only for local status/diff inspection. Publication failure must never cause scientific recomputation.

Do not merge to `main`.
