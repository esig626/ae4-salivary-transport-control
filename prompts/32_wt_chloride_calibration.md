# Task 32: WT chloride calibration with frozen mechanistic NHE1

Repository: `esig626/ae4-salivary-transport-control`
Branch: `codex/task-32-wt-chloride-calibration`
Starting scientific state: Task 31 commit `fb4bc7d9a569051f5e41cf84837046cba6952d4f`.

Read `AGENTS.md` first. This is an implementation/calibration task, not another transporter diagnosis.

## Central objective

Task 31 replaced the pathological generic NHE1 law with the Cha et al. eight-state/Mod2 mechanism. That repair removed the extreme intracellular sodium requirement, but the accepted R09 WT state still had intracellular chloride about 60.3 mM, well above the measured WT value 50.10 +/- 1.50 mM, and AE4 removal still did not produce an accepted exact-null REST.

The next step is to give the existing routed AE4 pathway the WT chloride-loading magnitude required by WT physiology, while keeping the Task 31 NHE1 mechanism and all other transporter laws unchanged.

Do NOT fit knockout chloride, knockout pH or secretion. The only calibration data in this task are WT R09 resting pH and WT R09 resting chloride.

## Frozen mechanism

Use the Task 31 Cha NHE1 implementation exactly as published in the branch:

- eight-state exchange cycle;
- Mod2 intracellular proton modifier;
- Table S1 kinetic constants unchanged;
- n_H = 3, m_H = 1;
- same law in WT and AE4-null;
- zero carbon source and zero net charge source;
- 1 Na inward per 1 H outward.

Do not introduce another NHE1 family and do not alter the Cha kinetic constants. The effective NHE1 carrier amount remains an uncertain WT abundance parameter and may be adjusted only as specified below.

Keep the existing routed AE4 topology and stoichiometry unchanged. Do not introduce opposing Na/K cycles, a new AE4 mechanism, signalling, hidden acid/base terms or genotype-specific regulation.

## Primary background

Use R09 as the calibration background because Task 31 obtained an admissible WT state there with physiological pH, Na and volume. Treat R10 only as an out-of-sample stress check after the R09 WT/null pair has been frozen and accepted.

Task 31 R09 WT values are a starting state, not targets. Reuse cached accepted state information rather than recomputing it unnecessarily.

## WT-only observables

Use the 2015 WT resting observations:

- pH_i = 6.91 +/- 0.07;
- Cl_i = 50.10 +/- 1.50 mM.

The purpose is not to force the exact central values if the final state is already inside both reported intervals. These are WT calibration/acceptance data only.

No AE4-null data may be accessed during calibration beyond knowing that the later intervention is exact AE4 loss.

## Exactly two uncertain abundance/activity scalars

Only these two scalar quantities may be altered:

1. **NHE1 effective carrier amount** used by the Task 31 Cha implementation.
2. **WT AE4 activity/capacity scalar** multiplying the existing routed AE4 rate law.

The AE4 scalar must be distinct from genotype expression. WT genotype expression remains 1.0; exact-null expression remains 0.0. Prefer scaling an existing multiplicative AE4 activity/capacity coefficient. If no clean scalar exists, add one explicit positive multiplicative activity field with units/documentation and tests. Do not change AE4 stoichiometry, routing fractions or state dependence.

All NKCC1, AE2, pump, K, CaCC, CO2, paracellular, water, bath, hydraulic and regulatory parameters remain frozen.

## Calibration procedure: bounded and sequential, not a search

No multidimensional optimiser, grid, pair sweep or random search is allowed.

Use this fixed sequence only:

1. Start from the Task 31 frozen R09 NHE1 carrier amount and accepted R09 WT state.
2. With NHE1 carrier amount fixed, use one bracketed scalar AE4-activity solve to bring stationary R09 WT Cl_i into 48.6-51.6 mM. Prefer the centre 50.10 mM only as the scalar root target.
3. Freeze that AE4 activity temporarily. If stationary R09 WT pH_i is outside 6.84-6.98, use one bracketed scalar NHE1-carrier solve to return pH_i to 6.91. If pH already lies in range, do not recalibrate NHE1.
4. Re-evaluate WT chloride once. If it left the WT chloride interval after the NHE1 update, permit ONE final bracketed AE4 scalar correction to 50.10 mM with NHE1 fixed. Do not iterate further.
5. Accept the resulting WT calibration only if full stationary closure passes and pH_i and Cl_i are both inside the reported WT intervals, Na_i <= 30 mM, cell volume < 3 pL, concentrations are positive and the existing charge/current/conservation gates pass.

If a required scalar root cannot be bracketed with finite positive activity, or the final state misses either WT interval after the permitted sequence, stop. Do not widen into another parameter family.

Do not impose an arbitrary AE4 percentage such as 10%, 20% or 30%. The resulting AE4 chloride-loading share is an output of the WT calibration, not an input.

Record the final AE4 positive chloride-loading share relative to:

- positive NKCC1 chloride loading;
- positive AE2 chloride loading;
- total positive basolateral chloride loading;
- net basolateral chloride loading.

Report signed AE2 and other counterfluxes separately.

## Freeze before genotype access

Once the R09 WT calibration passes:

- write and hash a compact frozen-parameter record;
- freeze NHE1 carrier amount and AE4 WT activity scalar;
- do not revise either after any AE4-null result is seen.

Then set AE4 expression exactly to 0.0 with everything else frozen.

## Exact AE4-null REST

Solve R09 exact-null REST from the accepted calibrated WT state using bounded continuation/direct stationary methods already present in the repository. This is a stationary solve, not genotype fitting.

Acceptance requires full stationary closure and the normal numerical/conservation/current gates, plus:

- finite positive concentrations;
- Na_i <= 30 mM;
- cell volume < 3 pL.

After the null result is frozen, compare against the experimental context:

- pH_i 6.89 +/- 0.02;
- Cl_i 36.50 +/- 1.60 mM.

These knockout measurements are validation only. Do not tune to them.

Report whether loss of AE4 is accommodated by NHE1, AE2, NKCC1, carbon storage, cation changes, paracellular return or water/osmotic changes. Do not call an unclosed candidate a phenotype.

## Dynamics

Only if both calibrated R09 WT and exact-null REST are admissible, run the central CCh + IPR / Ca = 0.25 uM 0-600 s matched trajectories.

Report:

- ductal cumulative secretion and null/WT ratio;
- intracellular chloride and pH trajectories;
- integrated positive AE4, NKCC1 and AE2 chloride transport;
- apical CaCC chloride export;
- paracellular chloride return;
- lumen osmolarity and q_out.

The 30-35% experimental secretion deficit is validation context only after results are frozen. Do not tune to it.

## R10 stress check

Only after the R09 calibration and exact-null evaluation are complete, and only if numerical budget remains:

- apply the frozen R09 NHE1 carrier amount and AE4 activity scalar to R10 with no recalibration;
- solve R10 WT REST;
- if and only if R10 WT is admissible, solve R10 exact-null REST;
- no R10 parameter adjustment is permitted.

R10 failure is an out-of-sample limitation, not permission to restart calibration.

## Numerical budget

Absolute limits across the whole task:

- at most 18 stationary solver calls;
- at most 30,000 stationary residual evaluations;
- at most 6 AE4 scalar evaluations;
- at most 4 NHE1 carrier scalar evaluations;
- at most 2 stimulated integrations total;
- one scientific worker;
- one BLAS thread;
- at most 30 minutes numerical execution.

Persist counters after each case. Stop at the first exhausted limit. Do not reset counters or rerun cached successful cases.

## Prohibited

- no combinatorial or factorial search;
- no parameter grid;
- no random/global/Bayesian/evolutionary optimisation;
- no automatic singles-to-pairs escalation;
- no knockout fitting;
- no secretion fitting;
- no new NHE1 family;
- no new transporter;
- no AE4 stoichiometry/topology change;
- no changes to NKCC1, AE2, pump, K, CaCC, CO2, paracellular or water laws;
- no arbitrary pH clamp, bicarbonate sink or acid source;
- no stimulated run before an admissible matched R09 WT/null resting pair exists.

## Required outputs

Write compact UTF-8 outputs under:

- `results/32_wt_chloride_calibration/`
- `analysis/32_wt_chloride_calibration/`

Include at minimum:

- `contract.json`;
- `budget.json`;
- `wt_calibration.json`;
- `frozen_parameters.json`;
- `resting_states.csv`;
- `flux_ledger.csv`;
- `chloride_allocation.csv`;
- `verification.json`;
- `final_answer.md`;
- trajectory summary text/CSV only if dynamics run.

The final answer must directly state:

1. Did WT-only calibration recover R09 pH and chloride simultaneously?
2. What NHE1 carrier amount and AE4 activity scalar were required?
3. What fraction of WT positive chloride loading is carried by AE4 after calibration?
4. Does exact AE4-null REST close with ordinary Na and volume?
5. Without tuning, where do null pH and chloride land relative to experiment?
6. If dynamics ran, what is the secretion direction and magnitude?
7. Does R10 survive as an out-of-sample check?
8. Which conclusion remains unproven?

## GitHub publication

Use the connected GitHub integration for all remote writes. Do not depend on shell Git authentication, `gh auth`, PATs or SSH credentials.

Publish completed source edits, tests and compact text results directly to `codex/task-32-wt-chloride-calibration`. Do not merge to `main`. Publication failure must never trigger scientific recomputation.
