# Task 34: stimulated NHE1 alkalinity supply and AE4 secretory role

Repository: `esig626/ae4-salivary-transport-control`
Branch: `codex/task-34-stimulated-nhe1-alkalinity-ae4`
Starting head: Task 33 final head `fe754f409972ac873c4afaf3bcbb94c90feaf435`.

Read `AGENTS.md` first.

## Scientific question

Task 33 showed that forcing an approximately 70/30 chloride-loading partition at REST is not a sensible way to make AE4 relevant. The accepted R09 resting state remains physiologically usable, while the existing steady acid/base balance throttles AE4 flux to a few percent even when AE4 capacity is made very large.

The experimental biology points to a different question: during CCh + isoproterenol stimulation, NHE-dependent alkalinization supplies intracellular bicarbonate/alkalinity, while beta/cAMP/PKA signalling increases AE4 activity. Test whether adding only the missing stimulus-dependent NHE1 activation to the already-built mechanistic model is enough to give AE4 a substantial role during secretion.

This task is a **single-mechanism test**, not an optimisation programme.

## Primary sources and evidence boundary

Use only these primary sources for the new mechanism and stimulus logic:

1. Peña-Münzenmayer et al. 2015, JBC, DOI `10.1074/jbc.M114.612895`, PMCID `PMC4409235`.
   - CCh + IPR stimulation causes intracellular alkalinization in mouse SMG acinar cells.
   - With EIPA, the same stimulation causes acidification, supporting stimulation-dependent NHE activation as the source of the alkalinization.
   - The NHE-dependent alkalinization rate is not detectably changed in Ae4-null or Ae2-null cells.
   - Resting pH is nearly unchanged between genotypes.
   - Use Fig. 2C/D only to obtain the WT/control initial alkalinization-rate target and its graphical uncertainty. Do not use knockout values for calibration.

2. Peña-Münzenmayer et al. 2021, AJP Gastrointestinal and Liver Physiology, PMCID `PMC8887885`.
   - isoproterenol increases native AE4 Cl/HCO3 exchange activity;
   - the increase is associated with cAMP-dependent PKA;
   - H89 prevents the IPR-mediated activation in the reported assay.

The repository already contains the beta/cAMP/PKA -> AE4 capacity-regulation machinery in `src/modern_full_model/camp_pka.py`, and `ae4_routing_only.py` passes `regulation_gain` into the routed AE4 evaluator. **Do not build a second AE4 regulatory mechanism.** Verify the existing path with a focused test and then leave it unchanged.

Do not introduce NBC, another bicarbonate transporter, a carbonic-anhydrase stimulation law, another NHE family, or any other rescue in this task.

## Fixed resting model

Use the accepted Task 31 R09 WT resting model and its calibrated mechanistic NHE1 carrier amount:

`2.339370005697548e-05 fmol`.

The Task 31 R09 WT resting state is the fixed initial state for this task:

- pH approximately 6.91;
- Na approximately 11.64 mM;
- K approximately 116.76 mM;
- Cl approximately 60.30 mM;
- cell volume approximately 1.448 pL;
- resting q_out approximately `0.0010757073853493032 pL/s`.

Do not re-optimise REST.
Do not use any Task 33 optimised candidate as a new baseline.
Do not force resting NKCC1 to 70%.
Do not alter any resting transporter capacity, water parameter, bath, geometry, AE4 stoichiometry/routing, membrane allocation fraction, or NHE1 carrier amount.

At rest the new stimulation factor must be exactly 1.0, so the Task 31 R09 WT state and fluxes are unchanged to numerical precision.

## Exactly one new mechanism

Add one effective **stimulus-dependent NHE1 activity multiplier** to the existing Cha et al. 2009 eight-state + Mod2 NHE1 mechanism.

Implement it as a multiplicative factor on the effective active NHE1 carrier amount / cycle capacity, not by changing any Cha kinetic constant, affinity, stoichiometry, modifier exponent, or equilibrium condition.

Requirements:

- rest multiplier = exactly `1.0`;
- stimulated multiplier = one scalar `G_NHE_stim >= 1`;
- no new NHE state variable;
- no new time constant;
- no new pH clamp or bicarbonate source;
- the existing NHE1 mechanism still determines direction and state dependence;
- the multiplier acts only during the existing CCh/Ca secretory stimulus.

Use the existing model stimulus representation. If the model exposes only resting versus stimulated calcium rather than an explicit muscarinic flag, use the pre-existing central protocol distinction: the multiplier is 1.0 at the resting calcium input and `G_NHE_stim` under the standard stimulated CCh/Ca input. Do not invent an additional continuous calcium-response curve.

## One-scalar WT calibration only

Before any genotype/expression perturbation, obtain the WT/control initial NHE-dependent alkalinization-rate target from Peña-Münzenmayer et al. 2015 Fig. 2C/D.

Record:

- the plotted WT/control central estimate;
- units;
- the plotted uncertainty if recoverable;
- whether the value was read directly from axes/bar height or digitized;
- a conservative numerical acceptance interval reflecting the graphical precision.

If the figure does not support a defensible numerical interval, STOP and report `NHE_STIM_TARGET_NOT_QUANTIFIABLE`. Do not substitute a value from another tissue or species.

Calibrate only `G_NHE_stim`.

At the fixed Task 31 R09 WT state, apply the standard CCh + IPR stimulus at `t=0+` and compare the model initial pH slope with the WT/control target. Use the cheapest exact calculation available; prefer the RHS/instantaneous derivative over repeated time integrations.

Use one bracketed scalar solve only:

- hard gain range: `1.0 <= G_NHE_stim <= 5.0`;
- maximum 6 distinct gain evaluations total;
- no grid;
- no random search;
- no second scalar;
- no tuning to chloride partition or secretion.

If the target interval is already met at `G_NHE_stim = 1`, freeze 1.0 and proceed.
If the target cannot be bracketed or met inside `[1,5]`, STOP and report `NHE_STIM_GAIN_NOT_IDENTIFIED_WITHIN_BOUND`.

After calibration, hash/freeze `G_NHE_stim` before any 5% or 0% AE4 result is evaluated.

## Mandatory code-path checks before dynamics

With focused unit tests, verify all of the following:

1. rest NHE1 flux is unchanged when the stimulation mechanism is present but unstimulated;
2. the calibrated stimulation multiplier scales the Cha NHE1 carrier activity without changing stoichiometry or the zero-flux equilibrium;
3. setting NHE1 expression to zero removes the stimulated NHE source exactly;
4. existing beta/cAMP/PKA regulation changes the routed AE4 `regulation_gain` under the standard IPR/beta stimulus;
5. no new AE4 regulatory law was introduced;
6. water laws and resting water parameters are byte-for-byte / value-for-value unchanged.

If the existing AE4 PKA path is found to be disconnected by an implementation defect, repair only that wiring defect and add a regression test. Do not change its family, gain, time constants, or fit parameters.

## WT stimulated test

Use exactly one standard central WT stimulation trajectory:

- initial state: accepted Task 31 R09 WT rest;
- calcium: `0.25 uM`;
- beta/IPR input: the existing standard CCh + IPR beta input already used by the model;
- time: `0--600 s`;
- all parameters frozen except the already calibrated stimulus multiplier acting as defined above.

No calcium sweep. No beta sweep. No alternative regulatory family. No second solver trajectory unless the first integration is numerically invalid; if one numerical retry is needed, it must use the same equations and parameters.

WT dynamic acceptance requires:

- finite positive concentrations and volumes throughout;
- `Na_i < 40 mM` throughout;
- `30 <= Cl_i <= 80 mM` throughout;
- `6.6 <= pH_i <= 7.3` throughout;
- cell volume `< 3 pL` throughout;
- calcium stimulation increases secretion: mean `q_out` over 60--600 s must exceed resting q_out by at least 10%;
- cumulative secretion over 600 s must exceed the resting-flow integral by at least 10%;
- `q_out(600 s) >= q_out_rest`.

If WT fails this secretion/admissibility gate, STOP and report `WT_STIMULATED_SECRETION_GATE_FAILED`. Do not alter gain, calcium, AE4, NKCC1, CaCC, water, or any other parameter.

## Stimulated chloride-loading allocation: validation, not calibration

For the accepted frozen WT trajectory, compute integrated positive basolateral chloride loading over the sustained interval 60--600 s:

- `I_NKCC1 = integral max(J_NKCC1_Cl, 0) dt`;
- `I_AE4 = integral max(J_AE4_Cl, 0) dt`;
- `I_AE2 = integral max(J_AE2_Cl, 0) dt`;
- include any other already-existing positive basolateral chloride loader if present, but do not add one.

Define shares using the integrated positive pool.

Report the NKCC1, AE4 and AE2 shares. Compare NKCC1 with the approximate experimental ~70% context. Use `0.65--0.75` as a declared comparison band, **not** as a fitting target. Do not change `G_NHE_stim` if the share misses the band.

Also report the signed integrated AE2 flux separately so counterflux is not hidden by the positive-pool definition.

## AE4 perturbations: exactly 5% and 0%, nothing else

This task is a controlled **stimulated-expression perturbation**, not yet a final constitutive-knockout equilibrium claim.

To isolate whether the new stimulated alkalinity mechanism gives AE4 dynamic leverage, use the same accepted Task 31 R09 WT resting state as the initial condition for all three matched stimulated trajectories. At `t=0`, apply the identical CCh + IPR stimulus and set only AE4 expression to:

- WT reference: `1.00` (already run);
- test 1: `0.05`;
- test 2: `0.00`.

Do not solve or optimise separate genotype resting states in this task. The known exact-null resting-root issue is outside this single-mechanism dynamic test.

Run exactly one 5% trajectory and one 0% trajectory, with all non-AE4 parameters and `G_NHE_stim` frozen.

For each report:

- cumulative secretion 0--600 s;
- ratio to the frozen WT trajectory;
- minute-wise q_out or a compact 60 s sampled summary;
- intracellular pH and Cl trajectories;
- integrated NHE1, NKCC1, AE4 and AE2 fluxes;
- CaCC chloride export;
- paracellular chloride return;
- lumen osmolarity and water outflow;
- whether the effect becomes larger after approximately 2--3 min, as qualitative context only.

Do not fit to the experimental 35% secretion deficit.
Do not tune after seeing 5% or 0%.
Do not test any intermediate expression, AE2 deletion, R10, another calcium level, another beta input, or another NHE gain.

Clearly label the 5% and 0% results as acute/matched-initial-state stimulated-expression predictions rather than constitutive knockout equilibria.

## Absolute anti-combinatorial rules

This task contains **no multidimensional optimisation**.

Hard prohibitions:

- no Cartesian grid;
- no parameter sweep;
- no pairwise or factorial search;
- no random search;
- no Latin hypercube;
- no multistart cloud;
- no global/evolutionary/Bayesian optimiser;
- no Shapley analysis;
- no automatic singles-to-pairs escalation;
- no trying alternative NHE mechanisms;
- no trying alternative AE4 regulatory families;
- no adding NBC/CA stimulation as a fallback;
- no widening the gain bound after failure;
- no changing the source-derived alkalinization target after seeing model output;
- no phenotype fitting.

Hard numerical budget:

- at most 6 scalar gain evaluations;
- at most 1 short numerical retry for the WT dynamic integration if the first fails for solver reasons;
- exactly 1 accepted WT trajectory plus at most 1 trajectory at AE4=0.05 and 1 at AE4=0.00;
- at most 4 dynamic integrations total including any allowed WT solver retry;
- zero new stationary optimisation runs;
- one worker, one BLAS thread;
- at most 15 minutes numerical execution.

STOP immediately on the first applicable stop condition. Do not diagnose another mechanism inside Task 34.

## Required outputs

Write compact UTF-8 outputs under:

- `analysis/34_stimulated_nhe1_alkalinity_ae4/`
- `results/34_stimulated_nhe1_alkalinity_ae4/`

At minimum include:

- `source_target.md`;
- `mechanism_contract.json`;
- `gain_calibration.csv`;
- `frozen_gain.json` if identified;
- `wt_dynamic_summary.csv` if WT runs;
- `chloride_allocation.csv` if WT passes;
- `expression_comparison.csv` if 5%/0% run;
- `flux_integrals.csv` if dynamics run;
- `budget.json`;
- `verification.json`;
- `final_answer.md`.

No binary arrays, no giant trajectory tables, and no search dumps.

The final answer must state directly:

1. Was a quantitative WT NHE-dependent alkalinization target recoverable from the 2015 source?
2. What single `G_NHE_stim` was required, if any?
3. Did the fixed resting state remain unchanged?
4. Did the WT model produce calcium-activated sustained secretion?
5. What were the sustained 60--600 s NKCC1/AE4/AE2 positive chloride-loading shares?
6. With everything frozen, what were the secretion ratios at AE4=5% and AE4=0%?
7. Did the effect strengthen after ~2--3 min?
8. Does this single missing stimulation-dependent NHE mechanism plausibly restore AE4 secretory leverage, or not?
9. What remains unproven because this task used a matched WT initial state rather than constitutive genotype-specific resting equilibria?

## GitHub publication

Use the connected GitHub integration for all remote writes. Do not depend on shell Git authentication, `gh auth`, PATs or SSH credentials. Shell Git is only for local status/diff inspection.

Publish completed source edits, tests and compact outputs directly to `codex/task-34-stimulated-nhe1-alkalinity-ae4`.

Publication failure must never trigger scientific recomputation.
Do not merge to `main`.
