# Active phase: Task 34

Work only on `codex/task-34-stimulated-nhe1-alkalinity-ae4`.
Read and execute `prompts/34_stimulated_nhe1_alkalinity_ae4.md`.

This is a single-mechanism dynamic test. Do not reopen the Task 33 REST optimisation or force a 70/30 chloride partition at rest.

## Fixed scientific decisions

- Use the accepted Task 31 R09 WT resting state as the fixed initial state.
- Freeze the Cha et al. 2009 eight-state + Mod2 NHE1 mechanism and its calibrated resting carrier amount `2.339370005697548e-05 fmol`.
- Freeze all resting transporter capacities, water parameters, bath, geometry, AE4 routing/stoichiometry and membrane allocation fractions.
- The new NHE stimulation factor must equal exactly 1.0 at rest, so Task 31 REST is unchanged.
- The repository already contains beta/cAMP/PKA regulation of AE4. Verify the existing code path but do not create another AE4 regulatory mechanism or change its family/gain/time constants.
- Add/test only one missing mechanism: a stimulus-dependent multiplicative increase in effective NHE1 carrier activity under the existing CCh/Ca secretory stimulus.
- Calibrate only one scalar, `G_NHE_stim`, using the WT/control NHE-dependent alkalinization rate from Peña-Münzenmayer et al. 2015 Fig. 2C/D.
- Do not use knockout pH, knockout Cl, knockout secretion, stimulated chloride partition or the 35% secretion deficit to calibrate `G_NHE_stim`.
- Do not add NBC, stimulated carbonic anhydrase, another NHE model or any other fallback mechanism in Task 34.

## Scalar calibration only

- Hard bound: `1.0 <= G_NHE_stim <= 5.0`.
- Maximum 6 distinct gain evaluations.
- Use one bracketed scalar solve only.
- Prefer instantaneous RHS/pH-slope evaluation. Do not run a dynamic trajectory for each gain.
- If the 2015 figure cannot support a defensible numerical target, stop with `NHE_STIM_TARGET_NOT_QUANTIFIABLE`.
- If the target cannot be met inside the gain bound, stop with `NHE_STIM_GAIN_NOT_IDENTIFIED_WITHIN_BOUND`.
- Once identified, freeze/hash the gain before evaluating AE4=5% or AE4=0%.

## Dynamic test

Use only the standard central CCh + IPR stimulation:

- Ca = 0.25 uM;
- existing standard beta/IPR input;
- 0--600 s.

WT must be physically admissible and calcium must activate secretion as specified in the Task 34 prompt. If it fails, stop with `WT_STIMULATED_SECRETION_GATE_FAILED`. Do not alter any parameter.

For the accepted WT trajectory, report the integrated positive NKCC1/AE4/AE2 chloride-loading shares over 60--600 s. The approximate 70% NKCC1 share is validation context only, not a fitting target.

Then test exactly two matched-initial-state stimulated AE4 expression perturbations:

- AE4 = 0.05;
- AE4 = 0.00.

Use the same frozen Task 31 WT initial state and identical stimulus. This is an acute stimulated-expression mechanism test, not a constitutive genotype-equilibrium claim.

Do not solve separate genotype REST states in this task.
Do not test another AE4 level, another Ca level, another beta input, AE2 loss or R10.
Do not refit after seeing 5% or 0%.

## Absolute anti-combinatorial rules

There is no multidimensional optimisation in Task 34.

- No Cartesian grids.
- No parameter sweeps.
- No pairwise or factorial searches.
- No random search.
- No Latin hypercube.
- No multistart cloud.
- No global/evolutionary/Bayesian optimisation.
- No Shapley analysis.
- No automatic singles-to-pairs escalation.
- No alternative NHE mechanisms.
- No alternative AE4 regulatory families.
- No adding another bicarbonate mechanism after failure.
- No widening the gain bound or changing source-derived targets after seeing output.

## Hard compute stop

- Maximum 6 scalar gain evaluations.
- Maximum 4 dynamic integrations total, including at most one WT solver retry.
- Zero new stationary optimisation runs.
- One worker, one BLAS thread.
- Maximum 15 minutes numerical execution.

STOP immediately on the first applicable stop condition. Do not diagnose or implement another mechanism inside Task 34.

## GitHub publication

Use the connected GitHub integration for all remote writes.
Do not run `gh auth`, request a PAT, configure SSH or depend on shell `git push` credentials.
Shell Git is only for local status/diff inspection.
Publish compact UTF-8 source/test/result files directly to `codex/task-34-stimulated-nhe1-alkalinity-ae4`.
Publication failure must never trigger scientific recomputation.
Do not merge to `main`.
