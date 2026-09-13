# Active phase: Task 33

Work only on `codex/task-33-fixed-water-nkcc70-ae4-optimisation`.
Read and execute `prompts/33_fixed_water_nkcc70_ae4_optimisation.md`.

This is the direct WT allocation/optimisation task. Do not re-diagnose AE4 routing, NHE1, the 50 mM chloride point target, or the previously established NKCC1-dominant allocation.

## Fixed scientific decisions

- Freeze the Task 31 Cha et al. 2009 eight-state + Mod2 NHE1 implementation and its R09 WT-calibrated carrier amount `2.339370005697548e-05 fmol`.
- Keep the existing routed AE4 topology and working stoichiometry unchanged.
- Keep all water-law parameters unchanged and fix the WT resting outflow target at `0.0010757073853493032 pL/s`.
- Constrain final WT NKCC1 to `0.70 +/- 0.01` of the positive basolateral chloride-loading pool.
- Use a hard WT intracellular chloride guardrail `45 <= Cl_i <= 65 mM`. This is a range, not a point target, and it must not be widened during the task.
- Require WT `6.84 <= pH_i <= 6.98`, `Na_i <= 30 mM`, positive concentrations and cell volume `< 3 pL`.
- A stationary WT fit is not sufficient. Before freeze, the same candidate must pass the WT Ca = 0.25 uM, 0--600 s secretion-activation gate specified in the Task 33 prompt.
- Dynamic activation requires mean q_out over 60--600 s >= 1.10 times resting q_out, cumulative 600 s secretion >= 1.10 times the resting-flow integral, and q_out(600 s) >= resting q_out.
- The calcium input, stimulus equations and activation thresholds are fixed and must not be altered after seeing results.
- Optimise WT parameters only. The 5% and 0% AE4 outcomes are held out until the WT parameter payload is frozen.
- After freeze, test exactly AE4 expression 0.05 and 0.00. No other expression levels, backgrounds or genotypes.

## Absolute anti-combinatorial rules

- One bounded low-dimensional local WT optimisation only.
- Before numerical optimisation, write and checkpoint the immutable optimisation contract with the exact free parameters, bounds, objective, algorithm, restart rule, stationary gate and calcium-activation gate.
- Maximum 20 distinct WT parameter vectors total.
- Maximum one deterministic restart, from the first run's own best point only.
- Maximum 3 pre-freeze WT dynamic activation integrations, and only for candidates that already pass stationary constraints.
- No Cartesian grids, factor sweeps, pairwise scans, factorial designs, random search, multistart clouds, Latin hypercube, global/evolutionary/Bayesian optimisation, Shapley analysis, or automatic singles-to-pairs escalation.
- Do not try multiple optimiser families.
- Do not loop over alternative parameter subsets.
- Do not change bounds, targets, objective weights, chloride range, calcium gate or free parameters after seeing a failure.
- Do not enlarge the search because 5% or 0% gives the wrong phenotype.
- Free only the predeclared existing uncertain magnitude parameters in the Task 33 prompt and eliminate variables algebraically where possible.
- Do not use knockout chloride, knockout pH, knockout secretion, or the 0.65 secretion ratio during WT optimisation.
- Do not alter NHE1 density, water parameters, CO2, paracellular parameters, stimulus/regulation, AE4 topology/stoichiometry, or membrane allocation fractions.
- Do not merge or modify `main`.

## Hard stop budget

- At most 30 stationary solves total.
- At most 40,000 stationary residual evaluations total.
- At most 20 distinct WT parameter vectors total.
- At most 3 pre-freeze WT dynamic activation integrations total.
- At most 20 minutes numerical execution for WT optimisation plus activation-gate checks.
- One worker and one BLAS thread.
- After freeze, at most one 5% and one 0% dynamic trajectory. The accepted pre-freeze WT activation trajectory is the frozen WT denominator and must not be rerun just to duplicate it.

STOP immediately when any limit is reached. If no WT candidate has passed both stationary and calcium-activation gates within these predeclared bounds, report `WT_OPTIMISATION_FAILED_WITHIN_PREDECLARED_BOUNDS` and end the task. Do not diagnose another mechanism and do not keep searching.

If the 5% or 0% resting state fails after freeze, report it. Do not insert intermediate AE4 expression levels to help the solver.

## GitHub publication in ChatGPT Work

Use the connected GitHub integration for all remote writes. Do not depend on shell Git authentication.

- Do not run `gh auth`, request a PAT, configure SSH, or depend on `git push` credentials.
- Shell Git is allowed only for local status/diff inspection.
- Publish completed source edits, focused tests, frozen WT checkpoint and compact UTF-8 outputs directly to `codex/task-33-fixed-water-nkcc70-ae4-optimisation`.
- Publication failure must never trigger scientific recomputation.

The task ends after reporting either the bounded WT failure or the frozen WT plus 5% and 0% outcomes. Do not launch another repair from inside Task 33.
