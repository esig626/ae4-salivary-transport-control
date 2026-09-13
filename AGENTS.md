# Active phase: Task 33

Work only on `codex/task-33-fixed-water-nkcc70-ae4-optimisation`.
Read and execute `prompts/33_fixed_water_nkcc70_ae4_optimisation.md`.

This is the direct WT allocation/optimisation task. Do not re-diagnose AE4 routing, NHE1, the 50 mM chloride target, or the previously established NKCC1-dominant allocation.

## Fixed scientific decisions

- Freeze the Task 31 Cha et al. 2009 eight-state + Mod2 NHE1 implementation and its R09 WT-calibrated carrier amount `2.339370005697548e-05 fmol`.
- Keep the existing routed AE4 topology and working stoichiometry unchanged.
- Keep all water-law parameters unchanged and fix the WT resting outflow target at the accepted Task 31 R09 value `0.0010757073853493032 pL/s`.
- At the final WT baseline, constrain NKCC1 to exactly 70% of the positive basolateral chloride-loading pool.
- Do not impose WT chloride = 50.1 mM. Approximately 60 mM is not itself a failure.
- Optimise WT parameters only. The 5% and 0% AE4 outcomes are held out until the WT parameter payload is frozen.
- After freeze, test exactly AE4 expression 0.05 and 0.00. No other expression levels or genotypes.

## Optimisation restrictions

- One bounded low-dimensional local WT optimisation is authorised.
- No Cartesian grids, exhaustive pair searches, random/global/evolutionary/Bayesian optimisation, Shapley analysis, or automatic expansion of the free parameter set.
- Free only the predeclared existing uncertain magnitude parameters listed in the Task 33 prompt, and eliminate variables algebraically where possible.
- Do not use knockout chloride, knockout pH, knockout secretion, or the 0.65 secretion ratio during optimisation.
- Do not alter NHE1 density, water parameters, CO2, paracellular parameters, stimulus/regulation, AE4 topology/stoichiometry, or membrane allocation fractions.
- Do not merge or modify `main`.

## Compute budget

- At most 40 stationary solves total.
- At most 60,000 stationary residual evaluations total.
- At most one deterministic continuation/restart from the optimiser's own last feasible point.
- At most 2 pre-freeze WT stimulated integrations if needed for target verification.
- After freeze, at most 4 total genotype integrations covering 5% and 0% only.
- One worker, one BLAS thread, 45 minutes numerical execution.
- Stop at the first exhausted limit. Do not enlarge the search.

## GitHub publication in ChatGPT Work

Use the connected GitHub integration for all remote writes. Do not depend on shell Git authentication.

- Do not run `gh auth`, request a PAT, configure SSH, or depend on `git push` credentials.
- Shell Git is allowed only for local status/diff inspection.
- Publish completed source edits, focused tests, frozen WT checkpoint, and compact UTF-8 outputs directly to `codex/task-33-fixed-water-nkcc70-ae4-optimisation`.
- Publication failure must never trigger scientific recomputation.

The task ends after reporting the frozen WT solution and the 5% and 0% outcomes. Do not launch another repair from inside Task 33.
