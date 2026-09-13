# Active phase: Task 32

Work only on `codex/task-32-wt-chloride-calibration`.
Read and execute `prompts/32_wt_chloride_calibration.md`.

This is a WT-only chloride calibration task using the Task 31 mechanistic Cha NHE1 implementation. Do not start another NHE1 family or return to diagnosis of the old generic tanh law.

## Scientific direction

Task 31 removed the extreme intracellular sodium pathology but R09 WT chloride remained too high (~60.3 mM versus measured 50.10 +/- 1.50 mM). The next repair is therefore the magnitude of productive WT AE4 chloride loading, not another acid-base mechanism.

Keep the Cha NHE1 kinetic mechanism and Table S1 constants fixed. The only uncertain scalars permitted are:

1. the effective NHE1 carrier amount;
2. one multiplicative WT AE4 activity/capacity scalar distinct from genotype expression.

Use WT R09 pH and chloride only. Freeze both scalars before exact AE4-null evaluation. Knockout pH, knockout chloride and secretion are validation data only and must never enter calibration.

## Hard restrictions

- No combinatorial, factorial, pairwise, grid, random, global, Bayesian, evolutionary or high-dimensional search.
- No multidimensional optimiser.
- No automatic escalation to another transporter or mechanism.
- No arbitrary target AE4 share such as 10%, 20% or 30%; AE4 share is an output.
- Do not alter NHE1 kinetic constants, AE4 topology/stoichiometry, AE2, NKCC1, pump, K, CaCC, CO2, paracellular or water laws.
- Do not fit the AE4-null chloride, pH or 30-35% secretion phenotype.
- R09 is the calibration background. R10 is an out-of-sample stress check only after R09 is frozen.
- Calibration is the fixed bounded sequential scalar procedure in the Task 32 prompt. Do not iterate beyond it.
- Absolute numerical limits: 18 stationary solver calls, 30,000 residual evaluations, 6 AE4 scalar evaluations, 4 NHE1 scalar evaluations, 2 stimulated integrations, one worker, one BLAS thread and 30 minutes numerical execution.
- No dynamics unless an admissible matched R09 WT/exact-null REST pair exists.
- Do not merge or modify `main`.

## GitHub publication in ChatGPT Work

Use the connected GitHub integration for all remote writes.

- Do not run `gh auth`, request a PAT, configure SSH or depend on shell `git push` credentials.
- Shell Git is allowed only for local status/diff inspection.
- Publish completed source edits, focused tests and compact UTF-8 outputs directly to `codex/task-32-wt-chloride-calibration` through the connected GitHub integration.
- Do not use manual Git objects, Base64 upload loops, archives or binary trajectory publication as a workaround.
- A publication problem must never trigger scientific recomputation.

Success means WT R09 simultaneously satisfies the measured WT pH and chloride intervals using only the two declared abundance/activity scalars, after which the frozen exact AE4-null intervention is evaluated without tuning. Matching the knockout phenotype is not required for calibration success.
