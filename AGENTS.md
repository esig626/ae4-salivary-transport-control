# Active phase: Task 37 WT-only NBC validation

Work only on `codex/task-37-wt-nbc-validation`.
Read and execute `prompts/37_wt_nbc_validation.md`.

This task validates the corrected Task 36 architecture in WT only. It is not a calibration, search, genotype or mechanism-discovery task.

## Fixed scientific decisions

- Start from the accepted Task 31 R09 WT resting state.
- Keep the Task 31 Cha et al. mechanistic NHE1 carrier amount `2.339370005697548e-05 fmol`.
- Keep routed AE4 with `1 Cl : 1 monovalent cation : 2 HCO3` working stoichiometry.
- Keep the Task 36 electrogenic NBC surrogate exactly as implemented:
  `Na_o + 2 HCO3_o <-> Na_i + 2 HCO3_i`.
- The NBC source signature is exactly `(+1 Na, 0 K, 0 Cl, +2 TIC, +2 TA)` per inward cycle.
- NBC carries net `-1` charge-equivalent into the cell and participates in basolateral current closure.
- Keep `G_NBC = 0.11570913197464398 fmol/s` fixed.
- Keep `u_sec = clip((Ca - 0.058)/(0.25 - 0.058), 0, 1)` fixed.
- Keep the NHE1 stimulation multiplier `1 + 1.3 u_sec`, giving exactly `1.0x` at REST and `2.3x` at `Ca=0.25 uM`.
- Preserve AE2, NKCC1, pump, K channels, CaCC, CO2, paracellular transport, water laws, bath, geometry and existing AE4 beta/PKA regulation.
- Do not force 70/30 at REST.
- Do not fit the WT stimulated chloride partition. Compare it only with the 65--75% NKCC1 context after the run.

## Phase A: implementation and REST nesting only

Run the focused Task 36 NBC tests and the minimum directly required Task 31 NHE regressions.

Verify stoichiometry, electrogenic current sign, voltage-dependent affinity, reversal, source bookkeeping and exact zero NBC recruitment at REST.

Evaluate the accepted Task 31 R09 resting state directly under both the Task 31 model and Task 36 wrapper. Do not solve another root.

If a literal implementation bug prevents the already-declared Task 36 equations from being represented, make only the smallest correction consistent with those equations. Do not alter any scientific choice or parameter.

If focused tests or REST nesting fail after that one correction opportunity, stop with:

`TASK37_IMPLEMENTATION_OR_REST_NESTING_FAILED`

## Phase B: exactly one WT stimulation

Only after Phase A passes, run one intended WT trajectory:

- CCh `0.3 uM`
- isoproterenol `5 uM`
- Ca `0.25 uM`
- standard beta/IPR input
- `0--600 s`
- accepted Task 31 R09 WT state as the initial condition
- AE4 expression exactly `1.0`

No stationary solve is permitted.

One retry is allowed only for a clearly numerical ODE solver failure and may change only to another already-declared production stiff solver configuration.

If WT violates physical/conservation gates, stop with:

`WT_NBC_TRAJECTORY_PHYSIOLOGY_FAILED`

If WT fails the secretion activation gate, stop with:

`WT_NBC_SECRETION_GATE_FAILED`

Do not rescue either failure.

## Absolute anti-combinatorial rules

There is no search in Task 37.

- No optimisation.
- No stationary solve.
- No grids or sweeps.
- No one-at-a-time parameter exploration.
- No pairwise/factorial search.
- No random search.
- No multistart.
- No global/evolutionary/Bayesian optimisation.
- No Shapley analysis.
- No NBC capacity change.
- No alternate NBC stoichiometry or law.
- No NHE1 gain change or alternate NHE model.
- No AE4/NKCC1/AE2/pump/channel/CO2/paracellular/water changes.
- No alternative Ca, CCh, beta or isoproterenol input.
- No alternate AE4 regulatory family.
- No R10.
- No genotype perturbation.
- No AE4 5%.
- No AE4 0%.
- No automatic escalation to another mechanism.

If WT fails, report and STOP.

## Hard compute stop

- zero stationary solves;
- zero optimisation calls;
- exactly one intended WT integration;
- at most one numerical WT retry;
- maximum two integration attempts total;
- no genotype integrations;
- one worker;
- one BLAS thread;
- maximum 10 minutes numerical execution.

## GitHub publication

Use the connected GitHub integration for all remote writes.
Do not run `gh auth`, request a PAT, configure SSH or depend on shell `git push` credentials.
Shell Git is only for local status/diff inspection.
Publish compact UTF-8 outputs directly to `codex/task-37-wt-nbc-validation`.
Do not publish binary trajectories or archives.
Do not merge or modify `main`.

The task ends after the WT report. Do not continue to AE4 perturbations inside Task 37.
