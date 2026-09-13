# Active phase: Task 35

Work only on `codex/task-35-source-fixed-nhe1-activation`.
Read and execute `prompts/35_source_fixed_nhe1_activation.md`.

This is a single source-fixed dynamic mechanism test. Do not reopen Task 33 optimisation and do not repeat the failed Task 34 fluorescence-to-pH calibration.

## Fixed scientific decisions

- Use the accepted Task 31 R09 WT resting state as the common fixed initial condition.
- Freeze the Cha et al. 2009 eight-state + Mod2 NHE1 mechanism and its calibrated resting carrier amount `2.339370005697548e-05 fmol`.
- Freeze AE4 routing/stoichiometry, resting transporter capacities, water parameters, bath, geometry and membrane allocation fractions.
- Keep the existing beta/cAMP/PKA AE4 regulation unchanged.
- The only new mechanism is a muscarinic/Ca-linked stimulated NHE1 activity multiplier fixed from Park et al. 1999.
- Verify the primary source first. If it supports the reported `2.3-fold` muscarinic NHE1 activation, use exactly `G_NHE=2.3` during the standard stimulated condition and exactly `G_NHE=1.0` at rest.
- Do not fit or optimise the NHE gain.
- Do not force NKCC1 to 70% at rest. Stimulated chloride shares are outputs only.

## Absolute anti-combinatorial rules

Task 35 contains no search.

- Zero optimisation calls.
- Zero stationary solves.
- No second NHE gain.
- No gain sweep.
- No grids, factor sweeps, pairwise scans or factorial designs.
- No random, multistart, Latin-hypercube, global, evolutionary or Bayesian search.
- No alternative parameter subsets.
- No alternative NHE models or AE4 regulatory families.
- No NBC, stimulated carbonic-anhydrase mechanism or other bicarbonate-source fallback.
- No changing NKCC1, AE4, AE2 or water capacities after seeing output.
- No intermediate AE4 expression levels.
- No R10 or AE2-null tests.
- No phenotype fitting.
- No automatic escalation after failure.

## Dynamic sequence

1. Verify source support and existing AE4 beta/PKA wiring.
2. Run exactly one intended WT stimulation at Ca=0.25 uM, existing beta/IPR input, 0--600 s, with fixed stimulated `G_NHE=2.3`.
3. WT must pass the physical and secretion-activation gates in the Task 35 prompt. If it fails, stop. Do not rescue it.
4. Only if WT passes, run exactly AE4 expression 0.05 and then 0.00 from the same Task 31 R09 WT initial state.
5. Report stimulated chloride-loading shares and secretion ratios. The approximate 70/30 loading context and 35% secretion deficit are validation only, never fitting targets.

## Hard compute stop

- Maximum 4 dynamic integration attempts total, including at most one WT numerical retry for a purely solver-related failure.
- The spare retry cannot be used for another biological condition.
- One worker, one BLAS thread.
- Maximum 10 minutes numerical execution.

STOP at the first applicable source, WT, physical, numerical or budget failure. Do not diagnose or implement another mechanism inside Task 35.

## GitHub publication

Use the connected GitHub integration for all remote writes.
Do not run `gh auth`, request a PAT, configure SSH or depend on shell `git push` credentials.
Shell Git is only for local status/diff inspection.
Publish compact UTF-8 source/test/result files directly to `codex/task-35-source-fixed-nhe1-activation`.
Publication failure must never trigger scientific recomputation.
Do not merge to `main`.
