# Active phase: Task 39 source-fixed Palk NKCC1 full validation

Work only on `codex/task-39-palk-nkcc1-full-validation`.
Read `analysis/39_palk_nkcc1_full_validation/source_contract.md` first, then execute `prompts/39_palk_nkcc1_full_validation.md`.

This task replaces exactly one scientific component: the generic NKCC1 concentration-response law. Everything else from Task 38 is frozen.

## Fixed source law

Use the Palk et al. / Benjamin-Johnson two-state steady NKCC1 reduction used in Vera-Sigüenza et al. 2018 Eq. 17:

`X = Na_i * K_i * Cl_i^2`

`J_NKCC1 = alpha_eff * activity_multiplier * (a1 - a2_mM*X)/(a3 + a4_mM*X)`

with

- `a1 = 157.5`
- `a2_mM = 2.0096e-5`
- `a3 = 1.0306`
- `a4_mM = 1.3852e-6`.

Positive flux is one inward `1 Na : 1 K : 2 Cl` cycle. Do not clip negative flux.

The whole-cell prefactor is determined once, algebraically, from the accepted Task 31 R09 resting NKCC1 cycle flux. Expected audit value:

`alpha_eff = 0.017334746894096052 fmol/s`.

Recompute it independently. No optimiser is permitted.

## Freeze everything else

Do not change:

- Task 31 R09 resting state;
- Cha NHE1 and its carrier amount;
- Task 36 electrogenic `1 Na : 2 HCO3` NBC and fixed capacity;
- NHE1 stimulation multiplier;
- routed AE4 stoichiometry/routing;
- AE4 beta/cAMP/PKA regulation;
- AE2;
- Na/K pump;
- K or Cl channels;
- CO2/acid-base closure;
- paracellular pathways;
- water;
- bath/geometry;
- CCh/IPR/Ca/beta stimulus;
- the existing NKCC1 N1 stimulus-dependent activity multiplier.

The old generic NKCC1 law may remain selectable in code but is not active in Task 39 production runs.

## Required execution order

1. Source-equation and unit-conversion tests.
2. Exact algebraic resting prefactor check.
3. Exact Task 31 REST nesting check, with no root solve.
4. One WT 600 s stimulated run.
5. If WT passes, one AE4=0.05 run.
6. If 5% passes, one AE4=0.00 run.
7. Compare new NKCC1 compensation against frozen Task 38 results and report the held-out secretion phenotype without fitting.

Stop on the first scientific failure. Do not rescue it.

## Absolute anti-combinatorial rule

No optimisation, root solving, grids, sweeps, one-at-a-time exploration, pairwise/factorial search, random search, Latin hypercube, multistart, global/evolutionary/Bayesian optimisation, Shapley analysis, alternative NKCC1 coefficients, alternative prefactors, alternate stimulus multiplier, alternate NBC/NHE/AE4 mechanisms, extra genotype conditions, R10, or phenotype fitting.

Maximum three intended integrations and one purely numerical retry across the whole task. One worker, one BLAS thread, maximum 15 minutes numerical execution.

## GitHub publication

Use the connected GitHub integration for all remote writes.
Do not use shell Git authentication, PATs, SSH setup or `gh auth`.
Do not merge or modify `main`.
Publish compact UTF-8 source/test/result files only.

Task 39 ends after the final report or the first declared stop condition.
