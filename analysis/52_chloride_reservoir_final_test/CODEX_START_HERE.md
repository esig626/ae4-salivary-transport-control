# Task 52 — CODEX START HERE

This is the **single remaining funded/scientific Codex shot**. Treat the repository as authoritative memory and do not improvise a different mechanism.

## Required branch

Work only on:

`analysis/task-52-chloride-reservoir-final-test`

The prepared scientific parent is canonical merged main commit:

`d2b14c2821dd7e38703620389097cc8f1052977a`

Read `AGENTS.md` first.

## Binding reading order

Before any scientific action read, in order:

1. `docs/MANDATORY_RESEARCH_LEDGER.md`
2. `docs/MANDATORY_RESEARCH_LEDGER_PRE_TASK52_APPENDIX.md`
3. `docs/MANDATORY_RESEARCH_LEDGER_R51H.md`
4. `docs/MANDATORY_RESEARCH_LEDGER_R51I.md`
5. `docs/MANDATORY_RESEARCH_LEDGER_TASK52_START.md`
6. `docs/PHENOTYPE_TARGET_CONVENTION.md`
7. `docs/TASK50_EFFECTIVE_COUPLING_BENCHMARK.md`
8. `analysis/pre_task52_algebraic_gate/PRE_TASK52_ALGEBRAIC_GATE.md`
9. `analysis/pre_task52_algebraic_gate/CHLORIDE_RESERVOIR_SOLUTION.md`
10. `analysis/pre_task52_algebraic_gate/REDUCED_RESERVOIR_SIMULATION.md`
11. `analysis/51_beta_nkcc_swelling_vrac_mechanistic_reconstruction/BETA_NKCC_VRAC_RECONSTRUCTION_REPORT.md`
12. this task's `DECISION_LOG.md` and `NOVELTY_CHECK.md`
13. `prompts/52_chloride_reservoir_final_test.md`

Refer back to the ledger stack before every scientific decision. Do not rely on conversational memory.

## Non-negotiable scientific scope

Task52 tests exactly one architecture:

**measured chronic chloride reservoir + ordinary WT AE4 supply + shared non-AE4 chloride supply + shared apical secretory demand.**

It does not search mechanisms.

### Forbidden

- Palk/Benjamin NKCC in the Task52 candidate;
- `src/modern_full_model/nkcc1_palk2010.py` in the candidate dependency chain;
- Task41/Task50 AE4-expression-dependent CaCC multiplier;
- any genotype-specific apical channel conductance;
- any phenotype-selected NKCC cap, multiplier or rescue;
- parameter sweeps/grids, random/global/evolutionary search;
- post-result retuning;
- new signalling ODEs;
- new transporter families;
- changing AE4 stoichiometry/routing;
- changing pump, NHE1, NBC, K, CO2, water, geometry or paracellular laws merely to improve the phenotype;
- inventing a chronic 5%-AE4 chloride state.

### Required scientific parent

Use the accepted **pre-Palk Task37** full conservation-explicit lineage as the transport chassis. Historical later branches may be read for diagnostics/utilities only.

### Required initial-state construction

Use the measurement-constrained projection declared in `docs/MANDATORY_RESEARCH_LEDGER_TASK52_START.md`.

Central projection: keep Task37 intracellular Na and cell volume, impose measured genotype Cl/pH, solve only K and TIC from exact electroneutrality + bath-isotonicity. No stationary solve and no secretion fitting.

Alternate projection sensitivity: keep Task37 intracellular Na and TIC, solve K and volume from the same exact constraints.

If the central projection has no positive physically admissible solution, publish that as a hard stop before production trajectories. Do not invent another projection.

### Required paired supply construction

The WT evaluates the inherited pre-Palk NKCC and AE2 laws on its own state.

The matched KO receives exactly those same instantaneous NKCC and AE2 cycle fluxes, with exact stoichiometric source vectors. KO NHE1/NBC/pump/K/water/current closure remain local and dynamic. WT has ordinary AE4; KO AE4 is exactly zero.

Predeclared KO shared-supply sensitivities are x1.05 and x1.10 only.

### Required shared auxiliary current

Use the identical law in both genotypes:

`I_aux = beta * G_aux * (V_a - E_Cl)`

with only:

- `0 S`
- `2.32e-9 S`
- `4.49e-9 S`

No swelling gate and no AE4 expression input.

## Publication checkpoints

Use checkpoints `52A` through `52F` exactly as declared in the task ledger/prompt.

At every dependency boundary:

1. update `CURRENT_STATUS.md`;
2. append the mandatory ledger only if a new scientific fact/failure/qualification is established;
3. write a machine-readable checkpoint receipt;
4. commit;
5. push;
6. independently verify the remote SHA/tree/receipt;
7. only then proceed.

The final report must clearly separate SOURCE FACT, MODEL IDEALISATION, FORMAL DEDUCTION and NUMERICAL RESULT.

Stop after remotely verified `52F`. Do not start another task.
