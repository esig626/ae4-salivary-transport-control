# Task 51 Codex start here

## Canonical branch

`analysis/task-51-beta-nkcc-swelling-vrac-mechanistic-reconstruction`

## Scientific parent

Clean merged `main` at:

`4bb2c89870fe3f4d4bb887cf9d40ba7de377a725`

This parent contains completed Tasks 49 and 50, the frozen Task 50 benchmark, the corrected two-shot budget, and the source-backed Task 51 literature handoff.

## Mandatory reading order

Read `AGENTS.md` first.

Then read, in this exact order:

1. `docs/MANDATORY_RESEARCH_LEDGER.md`
2. `docs/PHENOTYPE_TARGET_CONVENTION.md`
3. `docs/TASK50_EFFECTIVE_COUPLING_BENCHMARK.md`
4. `docs/TASK51_BETA_NKCC_SWELLING_VRAC_EVIDENCE.md`
5. `analysis/51_beta_nkcc_swelling_vrac_mechanistic_reconstruction/NOVELTY_CHECK.md`
6. `analysis/51_beta_nkcc_swelling_vrac_mechanistic_reconstruction/DECISION_LOG.md`
7. `analysis/49_camp_vrac_secretory_branch_reconstruction/CAMP_VRAC_RECONSTRUCTION_REPORT.md`
8. `analysis/49_camp_vrac_secretory_branch_reconstruction/VRAC_IMPLEMENTATION.md`
9. `analysis/48_joint_experimental_constraint_reconstruction/JOINT_EXPERIMENTAL_CONSTRAINT_REPORT.md`
10. `analysis/39_palk_nkcc1_full_validation/final_answer.md`
11. `analysis/40_ae4_equal_cation_routing/final_answer.md`
12. `analysis/41_ae4_loss_algebraic_design/final_answer.md`
13. active source under `src/modern_full_model/`, especially `nkcc_stimulation.py`, `nkcc1_palk2010.py`, `camp_pka.py`, `nbc_minimal.py`, `model.py`, membrane/current closure, volume/water equations and Task 49 VRAC code.
14. `prompts/51_beta_nkcc_swelling_vrac_mechanistic_reconstruction.md`

Before **every scientific decision**, refer back to `docs/MANDATORY_RESEARCH_LEDGER.md`. It is binding project memory, not just startup reading.

## Mission in one line

Test exactly one literature-backed mechanistic chain:

`beta/IPR -> beta/cAMP activation of NKCC1 -> solute loading -> cell swelling -> VRAC-like apical anion conductance -> fluid secretion`

in parallel with the already represented

`beta/PKA -> AE4 activation`.

## Critical exclusions

- Do **not** use the Task 50 effective-coupling multiplier in the Task 51 model.
- Do **not** tune anything to the approximately 35% AE4-null secretion phenotype before the immutable prediction checkpoint.
- Do **not** introduce a direct AE4-expression dependence into VRAC/TMEM16A/apical conductance.
- Do **not** replace the Palk/Benjamin NKCC core, alter NKCC stoichiometry, impose an NKCC cap, or search NKCC mechanism families.
- Do **not** change NHE1, NBC, pump, K channels, AE4 stoichiometry/routing, calcium amplitude, water laws, geometry or paracellular laws.
- Do **not** add signalling ODEs. Beta/cAMP regulation is represented by the smallest effective algebraic activity term supported by source data.
- Do **not** search mechanism combinations, transporter subsets, parameter subsets, grids, random/evolutionary models or alternative VRAC gate families.

Task 51 changes one upstream regulatory input relative to Task 49: beta/cAMP regulation of the existing NKCC1 activity multiplier. The downstream positive-swelling VRAC law is the already conservation-tested Task 49 construction unless an exact implementation defect is found.

## Execution structure

Use one orchestrator and up to five read-only/advisory workers. Only the orchestrator may modify the canonical branch, accept scientific conclusions, run authoritative inference/trajectories, commit or push.

Publish and remotely verify every dependency boundary before dependent work. Required checkpoints are `51A` through `51F` as defined in the prompt.

This is the penultimate scientific shot. Do not stop after rediscovering Task 49's self-start failure; the entire purpose of Task 51 is to test the independently supported upstream beta-NKCC input that Task 49 deliberately excluded.

Stop after verified `51F`. Do not start Task 52 automatically.
