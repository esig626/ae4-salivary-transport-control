# Active repository instructions

`main` is the canonical scientific lineage. It contains the completed scientific work through merged Tasks 47 and 48. Task 49 is complete on its own analysis branch and must be treated as part of the cumulative exclusion history even before any later merge.

## Mandatory scientific memory

Before **any** new scientific analysis, model modification, numerical solve, fit, mechanism proposal, or Codex research run, read in this order:

1. `docs/MANDATORY_RESEARCH_LEDGER.md`
2. `docs/PHENOTYPE_TARGET_CONVENTION.md`
3. `docs/RESULTS_LEDGER.md`
4. `analysis/README.md`
5. the final reports for every overlapping mechanism named by the mandatory ledger
6. any task-specific prompt only after the above files have been read

The repository is the project memory. Do not rely on conversational memory or a short task summary in place of the ledger.

`docs/PHENOTYPE_TARGET_CONVENTION.md` is binding and supersedes earlier task wording that treated exact experimental timing or detailed signalling dynamics as hard reconstruction targets.

No new scientific compute is authorised until a `NOVELTY_CHECK.md` has been written for the proposed task. It must list overlapping prior tasks/files, explain exactly what remains unanswered, identify the independent evidence for any new equation or parameter, and state the falsification criterion.

Before claiming a mechanism is new, search `prompts/`, `analysis/`, `results/`, and remote branches for the mechanism and synonyms. This requirement exists because several later tasks rediscovered earlier acid-base and VRAC conclusions.

## Hard research constraints

- Do not repeat a mechanism or parameter axis closed in `docs/MANDATORY_RESEARCH_LEDGER.md` unless genuinely new independent evidence changes the scientific question.
- No transporter-subset enumeration, all-pairs/all-triples searches, Cartesian mechanism grids, best-subset/L0 searches, stochastic/global/evolutionary optimisation, or broad parameter sweeps to force the AE4 phenotype.
- Do not fit the approximately 35% AE4-null secretion deficit in isolation.
- Do not require the model to reproduce exact experimental onset times, minute landmarks, half-times, or the first-2-to-3-minute trajectory. Experimental timing is qualitative context unless a future task explicitly provides sufficient data to justify a temporal constraint.
- The final secretion reduction must be a robust model consequence, not something that appears only at one finely chosen observation time.
- Reduce CCh/beta/cAMP/PKA/calcium signalling to the smallest effective algebraic parameter or multiplier needed for transporter/channel activity. Do not add signalling ODEs or cascades merely to reproduce timing or manufacture the phenotype.
- Preserve the distinction between measured source data, model assumptions, target-selected constructions, numerical predictions and formal algebraic deductions.
- New scientific calculations require a new numbered task and branch. Each result or failure must be added to the mandatory ledger when the task finishes.
- Crash-safe publication is required for substantial tasks: push dependency-boundary checkpoints before beginning dependent work.

## Current unresolved phenotype

Any future reconstruction must address the experimental fingerprint without turning approximate biological timing into an exact mathematical target:

- substantially lower resting intracellular Cl in established AE4 knockout with resting pH close to WT;
- near-preserved AE4-KO CCh-only chloride uptake;
- positive but reduced IPR-only uptake in AE4 knockout;
- strong impairment under beta-containing stimulation;
- a substantial sustained/cumulative AE4-null secretion deficit of the experimentally observed order, robust to reasonable observation time rather than fitted to an exact trace;
- no large genotype difference in the specific isolated NKCC and stimulated NHE assays;
- essentially secretion-neutral AE2 knockout.

The inherited acute model is known to fail both the chronic resting-state phenotype and the beta-containing knockout response. These failures have already been analysed extensively; see the mandatory ledger before proposing a repair.

## Manuscript development

The active writing branch remains:

`manuscript/bmb-ae4-decade-reassessment`

For manuscript-only work, additionally read:

1. `README.md`
2. `manuscript/README.md`
3. `manuscript/figure_plan.md`
4. `manuscript/claims_and_evidence.md`

Writing constraints:

- Use British English.
- Preserve the distinction between source-derived experimental findings, model results and target-selected constructions.
- Do not claim that Catalan et al. 2025 established the proposed whole-cycle stoichiometries experimentally.
- Do not describe Task 41 as independent validation.
- Do not identify the NBC-like structural pathway with a specific salivary isoform unless new evidence is added.
- Keep LaTeX build products out of version control.
- Manuscript editing alone must not alter frozen scientific result files.

## Repository integrity

- Do not rewrite or delete `archive/`, numbered analyses or frozen results.
- Do not overwrite historical negative results merely because a later model differs.
- Do not merge a manuscript change until the LaTeX source compiles, citations and references resolve, and figures render without clipping.
