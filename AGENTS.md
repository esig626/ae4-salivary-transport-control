# Active repository instructions

`main` is the canonical scientific lineage. It contains the completed scientific work through merged Tasks 47 and 48 plus the repository-wide cumulative research ledger and phenotype-target convention. Task 49 is complete on its own analysis branch and is part of the cumulative scientific memory even before any later merge.

## Mandatory scientific memory and decision discipline

Before **any** scientific analysis, model modification, numerical solve, fit, mechanism proposal, interpretation decision, or Codex research run, read in this order:

1. `docs/MANDATORY_RESEARCH_LEDGER.md`
2. `docs/PHENOTYPE_TARGET_CONVENTION.md`
3. `docs/RESULTS_LEDGER.md`
4. `analysis/README.md`
5. the final reports for every overlapping mechanism named by the mandatory ledger
6. the task-specific `NOVELTY_CHECK.md`
7. only then the task-specific execution prompt

The repository is the project memory. Do not rely on conversational memory, a short task summary, or an agent's recollection in place of the ledger.

`docs/MANDATORY_RESEARCH_LEDGER.md` is the authoritative cumulative scientific ledger. It contains the experimental constraints, model architecture, exact algebraic facts, task-by-task results, exclusions, unresolved questions, interpretation classes, and remaining-shot budget. **Codex must refer back to it before every scientific decision, not merely once at task start.**

`docs/PHENOTYPE_TARGET_CONVENTION.md` is binding and supersedes earlier task wording that treated exact experimental timing or detailed signalling dynamics as hard reconstruction targets.

No new scientific compute is authorised until a task-local `NOVELTY_CHECK.md` exists and explicitly lists:

- the scientific question;
- every overlapping ledger entry/task/file;
- why the proposed work is not already answered;
- the exact equation/parameter being changed, if any;
- whether any parameter is source-identified, assumed, or phenotype-calibrated;
- what outcome would falsify the proposal.

Before claiming a mechanism is new, search `prompts/`, `analysis/`, `results/`, source modules and remote branches for the mechanism and synonyms.

## Decision log requirement

Every substantial scientific decision during a task must be written to the task's `DECISION_LOG.md` **before dependent compute**. Each entry must cite the relevant section/result identifiers from `docs/MANDATORY_RESEARCH_LEDGER.md` and state why the decision does not reopen a closed axis.

If a proposed decision conflicts with the master ledger, the agent must stop and report the conflict rather than improvise around it.

## Ledger publication rule

A scientific checkpoint is incomplete unless every new result, failure, exclusion, qualification or superseding interpretation established by that checkpoint is also added to `docs/MANDATORY_RESEARCH_LEDGER.md` in the same dependency boundary.

Do not silently rewrite history. Preserve the earlier entry and append a superseding/qualifying entry with the relationship made explicit.

## Hard research constraints

- There are **three remaining funded/scientific Codex execution shots** before Task 50 begins. Do not consume one on rediscovery.
- Do not repeat a mechanism or parameter axis closed in `docs/MANDATORY_RESEARCH_LEDGER.md` unless genuinely new independent evidence changes the scientific question.
- No transporter-subset enumeration, all-pairs/all-triples searches, Cartesian mechanism grids, best-subset/L0 searches, stochastic/global/evolutionary optimisation, or broad parameter sweeps to force the AE4 phenotype.
- Do not fit the approximately 35% AE4-null secretion deficit in isolation and then call the result validation. If the phenotype is deliberately used to select a parameter, label the result `TARGET-CALIBRATED CONSTRUCTION` everywhere.
- Do not require exact experimental onset times, minute landmarks, half-times, or the first-2-to-3-minute trajectory. Experimental timing is qualitative context unless a future task provides sufficient data to justify a temporal constraint.
- The final secretion reduction must be robust rather than appearing only at one finely selected observation time.
- Reduce CCh/beta/cAMP/PKA/calcium signalling to the smallest effective algebraic parameter or multiplier needed for transporter/channel activity. Do not add signalling ODEs or cascades merely to reproduce timing or manufacture the phenotype.
- Preserve the distinction between source facts, formal deductions, assumptions/sensitivities, numerical predictions, target-calibrated constructions, and negative/exclusion results.
- New scientific calculations require a new numbered task and branch.
- Crash-safe publication is required for substantial tasks: push and remotely verify every dependency-boundary checkpoint before dependent work.

## Current modelling objective

The remaining problem is not an open mechanism search. The repository already establishes that conservation-only/state-compensating models largely erase AE4 loss, while the Task 41 target-calibrated one-parameter apical-export coupling generates the correct order of secretion deficit.

The immediate next task, if launched, is the narrow Task 50 reduced construction predeclared in the mandatory ledger:

`g_Cl,eff = g_parent * [1 - lambda * a_Ca * beta * (1 - e_AE4)]`

with fixed `lambda = 0.89488127156712`, inherited from Task 41's already phenotype-calibrated `b=0.10511872843288`.

Task 50 is not allowed to search or optimise `lambda`, add signalling states, change NKCC/NHE/NBC/pump/K/AE4 stoichiometry, or invent a molecular interpretation. Its purpose is to test/repackage the already constructive Task 41 direction as the smallest beta-conditioned effective coupling while preserving WT, CCh-only and AE2 nesting.

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
- Do not describe Task 41 or Task 50 as independent validation if their effective coupling parameter derives from the phenotype target.
- Do not identify the NBC-like structural pathway with a specific salivary isoform unless new evidence is added.
- Keep LaTeX build products out of version control.
- Manuscript editing alone must not alter frozen scientific result files.

## Repository integrity

- Do not rewrite or delete `archive/`, numbered analyses or frozen results.
- Do not overwrite historical negative results merely because a later model differs.
- Do not merge a manuscript change until the LaTeX source compiles, citations and references resolve, and figures render without clipping.
