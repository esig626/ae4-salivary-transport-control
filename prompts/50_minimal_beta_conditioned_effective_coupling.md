# Task 50 — minimal beta-conditioned effective secretory coupling

## Mission

This task has one narrow purpose:

> Recast the already successful but target-selected Task 41 apical-secretory construction as the smallest beta-conditioned effective coupling, so that the combined CCh+IPR phenotype is retained while WT, rest, AE4-KO CCh-only and AE2-KO behaviour nest the parent model exactly.

This is **not** a mechanism search, not a parameter search, not an independent validation task, and not a signalling reconstruction.

Only three funded/scientific Codex shots remain including this run. Do not spend this run rediscovering anything in the master ledger.

## Canonical branch

Work only on:

`analysis/task-50-minimal-beta-conditioned-effective-coupling`

Do not modify `main` directly.

## Binding reading order

Read `AGENTS.md` first. Then read, in order:

1. `docs/MANDATORY_RESEARCH_LEDGER.md`
2. `docs/PHENOTYPE_TARGET_CONVENTION.md`
3. `docs/RESULTS_LEDGER.md`
4. `analysis/50_minimal_beta_conditioned_effective_coupling/NOVELTY_CHECK.md`
5. `analysis/50_minimal_beta_conditioned_effective_coupling/DECISION_LOG.md`
6. `analysis/40_ae4_equal_cation_routing/final_answer.md`
7. `analysis/41_ae4_loss_algebraic_design/final_answer.md`
8. `analysis/41_ae4_loss_algebraic_design/attempt_ledger.md`
9. `analysis/42_catalan_2025_ae4_mechanism_classes/final_answer.md`
10. `analysis/48_joint_experimental_constraint_reconstruction/JOINT_EXPERIMENTAL_CONSTRAINT_REPORT.md`
11. Task 49 final report on `analysis/task-49-camp-vrac-secretory-branch-reconstruction`
12. the active Task 40/41 source modules, especially `src/modern_full_model/ae4_cacc_recruitment.py`, Task 41 selected driver, protocol/stimulus mapping and electrical closure.

Before every new scientific decision, re-check the master ledger and record the decision in `DECISION_LOG.md` before dependent compute.

## Scientific classification

The Task 41 selected construction is phenotype-calibrated. Therefore Task 50 is also a **TARGET-CALIBRATED CONSTRUCTION**.

Do not use language implying independent validation or a discovered molecular mechanism.

The goal is model reduction and correct protocol specificity.

## One permitted model change

Implement exactly:

`g_Cl,eff = g_parent * [1 - lambda * a_Ca * beta * (1 - e_AE4)]`

with fixed

`lambda = 0.89488127156712`.

This is exactly `1 - b` from the already selected Task 41 value

`b = 0.10511872843288`.

Definitions:

- `g_parent`: the same parent stimulated apical chloride conductance used by Task 41;
- `a_Ca`: the existing normalized calcium/channel activation already used in Task 41;
- `beta`: the existing algebraic protocol input. It is zero when the beta/IPR arm is absent and one in the standard beta-containing stimulated condition. It is **not** a new dynamic state;
- `e_AE4`: the existing AE4 expression fraction;
- `lambda`: fixed, not fitted or searched.

The factor must be applied at the same place as the Task 41 conductance modifier, before full electrical closure. It must not multiply secretion or water output directly.

## Exact algebraic equivalence to prove before code changes

Task 41 selected factor:

`(1-a) + a[b + (1-b)e]`

must be expanded explicitly to

`1 - (1-b) a (1-e)`.

Therefore, at `beta=1`, the proposed Task 50 factor must be mathematically identical to Task 41.

Prove all of the following before any production calculation:

1. **WT nesting:** `e_AE4=1` -> multiplier exactly 1 for every `a_Ca` and `beta`.
2. **REST nesting:** `a_Ca=0` -> multiplier exactly 1.
3. **CCh-only AE4-KO nesting:** `beta=0` -> multiplier exactly 1 even when `e_AE4=0` and calcium activation is positive.
4. **AE2-KO nesting:** AE4 remains present, `e_AE4=1` -> multiplier exactly 1 under every protocol.
5. **Combined-stimulus Task 41 equivalence:** `beta=1` -> exact selected Task 41 factor for WT, AE4=5%, and AE4=null.

If any statement fails, stop and publish the contradiction. Do not repair it by adding another parameter.

## Absolute prohibitions

Do not:

- search, fit, profile or alter `lambda`;
- add another coupling parameter;
- add cAMP, PKA, beta-receptor, phosphorylation or calcium signalling ODEs;
- fit onset times, half-times, minute landmarks or exact trajectory shape;
- change NKCC1 kinetics, its 1.75 stimulation gain, stoichiometry, capacity or add a beta-NKCC arm;
- change AE4 cycle kinetics, stoichiometry or Na/K routing;
- change NHE1, AE2, NBC, pump, K channels, Ca amplitude, paracellular transport, water, geometry or bath;
- add VRAC, swelling gates or another apical channel;
- add genotype-specific fudge factors;
- reopen Tasks 12-49 mechanism families;
- perform any parameter sweep, mechanism panel, optimiser, grid, multistart, stochastic search or best-subset calculation;
- retune after examining any result.

## Compute discipline

Reuse frozen results whenever exact equivalence proves they already answer the question.

Do not rerun Task 40 or Task 41 production trajectories simply to reproduce numbers already frozen in the repository.

Required hierarchy:

1. algebraic proof;
2. unit/source/current-closure equivalence tests;
3. frozen artifact/hash reuse when vector field + initial state + protocol are identical;
4. only if a required equivalence cannot be established from code/tests, run the minimum single case needed to resolve it.

No broad production panel is authorised.

## Task 50A — ledger and equivalence audit

Before modifying source:

- verify branch/head and required files;
- re-read the master ledger and novelty check;
- search the repository/remote branches for an already beta-conditioned form of the Task 41 coupling;
- write `EQUIVALENCE_PROOF.md` containing the algebra above and the precise code mapping of `a_Ca`, `beta`, `e_AE4`, and conductance insertion point;
- record all decisions in `DECISION_LOG.md`;
- update `CURRENT_STATUS.md`;
- update the master ledger only if this audit establishes a new fact;
- commit, push and verify the remote SHA before implementation.

No scientific trajectory or parameter calculation in 50A.

## Task 50B — minimal implementation and exact nesting tests

Implement the smallest possible wrapper/change around the existing Task 41 conductance-recruitment machinery.

Requirements:

- one fixed `lambda` constant only;
- no new model state;
- no new kinetic parameter;
- no new protocol timing;
- no direct secretion multiplier;
- unchanged parent conductance and electrical/current closure.

Tests must demonstrate, to exact or machine precision as appropriate:

- WT multiplier = 1;
- REST multiplier = 1;
- beta=0 multiplier = 1 at AE4 null;
- AE2 KO/vector field nests the parent because AE4 remains 1;
- beta=1 Task 50 multiplier equals the Task 41 selected multiplier over representative activation/expression values, including 1, 0.05 and 0;
- fixed-state full RHS/current/charge/chloride-source equivalence to Task 41 at beta=1;
- fixed-state full RHS equivalence to the Task 40/parent model at beta=0;
- no change to unrelated transporter/water/source evaluations.

If exact nesting fails, stop. No tuning.

Publish 50B immediately with ledger/status updates and remote verification.

## Task 50C — frozen-result reuse and minimal phenotype statement

Do not treat previously published trajectories as if newly predicted.

If exact equivalence is established:

- reuse/hash-verify Task 41's frozen combined-stimulus WT, AE4=5%, and AE4=null trajectories/results;
- report the already established cumulative deficits as the Task 50 combined-stimulus consequence of an algebraically equivalent vector field:
  - AE4=5%: approximately `23.1634%`;
  - AE4=null: approximately `30.2612%`;
- state explicitly that these magnitudes are inherited from a phenotype-calibrated construction, not new held-out validation;
- prove that CCh-only AE4-null is exactly the parent model by beta=0 vector-field nesting. Reuse an existing parent result if a directly matching frozen protocol exists; otherwise the proof of identical IVP is sufficient unless numerical regression is needed for software verification;
- prove AE2 KO is exactly the parent model because `e_AE4=1`;
- restate that WT is exactly unchanged.

Do **not** require exact early/late timing agreement. Do check that the inherited combined-stimulus deficit is substantial and not a one-instant artefact, using the already frozen Task 41 cumulative/endpoint summaries rather than a new timing fit.

Report explicitly what remains unresolved:

- chronic established-genotype resting Cl/pH adaptation;
- positive AE4-KO IPR-only uptake in the beta-only arm;
- microscopic identity of the effective network coupling.

Do not fix those in this task.

Publish 50C immediately.

## Task 50D — final reduced-model report and ledger update

Create:

`analysis/50_minimal_beta_conditioned_effective_coupling/MINIMAL_EFFECTIVE_COUPLING_REPORT.md`

The report must distinguish:

### What Task 50 establishes

- the Task 41 constructive direction can be written with one beta-conditioned algebraic scalar;
- the combined CCh+IPR behaviour is mathematically identical to the already frozen Task 41 selected construction when beta=1;
- WT, rest, AE4-KO CCh-only and AE2-KO parent behaviour nest exactly under the declared conditions;
- no signalling dynamics or mechanism search are required for this reduced construction.

### What Task 50 does not establish

- no proof that AE4 molecularly recruits TMEM16A;
- no independent estimate of `lambda`;
- no solution to chronic knockout resting adaptation;
- no solution to the IPR-only KO uptake discrepancy;
- no microscopic cAMP/PKA mechanism.

### Permitted article-facing interpretation

Use language equivalent to:

> Conservation-only compensation largely masks AE4 deletion in the reconstructed transport model. A single phenotype-calibrated, beta-conditioned effective coupling between AE4 availability and stimulated apical chloride secretory capacity is sufficient to recover a substantial AE4-loss secretion phenotype while leaving WT, rest, CCh-only and AE2-control behaviour nested in the parent model. The coupling is an effective network parameter, not an identified molecular AE4-TMEM16A interaction.

Update `docs/MANDATORY_RESEARCH_LEDGER.md` in the same final checkpoint with the actual Task 50 result and exact verified commit(s).

Push and verify the remote SHA. Stop. Do not automatically start Task 51.

## Success criterion

Task 50 succeeds if the one fixed beta-conditioned scalar is implemented and the exact nesting/equivalence statements are verified without any additional mechanism or parameter search.

Scientific success does **not** require discovering the microscopic origin of `lambda`.

If the algebraic/software equivalence fails, Task 50 is a negative result. Preserve it and stop rather than consuming the run on another idea.
