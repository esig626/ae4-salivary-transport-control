# Task 46I addendum: parallel-agent CBM execution without combinatorial search

This addendum is mandatory for Task 46H and supersedes any wording in Task 46D, 46E or 46H that would require independent structural tests to be executed serially when they can be evaluated safely in parallel from the same frozen CBM state.

The scientific search remains hypothesis-driven and non-combinatorial. The computational execution should be parallel wherever tasks are independent.

## Core distinction

Do not confuse:

- **parallel evaluation of independent predeclared hypotheses**, which is encouraged; with
- **combinatorial model search**, which remains forbidden.

Once a CBM network, evidence split, stoichiometric matrix and bound set are frozen at a published checkpoint, independent analyses that do not change that frozen base may be fanned out to separate agents and run concurrently.

Examples that should be parallelised include:

- independent provenance audits of different transporter bounds;
- WT FVA endpoint groups;
- AE4-null FVA endpoint groups;
- the predeclared targeted flux-coupling questions;
- the predeclared one-at-a-time structural constraint tests;
- independent verification of rank, conservation, charge and endpoint feasibility;
- sensitivity to individually justified alternative numerical bounds when those alternatives were predeclared before computation.

Parallel agents must not recursively generate new hypotheses or launch combinations of their assigned perturbations.

## One orchestrator, read-only workers

Use one lead/orchestrator agent responsible for the Task 46 branch, checkpoint state, scientific ledger and Git publication.

Worker agents should operate from the same immutable published base for the current milestone. They may use isolated scratch directories/worktrees if available, but they must not independently push competing histories to the Task 46 branch.

Worker agents must return compact machine-readable outputs plus a short scientific interpretation to the orchestrator.

Only the orchestrator may:

- change `CURRENT_STATUS.md`;
- create the official `CHECKPOINT_46H0X_*.md` file;
- integrate milestone results into the canonical Task 46 output tree;
- commit and push the milestone;
- advance the published base used by subsequent milestones.

This avoids Git races while still obtaining parallel scientific computation.

## Milestone parallelisation plan

### 46H00: recovery and evidence split

Parallelise evidence review by domain if useful:

- transporter stoichiometry/localisation;
- WT physiology/bounds;
- AE4 phenotype and holdouts;
- NKCC1 evidence;
- carbon/alkalinity/charge constraints.

The orchestrator reconciles the reviews into one frozen construction/holdout ledger before publication.

### 46H01: network construction

The lead agent constructs the single canonical network and S matrix.

Once a draft is frozen locally, launch independent verification agents in parallel for:

1. Na/K/Cl conservation and transporter stoichiometry;
2. TIC/TA bookkeeping;
3. charge/electroneutrality consistency;
4. rank/nullspace calculation;
5. bound/provenance audit.

These agents verify one network. They do not propose alternative network architectures unless they identify a concrete error.

The orchestrator resolves any errors, then publishes the one accepted network.

### 46H02: WT feasible-space analysis

Use the frozen 46H01 network.

FVA may be parallelised across reaction groups using the local VFFVA engine. Additional agents may independently inspect:

- chloride-loading routes;
- sodium/pump balance;
- alkalinity/carbon support;
- apical secretion/output routes.

All agents evaluate the same WT feasible region. No agent may alter bounds independently.

### 46H03: AE4-null feasible-space analysis

Use one frozen AE4-null region with shared non-genotype-specific capacities.

Parallelise endpoint groups and independent diagnostics. Agents may separately quantify:

- maximal secretion proxy;
- NKCC1 range;
- chloride loading/export balance;
- TIC/TA support;
- pump demand;
- complete feasible witnesses for key extrema.

Again, all agents use exactly the same KO region.

### 46H04: paired WT/KO NKCC structural discriminator

The orchestrator freezes the justified cross-genotype NKCC constraint before workers run.

Then parallelise analyses of the paired region, for example:

- maximum KO secretion;
- minimum NKCC change required for selected secretion levels;
- feasible WT/KO flux differences;
- chloride and alkalinity consequences.

Do not let separate agents invent different epsilon values after seeing results. Any numerical tolerance must be predeclared from evidence.

### 46H05: targeted flux coupling

Predeclare at most six biologically motivated coupling questions before launching workers.

Assign one coupling question per agent and evaluate them concurrently.

Each agent may compute exact/interval coupling information needed for its assigned pair or small block only.

Do not generate an all-pairs coupling matrix. Do not expand from an interesting result into an unapproved family of nearby pairings.

The orchestrator synthesises the six or fewer results into one structural-redundancy conclusion.

### 46H06: minimal structural discriminator / repair

Predeclare at most six candidate structural constraint groups based on 46H00-46H05.

The six one-at-a-time tests should be run **in parallel**, one worker per candidate constraint group, from the same frozen base.

Each worker tests only its assigned single structural change and reports:

- whether WT remains feasible;
- whether the AE4-null excessive-rescue problem changes;
- maximal KO secretion/output proxy;
- NKCC1 requirement;
- any violated hard physiological gates;
- whether the single change is sufficient, necessary-looking, irrelevant or invalid.

Do not have workers test pairs, triples or subsets.

After all single-change results return, the orchestrator ranks them.

If at least one single change gives a scientifically defensible structural resolution, stop there and do not launch joint combinations merely to seek a numerically prettier answer.

If no single change is sufficient but the returned results mathematically demonstrate that a joint change is necessary, permit **one** sparse joint analysis. Use an L1, lexicographic or minimal-cardinality formulation where practical. This is one optimisation problem, not enumeration of subsets.

Do not enumerate all pairs, triples or subsets even in parallel.

### 46H07: decision report

Parallel agents may independently review the assembled evidence for:

- mathematical consistency;
- physiological interpretation;
- evidence/holdout leakage;
- overclaiming;
- reproducibility.

The orchestrator writes one final `CBM_DECISION_REPORT.md`, publishes it, and stops.

## Concurrency budget

Use parallelism aggressively but coherently.

- Up to 6 scientific worker agents may run concurrently for independent analyses.
- A separate verification/reviewer agent may be used where useful.
- Do not create agents whose only purpose is to duplicate identical calculations unless an independent implementation is scientifically valuable.
- Do not create recursive agent trees.
- Worker agents do not spawn further agents unless the orchestrator explicitly assigns that subtask.

## Anti-combinatorial rules remain hard

Parallelism does not relax the Task 46D search limits.

Still forbidden:

- Cartesian products of transporter variants;
- all-subset structural searches;
- pair/triple enumeration after single-change tests;
- multidimensional parameter grids;
- arbitrary objective sweeps;
- exhaustive elementary-flux-mode enumeration;
- multiple competing CBM architectures;
- worker-generated model variants that were not predeclared by the orchestrator.

The rule is:

> **fan out calculations, not hypotheses.**

## Publication / crash-resilience rule

Parallel work does not change checkpoint discipline.

At the end of every 46H milestone:

1. collect all worker outputs;
2. independently check consistency and failures;
3. integrate only the accepted milestone products;
4. update `CURRENT_STATUS.md`;
5. create the immutable milestone checkpoint and machine-readable receipt;
6. commit;
7. push immediately;
8. verify the remote SHA and clean tree;
9. only then start work whose assumptions depend on that milestone.

Independent work within a milestone may run concurrently. Work that depends scientifically on an unpublished milestone must wait for publication.

## Efficiency principle

Do not rerun completed expensive analyses simply because another worker finished later. Cache and reuse verified outputs from the same frozen model identity.

Use VFFVA for broad endpoint calculations after the already completed Task 46G solver verification. Reference FVA should be used only for targeted audit checks, not duplicated for every milestone unless a discrepancy appears.

The objective is to finish the CBM structural diagnosis quickly without sacrificing scientific control.
