# Task 46H: milestoned CBM structural screening execution

This task is the next mandatory scientific stage of Task 46. It begins from the remotely published Task 46G checkpoint and must complete the constraint-based structural screening before any new kinetic or ODE model is constructed.

## Repository and starting state

Use repository:

`esig626/ae4-salivary-transport-control`

Work only on branch:

`analysis/task-46-physiology-constrained-model-reconstruction`

Resume from remotely published commit:

`9553d04442aeef5b0a8e0ba5f9df07f0986d7074`

Do not modify `main`.

Read `AGENTS.md` first. Then read and obey the Task 46 specifications, with particular priority to:

- `prompts/46d_anti_combinatorial_model_search_discipline.md`
- `prompts/46e_cbm_structural_screening_before_kinetic_model_discovery.md`
- `prompts/46g_use_local_cbm_core_and_isolate_carbonscope.md`
- this Task 46H prompt.

Where there is tension, this prompt controls execution order for the CBM stage; Task 46D controls search discipline; Task 46G controls repository isolation.

## Hard repository isolation

Use only the already verified local solver under:

`analysis/46_physiology_constrained_model_reconstruction/vendor/carbonscope_flux_core/`

Do not clone, modify, branch, commit, push, open pull requests in, or otherwise alter `esig626/CarbonScope`.

Do not access CarbonScope operationally unless a read-only source comparison is genuinely necessary. There should be no need to do so for this task.

Treat the local vendor core as frozen after Task 46G verification. Do not edit it for convenience. If a small Task-46-specific paired-LP or coupling utility is required, implement it above the vendor layer inside the Task 46 analysis directory and test it independently.

## Scientific objective

Build one compact, biologically curated constraint-based model of the salivary acinar transport architecture and use it to determine what sustained flux relationships are structurally possible in WT and AE4-null cells.

The CBM is a screening and diagnosis layer. It is not the final dynamic model and must not be tuned into one.

The specific questions are:

1. Is physiologically plausible WT sustained transport feasible under independently justified constraints?
2. With AE4 removed and all non-genotype-specific capacities shared, can the network still sustain near-WT apical chloride export / secretion proxy?
3. Is strong NKCC1 rescue structurally necessary, merely optional, or excluded once carbon/alkalinity, pump and charge constraints are retained?
4. Are AE4 and NKCC1 actually interchangeable flux routes once the full retained balance structure is represented?
5. Which one or two mechanistic constraints most strongly control the excessive rescue seen in the old dynamic model?
6. Which kinetic block or blocks should therefore be rebuilt next?

Do not use the CBM to fit the central ~35% secretion deficit. The exact null magnitude, 5% AE4 trajectory, timing phenotype, detailed sodium/chloride trajectories, initial CaCC behaviour and other dynamic observables remain later validation material unless an earlier Task 46 specification explicitly classifies one of them as construction evidence.

The qualitative fact that AE4 loss causes a substantial sustained secretion deficit may be used as a biological discriminator, but do not tune flux bounds until the CBM reproduces 35%.

## Evidence discipline

Before the first scientific optimisation, freeze a construction/holdout split.

Construction evidence may include:

- exact transporter stoichiometry;
- localisation;
- conserved-species bookkeeping;
- thermodynamically justified directions;
- WT physiological ranges or capacities with source provenance;
- AE4-null condition `v_AE4 = 0`;
- the independently motivated salivary NKCC1 evidence only to the extent needed for the later paired WT/KO structural discriminator;
- the qualitative fact that AE4 loss causes a substantial sustained secretion deficit.

Reserve the more detailed phenotype information for later validation wherever possible.

For every numerical bound classify it as one of:

- directly measured / experimentally bounded;
- literature-derived;
- inherited physical constant;
- algebraically derived;
- conservative screening bound;
- numerical big-M / computational bound.

A computational big-M must never be presented as physiology. Demonstrate that conclusions of interest are insensitive to it or replace it with a finite physiological bound.

## Model scope

Construct exactly one base CBM architecture during this task.

Represent at minimum the retained balances needed for:

- Na;
- K;
- Cl;
- total inorganic carbon or the minimal equivalent carbonate bookkeeping needed for steady flux structure;
- total alkalinity;
- charge/electroneutrality audit.

Charge may be algebraically redundant. Retain the full biological audit but solve only an independent equality system when rank analysis shows redundancy.

Represent the principal transport processes required by the existing AE4 model and evidence, including where justified:

- NKCC1;
- AE4;
- AE2;
- NHE1;
- NBC-like bicarbonate entry;
- Na/K ATPase;
- apical chloride exit / CaCC;
- effective K exit;
- CO2 / acid-base exchange terms required by the retained TIC/TA balances;
- paracellular / luminal bookkeeping required to close the sustained transport problem.

Do not infer stoichiometries from transporter names. Extract and document them from pinned model/source evidence and check charge consistency explicitly.

External reservoirs need not be balanced. Every balanced row and every reaction column must have a documented biological meaning.

The primary CBM output is sustained apical chloride export and related conserved fluxes. Do not call this literal saliva flow unless the model explicitly includes and justifies the required water/osmotic conversion.

## No combinatorial rampage

This is a hard constraint, not a suggestion.

### Single architecture rule

Use one curated base network. Do not create variants of the whole network.

If a conservation audit proves a missing reaction is required for feasibility, add only that specific reaction, document why, checkpoint the correction, and continue. Do not spawn architecture alternatives.

### No grids

Do not run multidimensional parameter sweeps, random search, global optimisation or Cartesian products of flux bounds.

Do not enumerate elementary flux modes.

Do not construct all subsets of candidate repairs.

Do not explore multiple objectives for curiosity.

### Targeted coupling budget

Predeclare no more than six coupling questions before the coupling milestone. Prefer the following unless the network audit proves one irrelevant:

1. AE4 versus NKCC1;
2. AE4 versus total basolateral chloride loading;
3. NKCC1 versus apical chloride export;
4. AE4 versus net alkalinity-support flux;
5. pump demand versus sustained chloride export;
6. total basolateral chloride loading versus apical chloride export.

Do not compute every pairwise reaction coupling in the network.

### Structural-repair budget

Predeclare no more than six biologically interpretable constraint groups that could plausibly control excessive KO rescue.

First test those groups one at a time.

Do not enumerate pairs, triples or arbitrary subsets.

If no single group is sufficient and a joint analysis is mathematically necessary, perform at most one sparse L1/lexicographic joint relaxation or tightening problem across the predeclared groups. Do not enumerate support patterns.

If that still does not identify a compact explanation, stop and report structural underdetermination or missing biology. Do not expand the candidate list automatically.

### Sampling discipline

FVA and exact LP extrema are preferred.

Do not sample large feasible-state ensembles merely because the local solver can. Complete jointly feasible witness states may be generated only when needed to demonstrate attainability or illustrate degeneracy, and keep the witness set compact.

## Milestone publication rule

This task is intentionally broken into small crash-safe milestones.

At the start of every milestone:

1. `git fetch origin`;
2. switch to `analysis/task-46-physiology-constrained-model-reconstruction`;
3. reconcile only by `git merge --ff-only origin/analysis/task-46-physiology-constrained-model-reconstruction`;
4. read the current remote `CURRENT_STATUS.md`;
5. verify the previous checkpoint is remotely present before doing new work.

At the end of every milestone:

1. verify all outputs for that milestone;
2. update `CURRENT_STATUS.md`;
3. create an immutable checkpoint file named `CHECKPOINT_46H0X_<NAME>.md`;
4. write a machine-readable receipt under `output/46h/` containing inputs, commands, solver versions, checks and principal results;
5. commit only that milestone;
6. push immediately;
7. verify the remote branch contains the new commit SHA;
8. verify the working tree is clean;
9. only then begin the next milestone.

Never force-push. Never accumulate two milestones locally before pushing.

If publication fails, stop scientific work until publication is recovered.

## Milestone 46H00: recovery, evidence split and protocol freeze

Do not rerun Task 44 or Task 46G.

Confirm the branch starts from `9553d04442aeef5b0a8e0ba5f9df07f0986d7074` or a fast-forward descendant containing it.

Create/update:

- `CBM_PROTOCOL.md`;
- `CBM_EVIDENCE_SPLIT.md`;
- `output/46h/constraint_budget.json`.

The constraint budget must list:

- one base network only;
- no more than six targeted coupling questions;
- no more than six candidate structural constraint groups;
- no exhaustive subsets;
- no multidimensional sweeps;
- no new kinetics in Task 46H.

Publish checkpoint 46H00 immediately.

## Milestone 46H01: network, S matrix and provenance freeze

Build the curated transporter ledger and exact steady-state constraint system.

Deliver at minimum:

- `cbm/reaction_ledger.csv`;
- `cbm/balance_ledger.csv`;
- `cbm/stoichiometric_matrix.csv`;
- `cbm/bounds_provenance.csv`;
- `cbm/build_cbm.py` or equivalent deterministic model constructor;
- a human-readable network audit.

Verify:

- every reaction stoichiometry;
- every balance sign;
- charge consistency of each transporter;
- rank of S;
- independent-row rank after any exact redundancy removal;
- nullspace dimension;
- no impossible source/sink created by bookkeeping;
- WT/KO shared-capacity definitions are explicit but phenotype has not yet been optimised.

Do not optimise KO phenotype in this milestone.

Publish checkpoint 46H01 immediately.

## Milestone 46H02: WT feasible region

Using the frozen network, determine whether WT sustained transport is feasible.

Use the local verified VFFVA/FBA core.

For full feasible-region FVA, do not retain an arbitrary biological objective. If the local API requires an FBA anchor, use the declared objective only to obtain a feasible witness and run the FVA region with no fraction-of-optimum retention. State this explicitly.

Compute FVA for the principal transporter and output fluxes.

Report:

- whether WT is feasible;
- one complete feasible WT witness;
- full-region ranges for principal fluxes;
- which fluxes are structurally fixed, tightly bounded or weakly determined;
- active bounds at the WT secretion maximum and minimum where scientifically useful;
- whether any conclusion depends on a computational big-M.

Do not fit anything to the AE4-null phenotype.

Do not rerun reference FVA across the full scientific network merely to revalidate the engine. Task 46G already verified the solver. Use reference FVA only for debugging a discrepancy.

Publish checkpoint 46H02 immediately.

## Milestone 46H03: AE4-null feasible region under shared capacities

Set AE4 flux to zero.

Keep all non-genotype-specific transporter capacities identical to WT unless an independently sourced genotype change was frozen in 46H00.

Do not enlarge NKCC1 capacity because the old kinetic model increased NKCC1 flux.

Compute the AE4-null full feasible region and quantify:

- maximum possible apical chloride export / secretion proxy;
- minimum and maximum NKCC1 flux;
- total basolateral chloride loading range;
- pump demand range;
- alkalinity-support relationships;
- whether near-WT output remains structurally feasible;
- a complete jointly feasible witness for any claimed maximal rescue.

Distinguish clearly between transporter capacity and realised flux.

Do not yet impose a protocol-mismatched NKCC flux-equality constraint.

Publish checkpoint 46H03 immediately.

## Milestone 46H04: paired WT/KO NKCC evidence constraint

Build the smallest Task-46-specific paired LP formulation needed to express cross-genotype constraints. Keep this code outside the frozen vendor directory.

First distinguish:

- shared NKCC1 capacity/abundance across genotypes;
- realised NKCC1 flux under normal stimulation;
- the isolated HCO3-free NKCC assay observable.

Do not equate these automatically.

Use the strongest numerical WT/KO NKCC relation actually supported by the preserved experimental evidence. If a numerical tolerance is not justified, do not invent epsilon.

If only qualitative equality/no-difference is supported under a distinct assay protocol, implement the corresponding assay-matched structural diagnostic separately and state what cannot be transferred to normal stimulation.

Then determine the maximum KO secretion proxy compatible with the defensible NKCC evidence.

Report whether the old-model style ~90% chloride replacement is:

- structurally required;
- merely allowed by loose constraints;
- or incompatible with the independent NKCC evidence.

Publish checkpoint 46H04 immediately.

## Milestone 46H05: targeted coupling and redundancy analysis

Use only the predeclared coupling questions from `output/46h/constraint_budget.json`.

For each, determine rigorous feasible ratio/range or linear coupling information using LP extrema and exact balance identities where available.

Do not calculate an all-pairs coupling matrix.

Specifically determine whether AE4 and NKCC1 remain substitutable once the retained Cl, Na, K, TIC/TA, pump and apical-export constraints are enforced.

Separate:

- exact stoichiometric identities;
- LP-derived interval relationships;
- conclusions that depend on numerical bounds.

Publish checkpoint 46H05 immediately.

## Milestone 46H06: minimal structural discriminator / repair analysis

Use only the predeclared maximum six structural constraint groups.

The goal is not to fit 35%.

The goal is to identify the smallest biologically interpretable constraint change that prevents implausibly strong KO rescue while preserving WT feasibility.

Proceed in this order:

1. test each predeclared constraint group individually;
2. quantify how much each one reduces the maximal KO secretion proxy while preserving WT feasibility;
3. compare any required threshold to independent evidence;
4. reject a proposed restriction if it requires an unsupported value;
5. only if no single group is sufficient and there is a clear mathematical reason, run one sparse L1/lexicographic joint analysis across the same predeclared groups.

Do not enumerate pairs or subsets.

Do not tune to the central 35% phenotype value.

If the feasible family remains broad enough to contain both near-WT and strongly deficient KO states, report structural non-identifiability rather than selecting a convenient point.

If no defensible compact restriction excludes excessive rescue, state that the steady single-cell CBM is insufficient and identify the missing class of biology/measurement needed.

Publish checkpoint 46H06 immediately.

## Milestone 46H07: CBM decision report and hard stop

Create `CBM_DECISION_REPORT.md` and machine-readable summary.

Answer explicitly:

1. Is WT sustained transport feasible?
2. Is near-WT AE4-null output feasible under shared capacities?
3. What is the maximal AE4-null output under the defensible NKCC evidence?
4. Is strong NKCC rescue required, optional or excluded?
5. Are AE4 and NKCC1 structurally interchangeable? Under exactly which assumptions?
6. Which carbon/alkalinity constraints materially limit AE4 flux?
7. Which one or two constraints dominate the excessive-rescue problem?
8. Is the CBM structurally predictive or underdetermined?
9. What is the smallest biologically defensible next kinetic reconstruction?
10. Which one or two kinetic blocks, and only those blocks, should Task 46 rebuild next?

Include a compact table separating:

- exact mathematical identities;
- LP/FVA numerical conclusions;
- evidence-dependent conclusions;
- unresolved biological questions.

Publish checkpoint 46H07 immediately.

Then STOP.

Do not start NKCC1 kinetic model discovery, parameter calibration or ODE reconstruction in this task, even if time remains. The next kinetic task must be written after the CBM decision report is inspected.

## Failure is an acceptable result

If the curated CBM cannot distinguish the mechanism, that is a valid scientific result.

Do not rescue the task by adding transporter variants, more free bounds or additional mechanisms.

Publish the negative result, explain exactly why the feasible region remains underdetermined, and identify the smallest additional experiment or biological mechanism required to resolve it.

## Final completion criteria

Task 46H is complete only when:

- checkpoints 46H00 through 46H07 are each present remotely as separate commits;
- every checkpoint has a reproducible machine-readable receipt;
- `CURRENT_STATUS.md` points to 46H07;
- `CBM_DECISION_REPORT.md` exists;
- the branch is clean;
- `esig626/CarbonScope` has not been modified;
- no new kinetic model has been started;
- no combinatorial model/parameter search was performed.
