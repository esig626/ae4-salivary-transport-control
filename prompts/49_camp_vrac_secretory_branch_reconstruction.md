# Task 49: source-backed cAMP/VRAC secretory branch reconstruction

## Mission

Task 49 is not another mechanism search. It tests one specific, independently documented piece of mouse submandibular-gland physiology that is absent from the current model and has not been tested in the repository:

> **Muscarinic/Ca2+-dependent secretion uses apical TMEM16A/CaCC, whereas beta-adrenergic/cAMP-dependent secretion is TMEM16A-independent and uses a distinct volume-regulated anion conductance, with beta-adrenergic stimulation also activating NKCC1.**

The current Task 48 model has only the Ca/TMEM16A apical chloride outlet. Beta input affects AE4, but after AE4 deletion it disappears from the core state equations. This is why the inherited AE4-null model forces CCh and CCh+IPR to be the same experiment. Task 49 asks whether restoring the independently established cAMP secretory branch removes that structural error and, without fitting the AE4-null saliva phenotype, produces the observed stimulus dependence and sustained loss.

This is one of four remaining scientific execution opportunities. Do not spend it rediscovering results already in the repository.

## Repository and branch

Scientific base: merged Tasks 47-48 on `main` at

`5ca70c61e8ccc49f084bac1b91ea05c612d197e8`

Work only on:

`analysis/task-49-camp-vrac-secretory-branch-reconstruction`

Read `AGENTS.md` first.

Then read, in this order:

1. `analysis/49_camp_vrac_secretory_branch_reconstruction/CODEX_START_HERE.md`
2. `analysis/49_camp_vrac_secretory_branch_reconstruction/EXCLUSION_LEDGER.md`
3. `analysis/49_camp_vrac_secretory_branch_reconstruction/SOURCE_SEED.md`
4. `analysis/48_joint_experimental_constraint_reconstruction/JOINT_EXPERIMENTAL_CONSTRAINT_REPORT.md`
5. `analysis/48_joint_experimental_constraint_reconstruction/PROTOCOL_MAP.md`
6. `analysis/48_joint_experimental_constraint_reconstruction/OBSERVATION_OPERATORS.md`
7. `analysis/48_joint_experimental_constraint_reconstruction/constraints.json`
8. `analysis/41_ae4_loss_algebraic_design/final_answer.md`
9. `analysis/42_catalan_2025_ae4_mechanism_classes/final_answer.md`
10. `analysis/39_palk_nkcc1_full_validation/final_answer.md`
11. `analysis/37_wt_nbc_validation/final_answer.md`
12. the active source in `src/modern_full_model/`, especially `nkcc_stimulation.py`, membrane/current closure, chloride balances, water/volume equations, and the Task 48 protocol layer.

Do not treat older historical equation skeletons as the active implementation.

## What has already been eliminated

The exclusion ledger is binding. In particular, do not revisit:

- AE4-only kinetic/stoichiometric families;
- alternative Na/K routing as the solution;
- the seven Catalan 2025 source classes;
- replacing NKCC1 with the Palk/Benjamin core law;
- a hard NKCC cap;
- Na/K-pump capacity/coupling as the missing phenotype mechanism;
- a K-channel/pump recycling ceiling;
- generic NHE1 replacement;
- calcium-amplitude escalation;
- AE4-expression-dependent TMEM16A/CaCC recruitment from Task 41;
- CBM resource caps;
- any gland-scale/spatial mechanism.

Task 41 is a target-selected constructive result only. It showed that changing apical chloride export can create the required magnitude, but its AE4-dependent CaCC coupling is not biologically validated and has the wrong timing/state phenotype. Do not resurrect it. Task 49 tests the independently documented **different cAMP apical channel pathway** instead.

## Primary independent evidence

Verify every numerical claim before inference.

### Catalan et al., PNAS 2015

`A fluid secretion pathway unmasked by acinar-specific Tmem16A gene ablation in the adult mouse salivary gland`

DOI `10.1073/pnas.1415739112`, PMID `25646474`, PMCID `PMC4343136`.

The source establishes in adult mouse SMG that:

- TMEM16A is obligatory for muscarinic/Ca2+-dependent secretion;
- beta-adrenergic IPR secretion persists in acinar-specific TMEM16A knockout;
- IPR secretion is also independent of CFTR and ClC-2;
- IPR secretion requires extracellular chloride and bicarbonate;
- DCPIB and NPPB strongly suppress IPR secretion, supporting a VRAC-like apical conductance;
- IPR produces a slowly developing acinar-cell swelling and a similarly delayed secretory response;
- IPR activates a DCPIB-sensitive chloride conductance.

Extract the quantitative source data, including figure-derived data only where necessary. If figure digitisation is required, save the figure, axis calibration, digitised points, uncertainty and provenance before fitting anything. Do not eyeball values into code.

### Pena-Munzenmayer et al., JBC 2015

DOI `10.1074/jbc.M114.612895`, PMCID `PMC4409235`.

This remains the AE4 phenotype dataset. For Task 49 it is a **prediction/validation dataset**, not the source used to choose the new cAMP-channel parameters. The critical held-out features are:

- AE4 KO secretion initially comparable for about 2-3 min, then sustained loss, total `35 +/- 4.7%` lower at 10 min under CCh+IPR;
- AE4 KO CCh-only chloride uptake `2.30 +/- 0.10 x10^-3 s^-1`;
- AE4 KO CCh+IPR uptake `0.90 +/- 0.09 x10^-3 s^-1`;
- AE4 KO IPR-only uptake `0.20 +/- 0.03 x10^-3 s^-1`;
- WT counterparts and all resting/AE2/NHE/bicarbonate constraints already frozen in Task 48.

Do not use the 35% deficit, the AE4-null combined-stimulus uptake, or the early/late AE4-null secretion trace to choose a new parameter.

### Pena-Munzenmayer et al., AJP GI 2021

DOI `10.1152/ajpgi.00145.2021`, PMCID `PMC8887885`.

This independently establishes beta/cAMP/PKA activation of AE4 and explicitly notes that TMEM16A is involved in Ca-dependent but not cAMP-dependent secretion and that the cAMP apical conductance is likely VRAC-like. Existing AE4 PKA regulation is already present. Do not add it again.

### NKCC beta regulation

Primary salivary-cell studies report strong beta-adrenergic/cAMP-dependent upregulation/phosphorylation of NKCC1, including Kurihara et al. 2002, DOI `10.1152/ajpcell.00352.2001`, PMID `11880270`. These studies include rat parotid data and are not automatically quantitative mouse-SMG constants. Use them to establish the **sign and existence** of the beta arm. Only use a numerical beta-NKCC gain if a source-compatible value or defensible cross-source envelope is established independently of the AE4-null phenotype.

## The precise architectural defect being tested

The active model currently implements essentially

`J_apical_Cl = J_TMEM16A(Ca, state)`

and beta input affects AE4 but not an independent apical cAMP channel. In AE4 knockout, beta therefore vanishes from the core vector field.

Task 49 must replace the apical chloride-output architecture by the source-backed parallel form

`J_apical_Cl = J_TMEM16A + J_VRAC`

where:

- `J_TMEM16A` is the existing Ca/CCh branch, unchanged;
- `J_VRAC` is a distinct beta/cAMP-associated, volume-sensitive apical chloride conductance;
- VRAC contributes to the intracellular Cl balance, luminal Cl balance and apical electrical closure exactly like a real chloride current, not as a secretion multiplier;
- VRAC has no direct AE4-expression multiplier;
- VRAC is inactive at the declared resting state;
- CCh-only must nest the Task 48 CCh-only model to numerical precision when the VRAC gate is inactive.

The simplest admissible volume law is preferred. Start from a one-scale local swelling law such as an effective conductance proportional to positive fractional swelling relative to the genotype-specific resting volume. Do not introduce a Hill coefficient, arbitrary delay, knockout multiplier or extra threshold unless the independent source data require it. The existing volume dynamics should supply the delay if possible.

The exact current sign and chloride source terms must be derived and independently checked against the repository conventions before integration.

## Beta-dependent NKCC1 arm

The current `N1_CCH_ALGEBRAIC` gain is driven only by the calcium/CCh arm and was a WT-only effective construction. Independent physiology supports beta/cAMP activation of NKCC1 as well.

Extend the regulator minimally so NKCC1 can respond to both Ca/CCh and beta/cAMP inputs while preserving:

- the source-fixed Palk/Benjamin concentration-response core;
- 1 Na : 1 K : 2 Cl stoichiometry;
- reversal and all existing conservation identities;
- the existing CCh/Ca arm at beta=0;
- the resting multiplier exactly one.

Use the smallest no-interaction form unless independent evidence requires an interaction. For example, an additive activation above baseline is preferable to introducing an unmeasured Ca×beta synergy parameter.

The beta gain must be fixed or constrained from independent beta/NKCC evidence and/or the independent WT IPR-only reconstruction. It must **not** be selected from AE4-null combined-stimulus secretion or uptake.

If the independent data cannot identify both the beta-NKCC gain and VRAC scale sufficiently to make a prediction, report non-identifiability and stop rather than using the AE4 phenotype to close the system.

## Why this mechanism is scientifically capable of addressing the Task 48 contradiction

Do not assume success, but test the actual causal prediction.

Under CCh alone, the new cAMP branch is off, so the Task 48 CCh trajectory should remain unchanged.

Under CCh+IPR, beta can increase basolateral chloride loading while the distinct VRAC-like apical conductance adds a sustained chloride exit. Net SPQ chloride reuptake can therefore be lower even when basolateral loading is activated, because the measured concentration slope is net cellular chloride change, not a transporter flux. In WT, AE4's beta-activated chloride loading may help offset the extra efflux. In AE4 KO that contribution is absent, so intracellular stores and the remaining loaders must support both apical branches. This provides a source-backed way to produce initially preserved secretion followed by a sustained deficit without imposing an AE4-dependent CaCC penalty.

This paragraph is a predeclared causal hypothesis, not a result.

## Task 49A: evidence freeze and no-repeat audit

Use up to five parallel workers, but only one orchestrator may modify the canonical branch or accept scientific conclusions.

Suggested workers:

1. **Repository-history worker:** independently verify the exclusion ledger and confirm VRAC/DCPIB/LRRC8 has not already been implemented/tested under another name.
2. **PNAS evidence worker:** extract/digitise the independent IPR/TMEM16A/VRAC/swelling/current constraints from DOI `10.1073/pnas.1415739112`.
3. **NKCC regulation worker:** reconstruct the beta/cAMP-NKCC evidence and determine what quantitative constraint, if any, is defensible for mouse SMG.
4. **Conservation/electrical worker:** derive the minimal VRAC source/current insertion and verify charge, chloride and water/electrical bookkeeping independently.
5. **Validation-design worker:** predeclare calibration data, held-out AE4 data, success/failure gates and identifiability checks.

Publish 49A before any production equation is changed.

Required outputs include a machine-readable source ledger, exclusion ledger verification, and a clear declaration of which quantities may be calibrated and from which non-AE4-null observations.

## Task 49B: minimal source-backed implementation

Implement only:

1. one apical VRAC-like chloride conductance driven by the minimum source-supported beta/volume logic;
2. one beta input to the existing NKCC regulatory multiplier if supported by 49A.

Do not change AE4 stoichiometry/kinetics, NHE1, NBC, pump, K channels, CaCC, paracellular transport, water permeabilities, geometry or the calcium stimulus.

Do not change the Task 48 genotype resting solver. Task 49 does **not** repair the chronic resting-state mismatch. That is a separate unresolved problem and must remain visible.

Unit tests must establish:

- exact/near-exact REST nesting;
- CCh-only nesting when beta/VRAC are off;
- zero VRAC under its off condition;
- correct Nernst/current reversal and chloride source sign;
- exact charge/current closure with VRAC included;
- beta-NKCC arm does not move NKCC thermodynamic reversal;
- disabling TMEM16A does not disable VRAC;
- disabling VRAC does not alter TMEM16A or AE4;
- no direct dependence of VRAC on AE4 expression.

Publish 49B immediately.

## Task 49C: calibrate/validate only against independent cAMP-pathway data

Before evaluating the AE4-null JBC 2015 secretion phenotype, use the PNAS 2015 cAMP-pathway data and independently supported NKCC beta evidence to fix the new parameter(s).

Use the smallest identifiable parameter set. Prefer analytic or one-/two-dimensional deterministic inference over a broad optimiser. Do not perform a large grid or mechanism family search.

At minimum test whether the fixed model reproduces the independent qualitative/quantitative facts:

- IPR-only secretion persists when TMEM16A is absent;
- VRAC blockade markedly reduces IPR-only secretion;
- IPR induces the correct direction and approximate time scale of cell swelling and delayed secretion;
- the new branch is chloride dependent where the protocol can be represented faithfully;
- CCh-only secretion remains TMEM16A dependent and unchanged from its parent model.

If the independent source data cannot support a unique or bounded reconstruction, stop here and report that before touching the AE4-null phenotype.

Publish 49C immediately with all chosen parameters, provenance and uncertainty/sensitivity envelope frozen.

## Task 49D: immutable pre-AE4 prediction checkpoint

With the cAMP branch fixed exclusively from 49C evidence, generate and publish the predictions for:

- WT: REST, CCh, IPR, CCh+IPR;
- AE4 KO: REST, CCh, IPR, CCh+IPR;
- AE2 KO/control where executable and relevant.

The AE4-null 35% secretion value and combined-stimulus uptake data must not select or alter a parameter before this commit is remote.

Record at least:

- TMEM16A and VRAC chloride flux/current separately;
- NKCC1, AE4 and AE2 chloride loading;
- intracellular Na/K/Cl, pH and volume;
- apical/basolateral voltage;
- water flow and cumulative secretion;
- VRAC activation variable or volume gate;
- early and late secretion windows matching the Task 48 definitions;
- exact CCh-versus-combined knockout trajectory difference.

Publish 49D and verify the remote SHA before reveal.

## Task 49E: held-out AE4/Ae2 validation

Only after 49D is remote, compare the fixed predictions against the Task 48 frozen JBC 2015 constraints.

The important questions are:

1. Is the exact AE4-null CCh = CCh+IPR structural invariance removed for the correct source-backed reason?
2. Is AE4-null CCh-only uptake/secretion preserved while beta-containing stimulation is selectively impaired?
3. Is IPR-only AE4-null uptake reduced in the correct direction?
4. Does combined AE4-null secretion remain near WT early and diverge later?
5. What is the 10-min cumulative AE4-null secretion deficit, without retuning?
6. Does AE2 KO remain essentially secretion neutral?
7. Does the model avoid the Task 41 pathology of an immediate large deficit and AE4-expression-dependent TMEM16A recruitment?
8. What resting Cl/pH mismatch remains unchanged from Task 48?

A successful stimulated reconstruction does not erase the resting-state failure. Report the two separately.

No post-reveal retuning.

Publish 49E immediately.

## Task 49F: article-facing decision

Create

`analysis/49_camp_vrac_secretory_branch_reconstruction/CAMP_VRAC_RECONSTRUCTION_REPORT.md`

The report must place Task 49 in the actual repository narrative:

- Tasks 12/38-42 established AE4-only and stoichiometric insufficiency and excessive NKCC rescue;
- Task 46 showed structural rescue is feasible;
- Task 47 rejected a binding K-recycling ceiling;
- Task 48 showed the inherited model mathematically erases beta input from the AE4-null core;
- independent mouse SMG physiology identifies a separate beta/cAMP apical chloride pathway absent from the model;
- Task 49 tests that pathway prospectively, with parameters fixed before the AE4-null phenotype comparison.

State whether the cAMP/VRAC reconstruction succeeds, partially succeeds, or fails, and exactly which experimental constraints remain wrong. Do not call a partial stimulated success a complete model while the chronic resting Cl/pH phenotype remains unresolved.

## Absolute prohibitions

Do not:

- fit the 35% AE4-null saliva deficit;
- add an AE4-expression multiplier to VRAC;
- revive Task 41 AE4-dependent CaCC recruitment;
- lower NKCC capacity to force a deficit;
- change K-channel or pump capacities;
- alter NHE1 or NBC in this task;
- change AE4 source stoichiometry/routing;
- alter calcium amplitude;
- introduce CFTR or ClC-2 as candidate branches when the primary PNAS experiment already excludes them from IPR secretion;
- perform a mechanism family search;
- fit protocol-specific SPQ scale factors to erase residuals;
- silently approximate blocked zero-substrate assays and call them exact;
- modify manuscript files or previous numbered analyses.

## Crash-safe publication

Every milestone 49A through 49F must be committed and pushed immediately. Each must contain:

- updated `CURRENT_STATUS.md`;
- a machine-readable receipt under `analysis/49_camp_vrac_secretory_branch_reconstruction/output/checkpoints/`;
- parameter/source hashes where relevant;
- remote commit SHA verification before dependent work starts.

Never force-push or rewrite a published milestone. If the workspace dies, resume from the latest verified remote checkpoint rather than rerunning completed work.

Only the orchestrator may integrate worker outputs, run authoritative inference/production trajectories, accept scientific corrections, commit or push to the canonical Task 49 branch.

Stop after 49F. Do not automatically begin a Task 50 mechanism search.
