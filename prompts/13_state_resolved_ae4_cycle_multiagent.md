# Task 13. State-resolved AE4 cycle and whole-cell phenotype reconstruction

Work in repository `esig626/ae4-salivary-transport-control` on branch `codex/state-resolved-ae4-cycle`.

Read `AGENTS.md` first. Then read all prior reconstruction outputs that materially constrain this task, including at minimum:

- `analysis/10_identifiability_discrimination/`
- `analysis/11_forensic_reconstruction/`
- `analysis/12_ae4_mechanism_reconstruction/`
- `model/SOURCE_MAP.md`
- `model/equations.md`
- `model/parameters.md`
- `model/fluxes.md`
- `docs/PROVENANCE.md`
- `docs/DECISIONS.md`
- `docs/RESULTS_LEDGER.md`

Do not edit anything under `archive/`.

## Core scientific problem

Task 12 tested a broad family of **coarse-grained AE4 flux laws** on the best internally coherent historical seven-state chassis and did not recover the held-out AE4-knockout secretion phenotype. That result must not be summarized as “AE4 itself is insufficient.” The biology already establishes that AE4 is functionally important for salivary secretion.

The unresolved question is instead:

> What mechanistic structure is missing between molecular AE4 transport and the whole-cell secretion phenotype?

There are two live possibilities that this task must distinguish.

1. The AE4 transporter was represented too coarsely. A state-resolved, cation-specific alternating-access cycle may produce network behavior that cannot be captured by pooled Na+K or parallel scalar flux laws.
2. The surrounding whole-cell transport model omits a source-supported cation, chloride, or acid-base homeostasis module required for loss of AE4 to propagate into the observed phenotype.

Do not choose between these possibilities in advance.

The goal is to **solve the reconstruction problem**, not to produce another negative scan.

## Primary experimental and structural sources

Before modeling, obtain and audit the primary sources. At minimum use:

1. Peña-Münzenmayer et al., 2015, JBC, DOI `10.1074/jbc.M114.612895`.
   - AE4-null salivary secretion phenotype.
   - Time dependence of the stimulated secretion deficit.
   - Resting intracellular chloride and related phenotypes.
   - AE2-null comparisons.
   - NKCC1 activity/expression observations relevant to compensation.

2. Peña-Münzenmayer et al., 2016, JGP, DOI `10.1085/jgp.201611571`.
   - Direct Na transport.
   - Direct K transport.
   - Electroneutrality and voltage independence.
   - EC50 and Hill fits for Na and K.
   - Alternative electroneutral stoichiometries considered in the paper.
   - Thermodynamic arguments showing why Na-only transport is problematic physiologically.
   - The proposal that K recycling through AE4 may support Na/K ATPase under some conditions.

3. AE4 PKA-regulation paper, DOI `10.1152/ajpgi.00145.2021`.
   - Beta/cAMP/PKA activation.
   - H89 sensitivity.
   - S173-associated regulation.
   - Quantitative fold changes and, if recoverable, activation kinetics or time scales.

4. Catalán et al., 2025, AJP Cell Physiology, DOI `10.1152/ajpcell.00346.2024`.
   - TM3-TM10 cation coordination site.
   - T448, T756, D709, T713 and related residues.
   - The T756A-T448I double mutant losing Na-dependent transport while retaining substantial K-dependent transport.
   - Evidence that Na and K coordination are not functionally identical.
   - The proposed state-level cycle hypotheses involving Na and K at different stages of the transport cycle.
   - Candidate HCO3 versus carbonate transport interpretations.

5. The 2018 apical Na/K ATPase and K-channel salivary paper already in the repository literature.
   - Apical versus basolateral pump distribution.
   - Apical versus basolateral K conductance.
   - Independent constraints on the cation-current architecture.

6. The 2018 AE4 mathematical paper and the archived code as historical model evidence only.

If the literature cited by these sources contains a necessary mechanistic constraint, obtain the primary source rather than relying on a secondary statement.

Create `analysis/13_state_resolved_ae4/evidence_freeze.md` with a table containing:

- quantity or qualitative constraint;
- exact source;
- direct measurement versus inference versus hypothesis;
- numerical value and uncertainty if available;
- whether it may be used for calibration;
- whether it is held out for validation.

Do not begin model fitting until this evidence table is frozen.

## Held-out validation hierarchy

The main secretion phenotype must remain an out-of-sample prediction.

### Stage A holdout

At the start of the task, withhold **all AE4-null whole-gland outcomes** from calibration. In particular do not fit to:

- total stimulated saliva reduction;
- the time course of the AE4-null secretion deficit;
- AE4-null intracellular chloride;
- AE4-null pH;
- AE4-null cell volume, if directly measured.

Build the transporter and whole-cell WT model without these targets.

### Stage B holdout

If every source-supported state-resolved AE4 model fails on the current chassis, then AE4-null chloride, pH, and volume may be released as **localization constraints** for reconstructing the missing whole-cell module.

Even in Stage B, the following remain strictly held out:

- total AE4-null secretion magnitude;
- early versus late stimulated secretion pattern;
- any integrated or endpoint secretion ratio.

A model that uses the target secretion ratio during optimization has not explained that phenotype.

## Multi-agent organization

Run this as a genuine multi-agent investigation. At minimum instantiate the following independent roles.

### Agent A. Experimental evidence curator

Responsibilities:

- build the evidence freeze;
- verify every phenotype number against the primary paper;
- extract time-resolved AE4-null secretion data, digitizing figures when necessary and documenting digitization uncertainty;
- separate WT, AE2-null, AE4-null, and double-null evidence;
- extract Na/K dose-response and mutant constraints from 2016 and 2025;
- identify which data correspond to salivary acinar cells versus heterologous systems;
- flag species, temperature, bath composition, and stimulation differences.

Agent A must not design the model.

### Agent B. AE4 structural-cycle modeler

Responsibilities:

- translate the 2016 and 2025 AE4 evidence into the smallest plausible state-resolved transport-cycle families;
- distinguish what is required by experiment from what is merely proposed;
- derive candidate alternating-access state graphs;
- incorporate Na- and K-specific coordination differences suggested by mutagenesis and MD;
- identify which microscopic parameters are actually identifiable from available transporter-level data;
- reduce unidentifiable microscopic parameters into defensible lumped combinations.

Agent B must not use the whole-gland knockout phenotype.

### Agent C. Thermodynamic and mathematical auditor

Responsibilities:

- enforce electroneutrality;
- derive cycle affinities from chemical potentials;
- enforce detailed balance or local detailed-balance constraints where applicable;
- reject cycles that violate thermodynamic direction under source conditions;
- derive exact invariants, limiting cases, and source vectors;
- identify when two apparently different Markov schemes are observationally equivalent at steady state;
- derive a minimal state reduction when possible without destroying cation-specific behavior.

### Agent D. Whole-cell chassis reconstruction agent

Responsibilities:

- reconstruct the most defensible WT salivary chassis from published and historical evidence;
- document all differences between the seven-state historical chassis and the best source-supported modern topology;
- preserve AE2 explicitly when required for the comparison;
- evaluate apical/basolateral Na/K ATPase and K-channel architecture;
- audit NHE1, CO2/HCO3 buffering, NKCC1, membrane potential, water flux, and luminal outflow;
- build modular code so each chassis addition can be switched on or off independently.

### Agent E. Numerical inference and continuation agent

Responsibilities:

- implement candidate AE4 cycles and chassis modules independently from the archived code;
- calibrate only to allowed evidence;
- solve resting states robustly using multistart and continuation;
- preserve branches instead of selecting convenient roots silently;
- integrate stimulation protocols with stiff solvers;
- perform nested model comparison with explicit penalties for unnecessary parameters;
- maintain machine-readable ledgers of all accepted and rejected models.

### Agent F. Residual localization agent

Responsibilities:

- if Stage A fails, use the AE4-null vector field and Stage B ionic phenotypes to calculate the **required missing balance/source direction**;
- express the discrepancy as a vector in Na, K, Cl, HCO3/H, membrane-current, and water/volume balance space;
- project that vector onto candidate transporter-module stoichiometric signatures;
- solve sparse or minimum-complexity source reconstruction problems;
- distinguish exact span results from approximate numerical alignments;
- identify the smallest known module or module pair capable of supplying the missing direction.

### Agent G. Adversarial scientific auditor

Responsibilities:

- attempt to falsify every claimed successful mechanism;
- search for target leakage;
- test whether success depends on one arbitrary bound, hidden scaling factor, or solver choice;
- check whether a supposedly state-resolved effect could be reproduced by a simpler scalar flux law;
- check whether any added chassis module is independently justified rather than introduced because it fixes the phenotype;
- verify that claims do not discredit the historical article beyond the evidence.

### Agent H. Independent numerical reproducer

Responsibilities:

- reproduce final accepted calculations from the saved equations and parameter files without reusing Agent E's intermediate state;
- use an independent solver stack or implementation where feasible;
- compare steady states, trajectories, conservation errors, and held-out ratios;
- report discrepancies before synthesis.

### Lead synthesis agent

The lead agent must reconcile disagreements explicitly. It may not simply choose the most optimistic agent output.

## Round 0. Evidence freeze and exact question decomposition

Create:

- `analysis/13_state_resolved_ae4/evidence_freeze.md`
- `analysis/13_state_resolved_ae4/question_tree.md`

The question tree must separate:

1. transporter-level adequacy;
2. WT whole-cell adequacy;
3. AE2-null adequacy;
4. AE4-null ionic adequacy;
5. AE4-null secretion adequacy;
6. time-resolved stimulation adequacy.

A candidate must not jump directly from 1 to 5.

Verify the often-quoted “30%” versus “35 +/- 4.7%” secretion reduction from the primary 2015 data and record exactly which quantity each number refers to. Do not propagate an approximate value without source context.

## Round 1. Build state-resolved AE4 candidate families

The previous task already tested scalar and coarse-grained families C1-C8. Do not merely repeat them.

Construct new candidate families that can express **different Na and K roles within one alternating-access cycle**.

At minimum test the following mechanistic classes, with exact state graphs documented.

### SR1. Minimal alternating-access 1:1:2 cycle

Represent outward- and inward-facing states explicitly. Allow one cation species per cycle but separate Na and K microscopic binding affinities/rates. Do not pool concentrations as `[Na]+[K]` unless this reduction is derived from the state model.

### SR2. Parallel Na and K alternating-access cycles sharing conformational states

Na and K may use the same transporter but with cation-specific binding/release constants and shared conformational transitions. This is different from adding two independent scalar flux laws.

### SR3. Asymmetric sequential Na/K cycle

Test the 2025-motivated possibility that extracellular Na and intracellular K participate at different stages of the transport cycle. Construct the smallest thermodynamically consistent scheme compatible with the proposed `Cl:Na:HCO3:K` style cycle. Explicitly determine whether a `1:1:1:1` cycle is electroneutral and thermodynamically feasible under salivary concentrations.

### SR4. Carbonate variant

If primary evidence supports retaining it, implement the smallest explicit carbonate-containing state model needed to test the 2025 proposal. Add carbonate chemistry explicitly. Do not smuggle carbonate into a bicarbonate variable.

### SR5. Cation coordination plus regulation

Allow cation-specific binding to alter conformational transition rates, consistent with cooperative Hill behavior or allosteric regulation, while retaining direct cation transport because the 2016 experiments require it.

### SR6. PKA-regulated state transitions

Instead of applying one scalar activity multiplier, test whether source-supported phosphorylation can alter one or more state-transition rates. Keep the number of new regulatory parameters minimal and constrained by the 2021 data.

For every candidate family:

- derive the state graph;
- specify transported stoichiometry per completed cycle;
- write the cycle affinity;
- identify the equilibrium/reversal condition;
- identify cation-specific parameters;
- identify parameters constrained by 2016 or 2025 experiments;
- identify parameters that remain free;
- identify whether the scheme reduces to any Task 12 coarse-grained law in a limit.

Save these in `analysis/13_state_resolved_ae4/state_cycle_models.md`.

## Round 2. Transporter-level calibration without whole-gland knockout data

Fit or constrain the SR families using transporter-level evidence only.

At minimum require consistency with:

- Na dependence;
- K dependence;
- direct Na transport;
- direct K transport;
- negligible voltage dependence/electroneutrality;
- Na and K EC50 values;
- Na and K Hill slopes;
- broad monovalent cation permissiveness where relevant;
- 2025 mutant hierarchy, especially the Na-versus-K behavior of T756A-T448I;
- transporter reversal or thermodynamic direction under reported ionic conditions;
- PKA regulation when included.

Do not overfit every curve with a microscopic parameter. Where the data identify only effective rates or affinities, preserve that nonidentifiability explicitly.

Perform model comparison among SR1-SR6 using predictive adequacy and complexity penalties. A more complicated scheme must earn its parameters by explaining a constraint that a simpler scheme cannot.

Create:

- `analysis/13_state_resolved_ae4/transporter_fit.md`
- `results/13_state_resolved_ae4/transporter_fit_summary.csv`

## Round 3. Integrate into the fixed best current chassis

Take every transporter-level-surviving SR family and embed it into the best current internally coherent seven-state chassis **without changing non-AE4 components yet**.

Recalibrate only allowed WT/rest quantities. Preserve the held-out AE4-null data.

For each candidate:

- locate all WT resting roots on a documented domain;
- enforce source-supported resting Cl, Na, K, pH, and volume ranges;
- preserve a negligible AE2-null secretion phenotype;
- verify mass balance, electroneutrality, membrane-current closure, positivity, and thermodynamic signs;
- implement the source-supported stimulation protocol as closely as possible;
- if PKA dynamics are source-constrained, use them without knockout tuning.

Only after freezing the candidate reveal the Stage A AE4-null data.

Record:

- AE4-null ionic predictions;
- AE4-null secretion magnitude;
- early versus late secretion behavior;
- whether the observed delayed divergence is captured;
- which state-resolved feature caused any improvement relative to Task 12.

Create `analysis/13_state_resolved_ae4/stage_a_results.md`.

### Stage A success criterion

A state-resolved AE4 mechanism is provisionally successful only if it predicts the secretion phenotype without tuning to AE4-null outputs and does not destroy independent WT and AE2 constraints.

If one or more SR mechanisms succeeds, proceed to Rounds 7 and 8. Do not unnecessarily reconstruct the chassis.

If all fail, proceed automatically to Round 4. **Do not stop.**

## Round 4. Release AE4-null ionic data and localize the missing source direction

If Stage A fails, release AE4-null Cl, pH, and volume as localization constraints while keeping AE4-null secretion strictly held out.

Construct a mathematically explicit discrepancy analysis.

For the AE4-null model, determine what additional steady or stimulus-dependent source terms would be required to move the model toward the observed AE4-null ionic state.

Represent the required correction in balance coordinates such as:

`(Na, K, Cl, HCO3/H, membrane current, volume/water)`.

Then form stoichiometric/signature vectors for independently supported modules:

- Na/K ATPase;
- apical Na/K ATPase contribution;
- basolateral Na/K ATPase contribution;
- apical K conductance;
- basolateral K conductance;
- NHE1;
- NKCC1;
- CO2/HCO3 buffer/transport;
- paracellular Na/K pathways;
- water/luminal outflow coupling.

Use exact linear span/rank arguments where possible. When nonlinear state feedback matters, use local Jacobian/source-response directions and verify with re-solved equilibria.

Solve sparse/minimum-complexity reconstruction problems such as:

- one-module correction;
- two-module correction;
- source-constrained module pair;
- minimum norm correction subject to physiological sign restrictions.

The goal is not to fit the secretion ratio. The goal is to identify what network direction is missing from the AE4-null ionic state.

Create:

- `analysis/13_state_resolved_ae4/missing_direction.md`
- `results/13_state_resolved_ae4/module_projection.csv`

## Round 5. Build the smallest independently supported chassis extension

Starting from the Round 4 localization, add modules in order of evidential support and minimality.

### Priority module family A. Cation-current architecture

Use the experimental evidence for apical and basolateral Na/K ATPase and K-channel localization/activity. Test whether distributing pump and K conductance across membranes changes the response to AE4 loss while preserving WT physiology.

Do not assume pump redistribution changes total secretion simply because it changes current balance. Demonstrate the mechanism.

### Priority module family B. Acid-base coupling

Audit NHE1, CO2/HCO3 chemistry, buffering, and any fast-equilibrium reductions. Determine whether the AE4-null pH error implies missing acid-base flux rather than cation current alone.

### Priority module family C. Coupled cation plus acid-base reconstruction

If neither A nor B alone spans the required correction, test the smallest pair or coupled module predicted by Round 4.

### NKCC1 restriction

Do not use arbitrary NKCC1 upregulation as a rescue if the primary knockout study reports no corresponding increase. NKCC1 may vary only within independently supported physiological or stimulation ranges.

For each added module:

1. cite the independent source;
2. state the new equations;
3. state every new parameter;
4. identify what data calibrate those parameters;
5. show which previous residual or phenotype mismatch is repaired;
6. show whether the improvement survives removal of unnecessary parameters.

Create `analysis/13_state_resolved_ae4/chassis_reconstruction.md`.

## Round 6. Dynamic stimulation reconstruction

The AE4-null secretion phenotype is not merely one endpoint. The primary experiment reports a temporal structure in which secretion is initially comparable and the deficit becomes sustained later.

Treat this as an important mechanistic constraint.

Reconstruct the stimulation protocol using source-supported Ca and beta/cAMP/PKA inputs. Investigate whether the delayed separation can arise from:

- delayed PKA activation of AE4;
- intracellular chloride depletion dynamics;
- changing cation gradients;
- Na/K ATPase and K-channel recruitment;
- acid-base adaptation;
- lumen/osmotic feedback.

Do not fit a free delay constant to the knockout curve without independent justification.

Use WT time-course information, transporter-level regulation, or other independent data to constrain the stimulus dynamics first.

Then freeze the model and reveal the held-out AE4-null secretion time course.

Create:

- `analysis/13_state_resolved_ae4/dynamic_reconstruction.md`
- reproducible WT and knockout trajectory data under `results/13_state_resolved_ae4/`.

## Round 7. Held-out secretion test and model selection

At this stage the model must be frozen with respect to the AE4-null secretion outcome.

Evaluate:

- total secretion ratio;
- endpoint ratio;
- cumulative post-stimulus secretion ratio;
- first 2-3 minute behavior;
- sustained late deficit;
- WT absolute/relative secretion where units are certified;
- AE2-null secretion;
- ionic phenotypes.

A candidate is not successful merely because one endpoint lies near the target.

Compare successful candidates by:

1. number of new mechanisms;
2. number of free parameters;
3. number of independent data sets explained;
4. robustness to source-supported parameter uncertainty;
5. out-of-sample secretion prediction;
6. thermodynamic and conservation consistency.

Prefer the minimal model that explains the broadest evidence.

Create `analysis/13_state_resolved_ae4/heldout_validation.md`.

## Round 8. Structural explanation

If a model succeeds, do not stop at “the simulation matches.” Determine **why**.

Use exact balances, Jacobian/source decomposition, limiting cases, and knockouts of individual model features to identify the necessary mechanism.

Examples of acceptable conclusions include:

- cation-specific sequential AE4 cycling is necessary;
- dynamic PKA regulation creates the delayed WT/KO divergence;
- AE4 requires a specific pump/K-channel topology to translate chloride loading into secretion;
- acid-base support is jointly necessary with cation current;
- several distinct mechanisms are observationally equivalent.

Examples of unacceptable conclusions include:

- “parameter X is sensitive”;
- “the optimizer preferred model Y” without structural explanation;
- “AE4 is insufficient.”

Create `analysis/13_state_resolved_ae4/structural_explanation.md`.

## Round 9. Robustness and adversarial audit

The adversarial agent must attempt to break the result.

At minimum test:

- alternate steady roots;
- solver/tolerance dependence;
- reasonable WT parameter uncertainty;
- uncertainty in Na/K affinity and Hill coefficients;
- uncertainty in PKA fold changes;
- uncertainty in pump/K-channel localization fractions;
- alternative HCO3 versus carbonate assumptions where relevant;
- target leakage;
- whether success disappears if one unsupported parameter is fixed to a neutral value;
- whether a simpler nested model performs equivalently.

Use at least two independent numerical solvers or implementations for the final accepted model.

Create `analysis/13_state_resolved_ae4/adversarial_audit.md`.

## Round 10. Conditional modern inference scan

Only if a physiologically validated model or model family survives Rounds 7-9, perform a short modern inference scan.

Ask:

- Are multiple AE4 mechanisms observationally equivalent under secretion alone?
- Which additional ion or perturbation separates them?
- Is there a clean composite model-class discrimination problem?
- Can a minimal experimental design be stated analytically?

Do not force information theory into the paper. This round is only to determine whether the validated mechanism naturally creates an identifiability or discrimination result.

Create `analysis/13_state_resolved_ae4/inference_scan.md`.

## Required code architecture

Create new code outside `archive/` with modular separation between:

- state-resolved AE4 cycles;
- whole-cell chassis;
- experimental protocols;
- parameter/evidence definitions;
- solvers;
- validation metrics.

Suggested namespace:

- `src/state_resolved_ae4/`

Do not copy historical MATLAB source line for line. Historical code may be used as implementation evidence, but the new code must be independently written and provenance-labeled.

Add tests for:

- thermodynamic cycle consistency;
- electroneutrality;
- detailed-balance identities;
- Na/K-specific mutant limits;
- WT baseline closure;
- AE2 knockout behavior;
- AE4-null candidate independence where appropriate;
- conservation;
- solver reproducibility;
- target holdout discipline.

## Machine-readable ledgers

Maintain at least:

- `results/13_state_resolved_ae4/model_registry.csv`
- `results/13_state_resolved_ae4/parameter_provenance.csv`
- `results/13_state_resolved_ae4/calibration_targets.csv`
- `results/13_state_resolved_ae4/heldout_targets.csv`
- `results/13_state_resolved_ae4/model_scores.csv`
- `results/13_state_resolved_ae4/rejected_models.csv`

Every model row must record whether AE4-null secretion data were ever exposed before freezing. Any row with leakage is inadmissible as a predictive result.

## Required final answer

Create `analysis/13_state_resolved_ae4/final_answer.md` and answer all of the following directly.

1. What is the smallest AE4 state-resolved mechanism compatible with the transporter-level evidence?
2. Does cation-specific state structure explain anything that the Task 12 coarse-grained laws could not?
3. Can a source-supported state-resolved AE4 model reproduce the AE4-null secretion phenotype on the original best chassis?
4. If not, what exact balance/source direction is missing in the AE4-null state?
5. Which known module or smallest module combination supplies that direction?
6. Does adding that module, calibrated without AE4-null secretion, reproduce the secretion magnitude?
7. Does it reproduce the early-comparable then sustained-deficit time pattern?
8. Which features are necessary, which are sufficient, and which are merely one successful parameterization?
9. Are multiple successful mechanisms observationally equivalent?
10. What single experiment would best discriminate them?
11. Which conclusions are direct experimental facts, model-derived results, or new assumptions?
12. Does the result justify a new revisit paper, and what is the precise scientific contribution without saying the old article was wrong?

The final classification must be exactly one of:

- `PHENOTYPE RECONSTRUCTED — STATE-RESOLVED AE4 SUFFICIENT ON VALIDATED CHASSIS`
- `PHENOTYPE RECONSTRUCTED — AE4 PLUS MINIMAL CHASSIS COUPLING REQUIRED`
- `MULTIPLE SOURCE-CONSISTENT MECHANISMS RECONSTRUCT PHENOTYPE`
- `PHENOTYPE UNRESOLVED — SPECIFIC MISSING EXPERIMENT OR FLUX IDENTIFIED`

A generic `STOP` is not permitted.

For the final category, the missing experiment or flux must be specific enough that a researcher could act on it. “More data are needed” is not sufficient.

## Commit and branch discipline

Commit all work to `codex/state-resolved-ae4-cycle`.

Do not merge to `main`.

Do not modify `archive/`.

Do not draft manuscript prose.

## Acceptance criteria

Task 13 is complete only when:

- the primary evidence has been re-audited and frozen;
- at least one genuinely state-resolved AE4 family has been derived and calibrated to transporter-level evidence;
- the fixed-chassis Stage A test has been completed;
- if Stage A fails, the AE4-null missing balance direction has been mathematically localized;
- the smallest source-supported chassis extension has been tested without fitting AE4-null secretion;
- the time-resolved secretion phenotype has been evaluated;
- at least two independent numerical routes reproduce the final accepted result;
- an adversarial audit has been completed;
- the final answer identifies either a successful minimal mechanism or a specific irreducible missing experiment/flux.

Do not stop merely because one round fails. Continue through the hierarchy until the scientific bottleneck is genuinely resolved.