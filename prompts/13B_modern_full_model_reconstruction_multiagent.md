# Task 13B — Modern full salivary AE4 model reconstruction

## Mission

Task 13 did not solve the whole-cell reconstruction problem. It established that modern state-resolved AE4 transport models can encode real transporter-level structure, but no tested candidate produced a valid WT-calibrated whole-cell model on the inherited seven-state chassis. It also showed that AE4-null Cl and pH alone do not determine a unique missing source correction.

Task 13B has one objective.

> Build a physiologically credible modern full salivary acinar-cell model that satisfies independent WT and transporter constraints and then predicts the AE4-null secretion phenotype without fitting that phenotype.

This is not a reduction task. Do not attempt a three-ODE model, GSPT, bifurcation analysis, identifiability theory, or manuscript drafting until the full model works.

Read `AGENTS.md` first and obey it throughout.

Do not edit anything under `archive/`.

Do not merge to `main`.

Do not stop after one failed model generation.

---

# 1. Scientific target

The historical 2018 AE4 model was constructed before several later mechanistic results were available. The modern reconstruction must use later biology without treating any later molecular hypothesis as established fact.

The core problem is not to force one transporter formula to give a 35 percent knockout effect. The core problem is to reconstruct the full balance architecture through which experimentally established AE4 transport and regulation contribute to salivary secretion.

A successful model must explain, at minimum, why removal of AE4 changes intracellular chloride and pH and eventually reduces stimulated secretion while preserving independently measured WT physiology and the small functional role of AE2 in saliva output.

The held-out secretion result must emerge from the model after calibration is frozen.

---

# 2. Starting point from Task 13

Treat the following Task 13 findings as established project results that must be reproduced before extension.

1. State-resolved AE4 models can encode shared-carrier Na/K competition and source structure absent from a single scalar event-rate law.
2. The smallest source-supported shared-pool family was represented by SR2-type models with separate Na and K branches sharing carrier states.
3. Pure-cation dose data do not uniquely distinguish alternative microscopic mechanisms such as transported-cation cooperativity versus transported plus catalytic cation models.
4. On the inherited seven-state whole-cell chassis, zero tested state-resolved transporter candidates gave a valid WT physiological calibration.
5. After AE4 deletion, all candidates embedded in the same unchanged chassis share the same AE4-null vector field.
6. The bounded AE4-null root on that chassis substantially misses measured knockout Cl and pH.
7. Adding a simple apical/basolateral pump and K-channel split did not by itself produce a valid WT model.
8. Released knockout Cl and pH leave several unmeasured state coordinates free, so they do not define a unique one-protein correction.
9. Task 13 did not establish that AE4 itself is insufficient.
10. Task 13 did not establish that the 2018 paper was wrong.

Before changing the model, reproduce the decisive Task 13 numerical checks from saved artifacts. If any do not reproduce, resolve that discrepancy first.

---

# 3. Primary evidence base

Build a source ledger before coding new equations. Separate direct measurement from inference, model hypothesis, and historical assumption.

At minimum audit the following primary sources.

## 3.1 AE4 knockout physiology

Peña-Münzenmayer et al. 2015, JBC, DOI `10.1074/jbc.M114.612895`.

Extract exactly

- WT and AE4-null total stimulated saliva over the experimental protocol;
- the reported approximately `35 +/- 4.7%` reduction and sample sizes;
- the statement that early secretion is relatively comparable before a later sustained deficit;
- resting intracellular Cl in WT and AE4-null acinar cells;
- resting pH in WT and AE4-null cells;
- AE2-null and AE2/AE4-related phenotypes relevant to model validation;
- NKCC1 and NHE-related findings that constrain compensatory explanations;
- exact stimulus, tissue, timing, genotype, and preparation.

Do not convert descriptive statements into exact confidence intervals or equality constraints unless the paper provides them.

## 3.2 AE4 Na/K transport mechanism

Peña-Münzenmayer et al. 2016, JGP, DOI `10.1085/jgp.201611571`.

Extract

- direct Na transport evidence;
- direct K transport evidence;
- electroneutrality and voltage-independence evidence;
- Na and K EC50 values and Hill summaries with assay context;
- reversal and thermodynamic arguments;
- cation nonselectivity evidence;
- the candidate stoichiometries discussed by the authors;
- the authors' explicit statement that Hill coefficients may reflect more than one transported cation or allosteric regulation;
- HCO3/HCO3 exchange and any carbonate alternatives that remain unresolved.

Keep pure-cation heterologous assay conditions separate from native salivary physiology.

## 3.3 Dynamic beta/cAMP/PKA regulation of AE4

Peña-Münzenmayer et al. 2021, AJP Gastrointestinal and Liver Physiology, DOI `10.1152/ajpgi.00145.2021`.

This source is mandatory for Task 13B.

Extract

- beta-adrenergic stimulation protocol and measured AE4 response in native salivary acinar cells;
- H89 inhibition evidence;
- forskolin or adenylate-cyclase-related evidence if present;
- constitutively active PKAc evidence in heterologous cells;
- S173A and S273A results;
- measured response magnitudes, time traces, and observation windows where recoverable;
- any information about onset, persistence, reversibility, or washout;
- whether the data support direct phosphorylation of AE4 or only a PKA-dependent mechanism;
- whether a quantitative cAMP concentration or PKA activity trace was actually measured.

Do not assume the paper identifies a phosphorylation rate, dephosphorylation rate, cAMP time constant, or direct S173 phosphorylation if it does not.

## 3.4 2025 structural and mutagenesis evidence

Catalán et al. 2025, AJP Cell Physiology, DOI `10.1152/ajpcell.00346.2024`.

Extract

- T448, T756, D709, T713, G449, K879 evidence;
- Na versus K differences;
- the T756A-T448I double-mutant behavior;
- direct experimental findings versus MD hypotheses;
- proposed state-level transport cycles and alternative stoichiometries;
- HCO3 versus carbonate possibilities;
- which claims are measurements and which are speculation.

Do not hard-code the proposed 2025 sequential cycle as truth.

## 3.5 Apical cation topology

Almássy et al. 2018, Pflügers Archiv, DOI `10.1007/s00424-018-2109-0`.

Extract

- evidence for apical Ca-activated K current;
- evidence that Na/K-ATPase is distributed beyond the basolateral membrane and includes the apical pole;
- what was experimentally measured versus assumed in the mathematical model;
- any quantitative pump or conductance fractions that are directly supported versus merely model choices;
- consequences for luminal Na/K balance and membrane-current closure.

Do not import a localization fraction as measured unless the paper actually measured it.

## 3.6 Historical mathematical lineage

Vera-Sigüenza et al. 2018, Bull Math Biol, DOI `10.1007/s11538-017-0370-6`.

Use for

- equations and historical topology;
- original calibration targets;
- old water and lumen structure;
- old AE4 Markov formulation;
- reported model phenotype;
- stated discrepancies and limitations.

Treat it as a historical model, not an executable ground truth.

---

# 4. Multi-agent organization

Run a genuine multi-agent investigation. The lead agent must assign independent roles and reconcile disagreements explicitly.

At minimum create these roles.

## Agent A — experimental evidence and protocol reconstruction

Build `analysis/13B_modern_full_model/evidence_ledger.md` and machine-readable target tables.

Tasks

- extract all primary measurements and assay contexts;
- distinguish SMG versus parotid data;
- distinguish native acinar, whole-gland, CHO, and HEK experiments;
- record genotype, stimulus, bath, temperature, duration, and measured variable;
- identify which measurements may legally calibrate the model and which remain held out;
- identify missing numerical values rather than inventing them.

## Agent B — cAMP/PKA regulatory modeling

Build `analysis/13B_modern_full_model/camp_pka_model.md`.

Tasks

- reconstruct the smallest dynamic regulatory schemes consistent with 2021 evidence;
- consider nested models rather than one arbitrary pathway;
- determine which variables are actually required to represent beta-adrenergic activation of AE4;
- derive positivity and steady-state conditions;
- identify which kinetic parameters can be constrained from 2021 traces and which cannot;
- test whether an effective one-state regulatory variable is observationally equivalent to a two-state cAMP plus PKA model under available data;
- test common AE4 capacity modulation versus state-specific modulation of AE4 cycle transitions;
- include S173A/S273A logic as qualitative discrimination constraints where justified;
- never choose a time constant from the AE4-null secretion curve.

Candidate nested regulatory families should include at least

`R0` no dynamic regulation, static source-supported fold only;

`R1` one effective PKA/AE4 activation state;

`R2` effective cAMP plus downstream AE4 activation state;

`R3` phosphorylation/dephosphorylation state of AE4 with upstream PKA activation;

`R4` state-specific modulation of selected AE4 transition rates if supported by transporter data.

Do not retain R4 merely because it is flexible.

## Agent C — AE4 transport cycle and thermodynamics

Start from Task 13 state-resolved families.

Tasks

- reproduce SR2/SR5 transporter-level behavior;
- simplify any state graph that contains unidentifiable microscopic rates without losing experimentally supported behavior;
- preserve direct Na and K transport and electroneutrality;
- maintain local detailed balance;
- preserve mutation logic where source-supported;
- derive how the regulatory variable couples to the cycle;
- determine whether explicit carrier-state dynamics are required or whether QSS carrier occupancy is justified relative to whole-cell time scales;
- if carrier dynamics are retained, quantify the time-scale evidence independently of the knockout secretion curve.

The objective is not maximal molecular detail. The objective is the smallest source-consistent AE4 module usable inside a full physiological model.

## Agent D — whole-cell epithelial and acid-base reconstruction

Build the modern full model from first principles and source evidence rather than patching residuals blindly.

Audit every state, flux, membrane orientation, and conserved quantity.

Explicitly decide whether the model requires states for

- intracellular Na;
- intracellular K;
- intracellular Cl;
- intracellular HCO3;
- intracellular H or pH;
- dissolved CO2 or total inorganic carbon;
- cell volume;
- luminal Na;
- luminal K;
- luminal Cl;
- luminal bicarbonate or total inorganic carbon;
- apical and basolateral membrane potential;
- cAMP/PKA/AE4 regulatory state;
- explicit AE4 carrier occupancies if not QSS.

Derive or audit

- NKCC1;
- NHE1;
- AE2;
- AE4;
- Na/K-ATPase at both membranes where supported;
- K conductances at both membranes where supported;
- apical Cl conductance;
- paracellular Na and K pathways;
- water fluxes;
- luminal outflow;
- buffer chemistry;
- CO2 hydration/dehydration and open-bath coupling;
- charge balance and current closure.

A larger ODE system is allowed. Do not preserve seven states merely for convenience.

## Agent E — dimensional and time-scale audit

Build a complete dimensional ledger.

Tasks

- establish physical units for every state, parameter, rate, current, flux, and water term;
- resolve current-to-molar-flux conversions;
- resolve membrane area and compartment-volume scaling;
- verify Faraday factors and charge stoichiometry;
- resolve code time into physical time using independent source information where possible;
- identify parameters whose units remain uncertified;
- refuse to label code units as minutes without a valid map.

Task 13 failed to license the historical physical time scale. Task 13B must address this directly.

## Agent F — WT calibration and root geometry

Calibrate the modern full model using only allowed evidence.

Use deterministic multistart, continuation, or constrained root methods as appropriate.

Fit the smallest number of parameters required to satisfy WT gates.

Do not fit knockout secretion.

Create a hierarchy of calibration tiers

Tier 0 exact source-known parameters;

Tier 1 parameters estimated directly from primary assays;

Tier 2 parameters calibrated to WT resting physiology;

Tier 3 weakly constrained nuisance parameters with declared bounds;

Tier 4 unsupported parameters that should be removed or fixed unless absolutely necessary.

The accepted model should minimize Tier 3 and avoid Tier 4.

## Agent G — dynamic simulation and stimulus reconstruction

Establish a physically defensible stimulated WT protocol before any knockout comparison.

Tasks

- reconstruct the Ca stimulus used by the historical model and compare it with experimental stimulation data;
- determine whether carbachol and isoproterenol should enter through distinct pathways;
- implement dynamic beta/cAMP/PKA regulation driven by the isoproterenol arm;
- validate WT ionic and flow dynamics wherever data exist;
- calibrate physical time only from WT or regulatory data;
- freeze all dynamic parameters before knockout secretion is revealed.

## Agent H — adversarial scientific audit

Freeze an audit checklist before seeing final held-out results.

Attempt to falsify

- holdout leakage;
- hidden compensation;
- unit inconsistencies;
- arbitrary source directions;
- invalid membrane signs;
- false conservation;
- overfit regulatory dynamics;
- unjustified apical fractions;
- branch selection;
- remote roots;
- solver dependence;
- parameter nonidentifiability disguised as mechanism identification;
- historical claims that exceed evidence.

## Agent I — independent numerical reproduction

Reimplement the accepted equations or decisive failed result independently enough to detect implementation mistakes.

At minimum independently reproduce

- WT root;
- AE4 and AE2 genotype perturbations;
- charge/current closure;
- carbon and buffer conservation;
- regulatory state dynamics;
- held-out KO secretion ratio if a valid WT denominator exists;
- any final source-rank or impossibility claim.

## Lead agent

The lead agent owns the iteration loop, evidence firewall, and final classification. It must record disagreements among agents.

---

# 5. Model-generation protocol

Do not jump directly to a large final model. Build generations and retain only justified changes.

Create `analysis/13B_modern_full_model/generation_log.md`.

Each generation must document

- parent generation;
- specific defect being addressed;
- changed equations;
- added or removed states;
- added or removed parameters;
- source justification;
- calibration data used;
- validation data used;
- WT gate results;
- genotype gate results if legally released;
- numerical stability;
- decision to retain or revert.

Use the following minimum generation sequence unless evidence supports a better ordering.

## Generation G0 — reproduce Task 13 baseline

Reproduce the best Task 13 state-resolved model and inherited chassis exactly.

This is a negative control.

## Generation G1 — repair full acid-base and carbon chemistry

The inherited chassis has persistent pH problems. Build a conservation-correct H/HCO3/CO2 or total-carbon block.

Requirements

- explicit carbon conservation accounting;
- explicit buffer pool or justified effective buffer law;
- correct bath coupling;
- NHE1 coupling;
- AE2 and AE4 bicarbonate source terms;
- charge consistency;
- no infinite-buffer shortcut unless clearly treated as a limit and tested.

Calibrate with WT pH and source-supported acid-base data only.

## Generation G2 — modern cation topology

Introduce apical and basolateral Na/K-ATPase and K channels using Almássy 2018 topology.

Requirements

- membrane-specific electrochemical driving forces;
- luminal Na/K source terms;
- conserved total pump or conductance when only redistribution is supported;
- source-supported fraction bounds or explicit uncertainty if fractions are not measured;
- nested recovery of the baseline topology.

Do not invent genotype-specific pump recruitment unless independently supported.

## Generation G3 — dynamic cAMP/PKA regulation

Introduce the smallest supported dynamic regulatory model from Agent B.

Requirements

- isoproterenol or beta-adrenergic input drives the regulatory subsystem;
- WT regulatory kinetics calibrated only from 2021 data or other independent WT data;
- S173A/S273A qualitative behavior used as a discriminator where possible;
- no use of the knockout secretion timing or magnitude;
- test static versus dynamic regulation as nested models;
- test whether the regulatory state changes AE4 capacity only or specific cycle transitions.

## Generation G4 — integrated AE4 transport plus regulation

Embed the retained state-resolved AE4 model in the repaired whole-cell chassis.

Check whether QSS carrier occupancy remains valid once dynamic regulation is present.

If explicit carrier dynamics are required, add them and establish their time scale independently.

## Generation G5 — full WT stimulated model

Before any held-out knockout prediction, require a valid WT model under the combined carbachol plus isoproterenol protocol.

Freeze

- topology;
- all parameter values;
- all regulatory time constants;
- input waveforms;
- physical time map;
- water/flow scaling;
- root selection rules;
- solver tolerances used for production.

Only after G5 passes all WT gates may the AE4-null secretion outcome be opened.

## Generation G6 and later — targeted repair only

If the held-out phenotype fails after a valid WT model exists, do not globally retune.

Use independent non-secretion phenotype data to localize the discrepancy.

Any further model change must be justified by an independently measured missing process. Do not use the direction of the secretion error itself as the parameter target.

---

# 6. Physiological gates

Define machine-readable gates before final calibration. The lead agent may refine numerical tolerances only from source uncertainty and must record changes before knockout reveal.

At minimum include the following.

## Gate A — mathematical integrity

- positive concentrations and volumes;
- electroneutrality where imposed;
- charge-correct transporter cycles;
- current closure at both membranes;
- water balance;
- carbon and buffer conservation;
- zero net thermodynamic cycle affinity at equilibrium;
- correct reversal signs.

## Gate B — WT resting physiology

Use independent source-supported ranges for

- intracellular Cl;
- intracellular pH;
- intracellular Na if available;
- intracellular K if available;
- cell volume;
- membrane potential if available;
- luminal/resting states where relevant.

No boundary-touching solution counts as a valid calibration unless that boundary is itself physiological and source-supported.

## Gate C — transporter assays

The selected AE4 model must remain compatible with

- direct Na transport;
- direct K transport;
- electroneutrality;
- cation dose summaries;
- reversal;
- 2025 Na/K mutation asymmetry;
- 2021 PKA regulation hierarchy.

Do not demand precise simultaneous fitting of data sets that were measured under different assay contexts unless an explicit observation model maps between them.

## Gate D — AE2 validation

The model must not create a large secretion role for AE2 inconsistent with the biological evidence.

Where direct AE2-null Cl or pH data exist, evaluate them separately from secretion.

## Gate E — WT stimulated physiology

Require a physically calibrated time course and sensible stimulated secretion before forming any knockout ratio.

Normalized KO/WT behavior cannot rescue an invalid WT trajectory.

## Gate F — numerical robustness

- independent stiff solvers;
- tighter and looser tolerances;
- nearby initial states;
- deterministic root multistart;
- branch tracking where needed;
- conservation residuals reported explicitly.

---

# 7. Calibration firewall and staged data release

Create machine-readable files

- `results/13B_modern_full_model/calibration_targets.csv`
- `results/13B_modern_full_model/validation_targets.csv`
- `results/13B_modern_full_model/heldout_targets.csv`
- `results/13B_modern_full_model/reveal_log.json`

Hash or freeze them before optimization.

## Before held-out reveal

Allowed calibration evidence may include

- WT resting ions and pH;
- WT transporter-level assay data;
- WT cAMP/PKA regulatory data;
- source-supported membrane topology;
- WT stimulated data unrelated to the AE4-null comparison where independently available;
- AE2 data as an independent model validation gate.

The AE4-null saliva magnitude and time course are prohibited.

## Conditional localization evidence

If a valid WT model exists but the held-out phenotype fails, AE4-null resting Cl and pH may be used to diagnose the missing balance. If they were already exposed in Task 13, keep their role explicit and do not silently treat them as primary calibration data for the new model.

The key out-of-sample target remains secretion.

---

# 8. Dynamic cAMP/PKA modeling requirements

This is a mandatory part of Task 13B, but it must remain source disciplined.

At minimum analyze these nested possibilities.

## R0 static regulation

A source-supported static AE4 gain. Retain only as a negative/control model.

## R1 one effective regulatory state

For example

`dA/dt = activation(beta_input, A) - deactivation(A)`

with AE4 flux depending on `A`.

The exact functional form must be derived or justified. This model represents an effective beta/PKA/AE4 process and should be preferred if the data cannot resolve upstream cAMP separately.

## R2 cAMP plus AE4 activation

Represent an effective cAMP state and a downstream AE4 regulatory state.

For example, a saturating adenylate-cyclase/PDE balance for cAMP followed by a PKA-dependent activation/deactivation process.

Do not invent biochemical detail such as exact PDE isoforms unless relevant data require it.

## R3 phosphorylation-state model

Represent an active/inactive AE4 regulatory fraction with phosphorylation/dephosphorylation kinetics driven by PKA.

Use S173A as a qualitative constraint. Do not claim direct phosphorylation unless experimentally established.

## R4 state-specific transport-cycle regulation

Allow PKA regulation to modify selected microscopic AE4 transitions only if transporter/regulatory data materially require it.

Compare R4 against a common capacity modulation by nested model criteria. If R4 is not distinguishable, reject the extra detail.

For every retained regulatory model report

- parameter identifiability or profile uncertainty;
- WT activation time scale;
- whether the model predicts delayed AE4 engagement under the physiological stimulus;
- whether that delay was learned entirely from WT/regulatory data;
- sensitivity of final knockout predictions to allowed regulatory uncertainty.

---

# 9. Whole-cell equations must be rebuilt transparently

Do not make the new model a sequence of undocumented patches to `src/state_resolved_ae4/chassis.py`.

Create a new implementation under

`src/modern_full_model/`

with modules such as

- `states.py`
- `parameters.py`
- `transporters.py`
- `acid_base.py`
- `camp_pka.py`
- `membranes.py`
- `water.py`
- `model.py`
- `calibration.py`
- `validation.py`
- `run_analysis.py`

Exact filenames may differ, but the scientific modules must be separable and testable.

Every equation should have a provenance tag in code comments or a machine-readable registry using categories

- `PRIMARY_MEASUREMENT`
- `PUBLISHED_MODEL`
- `HISTORICAL_IMPLEMENTATION`
- `DERIVED_CONSTRAINT`
- `NEW_MODELING_DECISION`

Update `docs/PROVENANCE.md` and `docs/DECISIONS.md`.

---

# 10. Required analytical work before large numerical search

Do not rely only on optimization.

Derive as much as possible analytically.

At minimum derive

1. full species source matrix for all transporters and channels;
2. charge balance of every transport cycle;
3. total intracellular cation balance;
4. total intracellular anion/carbon balance;
5. membrane-current identities;
6. water/osmotic balance;
7. carbon conservation with bath exchange;
8. exact nesting limits that recover historical topologies;
9. equilibrium conditions for the cAMP/PKA subsystem;
10. limiting behavior as AE4 activity goes to zero;
11. limiting behavior under zero beta-adrenergic input;
12. structural rank of new pump/K redistribution directions;
13. which residual components each added module can and cannot alter.

Use these results to prevent impossible parameter searches.

---

# 11. Parameter fitting strategy

Do not fit everything simultaneously.

Use staged fitting.

## Fit 1 — transporter-level AE4

Fit only transporter-specific quantities to transporter assays.

No whole-cell phenotype data.

## Fit 2 — regulatory subsystem

Fit cAMP/PKA/AE4 regulatory parameters to 2021 WT regulatory experiments only.

If only an effective time scale is identifiable, collapse the model accordingly.

## Fit 3 — WT resting whole-cell state

Fit the minimal remaining whole-cell nuisance parameters required for WT closure.

Do not alter transporter assay parameters unless the assay-to-cell observation map explicitly requires context conversion.

## Fit 4 — WT stimulation/time map

Use WT dynamic data and source-supported timing information to calibrate physical time, calcium input, and any regulatory onset parameters.

Do not use knockout timing.

## Freeze

Freeze all parameters and generate immutable manifests before knockout reveal.

---

# 12. Held-out AE4-null test

Only after a valid WT denominator exists.

Set AE4 expression/activity to zero in a way consistent with genotype deletion. Remove all AE4 transport and regulatory contributions without changing unrelated parameters.

Evaluate separately

- 10-minute cumulative saliva KO/WT;
- late/endpoint flow ratio if physically meaningful;
- first 2 to 3 minute comparability;
- sustained later divergence;
- knockout intracellular Cl;
- knockout pH;
- knockout volume if an independent target exists;
- Na/K trajectories if predictions are available;
- conservation and solver robustness.

The primary secretion target is approximately `35 +/- 4.7%` lower total saliva over the experimental 10-minute protocol. Use the exact primary source definition rather than a rounded 30 percent narrative value.

Do not declare success from one scalar ratio if timing or ionic phenotypes are grossly wrong.

---

# 13. If the phenotype still fails

Do not globally retune.

Classify the residual into one of the following categories.

- transporter mechanism error;
- cAMP/PKA regulatory timing or gain uncertainty;
- cation-homeostasis topology;
- acid-base/carbon chemistry;
- calcium stimulus representation;
- water/flow conversion;
- missing measured transport pathway;
- missing genotype adaptation;
- uncertified physical unit mapping.

For the leading category, prove that the missing quantity is decision-critical.

This requires at least two source-consistent values or models that both satisfy all non-secretion evidence but produce materially different knockout secretion predictions.

If such a degree of freedom cannot be shown to change the decision, it is not the decisive missing measurement.

If outcome 4 is reached, specify one concrete experiment including

- genotype;
- tissue/preparation;
- stimulus;
- measured variable;
- time resolution;
- measurement compartment;
- expected discriminating outcomes;
- how the result would choose between surviving model classes.

Do not end with “measure pump activity” or “collect more data.”

---

# 14. Required outputs

Create

`analysis/13B_modern_full_model/`

containing at least

- `evidence_ledger.md`
- `task13_reproduction.md`
- `model_architecture.md`
- `camp_pka_model.md`
- `ae4_module.md`
- `acid_base_reconstruction.md`
- `cation_topology.md`
- `dimensional_ledger.md`
- `generation_log.md`
- `wt_calibration.md`
- `dynamic_validation.md`
- `heldout_validation.md`
- `adversarial_audit.md`
- `independent_reproduction.md`
- `final_answer.md`

Create machine-readable results under

`results/13B_modern_full_model/`

including at least

- calibration and holdout ledgers;
- parameter provenance;
- model-generation registry;
- accepted/rejected model table;
- root table;
- WT trajectories;
- genotype trajectories if licensed;
- solver crosschecks;
- conservation residuals;
- regulatory fit summaries;
- final frozen parameter manifest;
- artifact hashes where practical.

Add tests under `tests/` for every accepted scientific invariant and decisive numerical result.

---

# 15. Adversarial audit checklist

The adversarial agent must freeze a checklist before final held-out evaluation.

It must challenge at least

1. whether the knockout secretion target leaked into any objective, parameter bound, branch choice, time constant, delay, or model selection decision;
2. whether cAMP/PKA time constants were effectively chosen from the knockout curve;
3. whether S173A was overinterpreted as proof of direct phosphorylation;
4. whether static and dynamic regulation are actually distinguishable from the available WT data;
5. whether added acid-base states conserve total carbon and charge;
6. whether apical pump/K fractions are measured or merely assumed;
7. whether the accepted WT root is unique over the declared physiological domain or only one convenient branch;
8. whether all current-to-flux and area/volume conversions are dimensionally correct;
9. whether the physical time map is independent of the knockout phenotype;
10. whether the model reproduces WT physiology before normalized genotype ratios are computed;
11. whether parameter compensation hides a wrong source topology;
12. whether knockout deletion removes every AE4 contribution consistently;
13. whether AE2 behavior remains acceptable;
14. whether two independent solvers reproduce the decisive trajectory;
15. whether any accepted mechanism depends on parameters sitting on arbitrary search bounds;
16. whether a simpler nested model performs equally well;
17. whether model predictions remain qualitatively stable under source-supported uncertainty;
18. whether claims about the historical 2018 paper exceed what the new reconstruction demonstrates.

---

# 16. Completion criteria

A generic `STOP`, “more data are needed,” or “model is too complex” is not an acceptable final answer.

Finish only with one of the following exact classifications.

## Outcome 1

`FULL MODEL VALIDATED — HELD-OUT AE4 PHENOTYPE RECONSTRUCTED`

Requirements

- valid WT resting model;
- valid WT stimulated model;
- independent physical time map;
- modern AE4 transport and regulation consistent with source evidence;
- held-out AE4-null secretion magnitude and qualitative timing reproduced without retuning;
- knockout ionic phenotypes not grossly contradicted;
- independent numerical reproduction;
- adversarial audit passes.

## Outcome 2

`FULL MODEL VALIDATED — PHENOTYPE PARTIALLY RECONSTRUCTED; ONE SPECIFIC DYNAMIC OR FLUX DEFICIT REMAINS`

Requirements

- valid WT full model exists;
- knockout can be legitimately predicted;
- substantial part of phenotype is explained;
- exactly one remaining source-supported deficit is localized;
- its effect is demonstrated rather than guessed.

## Outcome 3

`MULTIPLE FULL MODELS VALIDATED — PHENOTYPE MECHANISM NONUNIQUE`

Requirements

- at least two source-consistent full models independently pass WT gates and held-out phenotype;
- they differ mechanistically in a scientifically meaningful way;
- the available evidence cannot discriminate them;
- one decisive experiment is specified.

## Outcome 4

`FULL MODEL NOT IDENTIFIABLE FROM EXISTING DATA — ONE DECISION-CRITICAL MEASUREMENT SPECIFIED`

Requirements

- exhaustive source-supported reconstruction attempts completed;
- failure is not due to solver or implementation problems;
- at least two non-secretion-compatible possibilities imply materially different held-out predictions;
- one missing measurement is shown mathematically/numerically to decide between them;
- the experiment is concrete and executable.

---

# 17. Final answer requirements

`analysis/13B_modern_full_model/final_answer.md` must answer directly

1. What is the final full-model state vector?
2. Which equations differ from the 2018 model and why?
3. Which parts of Task 13 were retained, simplified, or rejected?
4. What is the smallest defensible cAMP/PKA regulatory subsystem?
5. What does the 2021 evidence identify and what remains unmeasured?
6. Does cAMP/PKA regulation introduce a source-supported slow time scale?
7. Is AE4 regulation best represented as common capacity modulation or state-specific transport-cycle regulation?
8. Which acid-base/carbon formulation is required for WT pH closure?
9. Which apical/basolateral pump and K-channel topology is required or allowed?
10. Does the final model reproduce WT resting Cl, pH, volume, Na/K, and membrane constraints?
11. Does it reproduce a valid WT stimulated secretion trajectory on a physical time axis?
12. What does it predict for AE2 deletion?
13. What does it predict for AE4 deletion before any retuning?
14. Does it reproduce the 10-minute AE4-null saliva deficit?
15. Does it reproduce the early-comparable then sustained-deficit pattern?
16. Which conclusions are exact, numerical, source-measured, or assumption-dependent?
17. Which parameters remain weakly identified?
18. What is the minimal mechanism required by the successful reconstruction, if one exists?
19. If unsuccessful, what exact measurement is decision-critical and why?
20. Is the full model now sufficiently validated to justify a later reduction/GSPT task?

The answer to question 20 must be `NO` unless the full model has passed the WT and numerical gates. A future reduction task is not a substitute for fixing the full model.

---

# 18. Scientific language constraints

Do not write

- “the 2018 paper was wrong”;
- “AE4 itself is insufficient”;
- “the old model is invalid” without precise scope;
- “PKA phosphorylates S173” as an established fact unless the source directly establishes phosphorylation rather than dependence;
- “the 2025 cycle is the AE4 mechanism” unless the data actually identify it;
- “the phenotype is explained” if the model was fitted to the phenotype.

Preferred formulations are scoped statements such as

- “the inherited historical-lineage chassis does not satisfy the modern WT physiological gates under the tested reconstruction”;
- “the 2021 experiments require PKA-dependent AE4 regulation but do not uniquely identify its kinetic implementation”;
- “the 2025 data constrain cation coordination but leave the full transport sequence nonunique”;
- “the held-out knockout phenotype is predicted after all WT calibration and dynamic inputs are frozen.”

---

# 19. Execution rule

Continue through the full iteration protocol. Do not terminate when one module fails. Do not skip the cAMP/PKA reconstruction. Do not switch to reduction or bifurcation work.

The purpose of Task 13B is to get the full model working, or to prove exactly why the currently available evidence cannot determine one.