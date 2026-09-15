# Task 48: joint experimental-constraint reconstruction of AE4-dependent secretion

## Mission

This is a **data-constrained reconstruction**, not another mechanism search.

Task 47 has already rejected the hypothesis that the inherited K-recycling / pump / voltage block supplies the missing AE4-null ceiling. Do not repeat Task 47. Start from its final published commit:

`0a6adda00f0bb57d527f341412ceffa348b67a76`

Work only on:

`analysis/task-48-joint-experimental-constraint-reconstruction`

The scientific objective is now:

> Determine whether one physiologically coherent dynamic model can reproduce the full experimentally observed AE4/Ae2 phenotype across resting, stimulated, isolated-transporter and secretion protocols using one shared parameterisation, and if the inherited architecture cannot do so, make the smallest **directly data-required** correction.

There are no remaining scientific resources for broad exploration. Treat this run like a constrained experiment. Do not spend it testing speculative mechanism families.

## Read first

1. `AGENTS.md`
2. `analysis/47_dynamic_potassium_recycling_reconstruction/DYNAMIC_POTASSIUM_RECYCLING_REPORT.md`
3. Task 43 parameter-provenance outputs
4. Task 44 mathematical-analysis outputs
5. Task 46 CBM decision/handoff outputs
6. `model/equations.md`
7. `model/fluxes.md`
8. `model/parameters.md`
9. `analysis/48_joint_experimental_constraint_reconstruction/constraints_seed.json`

Use the seed JSON as a starting ledger, not as permission to ignore repository evidence. Verify every entry against the primary source before using it.

Primary source:
Peña-Münzenmayer et al., J Biol Chem 2015, DOI `10.1074/jbc.M114.612895`, PMCID `PMC4409235`.

Secondary mechanistic sources already identified in the seed:
- Peña-Münzenmayer et al., J Gen Physiol 2016, DOI `10.1085/jgp.201511529`
- Catalán/Peña-Münzenmayer et al., AJP GI 2021/2022, DOI `10.1152/ajpgi.00145.2021`
- Evans et al., J Biol Chem 2000, DOI `10.1074/jbc.M003753200`
- Park et al., J Biol Chem 2001, DOI `10.1074/jbc.M102901200`

## What changed relative to Tasks 41–47

The approximately 35% AE4-null secretion phenotype is **no longer a held-out scalar target**. It is one member of a much larger joint data vector.

You are explicitly authorised to infer parameters that Task 43 classifies as fitted, model-determined, placeholder or otherwise genuinely uncertain, but only by fitting **all mapped experimental constraints simultaneously**.

You are not authorised to:
- tune one parameter to the 35% endpoint;
- create genotype-specific multipliers;
- create protocol-specific fudge factors;
- introduce a phenomenological secretion penalty;
- invent a numerical tolerance from “not significant”;
- alter measured/literature-fixed parameters outside their documented uncertainty;
- perform a Cartesian-product sweep;
- enumerate mechanism families;
- rerun Tasks 46 or 47 to reproduce known results.

One parameter vector must serve all genotypes and protocols. A genotype changes only the experimentally deleted transporter. A protocol changes only its documented agonists, extracellular solution and inhibitors, plus a source-supported signalling input.

## Central experimental fingerprint

The model must face all of the following at once.

### AE4 knockout, physiological CCh + IPR

- Total 10-minute saliva: AE4 KO is `35 ± 4.7%` lower than WT (`n=6` each).
- Time course: KO secretion is comparable to WT for the initial `2–3 min`, then has a sustained reduction after approximately 3 min.
- Resting `[Cl-]i`: WT `50.10 ± 1.50 mM` (`n=7`), KO `36.50 ± 1.60 mM` (`n=6`).
- Resting pHi: WT `6.91 ± 0.07` (`n=4`), KO `6.89 ± 0.02` (`n=4`).
- Initial Cl exit, CCh + IPR, SPQ ΔF: WT `-0.21 ± 0.01` (`n=10`), KO `-0.18 ± 0.02` (`n=9`).
- Initial Cl uptake, CCh + IPR: WT `2.02 ± 0.10`, KO `0.90 ± 0.09` in `10^-3 s^-1`.
- Initial Cl uptake, CCh only: WT `2.18 ± 0.20`, KO `2.30 ± 0.10`.
- Initial Cl uptake, IPR only: WT `0.40 ± 0.07`, KO `0.20 ± 0.03`.
- Salivary HCO3 amount over 10 min: WT `5.3 ± 1.4`, KO `7.1 ± 0.8` microequiv, `p=0.28`.

This means the desired behaviour is not “35% less output”. It is:
1. a substantially lower KO resting chloride state;
2. essentially preserved resting pH;
3. essentially preserved initial chloride exit;
4. preserved muscarinic-only uptake;
5. impaired beta-adrenergic uptake;
6. strongly impaired combined uptake;
7. preserved NHE activity;
8. no detectable increase in isolated NKCC activity;
9. preserved bicarbonate output;
10. a secretion deficit that emerges mainly during sustained stimulation.

### AE2 knockout controls

Use all Table 1 AE2 control/KO values in `constraints_seed.json`.

In addition:
- AE2 KO saliva kinetics and total output are essentially unchanged from its own strain-matched control (`n=6` each).
- HCO3 output: control `4.0 ± 0.4` (`n=5`), AE2 KO `6.3 ± 1.8` (`n=6`), `p=0.28`.

This is a strong negative control. The reconstruction must not “solve” AE4 loss by turning AE2 into a transporter whose deletion would itself cause a large secretory phenotype.

### Isolated NKCC1 assay

Protocol:
- HCO3-free;
- `30 uM` ethoxyzolamide;
- `50 uM` T16Ainh-A01;
- Cl depletion followed by readdition.

Observed:
- no statistically detectable WT/AE4-KO difference in NKCC1-dependent initial Cl uptake;
- no statistically detectable control/AE2-KO difference;
- `50 uM` bumetanide reduces uptake by `95.6 ± 0.9%`.

This is **not** permission to impose exact equality in physiological stimulation. It is a separate assay that the model must simulate as a separate protocol.

### NHE assay

Under `0.3 uM CCh + 5 uM IPR`, stimulation-induced NHE-mediated alkalinisation is not detectably different between control, AE4 KO and AE2 KO. EIPA reverses the alkalinisation to acidification.

Task 47's compensation involved reduced NHE and NBC burden. A reconstruction that relies on strongly suppressing NHE in AE4 KO is therefore contradicted by the experiment.

### Isolated HCO3-dependent exchanger assays

With NKCC1 and Tmem16A inhibited:
- AE4 KO reduces HCO3-dependent Cl uptake by `38.2 ± 6.6%`.
- AE2 KO reduces it by `55.1 ± 3.4%`.
- AE4/AE2 double KO reduces it by `87.7 ± 3.1%`.

Beta-adrenergic activation:
- control: `1.56 ± 0.19` without IPR to `4.33 ± 0.91` with IPR;
- AE4 KO: `0.95 ± 0.10` without IPR to `1.08 ± 0.09` with IPR.

Na dependence in AE2-KO cells, where the remaining exchanger is AE4-dominated:
- Na present: `1.35 ± 0.16`;
- Na-free: `0.51 ± 0.09`;
- `p=0.001`.

### Independent AE4 biophysics and regulation

The 2016 JGP study shows:
- AE4 transport is electroneutral;
- transport is associated with Cl, HCO3, Na and K;
- external Na EC50 approximately `49 mM`, Hill coefficient `2.0`;
- external K EC50 approximately `62 mM`, Hill coefficient `1.8`;
- both Na and K can support/participate in AE4 transport.

Do **not** treat a fixed realised Na/K routing fraction as experimentally measured merely because the current model uses one.

The 2021/2022 PKA study shows:
- beta-adrenergic/cAMP stimulation increases AE4 activity;
- H89 prevents that activation;
- constitutively active PKA increases activity;
- the S173A mutation abolishes PKA activation.

Therefore, if the inherited model lacks an explicit cAMP/PKA-dependent AE4 activation mechanism, that is a **documented missing mechanism**, not a speculative class.

## Mandatory protocol reconstruction

Before optimisation, implement a protocol layer that can simulate at least:

1. WT rest.
2. AE4-KO rest.
3. AE2-KO rest.
4. WT CCh only.
5. AE4-KO CCh only.
6. WT IPR only.
7. AE4-KO IPR only.
8. WT CCh + IPR.
9. AE4-KO CCh + IPR.
10. AE2 control and AE2 KO CCh + IPR.
11. NKCC isolation in WT/control, AE4 KO and AE2 KO.
12. NKCC isolation + bumetanide.
13. NHE assay ± EIPA where representable.
14. HCO3-dependent exchanger isolation for WT/control, AE4 KO, AE2 KO and double KO.
15. Isolated exchanger ± IPR.
16. Isolated exchanger Na-present versus Na-free in AE2 KO.

Do not emulate an inhibitor by changing an unrelated parameter. Disable the experimentally inhibited pathway itself.

## Genotype-specific resting states are mandatory

The experimental AE4 KO is an established knockout animal, not an acute deletion from a shared WT resting state.

Task 47 explicitly noted this limitation.

Therefore:
- solve/relax WT to its own resting state;
- solve/relax AE4 KO with AE4 deleted to its own resting state;
- solve/relax AE2 KO with AE2 deleted to its own resting state;
- use those genotype-specific resting states as the starting conditions for their stimulation protocols.

Keep an acute-deletion-from-WT-state simulation only as a labelled counterfactual diagnostic if useful. It is not the experimental comparator.

A model that cannot generate the observed lower AE4-KO resting `[Cl-]i` while preserving resting pH has already failed an important part of the phenotype.

## Observation operators: do not compare unlike quantities

The experiment contains whole-gland outputs, concentrations and fluorescence-derived quantities. Build explicit observation operators.

Rules:

- `[Cl-]i` and pHi may be compared directly when the model has the same quantity.
- Whole-gland saliva and bicarbonate output must be compared by genotype ratios and time-course shape unless a separately justified acinus-to-gland scaling exists.
- SPQ `ΔF` is not automatically a molar Cl flux.
- SPQ-derived uptake rate in `10^-3 s^-1` is not automatically the model's transporter flux in fmol/s.
- If the fluorescence calibration/normalisation can be reconstructed from the paper, implement it.
- Otherwise compare dimensionless genotype/stimulus ratios or a normalized observation operator and state exactly what information is retained.
- “Not significant” is not exact equality.
- Never invent ±5%, ±10%, or any other tolerance.

If you digitize a figure, save:
- source figure;
- axis calibration;
- digitized points;
- digitization script;
- a documented digitization-error estimate.

## Parameter policy

Create `PARAMETER_ADMISSIBILITY.md` before fitting.

For every parameter touched by the protocols or candidate reconstruction, classify it using Task 43 as:

1. measured;
2. literature-derived;
3. fitted/calibrated;
4. model-determined;
5. placeholder;
6. unconstrained effective parameter.

Rules:
- measured: fixed or bounded by the source uncertainty;
- literature-derived: fixed or bounded by the source-supported range;
- fitted/model-determined/placeholder/unconstrained: may be inferred from the joint experimental data, but use one shared value across genotypes;
- do not allow a parameter to vary independently by genotype unless gene deletion mathematically removes that protein;
- do not allow independent CCh/IPR-specific transporter multipliers except where a documented signalling mechanism requires one;
- any new free parameter must correspond to a specific source-supported biochemical mechanism.

The existing AE4 constitutive parameters are model-determined according to the repository provenance. They are legitimate candidates for re-identification against the isolated AE4 data. That does **not** justify changing unrelated channel/pump parameters to make secretion land at 35%.

## Data-driven correction order

Do not open a mechanism search.

Evaluate corrections in this fixed order:

### R0: protocol-correct inherited architecture

First test the inherited equations with:
- genotype-specific rest;
- exact assay/inhibitor protocols;
- one shared parameter vector.

No equation changes.

### R1: AE4 cAMP/PKA regulation, if absent

Audit the AE4 law.

If it does not already contain experimentally supported beta-adrenergic/cAMP/PKA activation, add the **smallest explicit activation factor/state** capable of reproducing the isolated exchanger ±IPR data.

Requirements:
- the activation must affect AE4, not generic secretion;
- it must be absent when AE4 is deleted;
- it must be driven by the beta-adrenergic/cAMP input, not by genotype;
- fit its magnitude/kinetics to the exchanger assay data, not to whole-gland secretion;
- preserve electroneutrality and the source-supported AE4 ion coupling;
- do not add a secretion target term.

This is the primary reconstruction candidate because it is independently observed.

### R2: only discrepancy-forced equation repair

Only if R1 still cannot fit the joint data, inspect the standardized residual pattern and local sensitivity/Jacobian.

You may make **one** additional minimal equation-level correction only if:
- a specific observed constraint cannot be represented by the equation as written;
- the correction has direct independent physiological support;
- you can state which residual(s) require it before coding.

Do not test alternatives in parallel. Do not try a list of mechanisms. If no unique data-required correction exists, stop and report architectural underdetermination/failure.

## Optimisation / inference

This is a joint inverse problem.

Use all quantitatively mappable observations simultaneously.

Recommended objective:
- Gaussian observation term using reported means and SEs for model-predicted population means;
- source-based parameter prior/bound term where provenance supplies uncertainty/ranges;
- separate qualitative gates for observations without numeric data.

Do not silently treat SE as biological SD. Record exactly what likelihood interpretation you use.

Before high-dimensional inference:
1. compute local sensitivities of every observable to admissible parameters at the protocol-correct baseline;
2. remove parameters with essentially zero sensitivity across the full data vector;
3. identify severe collinearities;
4. fit only the identifiable/admissible subspace.

Prefer constrained least squares / trust-region methods and a small predeclared multi-start set if needed. Do not use a giant global sweep.

The objective is **joint fit quality**, not a weighted hack that lets saliva dominate everything else.

Report for every quantitative observation:
- experimental mean;
- SE and n where known;
- model prediction;
- standardized residual;
- protocol;
- observation operator used.

Report the global residual vector before and after any correction.

## Hard rejection conditions

A candidate reconstruction fails if it achieves the secretion deficit by any of the following:

- raising isolated AE4-KO NKCC activity substantially above WT despite the assay;
- suppressing AE4-KO NHE activity despite the NHE assay;
- producing a large resting pH shift;
- abolishing the experimentally preserved CCh-only Cl uptake;
- making AE2 deletion produce a large secretion phenotype;
- collapsing bicarbonate output;
- altering initial Cl exit enough to contradict the SPQ data;
- using a shared WT resting state for the established KO;
- requiring genotype-specific fitted capacity multipliers;
- matching the 35% endpoint while missing the sustained timing.

## Secondary cross-validation

After fitting the primary SMG AE4/Ae2 dataset, test without retuning:

- NKCC1 KO: >60% reduction in muscarinic-stimulated saliva in mouse parotid and loss of bumetanide-sensitive Cl influx.
- NHE1 KO: >95% loss of acinar Na/H exchange and impaired prolonged secretion, with compensatory NKCC expression.

These are different gland/stimulus contexts. Use them as direction/architecture validation, not as identical-protocol numerical constraints unless a justified translation is provided.

## Required outputs

Create under:

`analysis/48_joint_experimental_constraint_reconstruction/`

At minimum:

- `EXPERIMENTAL_CONSTRAINT_LEDGER.md`
- `constraints.json`  (verified replacement for the seed)
- `PROTOCOL_MAP.md`
- `PARAMETER_ADMISSIBILITY.md`
- `OBSERVATION_OPERATORS.md`
- `BASELINE_RESIDUALS.md`
- `output/baseline_residuals.csv`
- `output/final_residuals.csv`
- `output/parameter_vector.json`
- `output/protocol_predictions/`
- `RECONSTRUCTION_DECISION.md`
- `FINAL_REPORT.md`

`RECONSTRUCTION_DECISION.md` must answer explicitly:

1. Can the inherited equations, with genotype-specific rest and correct assay protocols, fit the joint data?
2. Does the model contain the experimentally observed cAMP/PKA activation of AE4?
3. If added, is that single documented correction sufficient?
4. Which constraints remain incompatible after the best admissible fit?
5. Is NKCC compensation still excessive under the isolated-NKCC constraint?
6. Is NHE preserved as experimentally observed?
7. Does AE2 KO remain essentially secretion-neutral?
8. Does AE4 KO reproduce the lower resting chloride without a resting-pH artefact?
9. Does the model reproduce the early-near-WT / late-deficit secretion timing?
10. Is there a defensible single architecture left, or does the current model class fail?

## Publication discipline

This task is expensive. Publish recoverable checkpoints.

Commit and push after:

1. verified constraint ledger + protocol map;
2. parameter admissibility + observation operators;
3. protocol-correct inherited baseline residual vector;
4. R1 cAMP/PKA AE4 reconstruction and fitted joint residual vector;
5. final decision/report.

If the inherited architecture already succeeds, do not add R1 unnecessarily.
If R1 fails and there is no uniquely justified R2, stop. A clean negative result is better than spending the remaining runs teaching the model new ways to cheat.

Do not merge to `main`.
Do not modify frozen Task 41–47 outputs.
Do not begin Task 49.
