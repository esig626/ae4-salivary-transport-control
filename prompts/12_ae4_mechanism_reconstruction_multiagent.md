# Task 12. Multi agent AE4 mechanism reconstruction

Work in repository `esig626/ae4-salivary-transport-control` on branch `codex/ae4-mechanism-reconstruction`.

Read `AGENTS.md` first. Then read the completed evidence from

- `analysis/00_inventory/`
- `analysis/10_identifiability_discrimination/`
- `analysis/11_forensic_reconstruction/`
- `model/`
- `docs/PROVENANCE.md`
- `docs/DECISIONS.md`
- `docs/RESULTS_LEDGER.md`

Do not edit anything under `archive/`.

## Central question

The experimental AE4 knockout phenotype is an approximately 30 percent reduction in stimulated salivary secretion. The archived executable model does not reproduce it, and the exact 2018 implementation is missing.

Do not assume that the archived AE4 law is correct. Do not assume that the published 2018 AE4 law is correct.

Determine:

> **What structural property of the AE4 transport mechanism is required for an otherwise minimally changed salivary secretion model to reproduce the experimentally observed AE4 knockout phenotype?**

The answer may be that a particular AE4 mechanism succeeds, that several mechanisms remain equivalent, or that changing AE4 alone is insufficient and another missing component can be localized. You must reach one of those answers.

A generic `STOP` is not an acceptable outcome.

---

# Multi agent organization

Run this as a genuine multi agent investigation. Use at least the following independent roles. Agents should work in parallel where possible and challenge one another's assumptions.

## Agent A. Experimental evidence and AE4 biology

Build the primary evidence ledger from the original experimental literature.

Tasks:

1. Find and read the primary AE4 transport study or studies underlying the 2018 model.
2. Find and read the primary salivary knockout study underlying the approximately 30 percent secretion reduction.
3. Extract the experimental facts needed to constrain an AE4 mechanism, including where available
   - exact knockout flow reduction and uncertainty
   - stimulation protocol
   - intracellular chloride phenotype
   - pH phenotype
   - Na versus K dependence
   - cation affinities or EC50 values
   - Hill coefficients
   - reversal or directionality evidence
   - stoichiometries explicitly established, excluded, or merely proposed
   - whether cations are demonstrated to be transported or may instead act allosterically
   - AE2 knockout phenotype
   - evidence for compensatory transporter changes
4. Separate direct measurements from mechanistic interpretation.
5. Record exact source locations and citations.

Output:
`analysis/12_ae4_mechanism_reconstruction/experimental_evidence.md`

Also create a machine readable phenotype ledger in
`results/12_ae4_mechanism_reconstruction/phenotype_targets.csv`.

Important: do not rely on the 2018 modelling paper to define the experimental AE4 mechanism when the primary experimental paper is available.

## Agent B. Chassis and historical model audit

Define the smallest trustworthy non AE4 chassis on which AE4 variants can be compared.

Tasks:

1. Use Task 11 to identify the archived model variant with the best internally closed WT baseline.
2. Reproduce its baseline independently outside `archive/`.
3. Document every non AE4 equation and parameter retained.
4. Identify all places where the archived AE4 mechanism affects balances directly or indirectly.
5. Identify which differences between archived code and 2018 publication are AE4 specific and which are chassis differences.
6. Establish conservation, electroneutrality, water balance, current balance, and dimensional checks for the chassis.
7. Create a clean test harness that allows AE4 law and stoichiometric balance insertion to be swapped without altering the rest of the model.

Primary comparison rule:
**Non AE4 equations and parameters remain fixed across AE4 candidates.**

Output:
`analysis/12_ae4_mechanism_reconstruction/chassis.md`

Implementation should live in a new namespace such as
`src/ae4_mechanism_reconstruction/`.

Do not port historical code verbatim. Reimplement the documented equations independently.

## Agent C. AE4 mechanism derivation

Construct the complete source justified family of AE4 mechanisms.

Start from the experimental evidence, not from one preferred mathematical form.

At minimum examine the following families if compatible with primary evidence.

### C1. Historical Na only model

Use as a negative control and lineage reference.

### C2. Published 2018 Na plus K mechanism

Reconstruct the exact published mechanism, including the way Na and K contributions enter the cycle and the cell balances.

### C3. Explicit parallel Na and K branches

Represent Na coupled and K coupled cycles separately rather than replacing them by a pooled `[Na]+[K]` term. Allow source justified relative affinities or weights.

### C4. Cooperative cation dependence

Use the experimentally reported cation dependence, including measured Hill coefficients and affinity parameters where available. Test whether cooperativity belongs in transport stoichiometry, occupancy, or a gating factor.

### C5. Cation dependent allosteric gating

Construct a mechanism in which cations regulate AE4 turnover without necessarily being transported in the same stoichiometric cycle, if this remains compatible with the experimental paper.

### C6. Alternative electroneutral stoichiometries

Extract every stoichiometry explicitly considered or left viable by the primary experimental study. Do not rely on remembered ratios. Derive the balance signatures and thermodynamic driving force for each candidate.

### C7. Thermodynamically consistent reversible formulations

For each viable stoichiometry, construct at least one reversible law whose sign is controlled by the corresponding electrochemical affinity and which obeys the correct equilibrium condition. Compare it with the polynomial or empirical historical laws.

### C8. Saturating versus mass action transport

Where data permit, compare mass action, occupancy based, and saturating formulations while preserving the same stoichiometry.

For every candidate state explicitly

- transported species and stoichiometry
- cation selectivity
- equilibrium/reversal condition
- number of free AE4 specific parameters
- how parameters are constrained from independent data
- which parts are experimentally measured versus newly assumed

Use dimensional analysis and thermodynamics to reject invalid variants before simulation.

Output:
`analysis/12_ae4_mechanism_reconstruction/ae4_candidate_family.md`

## Agent D. Calibration and predictive computation

Implement and test the candidate AE4 family on the fixed chassis.

### Critical anti leakage rule

The approximately 30 percent knockout phenotype is **held out validation**.

Do not fit any parameter to the knockout flow reduction in the primary comparison.

For each mechanism candidate:

1. Fit or normalize only AE4 specific parameters.
2. Use only
   - WT baseline closure
   - independent AE4 transport assay data
   - source supported resting concentrations and affinities
   as calibration information.
3. Freeze the resulting model.
4. Simulate WT, AE4 knockout, and AE2 knockout.
5. Compare the held out predictions with all available phenotype targets.
6. Record full residuals and not only secretion flow.

For fairness, use a common calibration protocol across candidates. If candidates have different numbers of free parameters, record that and use a complexity penalty or minimal parameterization preference.

Generate a table with, at minimum

- mechanism ID
- free AE4 parameters
- calibration targets used
- calibration error
- WT baseline closure error
- predicted AE4 KO flow ratio
- predicted AE2 KO flow ratio
- predicted intracellular chloride change
- predicted pH change if available
- conservation residuals
- thermodynamic validity

Output:
`analysis/12_ae4_mechanism_reconstruction/model_comparison.md`

Machine readable outputs go in
`results/12_ae4_mechanism_reconstruction/`.

## Agent E. Adversarial audit

Independently audit Agents A to D.

Look specifically for

- leakage of knockout data into calibration
- arbitrary rescaling of AE4 density
- hidden non AE4 retuning
- unit compensation masquerading as mechanism
- nonphysical reversal behaviour
- incorrect Na/K partition in mass balances
- overfitting due to extra mechanism parameters
- dependence of the conclusion on one arbitrary convention
- use of publication values contradicted by primary experimental data

Repeat key calculations independently where possible.

Output:
`analysis/12_ae4_mechanism_reconstruction/adversarial_audit.md`

No candidate may be declared successful until it survives this audit.

## Agent F. Structural mathematics and mechanism localization

Use the candidate results to explain *why* a mechanism succeeds or fails.

Do not merely report simulations.

Investigate, where applicable

- exact steady state balance identities
- how AE4 cation stoichiometry changes Na, K, Cl, and HCO3 balance signatures
- whether knockout flow depends primarily on chloride loading, membrane potential support, osmotic driving force, or another pathway
- limiting cases for Na only, K only, pooled cation, and cooperative gating
- whether a threshold or monotonicity result explains the required phenotype
- whether a minimal stoichiometric or kinetic feature can be identified analytically

Seek a statement of the form

> The knockout phenotype cannot occur unless AE4 has property X under assumptions Y.

or

> Property X is sufficient over a source supported parameter range.

Output:
`analysis/12_ae4_mechanism_reconstruction/structural_explanation.md`

---

# Iterative search protocol

The lead agent must not terminate after one round.

Perform the following loop until the answer is stable.

## Round 1. Source justified mechanisms only

Test all mechanisms directly supported or explicitly left viable by the primary experimental evidence.

## Round 2. Minimal structural interpolation

If no Round 1 candidate reproduces the held out phenotype, determine which balance or dynamical response is quantitatively deficient.

Introduce only minimal AE4 structural interpolations between source justified mechanisms, for example

- Na versus K weighting
- cation cooperativity strength
- transported versus regulatory cation fraction
- stoichiometric interpolation only when physically meaningful
- saturation strength

Do not vary unrelated channel, pump, water, or tight junction parameters in this round.

The purpose is to localize the missing AE4 property, not to optimize the phenotype blindly.

## Round 3. Chassis falsification test

If changing AE4 alone still cannot reach the experimental phenotype while maintaining WT closure and other data, quantify the best achievable phenotype over the entire defensible AE4 family.

Then identify the first non AE4 model component whose modification could plausibly account for the remaining discrepancy. Rank candidates mechanistically using exact balance pathways and source evidence.

Examples may include NKCC1 coupling, water permeability calibration, membrane potential closure, channel gating, or acid base coupling, but do not modify them until the AE4 only hypothesis has been exhausted.

At this stage the required answer is no longer `STOP`. It must be a localization result such as

> AE4 structure alone can account for at most X percent reduction on this chassis; the missing effect must enter through pathway Y.

## Round 4. Minimal joint correction only if required

If Round 3 localizes one additional component with strong evidence of mismatch between the historical and published models, test the smallest joint correction involving AE4 plus that one component.

Do not open a broad parameter search.

Repeat held out validation and adversarial audit.

---

# Decision criteria

A mechanism is considered a successful reconstruction only if all of the following hold.

1. The WT baseline remains internally consistent.
2. No non AE4 parameter was tuned to the knockout target in the primary comparison.
3. The AE4 knockout prediction is quantitatively compatible with the primary experimental reduction, using the experimental uncertainty when available.
4. AE2 knockout remains compatible with its experimental phenotype.
5. Other available phenotypes, especially intracellular chloride and pH, are not contradicted.
6. Conservation, electroneutrality, and thermodynamic checks pass.
7. The mechanism survives the adversarial audit.
8. A more complex mechanism is not preferred when a simpler mechanism performs equivalently.

Do not define success as reproducing exactly 30.000 percent.

---

# Required final answer

Create
`analysis/12_ae4_mechanism_reconstruction/final_answer.md`.

It must answer all of the following directly.

1. Can the experimentally observed AE4 knockout phenotype be reproduced by changing AE4 alone on the best available internally coherent chassis?
2. If yes, what is the **minimal structural AE4 property** required?
3. Which AE4 mechanisms fail, and why?
4. Is the published 2018 AE4 formulation sufficient, insufficient, or impossible to assess under a coherent parameterization?
5. Is the historical Na only formulation sufficient?
6. What role do Na and K play separately?
7. Is cation cooperativity required?
8. Is an alternative transport stoichiometry required?
9. Is cation transport itself required, or can cation dependent gating explain the phenotype?
10. If AE4 alone is insufficient, what is the quantified maximum effect attainable and which additional pathway is localized as missing?
11. Which conclusions are exact/analytical, which are numerical, and which depend on new modelling assumptions?
12. Does the result justify reopening the identifiability/hypothesis testing project?

The final answer must select one of these substantive conclusions.

- `AE4 MECHANISM IDENTIFIED`
- `AE4 FAMILY NONUNIQUE`
- `AE4 ALONE INSUFFICIENT; MISSING PATHWAY LOCALIZED`

Do not end with `STOP`.

---

# Outputs

Create at least

- `analysis/12_ae4_mechanism_reconstruction/experimental_evidence.md`
- `analysis/12_ae4_mechanism_reconstruction/chassis.md`
- `analysis/12_ae4_mechanism_reconstruction/ae4_candidate_family.md`
- `analysis/12_ae4_mechanism_reconstruction/model_comparison.md`
- `analysis/12_ae4_mechanism_reconstruction/adversarial_audit.md`
- `analysis/12_ae4_mechanism_reconstruction/structural_explanation.md`
- `analysis/12_ae4_mechanism_reconstruction/final_answer.md`
- new independent implementation under `src/ae4_mechanism_reconstruction/`
- comprehensive tests under `tests/`
- machine readable outputs under `results/12_ae4_mechanism_reconstruction/`
- reproducible figures under `figures/12_ae4_mechanism_reconstruction/`
- updated `docs/PROVENANCE.md`
- updated `docs/DECISIONS.md`
- updated `docs/RESULTS_LEDGER.md` only for independently checked results

---

# Prohibitions

- Do not edit `archive/`.
- Do not tune arbitrary whole model parameters until the knockout figure appears.
- Do not use the knockout phenotype in primary AE4 parameter fitting.
- Do not assume the 2018 paper is correct when it conflicts with conservation, units, primary biology, or reproducible implementation evidence.
- Do not assume the historical code is correct because it closes its fitted baseline.
- Do not report a parameter sweep as a mechanistic explanation.
- Do not hide failed candidates.
- Do not stop after the first negative result.
- Do not draft manuscript prose.

Commit all work to `codex/ae4-mechanism-reconstruction`. Do not merge to `main`.