# Codex task. Quick identifiability and discrimination paper

Work in repository `esig626/ae4-salivary-transport-control` on branch `codex/identifiability-discrimination`.

This branch was created from completed Phase 00 commit `2f71f23098037037836853fd96814a0e286cc0af`. Treat the Phase 00 inventory, source map, reconstruction blockers, provenance log, decision log, and published inconsistencies as established inputs. Do not redo Phase 00 unless a specific error is found.

Read `AGENTS.md`, `README.md`, `model/SOURCE_MAP.md`, `model/states.md`, `model/parameters.md`, `model/fluxes.md`, `model/equations.md`, `model/observables.md`, `analysis/00_inventory/inventory.md`, `analysis/00_inventory/open_questions.md`, `docs/PROVENANCE.md`, and `docs/DECISIONS.md` before doing anything else.

Do not edit any file under `archive/`.

## Goal

Develop the shortest defensible route to a publishable mathematical result on identifiability and discrimination of transporter mechanisms in the salivary pump leak model.

The default target is a compact paper. However, if the analysis reveals a substantially stronger mathematical structure, do not stop at the routine local identifiability result. Perform the mandatory escalation scan below and push the stronger direction far enough to determine whether it is real.

This is not a sensitivity analysis task. Do not produce a parameter sensitivity paper. Do not make AUTO, MATCONT, or bifurcation diagrams the main result. Continuation may be used only as a diagnostic if required to verify a steady branch or an inference-relevant singularity.

## Scientific starting point

Use the published 2018 AE4 model as the primary scientific source. Independently derive from the published equations a reduced steady state system in which AE2 and AE4 activity parameters are explicit unknowns. Historical unpublished material may be inspected only to understand what was previously attempted and must never be copied into the new derivation, source code, proofs, prose, or figures.

Phase 00 documented published inconsistencies, including the volume-balance sign conflict in Eq. (9). Before using any disputed equation, resolve the issue from conservation laws and published context if possible. Record the resolution and its justification. If it cannot be resolved unambiguously, carry both admissible conventions far enough to determine whether the main identifiability conclusions depend on the choice.

Let the independently derived steady state problem have the form

`F(u; theta) = 0`, with `theta = (G2, G4)`,

and define candidate observation maps using physiologically meaningful steady state outputs such as secretion rate and intracellular ion concentrations.

## Required work

### 1. Independent reduction and baseline check

- Derive the smallest useful steady state system from the published model.
- Document every algebraic elimination and assumption.
- Resolve or explicitly bracket every Phase 00 reconstruction blocker that affects the reduction.
- Verify the reduction numerically against the full published model at the baseline state.
- Put derivation notes in `analysis/10_identifiability_discrimination/reduction.md`.
- Put new implementation in `src/` and tests in `tests/`.

### 2. Identifiability mathematics

For the parameter to observation map

`theta -> H(theta) = h(u*(theta))`,

use the implicit function theorem to derive

`D_theta u* = -F_u^{-1} F_theta`

where valid, and hence

`D_theta H = h_u D_theta u* + h_theta`.

Then establish rigorous local identifiability statements.

At minimum, determine

- why one scalar steady state observable cannot generically locally identify two independent transporter activity parameters
- the local dimension of observational equivalence sets when the observation Jacobian is rank deficient
- which candidate two or more observable panels give full local rank at the physiological baseline
- whether any biologically plausible panel remains rank deficient for structural rather than accidental numerical reasons

State analytical propositions or theorems precisely and prove them where possible. Do not substitute numerical rank checks for proofs when a general rank or dimension argument is available.

### 3. Minimal measurement problem

Systematically examine a small biologically plausible set of observables, for example

- secretion rate `Q`
- intracellular Na
- intracellular K
- intracellular Cl
- cell volume if available in the chosen reduction
- luminal ion concentrations only if they are retained by the independently reconstructed model

Identify minimal observation panels that locally distinguish `G2` and `G4`.

Separate

- generic mathematical statements
- model specific analytical statements
- purely numerical findings

### 4. Observation geometry

Compute and plot the image of a physiologically relevant `(G2,G4)` parameter domain under selected observation maps.

Look specifically for

- self intersections
- near overlaps
- folds in the observation image
- one dimensional equivalence curves under scalar observation
- separation induced by adding a second observable
- boundaries or singular sets where the rank changes

Save machine readable data in `results/10_identifiability_discrimination/` and figures in `figures/`.

### 5. Mechanism discrimination only if justified

If the geometry supports it, formulate a small discrimination problem between two explicit transporter mechanism classes. Keep this section narrow unless the escalation scan demonstrates that the testing problem itself contains the stronger result.

Permitted route

- define two model classes by a mechanistic difference supported by published biology
- map each class into observation space
- determine whether the two image sets intersect
- under a simple explicit noise model, derive or compute the error of discriminating the classes from finite noisy measurements

Do not import the full composite Renyi finite blocklength theory unless the simple geometry is inadequate and there is a clear mathematical gain.

## Mandatory escalation scan

Before concluding that the project is only a short local identifiability paper, actively search for stronger structure. This is a required part of the task, not an optional discussion section.

Create `analysis/10_identifiability_discrimination/escalation_scan.md` and investigate the following possibilities.

### A. Global identifiability or exact equivalence

Ask whether any selected observation map is globally injective on the physiologically admissible parameter region, or whether exact nontrivial equivalence sets exist globally. Search for monotonicity, algebraic invariants, exact elimination arguments, topology, sign structure, or other analytical devices that can turn a local rank statement into a global result.

If a global theorem looks plausible, attempt it. Do not claim global identifiability from numerical sampling alone.

### B. Perturbation-induced distinguishability

Treat experimentally controllable quantities as possible inputs, such as extracellular ion concentrations, agonist level, or transporter inhibition, but only where the published model supports such variation.

Ask whether two parameters or mechanism classes that are observationally equivalent at baseline become distinguishable under one additional controlled perturbation. Search for a mathematically sharp statement of the form

- unidentifiable at one condition
- identifiable under two conditions

or a minimal number of perturbations required for discrimination.

This is higher priority than producing more baseline plots.

### C. Optimal experimental design

If perturbations change identifiability, determine whether one can formulate and solve a small design problem that chooses the perturbation or observation panel to maximise separation, rank robustness, determinant, smallest singular value, or another defensible discrimination criterion.

Do not turn this into a generic numerical optimisation exercise. Look for analytical structure or a low dimensional design result.

### D. General pump leak theorem

Check whether the reduced equations expose a statement that is not specific to AE2 and AE4 but follows from pump leak network structure, stoichiometry, conservation, or membrane orientation.

Examples worth testing include

- conditions under which a scalar secretion observable generically leaves a positive dimensional mechanism equivalence set
- conditions under which adding one ion concentration breaks that equivalence
- structural rank criteria for transporter activities in a class of pump leak networks
- exact relations between stoichiometric signatures and identifiable output panels

Do not force a general theorem. Escalate only if the algebra clearly supports one.

### E. Hypothesis testing

If two mechanism classes generate separated or partially overlapping observation sets, determine whether this gives a genuinely nontrivial testing problem.

Search for

- a positive minimum separation under a useful metric or divergence
- least favourable parameter pairs
- perturbations that change zero separation into positive separation
- a finite sample error bound that is materially stronger than a trivial Gaussian classification calculation

Only then consider importing composite hypothesis testing machinery. The testing theory must solve a problem created by the pump leak geometry, not be attached artificially.

### F. Inference-relevant singularities

If folds, branch collisions, or bifurcations appear, ask whether they create or destroy identifiability or mechanism separation. Use AUTO or continuation only if such a singularity is directly relevant to inference. Do not pursue a bifurcation catalogue for its own sake.

### Escalation decision rule

At the end of the scan, classify the project as one of

1. `QUICK PAPER` — local identifiability plus minimal measurement result is the strongest defensible contribution.
2. `PUSH HARDER` — a stronger theorem, perturbation result, general network statement, or nontrivial testing problem has concrete evidence and should become the main project.
3. `STOP` — the result is routine, already known, or too weak to justify publication.

For every `PUSH HARDER` classification, state exactly what was found, what remains unproved, and the next three technical steps. Do not use speculative language as evidence.

## Targeted novelty audit

Once candidate analytical results are known, perform a focused literature search on the exact claims. Search identifiability of pump leak and epithelial transport models, structural identifiability of transporter networks, experiment design for ion transport models, and mechanism discrimination in physiological ODE models.

Record the search strategy and closest prior results in `analysis/10_identifiability_discrimination/novelty_audit.md`.

Do not claim novelty merely because the exact salivary model has not been analysed before. The mathematical claim itself must survive the audit.

## Paper viability audit

Create `analysis/10_identifiability_discrimination/viability.md` answering, without advocacy

- What is the strongest new theorem or proposition?
- Is it generic or specific to this model?
- What is the strongest new physiological conclusion?
- What is merely a numerical illustration?
- Did the escalation scan uncover anything stronger than local identifiability?
- Does the novelty audit support publication?
- Is the result enough for a short mathematical biology paper?
- Should the project be classified as `QUICK PAPER`, `PUSH HARDER`, or `STOP`?
- What additional work, if any, is strictly necessary before writing?

If the result is too weak, say so explicitly rather than padding the analysis.

## Outputs

Create or update

- `analysis/10_identifiability_discrimination/reduction.md`
- `analysis/10_identifiability_discrimination/theory.md`
- `analysis/10_identifiability_discrimination/observables.md`
- `analysis/10_identifiability_discrimination/escalation_scan.md`
- `analysis/10_identifiability_discrimination/novelty_audit.md`
- `analysis/10_identifiability_discrimination/viability.md`
- appropriate new files in `src/`
- appropriate regression and rank tests in `tests/`
- machine readable outputs in `results/10_identifiability_discrimination/`
- reproducible figures in `figures/`
- `docs/RESULTS_LEDGER.md` only for results that have been checked
- `docs/DECISIONS.md` for material methodological decisions

## Prohibitions

- No conventional sensitivity study as a headline result.
- No parameter sweep presented as a theorem.
- No AUTO or MATCONT paper.
- No copying unpublished historical derivations, code, text, tables, or figures.
- No manuscript drafting until the viability audit is complete and at least one analytical result is validated.
- No claim of global identifiability from a local Jacobian rank calculation.
- No claim of mechanism discrimination if the corresponding observation sets overlap under the stated parameter domains.
- No forced hypothesis testing section if the geometry does not justify it.
- No forced general pump leak theorem if the model gives no structural evidence for one.

## Acceptance criteria

The task is complete when the repository contains

- an independently reconstructed reduced model based on published sources and Phase 00 inputs
- a mathematically correct local identifiability analysis
- a minimal observation result
- reproducible observation geometry
- a completed escalation scan covering global identifiability, perturbations, design, general pump leak structure, testing, and inference-relevant singularities
- a targeted novelty audit
- an objective viability classification of `QUICK PAPER`, `PUSH HARDER`, or `STOP`

Commit all work to `codex/identifiability-discrimination`. Do not merge to `main`.