# Codex task. Quick identifiability and discrimination paper

Work in repository `esig626/ae4-salivary-transport-control` on a new branch named `codex/identifiability-discrimination`.

Read `AGENTS.md`, `README.md`, `docs/QUICK_PAPER_PLAN.md`, `docs/PROVENANCE.md`, and `model/SOURCE_MAP.md` before doing anything else.

Do not edit any file under `archive/`.

## Goal

Develop the shortest defensible route to a publishable mathematical result on identifiability and discrimination of transporter mechanisms in the salivary pump leak model.

This is not a sensitivity analysis task. Do not produce a parameter sensitivity paper. Do not make AUTO, MATCONT, or bifurcation diagrams the main result. Continuation may be used only as a diagnostic if required to verify the steady branch.

## Scientific starting point

Use the published 2018 AE4 model as the primary scientific source. Independently derive from the published equations a reduced steady state system in which AE2 and AE4 activity parameters are explicit unknowns. Historical unpublished material may be inspected only to understand what was previously attempted and must never be copied into the new derivation, source code, proofs, prose, or figures.

Let the independently derived steady state problem have the form

`F(u; theta) = 0`, with `theta = (G2, G4)`,

and define candidate observation maps using physiologically meaningful steady state outputs such as secretion rate and intracellular ion concentrations.

## Required work

### 1. Independent reduction

- Derive the smallest useful steady state system from the published model.
- Document every algebraic elimination and assumption.
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

Save machine readable data in `results/10_identifiability_discrimination/` and figures in `figures/`.

### 5. Mechanism discrimination only if justified

If the geometry supports it, formulate a small discrimination problem between two explicit transporter mechanism classes. Keep this section narrow.

Permitted route

- define two model classes by a mechanistic difference supported by published biology
- map each class into observation space
- determine whether the two image sets intersect
- under a simple explicit noise model, derive or compute the error of discriminating the classes from finite noisy measurements

Do not import the full composite Renyi finite blocklength theory unless the simple geometry is inadequate and there is a clear mathematical gain.

### 6. Paper viability audit

Create `analysis/10_identifiability_discrimination/viability.md` answering, without advocacy

- What is the strongest new theorem or proposition?
- Is it generic or specific to this model?
- What is the strongest new physiological conclusion?
- What is merely a numerical illustration?
- Is the result enough for a short mathematical biology paper?
- What additional work, if any, is strictly necessary before writing?

If the result is too weak, say so explicitly rather than padding the analysis.

## Outputs

Create or update

- `analysis/10_identifiability_discrimination/reduction.md`
- `analysis/10_identifiability_discrimination/theory.md`
- `analysis/10_identifiability_discrimination/observables.md`
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

## Acceptance criteria

The task is complete when the repository contains an independently reconstructed reduced model, a mathematically correct local identifiability analysis, a minimal observation result, reproducible numerical geometry, and an objective paper viability assessment.

Commit all work to `codex/identifiability-discrimination`. Do not merge to `main`.
