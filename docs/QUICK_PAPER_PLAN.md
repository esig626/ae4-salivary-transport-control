# Quick paper pivot

## Objective

Produce a compact mathematical paper on identifiability and discrimination of transporter mechanisms in a reduced pump leak model of epithelial secretion.

The compact paper is the default route because it can be completed quickly. It is not a cap on the project. If the analysis exposes a genuinely stronger theorem or inference structure, the project should be escalated rather than artificially kept small.

This is not a sensitivity analysis paper and should not be organised around numerical continuation or AUTO. Continuation may be used only as a diagnostic tool when needed to verify the relevant steady branch or an inference-relevant singularity.

## Core question

What can be inferred about transporter activities or competing transporter mechanisms from steady state measurements of a pump leak secretion model?

Starting from the published 2018 AE4 secretion model and the completed Phase 00 source map, derive independently a reduced steady state representation with unknown AE2 and AE4 activity parameters. Let the steady state satisfy

`F(u; G2, G4) = 0`

with state vector `u`, and let measurable outputs be collected in

`h(u) = (Q, Na_i, K_i, Cl_i, ...)`.

Study the parameter to observation map

`(G2, G4) -> h(u*(G2,G4))`.

## Quick paper mathematical targets

1. Establish local identifiability conditions using the implicit function theorem and the rank of the observation Jacobian.
2. Prove that a single scalar steady state observable such as secretion rate cannot generically identify two independent transporter activity parameters.
3. Characterise local observational equivalence sets for underdetermined measurement panels.
4. Determine minimal biologically plausible measurement panels that locally identify `G2` and `G4`.
5. Compare measurement panels using rank and geometry of the image of the parameter map. Conditioning is a numerical diagnostic, not the main result.
6. If justified, formulate two competing transporter mechanism classes and study whether their predicted observation sets intersect or separate.
7. Add a statistical discrimination result only if it follows naturally from the geometry.

## Mandatory escalation targets

Before stopping at the quick paper, actively test whether the reduced model supports something stronger.

### Global identifiability

Can any useful observation panel be proved globally injective on the physiological parameter region? Conversely, can exact global equivalence sets be characterised?

### Perturbation-induced distinguishability

Can an experimentally controllable second condition convert a nonidentifiable baseline problem into an identifiable one? A clean result of the form `one condition insufficient, two conditions sufficient` would supersede a routine local rank paper.

### Optimal experimental design

If perturbations matter, is there a small analytical design problem for choosing the perturbation or measurement panel that maximises separation or rank robustness?

### General pump leak structure

Does the algebra reveal a rank, stoichiometric, conservation, or membrane-orientation result that extends beyond AE2 and AE4 to a class of pump leak networks?

### Hypothesis testing

Do competing mechanism classes generate separated observation sets, least favourable pairs, or perturbation-dependent transitions from overlap to separation? Only if this happens should finite sample testing theory be pushed harder.

### Inference-relevant singularities

If folds or bifurcations occur, do they create or destroy identifiability? If not, do not spend time cataloguing them.

## Decision rule

At the end classify the project as

- `QUICK PAPER` if the strongest result is the compact identifiability and minimal measurement theorem
- `PUSH HARDER` if there is concrete evidence for a stronger global, perturbational, general-network, design, or testing result
- `STOP` if the results are routine or already known

## Novelty requirement

Once the strongest candidate mathematical result is identified, perform a focused literature search on that exact claim. Do not infer novelty merely because the salivary example is new.

## Scientific discipline

The historical unpublished reduction may be inspected only to understand what was attempted. The new reduction, equations, code, proofs, and figures must be derived independently from published sources. Do not copy historical unpublished derivations, prose, code, or figures.

Phase 00 identified published inconsistencies, including the Eq. (9) volume-balance sign conflict. These must be resolved from conservation laws and published context where possible. If an ambiguity remains, test whether the analytical conclusion depends on the convention rather than silently selecting one.

## Intended quick paper shape

A compact paper should ideally contain

1. published model and independent steady state reduction
2. identifiability theorem or proposition
3. observational equivalence result
4. minimal measurement result
5. numerical geometry for the salivary example
6. optional mechanism discrimination section
7. a clear statement that the mandatory escalation scan found no stronger route, if that is the outcome

## Stop conditions

Do not force the project if

- the published model cannot be independently reduced or reproduced
- every biologically plausible observation panel is trivially full rank and yields no useful structure
- the only result is a numerical condition number or conventional sensitivity calculation
- the discrimination section requires a much larger new statistical theory without a geometry-driven reason
- the novelty audit shows the main mathematical result is already standard

The target is a small mathematical paper with at least one clean analytical result and a concrete physiological example, unless the evidence justifies pushing harder.
