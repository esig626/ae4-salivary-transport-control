# Quick paper pivot

## Objective

Produce a compact mathematical paper on identifiability and discrimination of transporter mechanisms in a reduced pump leak model of epithelial secretion.

The paper should not be a sensitivity analysis paper and should not be organised around numerical continuation or AUTO. Continuation may be used only as a diagnostic tool when needed to verify the relevant steady branch.

## Core question

What can be inferred about transporter activities or competing transporter mechanisms from steady state measurements of a pump leak secretion model?

Starting from the published 2018 AE4 secretion model, derive independently a reduced steady state representation with unknown AE2 and AE4 activity parameters. Let the steady state satisfy

`F(u; G2, G4) = 0`

with state vector `u`, and let measurable outputs be collected in

`h(u) = (Q, Na_i, K_i, Cl_i, ...)`.

Study the parameter to observation map

`(G2, G4) -> h(u*(G2,G4))`.

## Mathematical targets

1. Establish local identifiability conditions using the implicit function theorem and the rank of the observation Jacobian.
2. Prove that a single scalar steady state observable such as secretion rate cannot generically identify two independent transporter activity parameters.
3. Characterise local observational equivalence sets for underdetermined measurement panels.
4. Determine minimal biologically plausible measurement panels that locally identify `G2` and `G4`.
5. Compare measurement panels using rank, conditioning, and geometry of the image of the parameter map. Conditioning is a numerical diagnostic, not the main result.
6. If justified, formulate two competing transporter mechanism classes and study whether their predicted observation sets intersect or separate.
7. Add a small statistical discrimination result only if it follows naturally from the geometry. Under an explicitly stated noise model, quantify how separation of the mechanism sets controls finite sample discrimination. Do not import the full composite Renyi programme unless it is genuinely required.

## Scientific discipline

The historical unpublished reduction may be inspected only to understand what was attempted. The new reduction, equations, code, proofs, and figures must be derived independently from published sources. Do not copy historical unpublished derivations, prose, code, or figures.

## Intended paper shape

A compact paper should ideally contain

1. published model and independent steady state reduction
2. identifiability theorem or proposition
3. observational equivalence result
4. minimal measurement result
5. numerical geometry for the salivary example
6. optional small mechanism discrimination section

## Stop conditions

Do not force the project if

- the published model cannot be independently reduced or reproduced
- every biologically plausible observation panel is trivially full rank and yields no useful structure
- the only result is a numerical condition number or conventional sensitivity calculation
- the discrimination section requires a much larger new statistical theory than the identifiability result warrants

The target is a small mathematical paper with at least one clean analytical result and a concrete physiological example.
