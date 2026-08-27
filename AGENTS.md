# AGENTS.md

These instructions apply to all work in this repository.

## Scientific objective

The current priority is a compact mathematical paper on identifiability and discrimination of transporter mechanisms in a reduced salivary pump leak model. The central question is what transporter information can and cannot be recovered from steady state physiological measurements.

This is not a sensitivity analysis project. Numerical continuation, AUTO, MATCONT, and parameter sweeps may be used as diagnostics but must not constitute the main result.

## Non negotiable provenance rules

1. Treat historical unpublished material as immutable evidence.
2. Do not copy unpublished prose, figures, derivations, or code from the historical material into the new manuscript or new source code.
3. The published 2018 AE4 model and other published literature may be used as scientific sources and must be cited in the relevant source map.
4. Newly implemented code and mathematical reductions must be derived independently from published equations and documented source by source.
5. If historical material is consulted to understand what was attempted, record that fact in `docs/PROVENANCE.md` or the relevant analysis README.
6. Do not attribute sole authorship to historical collaborative work. The new project must stand on newly derived results.

## Current research workflow

The quick paper plan in `docs/QUICK_PAPER_PLAN.md` supersedes the earlier sensitivity centred workflow.

The priority sequence is

1. independently reconstruct the published baseline model
2. independently derive the smallest useful steady state reduction
3. define parameter to observation maps for AE2 and AE4 activities
4. establish local identifiability and observational equivalence results
5. identify minimal biologically plausible measurement panels
6. compute observation geometry as a numerical illustration
7. add a small mechanism discrimination section only if the geometry supports it
8. conduct an objective paper viability audit before manuscript drafting

Before accepting a result

- save machine readable outputs in `results/`
- add or update tests where possible
- record accepted results in `docs/RESULTS_LEDGER.md`
- record important methodological choices in `docs/DECISIONS.md`
- distinguish reproduced results from new results

## Model implementation rules

- Prefer explicit units for every state, parameter, and flux.
- Preserve mass balance, charge conventions, and stoichiometry exactly as specified by the chosen published model.
- Keep parameter definitions separate from solver code.
- Keep observables such as fluid flow separate from state equations.
- Add regression tests for the baseline steady state before new analysis.
- Do not use raw dimensional sensitivities as scientific conclusions.
- For the steady state map `F(u;theta)=0`, verify implicit derivatives against numerical derivatives whenever feasible.
- Do not claim global identifiability from local rank calculations.
- Do not claim mechanism discrimination when the predicted observation sets overlap over the stated parameter classes.

## Mathematical standard

Prefer analytical results over numerical diagnostics. A parameter sweep, condition number, or continuation plot is not a theorem. Use the implicit function theorem, rank arguments, level set geometry, and related analytical tools where they genuinely apply. State clearly which conclusions are generic, which are model specific, and which are numerical only.

## Claims discipline

Do not promote an observation to a biological mechanism without testing plausible alternatives. Do not pad a weak result. If the identifiability structure is trivial or the proposed discrimination problem collapses, record that result and reassess the paper.

## Manuscript discipline

Do not draft substantive manuscript sections until the viability audit is complete and the main analytical result has been checked and frozen in `docs/RESULTS_LEDGER.md`. The manuscript must clearly distinguish reproduction of prior published results from new findings.
