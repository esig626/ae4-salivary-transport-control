# AGENTS.md

These instructions apply to all work in this repository.

## Scientific objective

The current priority is a compact mathematical paper on identifiability and discrimination of transporter mechanisms in a reduced salivary pump leak model. The central question is what transporter information can and cannot be recovered from steady state physiological measurements.

This is not a sensitivity analysis project. Numerical continuation, AUTO, MATCONT, and parameter sweeps may be used as diagnostics but must not constitute the main result.

The compact paper is the default target, not a ceiling. If the reduced model reveals a stronger analytical structure, the work must explicitly test whether it supports global identifiability, perturbation-induced distinguishability, a general pump leak network theorem, optimal experimental design, a nontrivial hypothesis testing problem, or an inference-relevant singularity.

## Non negotiable provenance rules

1. Treat historical unpublished material as immutable evidence.
2. Do not copy unpublished prose, figures, derivations, or code from the historical material into the new manuscript or new source code.
3. The published 2018 AE4 model and other published literature may be used as scientific sources and must be cited in the relevant source map.
4. Newly implemented code and mathematical reductions must be derived independently from published equations and documented source by source.
5. If historical material is consulted to understand what was attempted, record that fact in `docs/PROVENANCE.md` or the relevant analysis README.
6. Do not attribute sole authorship to historical collaborative work. The new project must stand on newly derived results.

## Phase 00 status

Phase 00 is complete at commit `2f71f23098037037836853fd96814a0e286cc0af` and is part of the history of this branch. Treat its inventory, source map, reconstruction blockers, provenance record, decision log, and published inconsistencies as inputs. Do not redo the inventory unless a concrete error is found.

## Current research workflow

The quick paper plan in `docs/QUICK_PAPER_PLAN.md` and the executable prompt in `prompts/10_identifiability_discrimination_quickpaper.md` supersede the earlier sensitivity centred workflow.

The priority sequence is

1. independently reconstruct the published baseline model using the completed Phase 00 source map
2. resolve or bracket published inconsistencies that affect the reduction
3. independently derive the smallest useful steady state reduction
4. define parameter to observation maps for AE2 and AE4 activities
5. establish local identifiability and observational equivalence results
6. identify minimal biologically plausible measurement panels
7. compute observation geometry as numerical support
8. perform the mandatory escalation scan
9. perform a targeted novelty audit on the strongest candidate results
10. classify the project as `QUICK PAPER`, `PUSH HARDER`, or `STOP`
11. draft only after the main analytical result is checked

Before accepting a result

- save machine readable outputs in `results/`
- add or update tests where possible
- record accepted results in `docs/RESULTS_LEDGER.md`
- record important methodological choices in `docs/DECISIONS.md`
- distinguish reproduced results from new results

## Model implementation rules

- Prefer explicit units for every state, parameter, and flux.
- Preserve mass balance, charge conventions, and stoichiometry exactly as specified by the chosen published model unless a published inconsistency is being explicitly resolved.
- Keep parameter definitions separate from solver code.
- Keep observables such as fluid flow separate from state equations.
- Add regression tests for the baseline steady state before new analysis.
- Do not use raw dimensional sensitivities as scientific conclusions.
- For the steady state map `F(u;theta)=0`, verify implicit derivatives against numerical derivatives whenever feasible.
- Do not claim global identifiability from local rank calculations.
- Do not claim mechanism discrimination when the predicted observation sets overlap over the stated parameter classes.
- If a published equation is internally inconsistent, do not silently choose a convention. Document the conservation argument and test whether the result depends on the choice.

## Mathematical standard

Prefer analytical results over numerical diagnostics. A parameter sweep, condition number, or continuation plot is not a theorem. Use the implicit function theorem, rank arguments, level set geometry, elimination, monotonicity, structural rank, topology, and related analytical tools where they genuinely apply.

Do not stop at a routine local rank calculation if the algebra indicates stronger structure. In particular, actively test whether

- local equivalence extends to an exact or global equivalence manifold
- a second experimental condition breaks an otherwise structural nonidentifiability
- transporter stoichiometry gives a general rank criterion for a class of pump leak networks
- a small experimental design problem has an analytical solution
- separated mechanism image sets create a genuinely nontrivial testing problem
- a fold or bifurcation changes inferential distinguishability

State clearly which conclusions are generic, which are model specific, and which are numerical only.

## Escalation discipline

The escalation scan is mandatory, but escalation itself is evidence based.

Use `PUSH HARDER` only when there is concrete mathematical structure beyond the quick paper and a plausible route to proof. Do not invent a general theory because it sounds interesting. Conversely, do not bury a stronger result merely to finish quickly.

For any proposed stronger result, perform a focused literature audit of the mathematical claim before treating it as novel.

## Claims discipline

Do not promote an observation to a biological mechanism without testing plausible alternatives. Do not pad a weak result. If the identifiability structure is trivial or the proposed discrimination problem collapses, record that result and reassess the paper.

## Manuscript discipline

Do not draft substantive manuscript sections until the viability audit is complete and the main analytical result has been checked and frozen in `docs/RESULTS_LEDGER.md`. The manuscript must clearly distinguish reproduction of prior published results from new findings.
