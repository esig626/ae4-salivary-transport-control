# AGENTS.md

These instructions apply to all work in this repository.

## Scientific objective

Develop a new, reproducible analysis of AE4 and AE2 control of salivary acinar cell fluid secretion. The work must distinguish genuine network structure from parameter scaling and must not assume that the explanation in the historical manuscript is correct.

## Non negotiable provenance rules

1. Treat historical unpublished material as immutable evidence.
2. Do not copy unpublished prose, figures, derivations, or code from the historical material into the new manuscript or new source code.
3. The published 2018 AE4 model and other published literature may be used as scientific sources and must be cited in the relevant source map.
4. Newly implemented code must be written independently from the published equations and documented source by source.
5. If historical material is consulted to understand what was attempted, record that fact in `docs/PROVENANCE.md` or the relevant analysis README.
6. Do not attribute sole authorship to historical collaborative work. The new project must stand on newly derived results.

## Research workflow

Work through the numbered `analysis/` phases in order unless a documented reason requires otherwise.

Before advancing from one phase to the next

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
- Add regression tests for the baseline steady state before perturbation studies.
- Never interpret a raw derivative with respect to a dimensional parameter as a cross parameter measure of biological control.
- Report dimensionless control coefficients when comparing AE2 and AE4.
- For steady state sensitivities, verify numerical derivatives against implicit differentiation whenever feasible.
- Record solver tolerances and continuation settings.

## Claims discipline

Do not promote an observation to a biological mechanism without testing plausible alternatives. In particular, the old claim that the large AE4 derivative explains the knockout difference must be treated as a hypothesis to audit, not as an established fact.

## Manuscript discipline

Do not draft substantive manuscript sections until the relevant result is frozen in `docs/RESULTS_LEDGER.md`. The manuscript must clearly distinguish reproduction of prior published results from new findings.
