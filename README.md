# AE4 Salivary Transport Control

This repository develops a new analysis of salivary acinar cell transport control. The published AE4 secretion model is the scientific starting point. Older unpublished material is retained only as a historical archive and provenance record. https://link.springer.com/article/10.1007/s11538-017-0370-6

## Central question

Why does loss of AE4 produce a substantial secretory phenotype while loss of AE2 produces little or no effect, and does that distinction survive a scale invariant analysis of transporter control?

The project will separate three possibilities.

1. The apparent AE4 dominance is a genuine structural property of the transport network.
2. It is partly or largely an artefact of parameter scaling in the older sensitivity calculation.
3. The phenotype depends specifically on AE4 monovalent cation coupling and can therefore discriminate among alternative microscopic transport mechanisms.

## Repository map

- `archive/` immutable historical material and provenance notes
- `literature/` published sources, evidence tables, and literature audit
- `model/` canonical equations, parameters, observables, and source mapping
- `src/` newly written reusable model and analysis code
- `analysis/` numbered research analyses in execution order
- `tests/` regression, conservation, and reproducibility tests
- `data/` curated or derived data that may legally be version controlled
- `results/` machine readable outputs used for figures and tables
- `figures/` publication figures and their source notes
- `manuscript/` the new manuscript only
- `docs/` research plan, provenance, decisions, and result ledgers
- `prompts/` reproducible Codex task prompts

## Execution order

1. Inventory the historical material and map every equation and parameter to a published source.
2. Reconstruct and verify the published baseline model with newly written code.
3. Reproduce the AE2 and AE4 knockout curves.
4. Replace raw parameter derivatives with dimensionless control coefficients.
5. Decompose secretion control through the steady state Jacobian and transport network.
6. Test alternative AE4 sodium and potassium coupling mechanisms.
7. Determine which mechanisms are identifiable from available observables.
8. Generate prospective whole cell predictions for experimentally characterised AE4 variants.
9. Perform robustness and uncertainty analyses.
10. Freeze the validated result set before manuscript drafting.

## Scientific rule

No result from the historical material is automatically a result of this project. New claims must be independently derived from published sources and newly written code, then entered in `docs/RESULTS_LEDGER.md` with a reproducibility pointer.

## Current status

Repository structure established. Scientific validation has not yet begun.
