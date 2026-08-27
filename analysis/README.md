# Analysis order

The numbered directories define the scientific execution order.

- `00_inventory` provenance and source inventory
- `01_baseline_reproduction` published baseline reconstruction
- `02_knockout_reproduction` AE2 and AE4 perturbation curves
- `03_dimensionless_control` scale invariant sensitivity analysis
- `04_jacobian_decomposition` structural decomposition of secretion control
- `05_mechanism_variants` alternative AE4 cation coupling models
- `06_identifiability` distinguishability of microscopic mechanisms
- `07_variant_predictions` prospective predictions for AE4 variants
- `08_robustness` uncertainty, parameter robustness, and final stress tests
- `10_identifiability_discrimination` independent steady-state reduction,
  identifiability theory, escalation, novelty, and viability audit
- `11_forensic_reconstruction` exhaustive historical-file audit, MATLAB
  lineage, paper/code concordance, translated reproduction, blocker
  reassessment, and conditional reconstruction decision
- `12_ae4_mechanism_reconstruction` multi-agent primary-evidence audit,
  fixed-chassis AE4 candidate reconstruction, iterative exclusion, and
  missing-pathway localization

A phase is complete only when its acceptance criteria are met and the corresponding results are entered in `docs/RESULTS_LEDGER.md`.

Phase 10 can also terminate with a documented `STOP` when the published record
does not support the requested model-specific numerical claims. In that case,
only checked structural identities and diagnostics enter the results ledger.

Phase 11 can retain that `STOP` after inspecting the historical implementation
when exact 2018-code identity, the central phenotype, or a physiological
parameter-to-observation map cannot be validated. Historical implementation
evidence remains distinct from published scientific evidence.

Phase 12 replaces the generic `STOP` with the substantive chassis-scoped
conclusion `AE4 ALONE INSUFFICIENT; MISSING PATHWAY LOCALIZED`. Its local README
defines the evidence boundary and reading order. This does not reopen the
physiological identifiability project because no jointly source-admissible WT
and knockout model has yet been obtained.
