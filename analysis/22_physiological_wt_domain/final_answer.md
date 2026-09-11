# Task 22 scientific result

**No admissible WT exists within the predeclared Task 22 domain and the
unchanged production chassis. The exact failure class is
`PROVEN_STRUCTURAL_CONTRADICTION`.** The specific conflict is physiological
and thermodynamic: the resting carbon/alkalinity balances and declared
luminal domain force an intracellular CO2 condition that makes productive
pooled AE4 chloride loading impossible.

This is a domain-limited result about the frozen chassis. It does not
disprove pooled/no-slip AE4 as a general mechanism.

## Why the exclusion is stronger than optimizer failure

The allowed luminal pH, TIC and volume, together with the inherited tiny
lumen buffer, imply `TA_l - TIC_l < -0.00645 mM`. The production carbon-minus-
alkalinity identities then require positive bath-to-cell CO2 flux even when
all four relevant resting residuals take their most permissive allowed
values. Thus intracellular CO2 is below bath CO2. With the predeclared
intracellular pH, chloride and Na/K limits, the dimensionless pooled AE4
affinity obeys

`A_AE4 < -0.0764726613570... < -0.0764`.

This contradicts the required strictly positive affinity and productive
chloride flux. All ten immutable root backgrounds satisfy the analytic
certificate. It does not depend on optimizer status, the capacity upper
limits, the osmotic screens or the anti-cancellation limit. The complete
derivation is in [structural_contradiction.md](structural_contradiction.md),
with machine-readable per-root bounds in
`results/22_physiological_wt_domain/structural_certificate.json`.

## Search and decision

| Item | Result |
|---|---:|
| Inherited WT root starts | 10 |
| Task 21 resting-witness starts | 7 |
| Additional WT-only starts | 10 |
| Completed and preserved bounded attempts | 27 |
| Admissible resting WT states | 0 |
| Complete admissible WT ensemble | Empty |
| Task 22 WT 600 s integrations | 0; no admissible resting initialization |
| Task 22 genotype calculations | 0 |
| Historical genotype regressions during Task 22 | 0 |

All 27 endpoints failed resting RHS closure and positive pooled AE4
affinity/loading; four also missed the intracellular Na band. All stayed
inside the declared capacity box. The best maximum scaled resting RHS was
approximately 0.00188291, against the unchanged 1e-7 limit. The raw
optimization failures are classified separately from the domain-wide proof.
No failed endpoint is presented as an accepted calibration or selected by
its later phenotype. The lexicographically ranked admissible set is empty.

The evidence includes every attempted full conserved state and parameter
vector, a complete decision table, intracellular/luminal physiology,
production osmolarity, all capacity folds, signed source ledgers and
anti-cancellation diagnostics. WT dynamic status records explicitly state
why the calcium panel was not run. There are no Task 22 AE4-low, AE2-loss
or matched genotype comparison files.

## Required checkpoints and preservation

| Gate | Pushed and remotely verified commit |
|---|---|
| Complete physiology/capacity contract, before any new optimization | [a1ec05b706dfaf4c1574abe3b33ea39f0195d472](https://github.com/esig626/ae4-salivary-transport-control/commit/a1ec05b706dfaf4c1574abe3b33ea39f0195d472) |
| Complete WT evidence and empty ensemble | [2b73ee293d135809904ea2a5aec08724a6520b9c](https://github.com/esig626/ae4-salivary-transport-control/commit/2b73ee293d135809904ea2a5aec08724a6520b9c) |

The remote WT manifest blob is
`0cb675ad1ae68136eef91cf6106d97af4069017a`. The complete remote tree matches
the staged local tree exactly. All 54 WT evidence files in the manifest
and all 79 inherited source/witness hashes were verified. The receipt and
verification transcript are `wt_checkpoint_receipt.json` and
`validation_checkpoint.txt` under the Task 22 results directory.

The genotype guard returns `STOP_NO_ADMISSIBLE_WT`. No genotype calculation,
historical genotype regression, phenotype comparison or genotype-specific
capacity refit was performed. The bounds and numerical tolerances were not
widened. `POOLED_CATION_112_NO_SLIP`, inherited scientific equations,
Task 21 evidence and `archive/` are unchanged. Archive preservation is
verified by its unchanged Git tree; the inherited sparse local store's
unavailable historical blobs are documented in [methods.md](methods.md).

## Validation and review disposition

**46 focused tests and 76 broad WT-only tests passed, plus 10 subtests.**
The focused suite includes production replay of all 27 saved endpoints,
every new physiology/capacity/osmotic/anti-cancellation screen, the proof
identities and conservative margins, and checkpoint-ordering protections.
The broader allowlist covers the inherited WT numerical, dimensional,
acid-base, transport and regulatory infrastructure. The full historical
suite was not invoked because it includes prohibited genotype regressions.
Exact commands, exit codes and transcripts are preserved in the results.

Task 22 is complete at the required empty-WT stop. The branch intentionally
contains Task 21's unmerged ancestry and is submitted to `main` for review.
Do not merge automatically. No claim about AE4-loss or AE2-loss secretion
follows from this empty WT ensemble.
