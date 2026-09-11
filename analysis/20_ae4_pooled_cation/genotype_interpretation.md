# Task 20 genotype interpretation

The pushed WT checkpoint `a2902286c148c886ee2d8188dabbd020d0a2d32f` contains
20 structural infeasibility decisions and zero accepted WT candidates. Its
branch, tree and manifest blob were remotely verified before this report
accessed legacy genotype results. Every frozen hash remains checked by the
post-freeze entry point.

There are no eligible pooled WT dynamics, AE4 near-loss continuations or AE2
deletion continuations. Each of the three requested dynamic result tables
contains all 60 root/condition/calcium rows, with
`NOT_RUN_STRUCTURALLY_INFEASIBLE_WT`. Totals, ratios, solver outcomes and
numerical gates are null. The tables record a feasibility exclusion before
integration; they do not present a failed ODE solve or a zero secretion value.

The required four-condition comparison contains 120 rows: 60 reused legacy
rows and 60 unavailable pooled comparisons. The Task 18 result tables are
checked against their inherited Git blobs and result digests. Their applicable
source trajectories and saved 10% payloads are verified too. Ratios are
checked against the corresponding saved WT and genotype integrals. No legacy
trajectory is recomputed and no 30% condition is newly evaluated.

| Architecture and WT allocation | Root/calcium pairs | AE4 5% / WT | AE2 loss / WT |
|---|---:|---:|---:|
| Legacy, inherited loading | 30 | 1.082286917–1.216630071 | 1.000317950–1.001600801 |
| Pooled, matched inherited loading (P0) | 0 eligible of 30 | unavailable | unavailable |
| Legacy, 10% loading | 30 | 1.186165043–1.489915201 | 1.000260099–1.000997292 |
| Pooled, 10% loading (P10) | 0 eligible of 30 | unavailable | unavailable |

Legacy AE4 near-loss increases secretion in every reused root/calcium pair.
Legacy AE2 effects are comparatively small (less than 0.161% in these
comparators). Neither fact establishes the response of the pooled model.

For pooled P0 and P10, AE4-loss secretion direction, resting chloride direction,
pH changes, Na/K changes, compensating NKCC1/NHE1/AE2/pump responses and AE2
neutrality are all unevaluable. `state_flux_summary.csv` preserves their null
fields alongside distinctly labeled inherited-target and reference-only
diagnostics. There are no accepted pooled resting or genotype states hidden
in that table.

The sign obstruction is consistent across all ten roots and both provenance
families. It occurs at the fixed resting state, before any of the three
stimulated calcium conditions can be evaluated. No pooled phenotype robustness
claim across calcium is possible. The approximately 35% experimental deficit
does not enter a fit, target, acceptance band, ranking or parameter revision.

Removing independent cation slip is established at the level of the declared
transport law. Whether it repairs secretion cannot be tested within this
fixed-state experiment. Further work would need a separately authorized
mechanism or fixed-state/chemistry premise change; no such change is made here.
