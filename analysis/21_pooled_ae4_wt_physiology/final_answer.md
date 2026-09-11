# Task 21 result for scientific review

**WT calibration remains numerically unresolved. No fully admissible WT
calibration was found, so no Task 21 AE4-loss or AE2-loss phenotype was evaluated.**
This is not a proof of structural or thermodynamic infeasibility.

| WT result | Count |
|---|---:|
| Independent inherited seeds | 10 |
| Saved solver outputs, including failed attempts | 170 |
| Initial resting witnesses passing production gates | 7 |
| Resting witnesses passing serialized-parameter replay | 6 |
| Completed 600 s WT trajectories | 21 / 21 |
| Trajectories passing all numerical/conservation gates | 0 / 21 |
| Fully admissible WT calibrations retained | 0 |
| Task 21 genotype calculations or genotype-specific refits | 0 |

The six replay-valid resting witnesses attain Cl_i = 50.10 mM and pH_i = 6.91,
with strictly positive pooled AE4 affinity and chloride loading, same-direction
Na/K transport, and the declared Na/K/volume screens. Thus releasing the exact
legacy state permits productive resting fits. The search does not certify
global capacity optimality: several constrained optimization and joint
tie-break attempts failed, and all search decisions are preserved.

The resulting states are extreme in quantities outside the declared
intracellular screens. Across the seven tested witnesses, lumen Na is about
1,637–1,771 mM, lumen Cl about 1,392–1,518 mM, and lumen volume about 63–69 pL.
The three replay-valid capacity-profile witnesses require maximum absolute
log folds of roughly 15.5–16.0, corresponding to departures of millions-fold.
Three other witnesses retain effectively zero capacities in some coordinates
and have much worse log-fold objectives. They were kept as feasibility
witnesses, not discarded using any later phenotype. The predeclared
physiological screens and inherited numerical tolerances were not altered.

Every WT trajectory reaches 600 s with all 602 requested samples, positive
core states and valid regulatory fractions. Nevertheless, every trajectory
fails both the inherited conservation gate and the 1e-9 fmol state-charge gate.
The maximum normalized conservation error per trajectory ranges from 6.32 to
414.96 times its allowed tolerance. The largest state-charge drift is
8.93e-4 fmol. One initial resting pass also fails conservation when its saved
parameters are replayed as Python floats; both diagnostics remain available.
Large balanced fluxes expose numerical sensitivity, but the evidence does not
establish that this is the only possible WT branch. Solver success and a good
Cl/pH fit therefore do not establish an admissible WT calibration.

## Checkpoints and held-out tests

The physiology contract was pushed before optimization at
[`9820a6b62b29cd9895b9020b85ae89ad7d44b3fb`](https://github.com/esig626/ae4-salivary-transport-control/commit/9820a6b62b29cd9895b9020b85ae89ad7d44b3fb).
The complete WT evidence, including all failed attempts and trajectories, was
pushed and remotely verified at
[`317b2781ac094a3e9d7001ab7da6cf1b7cd02e40`](https://github.com/esig626/ae4-salivary-transport-control/commit/317b2781ac094a3e9d7001ab7da6cf1b7cd02e40).
The receipt verifies the exact manifest blob and tree. The frozen admissible
ensemble is empty, and the genotype authorization function refuses execution.

AE4 5%, AE2-loss and matched-comparison tables explicitly say
`NOT_RUN_NO_ADMISSIBLE_WT`; their secretion totals and ratios are unavailable.
There is no conclusion about AE4-loss direction, AE2 specificity or agreement
with the experimental loss magnitude. The experimental magnitude was not
compared or used to tune WT.

## Validation and protocol deviation

The final focused Task 20/21 regressions pass: **36 passed**, including the
actual remote-checkpoint guard, production replay, no-slip signs, physiology
bounds, complete calcium panels, matched arithmetic and archive preservation.
No production mechanism, chemistry, regulation, solver tolerance or file in
`archive/` was changed.

**A broad regression run executed historical genotype tests before the Task 21
WT checkpoint push, contrary to the requested ordering for any phenotype
evaluation.** All Task 21 optimization had already completed; no Task 21
genotype was calculated and no WT bound, capacity, state or selection was
changed from those tests. This deviation is explicitly recorded in
[`validation_scope_deviation.json`](../../results/21_pooled_ae4_wt_physiology/validation_scope_deviation.json).
The broad suite is not passing: 576 passed, 21 failed and 9 errored, with one
skip. Its failures include unavailable inherited files and Git ancestors in
the partial-object checkout. The full transcript is retained; these failures
are not hidden by the focused passing result.

See [methods](methods.md), the [calibration table](../../results/21_pooled_ae4_wt_physiology/wt_calibration_table.csv),
[resting flux diagnostics](../../results/21_pooled_ae4_wt_physiology/resting_flux_diagnostics.csv),
[capacity folds](../../results/21_pooled_ae4_wt_physiology/capacity_folds.csv), and
[WT dynamic results](../../results/21_pooled_ae4_wt_physiology/wt_dynamic_results.csv).

Scientific review is required before merge. This branch records an unresolved
WT numerical result and stops the held-out phenotype experiment at its WT gate.
