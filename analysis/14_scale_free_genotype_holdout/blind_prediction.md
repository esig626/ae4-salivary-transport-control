# Blind Task 14 prediction

This report states the numerical predictions without an experimental comparison.

Exact resting continuation succeeded in 0/10 AE4 deletion branches and 10/10 AE2 deletion branches. All ten AE4 branches reached expression 0.05, then failed at exact zero under the inherited 3000 evaluation limit and numerical gates. All three local attempts at exact zero failed in each branch. The independent implementation reproduced the failure on all ten AE4 branches. Agreement between implementations at the retained positive expression steps does not establish that no exact null equilibrium exists; it establishes that the declared procedure did not obtain one.

No AE4 null stimulated trajectory is licensed. AE4 total, minute, early and sustained secretion ratios, the emergence statistic, ionic differences and the AE4 minus AE2 contrast are therefore unavailable for every calcium value. Failed exact zero iterates are retained as numerical attempts and are never treated as equilibria.

The full 90 row primary matrix contains 60 simulated WT and AE2 trajectories and 30 explicitly blocked AE4 entries. Production Radau completed 60 valid trajectories. BDF ran 12 representative WT and AE2 confirmations. All available solver comparisons pass: True. The independent resting audit ran 12 continuations, covering both representatives and every failed branch.

The complete valid AE2 predictions are summarised below. Reduction is defined as one minus the null to WT total secretion ratio; negative reduction means increased secretion.

| Calcium uM | Routing | Total ratio range | Median | Early ratio range | Sustained ratio range |
| :--- | :--- | :--- | :--- | :--- | :--- |
| 0.10 | AE4NA05 | 1.00125461 to 1.00139193 | 1.00132237 | 1.00121781 to 1.00136500 | 1.00127094 to 1.00140394 |
| 0.10 | AE4NA20 | 1.00031795 to 1.00038676 | 1.00035475 | 1.00031274 to 1.00038489 | 1.00032020 to 1.00038757 |
| 0.25 | AE4NA05 | 1.00143647 to 1.00155982 | 1.00150203 | 1.00145856 to 1.00160622 | 1.00142615 to 1.00153795 |
| 0.25 | AE4NA20 | 1.00039264 to 1.00046157 | 1.00042955 | 1.00039894 to 1.00047719 | 1.00038983 to 1.00045453 |
| 0.50 | AE4NA05 | 1.00148882 to 1.00160080 | 1.00154958 | 1.00154296 to 1.00168030 | 1.00146306 to 1.00156267 |
| 0.50 | AE4NA20 | 1.00041646 to 1.00048210 | 1.00045142 | 1.00042969 to 1.00050582 | 1.00041047 to 1.00047128 |

Every primary trajectory and every continuation attempt remains in the saved artifacts. The complete per root and calcium results are in the CSV files; full time grids, states, ionic quantities and accounting diagnostics are saved separately. Absolute cellular and mapped gland flows remain secondary diagnostics. The Task 13C files have not changed.

The declared numerical outcome before any experimental reveal is: **GENOTYPE CONTINUATION OR NUMERICS PREVENT A VALID HOLDOUT TEST**. This is a failure to obtain a valid AE4 prediction under the frozen rules, not an experimental rejection of a computed AE4 secretion curve.

Prediction contract SHA256: `e0fa0badc347024daf0fe7cac61ef4ab3e6e0e12ca44a60b37a598fa00031b0a`.
