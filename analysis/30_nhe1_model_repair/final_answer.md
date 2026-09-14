# Task 30 final answer

Run status: **STOPPED_AT_NUMERICAL_BUDGET** because the stationary residual evaluation limit was exhausted.

The published 2018 NHE1 law was implemented as a named evaluator and explicitly selected in both Task 30 model builders. The historical `tanh` evaluator remains selectable but was not used for the repaired resting states.

| Background | WT REST | Exact AE4 null REST | WT Na (mM) | Null Na (mM) | WT pH | Null pH | Null Cl (mM) |
|---|---|---|---:|---:|---:|---:|---:|
| R09 | admissible after WT scalar calibration | attempt incomplete at budget | 11.565526 |  | 6.872665 |  |  |
| R10 | not reached | not reached |  |  |  |  |  |

## Direct answers

1. **Implementation:** Yes. The evaluator reproduces three independently calculated Eq. 26 values to machine precision, retains both squared proton terms, and preserves electroneutral 1 Na inward to 1 H outward stoichiometry.
2. **WT REST:** R09 recovered an admissible stationary state only after the permitted WT only scalar calibration. The published `G_NHE1 = 0.0305 fmol/s` gave a numerically valid state at pH 6.766175, outside the declared WT pH band. The frozen calibrated value was `0.145419726007 fmol/s`; it gave pH 6.872665, Na 11.565526 mM, cell volume 1.447604 pL, full stationary closure, and no capacity violation. R10 was not reached.
3. **Exact AE4 null REST:** Not determined. The first R09 exact null solve was still in progress when the global residual evaluation ceiling was reached. There is no accepted null state and no basis for calling it admissible or failed.
4. **Null sodium:** Not determined. The incomplete solve must not be interpreted as either ordinary sodium or a recurrence of the Task 29 extreme values.
5. **Null pH:** Not determined. No accepted state is available for comparison with the 6.89 diagnostic centre.
6. **AE2 and NKCC1 after AE4 removal:** Not determined because the exact null stationary state did not complete.
7. **Dynamics:** None ran. The prerequisite admissible matched WT and exact null REST pair was unavailable, so no secretion direction is reported.
8. **Later chloride allocation:** Chloride allocation remains a separate later phase. No chloride allocation search or chloride fitting was performed here.

The complete available R09 WT state, transporter flux ledger, water fluxes, voltages, osmolarities, and numerical and conservation residuals are recorded in the CSV and JSON outputs. Missing exact null and R10 values are deliberately blank rather than converted to zero.

Budget used: 4 of 12 stationary calls, 10,000 of 10,000 stationary residual evaluations, 3 of 8 scalar G evaluations, and 0 of 4 stimulated integrations. One numerical worker and one numerical library thread were used.
