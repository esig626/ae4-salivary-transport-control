# Task 47 final checkpoint

Scientific execution is complete. No further mechanism search is opened.

The final article report is `DYNAMIC_POTASSIUM_RECYCLING_REPORT.md`. It explains the negative result for a binding K recycling ceiling in the inherited model, the complete dynamic K and Na balances, the effect of membrane driving voltage, and why uncertain physiological capacities do not justify an arbitrary correction. The fixed prediction fails the central experimental magnitude. The absence of a direct pump voltage law and the fixed bath NKCC reversal limitation remain explicit.

| Published milestone | Commit |
| --- | --- |
| Equation and provenance audit | `0df7a90a27a263757dbffc3917d00482749617b5` |
| Frozen baseline dynamics | `0aee25dcc3afd6945bd6203bc9fad9b00a3a6a1f` |
| Prediction freeze before comparison | `56ab6a00cb4b75c70e5857d2c58e7d46731d4d5f` |
| Fixed prediction validation | `56305f56c374e74a922e240016b1d4c349f00c96` |
| Article report and final verification | The commit containing this checkpoint |

No production equation, scientific parameter, initial state, main file, manuscript file or older analysis output changed. All work is published on `analysis/task-47-dynamic-potassium-recycling-reconstruction`, descending from `acd554eb430ae25968ef56e3a9acdd816ca4e688`.

Completed numerical work: three production integrations; no stationary or CBM solve; no optimisation, parameter sweep or altered parameter trajectory; 24 local derivative checks at saved states. Independent verification reconstructed channel, pump, NKCC, concentration and storage balances for all 1,806 saved states without importing the production model. Maximum absolute independent discrepancy was below `2.57e-13` in the respective reported units. The 88 prediction file hashes remain unchanged.

The WT summary serialization error was recovered from saved data without repeating its integration. Its unrecovered solver counters are not fabricated. The two other trajectories each passed 934 monitored state checks. All three pass the inherited 600 s gates, with the separate voltage and long time interpretation limits described in the report.

Two article figures were generated from saved outputs, exported as PNG and SVG, and visually inspected. The six lossless trajectory/state archives have verified compressed and uncompressed hashes. Report equations were checked for valid mathematical escapes and all figure links resolve. The optional trajectory replay wrapper was syntax checked but not run; it is an explicit reproduction tool, not additional scientific work.

Final result: null cumulative deficit 3.85748630%; 5% cumulative deficit 3.41757761%; sustained K channel efflux increase 18.06709986%; sustained pump turnover increase 1.85899125%. No fitting or retuning occurred before or after the prediction freeze. Whether real salivary K recycling has a tighter limit remains unresolved without independent functional kinetic data.
