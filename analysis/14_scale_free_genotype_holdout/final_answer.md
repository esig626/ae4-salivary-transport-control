# Task 14 final answer

**The frozen model did not produce a valid AE4 holdout prediction.** All ten AE4 resting continuation branches reached 5% expression but failed at exact deletion. The independent implementation reproduced every failure. All 60 exact zero attempts across the two implementations exhausted the declared 3000 evaluation budget; their scaled residuals were 0.0184222 to 0.0232276 against the required 1e-7. The failed iterates are not equilibria and were not used as initial states.

Every AE2 resting branch succeeded. The blind stage completed 30 WT and 30 AE2 production Radau trajectories over all ten roots and calcium 0.10, 0.25 and 0.50 uM. All 12 available representative BDF confirmations passed, with maximum relative discrepancy 4.65e-6 below the inherited 1e-4 tolerance. The largest dynamic accounting ratio was 0.00657 of its allowed tolerance; bulk charge error was at most 8.11e-13 fmol. Independent resting states agreed within 3.46e-11 in relative state difference wherever both implementations retained a root. All 105 blind tests passed.

The blind model predicts an AE2 null to WT total secretion ratio from **1.00031795 to 1.00160080**, a secretion increase of **0.0318% to 0.1601%**. Every root and calcium remains in the CSV files and routing family summaries. This is qualitatively consistent with the reported absence of a detected AE2 secretion difference. There is no reported equivalence margin, so it is not a statistical proof of equivalence.

After the checkpoint was pushed, the experiment revealed **35 ± 4.7% lower AE4 null total secretion (SEM)**, with comparable early flow followed by a sustained deficit. The model has no valid AE4 total ratio, minute ratio, early or sustained integrated ratio, emergence statistic, or exact null ionic prediction to compare with those data. Signed and absolute errors, interval membership and time course RMSE therefore remain explicitly unavailable in all 30 AE4 comparison rows. No root or calcium value was dropped or selected to improve agreement.

The small AE2 effect survives this task. **AE4 versus AE2 specificity is not established**, because the AE4 effect and the effect contrast are unavailable. The outcome is neither an AE4 magnitude match nor a demonstrated mismatch of a valid predicted curve. It is a numerical continuation failure under the frozen rules. This does not prove that an exact null equilibrium cannot exist.

The inherited absolute one SMG mapping discrepancy remains unresolved. The Task 13C best 0.50 uM result was 5.405993 uL/min against the inherited 9 to 10 uL/min mapping, and five roots failed the minute observation shape criterion. Those results are preserved. A common multiplicative scale would cancel from valid genotype ratios, but it cannot supply the missing exact AE4 resting states.

All model equations, parameters, original roots, regulation, calcium values and saved blind trajectories remained unchanged after reveal. Any investigation of why the exact deletion continuation fails belongs in a later task. This task makes no new calibration, mechanism or biological validation claim.



Blind checkpoint: `607ac21279d821b95d44ba63f35df27078503cf6`. Prediction contract SHA256: `e0fa0badc347024daf0fe7cac61ef4ab3e6e0e12ca44a60b37a598fa00031b0a`.

**GENOTYPE CONTINUATION OR NUMERICS PREVENT A VALID HOLDOUT TEST**
