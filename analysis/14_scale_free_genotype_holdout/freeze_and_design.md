# Task 14 freeze and design

The scientific state is commit `2f54e7c4f87b6da746d3a8427c7bdca40a77c3de`, tree `338e2891183546734ff9db81611a459136346dce`. Branch intake `9a1c3fb4e469fa016983b79b274981b1cef21b18` changes only AGENTS.md and the Task 14 prompt. Source files and all ten root payloads are hash verified. The Task 13B final manifest and Task 13C frozen manifest hashes are retained in the prediction contract.

The immutable prediction contract was written at 2026-09-09T11:39:34.306509+00:00, before any Task 14 genotype result existed. All ten roots, AE4 and AE2 exact deletion, calcium 0.10, 0.25 and 0.50 uM, the R1 member `R1_G125_P21_PROSE_TREFERENCE`, and the CCH_IPR 600 s protocol were declared together. Resting calcium is 0.058 uM. The NKCC1 calcium normalisation remains 0.10 uM with multiplier 1.75; changing the input amplitude does not change that parameter.

Expression continuation uses the inherited 21 points, three local starts, ten state coordinates, zero fitted parameters and fixed OTHER osmoles. All attempted and retained roots are saved. Selection uses only distance from the previous root divided by the declared coordinate span. The maximum allowed connected step is 0.25. No null dynamics may start from a WT resting state or an unconverged exact deletion attempt.

Production Radau uses relative tolerance 1e-7, amount absolute tolerance 1e-10 fmol, volume tolerance 1e-12 pL, regulatory tolerance 1e-10 and maximum step 2 s. BDF uses the same tolerances. Each pair uses the identical 602 point grid containing 0, the onset right limit 1e-6 s, and every integer second to 600. The inherited relative solver agreement tolerance is 1e-4. Both nulls and their WT denominator are checked for the lexicographic representative of each routing family at all three calcium values. A disagreement expands checks to the entire affected routing family and calcium. A missing exact resting state is recorded as unavailable, not as a Radau versus BDF disagreement.

The representatives are:

* AE4NA05: `N_ABS_NKCC_S4_PHIGH_KHIGH_AE4NA05_P1_H1:B00`.
* AE4NA20: `N_ABS_NKCC_S4_PHIGH_KHIGH_AE4NA20_P2_H1:B00`.

The independent resting implementation audits both representatives for both deletions and every anomalous branch. No independent resting audit changes a production initial state. Conservation, current closure, charge, positivity, exact deletion and full resting numerical rank remain required.

The paired metrics are total, instantaneous minute, cumulative minute, early 0 to 180 s and sustained 180 to 600 s ratios; late minus early ratio; AE4 and AE2 reduction fractions and their contrast; and resting and endpoint ion, pH and volume differences. Minute collection data, if available, are compared with integrated minute bins rather than instantaneous endpoint flow. Time course RMSE is unweighted across available matching relative points. No uncertainty interval, observation factor, time shift or other fitted transformation is introduced.

The inherited absolute whole gland discrepancy remains unchanged: the Task 13C best 0.50 uM result reached 5.405993 uL/min under the existing one SMG mapping. Five roots also failed its minute observation shape criterion. These diagnostics are preserved and are nonblocking for paired genotype quantities. A common multiplicative scale cancels from the ratios. No genotype specific effective cell count is estimated.

Every root and calcium remains in every summary, including failures. Exact experimental magnitude and time course remain sealed until the completed blind artifacts and tests have been committed and pushed. The later evaluator reads the saved predictions and writes only separately named comparison artifacts. No model or trajectory change is permitted after reveal.


Prediction contract SHA256: `e0fa0badc347024daf0fe7cac61ef4ab3e6e0e12ca44a60b37a598fa00031b0a`.
