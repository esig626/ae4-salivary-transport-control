# Task 26: exact AE4 null diagnostic

The formal exact zero AE4 solution increases secretion in every available matched comparison; residual AE4 activity is not required for this mathematical response.

Only 5 of ten null roots were numerically closed; all closed null roots lie outside the inherited resting domain. The other 5 roots remain unresolved. Secretion results are formal diagnostics for the closed subset, not an admissible full panel prediction.

## Design and provenance

The authorised preserved Task 25 workspace supplied the baseline. The requested historical commit was absent locally. All 89 source and input files recorded in the Task 24 contract matched their hashes. One saved WT JSON differed from the Task 25 input hash only in elapsed execution time; its full scientific content was verified against a local copy with the expected hash. No network, Git or publication operation was used during computation.

All ten retained roots were included. Each existing Task 24 WT core at Ca 0.10 was reused exactly for all three subsequent calcium levels and freshly audited without another WT root solve. Exact null REST was sought independently from matched Task 24 WT and 5% states. Expression was zero during every null solve. The 5% experiment was not rerun.

The twenty bounded direct null attempts all reached the inherited 3000 evaluation limit. Each was continued with a direct hybr solver using an invertible fixed initial Jacobian preconditioner, which changes no zero of the equations. For roots where neither conditioned solve closed, an additional direct attempt used reciprocal cell volume, with diagnostic ranges extending to 1e8 pL, pH 9 and TIC 200 mM. These numerical search ranges did not replace the inherited admissibility bounds. Every root was assessed with the original full RHS, current, charge, conservation, positivity, Jacobian rank and domain gates. All attempts are retained in warm_start_attempts.csv.

The two starting points agreed on a closed equilibrium for five roots in the AE4NA05 group. Those equilibria all exceeded the inherited volume and pH limits. The five AE4NA20 roots did not yield a closed equilibrium; their later attempts moved towards extremely large volumes and failed numerical closure. This is a numerical nonconvergence finding, not a proof that no finite or admissible equilibrium exists. No secretion result was used to select a root or change a model parameter.

Stimulation was run from all ten WT rests and the five numerically closed null rests, giving 45 actual trajectories. Fifteen requested null trajectories were not run because no closed null rest was available. Their cases remain explicitly marked as not run in trajectory_summary.csv. No approximate resting iterate was used for stimulation.

The equations, all fitted parameters, AE4 routing and stoichiometry, chloride allocation, NKCC1, NHE1, CO2 permeabilities, AE2, pumps, channels, acid base chemistry and water transport were unchanged. Stimulation retained CCh + IPR, R1 and N1 regulation, Radau, rtol 1e−7, amount and regulator atol 1e−10, volume atol 1e−12, maximum step 2 s, and all 602 inherited output times through 600 s. Secretion is the inherited trapezoidal integral of lumen outflow in pL per model cell, including basal flow. Full arrays were held in memory only.

## A and B. Does exact deletion reduce secretion?

The formal exact zero AE4 solution increases secretion in every available matched comparison; residual AE4 activity is not required for this mathematical response.

The following ranges apply only to the five closed null equilibria, all of which fail the inherited physiological domain checks. The remaining five roots have no reported null secretion value. These are formal model contrasts, not a valid full panel phenotype estimate.

| Ca (µM) | Available pairs | Null / WT ratio | Secretion change (%) | Median change (%) | Saved 5% change (%) on these roots |
| --- | --- | --- | --- | --- | --- |
| 0.10 | 5/10 | 1.035761 to 1.04297 | 3.57612 to 4.29702 | 3.97281 | 0.604721 to 0.685868 |
| 0.25 | 5/10 | 1.082329 to 1.1031 | 8.23293 to 10.31 | 9.2454 | 0.84388 to 0.946005 |
| 0.50 | 5/10 | 1.102411 to 1.127295 | 10.2411 to 12.7295 | 11.4084 | 0.935956 to 1.03688 |

Ranges describe retained roots, not confidence intervals or independently sampled animals. All thirty requested pairs are accounted for in comparison.csv; missing null quantities are left blank.

## C. NKCC1 replacement at REST

Across the five closed equilibria, lost AE4 chloride loading is 0.00453117 to 0.00454862 fmol/s. NKCC1 chloride loading changes by 0.00624868 to 0.00634789 fmol/s. The replacement fraction is 1.37766 to 1.39888, or 137.766 to 139.888% (median 138.8%).

The fraction is the signed increase in NKCC1 chloride influx divided by the positive AE4 chloride influx lost. NKCC1 carries two chloride ions per cycle, checked against the production source ledger. Values are not clipped to zero or one. This is a flux accounting result; other fluxes and intracellular conditions also change, and no NKCC1 intervention was performed to establish its isolated causal contribution.

## D. Carbon, pH and volume in the five closed null rests

| Quantity | Exact null range | Null minus matched WT range |
| --- | --- | --- |
| pH | 8.001803 to 8.013217 | 1.115765 to 1.123079 |
| HCO3 (mM) | 76.07526 to 78.15203 | 71.19469 to 73.22299 |
| TIC (mM) | 77.41189 to 79.5102 | 71.73064 to 73.78011 |
| Volume (pL) | 7.004518 to 8.556098 | 5.558474 to 7.103446 |

The inherited intracellular volume domain is 0.5 to 3 pL and pH domain is 6.2 to 8.0. Neither was raised to accept these equilibria. Changes are paired by root. All required resting fluxes and membrane potentials are retained in resting_states.csv, with numerical and physiological flags. Carbon and pH changes at 600 s are in comparison.csv and summary.json for the available pairs.

## E. Difference from the established 5% phenotype

The available secretion contrasts retain the increase seen in the saved Task 24 5% experiment. Exact deletion nevertheless causes a qualitatively different resting homeostasis problem: extreme swelling and alkalinisation in the five closed equilibria, and unresolved closure in the other five roots. A near neutral secretion ratio alone would conceal this failure. Only existing text outputs were read for the 5% comparison.

## F. Numerical and physiological admissibility

All inherited REST gates pass for 10/20 requested states: the ten WT controls. There are 15/20 numerically closed states. Of these, 5 null equilibria fail the inherited domain bounds. The other 5 null rests remain unresolved.

45/45 trajectories actually run pass the production numerical checks. 15 additional requested trajectories were not run. Only 30 trajectories start from rests passing all inherited gates, and 0/30 requested pairs have both WT and null rests admissible.

All seven AE4 post expression source quantities are exactly zero at every available null rest and all recorded samples of every completed null trajectory. Maximum absolute null source magnitude: 0 fmol/s. Agreement of the two warm starts for all closed null roots: True.

For the fifteen closed resting states, the largest scaled independent RHS is 1.76328e-10, amount RHS 9.7141e-14 fmol/s, volume RHS 1.76328e-14 pL/s, omitted RHS 9.7141e-14 fmol/s and current residual 3.55429e-26 A. The largest conservation residual divided by its inherited tolerance is 5.59552e-05 at REST and 0.000581313 during stimulation. All completed traces retain positive core quantities and nonnegative flow.

The raw failed solver coordinates, residuals and messages are retained; a solver termination message alone is never treated as steady state closure. No claim of global root uniqueness, nonexistence or stability follows from these numerical searches.

## Interpretation

This task does not establish an admissible exact null secretion phenotype across the ten retained roots. In the closed subset, the equations can sustain secretion without residual AE4 transport, so residual transporter self rescue alone does not explain that formal response. However, exact deletion exposes a much larger resting homeostasis failure. That failure must be resolved scientifically before treating the formal secretion ratios as a physiological prediction. The production model was not changed after observing these results.

## Reproduction and text outputs

Run `PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=src python -m modern_full_model.run_ae4_exact_null --workers 6` from the repository root. Run focused tests with `PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python -m unittest discover -s tests -p test_ae4_exact_null.py -v`. The runner refuses to overwrite existing results. The REST resume option continues preserved numerical attempts and refuses to rerun started stimulation.

All twenty requested REST cases, sixty stimulation cases and thirty comparisons have explicit records. Closed resting cores and complete RHS and conservation diagnostics are UTF8 JSON fields inside resting_states.csv. No state is fabricated for an unresolved rest. Source and input hashes, solver settings and numerical amendments are recorded in contract.json; validation and focused test outcomes are in verification.json. There are no trajectory binaries, archives or publication files.

Focused validation: all 10 Task 26 tests passed, with no failures or errors. These tests validate the accounting of unsuccessful rests and unrun cases as well as the completed numerical results; passing tests do not turn the inadmissible null rests into physiological predictions.
