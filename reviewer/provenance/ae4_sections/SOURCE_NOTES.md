# Sources and scope

## Repository snapshot

The repository was read at `c44160e1849828cde887a5e24618335c88f97af7`. The manuscript points to scientific source revision `764e32a648ee3b69265364788df5a59f76cb2843`.

The principal numerical source is `analysis/52_chloride_reservoir_final_test/PRODUCTION_SUMMARY_52D.md` together with its complete central case records under `output/reporter_recovery_52D/cases`. The source discipline is recorded in `manuscript/claims_and_evidence.md`, `manuscript/README.md` and the repository research ledgers.

No repository file was changed, and no new model trajectory or fit was run in preparing these sections.

## Figures

Figure 2 combines the retained central fluid flow plot with the recorded cumulative deficit at ten times. Its experimental point is the reported ten minute reduction of 35 ± 4.7 percent, with the error bar showing the standard error. There is no invented experimental trace or constant experimental band across time.

Figure 3 combines the retained intracellular chloride and outward driving force plots. Both belong to the central 2.32 nS paired simulation.

Figure 4 retains the exact chloride accounting plot. Its independent flux quadrature is distinct from the one second trapezoidal totals used for the headline secretion result.

Figure 5 is newly drawn from the three recorded cumulative secretion totals. It shows proportional changes relative to zero auxiliary conductance. The connected points are not a continuous parameter sweep or fitted response curve.

Figure 6 retains the seven class control pH plot from page 14 of the supplied revised manuscript. The vector crop and source file hash are recorded in `source_manifest.json`. Failed traces stop at their recorded failure times.

The five retained vector plots are associated with the repository figure names `fig05_fluid_flow`, `fig06_chloride`, `fig07_reservoir_budget`, `fig08_driving_force` and `fig03_class_ph`. Their numerical origins and historical versus central status are not interchanged.

## New deductions from the supplied equations

The donor allocation export identity, the reduced carrier scaling ambiguity and the resulting relation between local control coefficients are explicit mathematical deductions made during drafting. The first two were checked symbolically in `check_algebra.py`.

The full state Jacobian restriction follows from the exact constant intracellular charge balance. It is not a numerical eigenvalue computation or a claim of stability.

The logarithmic conductance contrasts are finite secants calculated from recorded secretion totals. They are not converged local sensitivity derivatives. No all transporter ranking is reported.

## Reported summaries

The 5 and 10 percent extra supply cases, alternate onset projection and single stimulus cases are preserved in `manuscript/reported_continuation_summary.md` and the supplied revised manuscript. Their underlying continuation trajectories were not available in this snapshot. The named local completion branch did not resolve through the repository connector.

Accordingly, the reported values retain their limited status in the text and tables. The IPR only pH failure near 488 seconds has no valid final secretion value.

## Wording corrected during drafting

The intracellular chloride heading no longer asserts that chloride alone determines the entire response. The initialisation also changes potassium and inorganic carbon.

The adrenergic heading describes regulation and stimulus dependence rather than claiming that the completed simulations isolate a causal increase in Ae4's contribution. Removing the beta input removes both Ae4 recruitment and the auxiliary apical conductance, and the isolated IPR record fails before the endpoint.

Neither the cumulative chloride accounting nor the difference between the experimental and simulated secretion percentages is interpreted as an additive allocation of causal responsibility.

The source class tests are separate constructions with local uptake feedback and class specific control resting states. They are not seven versions of the final matched supply pair and are not formal statistical hypothesis tests.

## Prior check records

`data/prior_results_document_checks.json` and `data/results_source_checks.json` are retained provenance records from the earlier Results package. They are not new validation of these sections. The present document check is `data/document_checks.json`, and the new static arithmetic and symbolic check is `data/algebra_and_arithmetic_checks.json`.
