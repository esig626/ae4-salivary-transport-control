# Reproducibility coverage

## Complete archived principal records

The three combined stimulation cases include the original model dependency snapshot, full precision onset inputs, complete stored state and diagnostic trajectories, flux quadrature, per case summaries and numerical checks. They were copied using their existing Git objects, not regenerated, refitted, interpolated or extracted from figures.

The source revision is `c44160e1849828cde887a5e24618335c88f97af7`. The scientific files originate in `analysis/52_chloride_reservoir_final_test/`. The complete central records are the reporter recovery 52D cases. Original historical files remain unchanged elsewhere in the repository.

The new plotting and portable launch scripts were checked for Python syntax during packaging. They have not undergone a fresh full numerical run in this publication step. The archived reports document the original computations, not a new independent execution by this packaging process.

## Results not yet deposited as complete reproducible calculations

The current manuscript gives a 41 point continuation over supply multipliers 1 to 2, a crossing at 1.4913598957, a local tangent, complex step differentiation and detailed numerical verification. The corresponding full source, trace and verification outputs were not located in the uploaded `ae4_sections.zip` or the inspected main repository snapshot. The retained paired constructor permits only 1, 1.05 and 1.10. It is not the later continuation implementation.

The uploaded package also gives reported CCh only and IPR only outcomes without the complete corresponding trajectory records. Those summaries are retained under `provenance/ae4_sections/data/reported_summaries.json` and are explicitly labelled as summaries.

Before claiming that all current manuscript figures and numerical checks are reproducible, deposit the original continuation code and its dependencies, the complete continuation trace and root/tangent verification records, the original script producing the current compensation figure, and the complete single stimulus output records. Do not reconstruct a smooth compensation curve from the handful of reported values and present it as the original numerical trace.

## Earlier archive and figure numbering

`ae4_sections.zip` contains an older manuscript package. This reviewer folder preserves its selected data, original summary plotting script and source notes. It does not present its old figure numbers as current manuscript figure numbers. In particular, the old Fig_5 concerns auxiliary conductance ratios, whereas the current manuscript Fig_5 is the nonlinear compensation response. The seven source class pH comparison is outside the present principal analysis.

The SHA256 digest of the original uploaded archive is `58d122f5c4acdd9bf8d9375f6c058fdf52ec557a856d091ad6792748af944286`.
