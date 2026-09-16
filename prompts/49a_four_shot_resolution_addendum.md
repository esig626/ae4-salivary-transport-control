# Task 49A addendum: four-shot bounded resolution waterfall

This addendum OVERRIDES any conflicting instruction in `prompts/49_camp_vrac_secretory_branch_reconstruction.md`.

## Hard budget

There are only four Codex scientific execution opportunities remaining, including this Task 49 run. Treat compute, context and human review as scarce. Do not spend a run on a one-question negative result if additional independent work can be completed safely in the same workspace after publishing that negative checkpoint.

The repository is the cumulative scientific memory. Before proposing or implementing anything, verify against the exclusion ledger and prior Tasks 12, 19, 31, 36-42, 46-48. Do not rediscover those results.

## Governing principle

Task 49 is a sequential, source-constrained reconstruction, not a mechanism search and not a multi-parameter fit to the AE4 saliva deficit.

Two experimentally demonstrated failures are already known and must be treated separately:

1. **Stimulated beta/cAMP architecture failure:** the inherited AE4-null equations erase beta input and therefore force CCh and CCh+IPR to be identical.
2. **Established-genotype resting-state failure:** the inherited shared equations do not make the measured WT, AE4-KO and AE2-KO resting Cl/pH states jointly stationary.

Use this single Codex run to obtain the maximum defensible information on BOTH failures, with an immutable checkpoint between them.

## Stage I: stimulated cAMP branch, VRAC ONLY

The first reconstruction is strictly the independently documented TMEM16A-independent cAMP/volume-regulated apical anion pathway.

### Binding restrictions

- Do **not** add or modify beta-dependent NKCC1 recruitment in Stage I.
- Keep the Palk/Benjamin NKCC1 core and the existing inherited NKCC stimulation exactly unchanged.
- Keep AE4, AE2, NHE1, NBC, pump, K channels, TMEM16A, paracellular transport, calcium input, geometry and water laws unchanged.
- Do not introduce any direct AE4-expression multiplier on the new conductance.
- Do not use AE4-KO secretion, AE4-KO combined-stimulus uptake, or the 35% deficit to calibrate the new branch.

### Source/calibration data

Use only independent cAMP-pathway data, primarily Catalan et al. PNAS 2015 (`10.1073/pnas.1415739112`), to identify the minimal conductance/gating law. Use WT/TMEM16A-KO IPR secretion, DCPIB/NPPB response, swelling/current timing and other directly relevant non-AE4-null observations. Digitise figures reproducibly if necessary and save calibration/error metadata.

Prefer the minimum law with one conductance scale and state dependence already present in the model (beta input and/or fractional swelling). Do not add a Hill coefficient, arbitrary delay, genotype multiplier or synergy unless demanded by independent source data. Do not create or compare a family of alternative VRAC laws. If the independent source cannot identify the single predeclared minimal law, report non-identifiability rather than branching into alternatives.

### Required prospective test

After the independent cAMP parameters are fixed, commit and push a PRE-AE4 prediction checkpoint. Only then compare with the frozen JBC 2015 AE4/AE2 constraints. No post-reveal retuning.

Publish the Stage-I conclusion even if negative.

## Stage II: measured resting-state inverse balance, in the SAME run

Regardless of whether Stage I succeeds, continue in the same Codex run after its result is safely pushed. Do not terminate merely because the VRAC hypothesis fails.

This stage must NOT guess a transporter first. Instead use the measured resting states as equations.

### Calibration data

Use the frozen Task 48 resting measurements:

- AE4 WT: Cl_i `50.10 +/- 1.50 mM`, pH_i `6.91 +/- 0.07`;
- AE4 KO: Cl_i `36.50 +/- 1.60 mM`, pH_i `6.89 +/- 0.02`;
- AE2 control: Cl_i `53.40 +/- 1.80 mM`, pH_i `6.87 +/- 0.01`;
- AE2 KO: Cl_i `54.50 +/- 1.80 mM`, pH_i `6.95 +/- 0.05`.

Use genotype-appropriate resting equilibria. Do not use stimulated secretion magnitude in this stage.

### Required calculation before any model edit

At each measured resting state, evaluate the inherited conserved-coordinate balance residuals for Na, K, Cl, TIC/alkalinity, charge/current and water/volume using the active source conventions. Quantify what steady source vector is missing for each genotype.

Then perform **one direct source-space calculation**, not a model-selection search:

1. assemble the matrix whose columns are the already represented pathway source signatures evaluated at the measured states;
2. compute its rank/SVD and project the measured residual vector onto its column space/cone;
3. if a continuous coefficient estimate is needed, solve one constrained least-squares problem using all predeclared columns simultaneously;
4. inspect the resulting identifiable direction and uncertainty.

Absolutely do **not** enumerate subsets of those columns. No all-pairs, all-triples, all-subsets, best-subset selection, L0 search, pathway on/off combinations, mechanism menus, candidate-family loops, genetic/evolutionary search, random search, or Cartesian parameter grids. Do not repeatedly refit after dropping/adding pathways. The purpose is geometric localization in the fixed source basis, not selection among combinations.

The question is:

> Does the measured residual lie in the span/cone of the fixed existing source basis, and if so which continuous source direction is identifiable; if not, what orthogonal residual source direction remains?

If the fixed basis is rank-deficient or the correction direction is not uniquely identifiable, publish that non-identifiability. Do not resolve it by discrete subset selection.

Explicitly include the already represented pathway signatures (NKCC1, NHE1, AE2, AE4, NBC, pump/K handling, CO2/carbon handling) at rest. If the residual lies in their span, identify the exact equation/gating assumption preventing the required shared flux relation. If it does not, report the residual source direction and the minimum additional conserved-coordinate signature implied by the projection residual. This is a vector-space statement, not permission to search through hypothetical transporters.

### Important prior-architecture warning

Task 36 set the stimulated NBC increment exactly to zero at rest solely to nest the prior WT state and explicitly stated that this was **not evidence that native basal NBC transport is absent**. Treat that as an architectural assumption to audit, not as the presumed solution. Do not simply switch on basal NBC unless the measured-state residual calculation supports its source direction.

Similarly, the source paper reports no detected genotype difference in stimulated NHE response; do not invent a genotype-specific NHE multiplier merely to close rest.

### Output

Publish:

- measured-state residual vectors;
- fixed source-signature matrix and rank/cone result;
- continuous projected/required correction direction with uncertainty;
- whether an existing shared pathway law can realise that direction;
- the exact next equation-level correction, only if uniquely identified;
- or a proof/diagnostic that the current pathway basis cannot span the measured resting phenotype.

Do not perform a large parameter search or any combinatorial subset search. If one uniquely indicated shared correction is identified and can be implemented/tested within this same run without contaminating Stage-I calibration, implement it and publish it as a separate checkpoint. Otherwise stop the model edit and preserve the mathematical localization for Task 50.

## Stage III: combine only independently fixed pieces

If Stage I yields a source-supported cAMP branch AND Stage II yields a uniquely identified shared resting correction, combine the two **without retuning either against AE4 saliva output**. Solve genotype-specific rests, run WT/AE4-KO/AE2-KO CCh, IPR and CCh+IPR protocols, freeze predictions, then compare with the complete Task 48 phenotype ledger.

If only one piece is identified, do not manufacture the other. Report exactly what remains.

## No automatic beta-NKCC modification

The original Task 49 prompt permitted adding beta-dependent NKCC regulation in parallel with VRAC. That permission is REVOKED for this run. A beta-NKCC change would confound the interpretation and is deferred unless a later measured-state/residual result uniquely requires it.

## Crash-safe publication

Only the orchestrator may modify the canonical branch or accept scientific corrections. Parallel workers are advisory.

Publish immediately after each dependency boundary:

- `49A`: repository exclusion verification + source/evidence freeze;
- `49B`: VRAC-only implementation and unit/conservation tests;
- `49C`: independent cAMP calibration and parameter freeze;
- `49D`: immutable pre-AE4 stimulated prediction;
- `49E`: held-out stimulated AE4/AE2 comparison;
- `49G`: measured resting-state residual/source-signature inverse analysis;
- `49H`: uniquely identified resting correction, if one exists and is implemented;
- `49I`: combined no-retune prediction, only if both independent pieces exist;
- `49F`: final article-facing decision report summarising all completed checkpoints.

Every checkpoint must update `CURRENT_STATUS.md`, include a machine-readable receipt, be committed and pushed, and have the remote SHA verified before dependent work begins.

## Definition of success for this run

Success does not mean forcing 35%.

This run is successful if it leaves the project with either:

1. a prospectively validated stimulated reconstruction plus a mathematically localized resting correction; or
2. a rigorous demonstration of which exact conserved-coordinate source/equation remains missing, sufficiently specific that Task 50 can implement one correction rather than search.

Do not spend this run producing another broad list of possibilities.