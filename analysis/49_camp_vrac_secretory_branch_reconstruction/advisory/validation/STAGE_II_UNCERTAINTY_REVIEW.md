# Stage II units and uncertainty review

Read-only advisory review of Task 49 `resting_geometry.py`, `resting_geometry.json`, `resting_uncertainty.json` and `source_svd_projection.json`. No model evaluation, finite difference, inference, SVD or trajectory was executed by this reviewer.

**Recommendation:** publish the conditional diagnosis and non-identification result, with the reporting clarifications below. No substantive formula/sign error was found by code and saved-output inspection. This is not acceptance of a unique correction.

The measured inputs retain their reported SE directly: chloride in mM and pH in pH units, without converting pH SE to proton-concentration SE or rescaling by sample size. The centred-difference columns therefore mean residual-rate change per mM or per pH unit, respectively. Each perturbation recomputes the declared algebraic lift, so these are derivatives along that conditional state surface; they are not derivatives holding all hidden physical state coordinates fixed.

The working covariance `J diag(SE²) Jᵀ` is correctly formed. The rowwise `sum(abs(J)*SE)` is the maximal linearised standard deviation over unknown correlations of the relevant input estimates. It is not a covariance matrix, confidence interval, nonlinear bound or proof that the cone certificate survives measurement/latent uncertainty. Cross-cohort independence is an explicit working assumption for the stored covariance. Nuisance uncertainties are unreported and are correctly left unassigned rather than silently set to measured error scales.

The lift imposes both charge constraints, `C=H` and zero cell-volume derivative, while retaining six cached nuisance coordinates per cohort. It preserves measured Cl/pH but is not an experimentally determined complete state. In particular, imposing `C=H` makes TIC and TA residuals agree; that agreement is a constructed necessary stationary identity, not an additional independent validation. The local nuisance Jacobian does not eliminate the unmeasured-state manifold.

The reported rank 48/nullity 36, approximately zero orthogonal residual and exact source identities do not identify pathway coefficients. Since the inherited residual itself is formed from the same sources, span membership follows from `F=Sj`. The block coefficients are a relaxed state-specific source-space representation, not independently fitted genotype multipliers or a solved common kinetic law. The existing warnings correctly say this.

The cone certificate has a clear limited meaning: at the declared AE4-KO state, native NHE1 and reverse AE2 both add intracellular alkalinity, while the required correction removes it. Nonnegative additions in those native directions cannot supply that correction. This local result, together with `H=E=C` and its scalar compatibility equation, identifies the obstructed balance at that lift. It does not prove the same affinity signs throughout the hidden-state manifold, select a shared equation correction, or justify basal NBC activation. Skipping 49H/49I is therefore appropriate under the correction gate.

Before publication, add these metadata/reporting clarifications without rerunning numerical work:

1. Specify units per raw residual row: Na/K/Cl/TIC in fmol/s; TA and charge/current audit rows in fmol charge-equivalents/s; water rows in pL/s. Current audit rows are not amperes; conversion is multiplication by Faraday's constant times `1e-15`. Jacobian entries use the corresponding row unit divided by their input-column unit; covariance entries use the product of their two row units.
2. Describe `min_norm_coefficient_representative` as minimum norm in the column-normalised coordinate `z = column_scale * delta`, not minimum unweighted norm of mixed-unit raw physical coefficients. It has no physiological preference or uniqueness implication.
3. Identify roundoff-level charge/current/imposed-water derivative and SE entries as numerical residue of exact/imposed identities. They are not estimates of biological measurement variation in those coordinates. Keep the raw audit output, but do not standardise those residuals into significance scores.

The openly recorded replay after the diagnostic-key export failure preserves the distinction between one distinct fixed basis and two actual SVD executions. Do not simplify that provenance to a claim that only one SVD call ever occurred. No further replay is needed for the reporting changes above.
