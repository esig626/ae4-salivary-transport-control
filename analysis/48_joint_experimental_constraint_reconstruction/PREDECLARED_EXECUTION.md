# Task 48 execution decisions before fitting

The orchestrator reviewed five independent advisory work streams at the common starting commit `5c945c4a7857657b5f7c76d5bf9e3d3abc7d5721`. All integration, scientific acceptance, inference and publication remain with the orchestrator. No worker branch is merged.

The selected architecture already has an AE4 beta driven activation state with gain 0.25 and time constant 30 s. R0 in Task 48 denotes this inherited architecture, including the older class named R1EffectiveActivation. Adding a second activation state is prohibited. Existing gain may be identified only if the isolated assay observation map is defensible.

## Protocol and observation decisions

Implement source recipe records, separate CCh and IPR arms, true transporter deletions and acute inhibitor masks. Keep an explicit unavailable prediction for each unrepresentable assay. Do not replace zero substrate by an arbitrary positive number. Carbonic anhydrase inhibition is not represented by instantaneous carbonate equilibrium. The selected NKCC law is conditioned on a fixed extracellular bath. These restrictions block a complete joint inverse problem unless a unique admissible correction can be established.

For supported positive bath diagnostics, use source Na 145, K 4.3, Cl 128.3 and HCO3 25 mM, pH 7.4 and temperature 310.15 K. Convert bicarbonate to TIC through inherited carbonate equilibrium. Add the documented Ca, Mg, glucose and HEPES particles as 17 mM external spectator osmoles. The net external spectator charge balances the explicit tracked ions; +4 mEq/L is supplied by the stated divalent salts, and the remaining carbonate/titrant difference is bookkeeping rather than a measured chemical species. HEPES is an externally clamped reservoir constituent, not a fitted cellular buffer. This is an effective representation of the positive B+ solution, not exact reconstruction of unreported pH titrant. No intracellular equation or transporter law changes. Publish every numerical result with this observation and bath qualification.

Resting states must be solved on the exact two charge constraints, using a shared parameter vector and the actual deletion. The saved WT rest is only a numerical seed. Check the full RHS and local stability. Never replace KO rest by WT rest.

SPQ is F0/F. Calibration values and the exact source slope windows are absent. Implement a conditional calibrated operator and report quantities that remain unavailable. Ratios across unequal chloride baselines do not automatically cancel calibration. The gland comparison uses integrated water and bicarbonate, not TIC. Reported SE describes uncertainty of a population mean, not biological SD.

## Inference gate and numerical budget

First test whether two supported source protocols are exactly indistinguishable in the knockout model. If a parameter independent contradiction survives a common observation operator, publish its analytic contrast and zero parameter sensitivity. Do not perform an expensive optimisation that cannot alter that contrast. This is an early infeasibility gate, not a claim to have fitted the remaining data.

If this gate fails the model, retain the complete residual ledger including unavailable entries; quantify resting states and supported stimulation diagnostics once, without fitting. The admissible parameter list remains a candidate list. No selected optimum or full numerical observation Jacobian may be claimed where protocols/operators are undefined. A sensitivity proof for the contradictory contrast is distinguished from a full Jacobian.

Use a single deterministic starting point for genotype resting roots and a small continuation in deletion only if the direct numerical root fails. This is a root solver, not parameter fitting. Numerical acceptance uses state scaled RHS and charge residuals, never biological percentage tolerances. A local stable root is not proof of global uniqueness.

R2 is accepted only if one uniquely supported equation repair represents the specific residual. Several possible biological explanations, or several missing assay dependencies, require a negative decision. Do not implement a list of alternatives. Publish 48D and 48E receipts even when the correction/inference is skipped or blocked, stating this explicitly. Final report must distinguish model rejection under a common operator from a theorem about every possible fluorescence calibration.
