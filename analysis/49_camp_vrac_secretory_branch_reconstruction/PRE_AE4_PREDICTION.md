# 49D immutable pre-AE4 prediction record

49C did not identify a source-calibrated conductance, so no complete numerical VRAC reconstruction is claimed. `output/pre_AE4_prediction.json` preserves every requested cohort/protocol with its available or unavailable status. No arbitrary parameter, new trajectory or held-out comparison is used to fill missing entries.

REST and CCh-only exactly nest the already published Task48 predictions for every conductance because beta is zero. Their artefacts are referenced by hash and are not recomputed. This does not validate their physiological resting values.

There is one stronger, parameter-independent prediction. Let `x*` be the exact AE4-KO resting core and set the gate reference to its own resting volume. In the inherited IPR-only equations, all beta action on the core has disappeared with AE4. At `x*`, the added positive-swelling gate is also zero. Therefore `x(t)=x*`, accompanied by the changing but disconnected AE4 regulatory coordinate, is a solution of the augmented ODE for every finite `g>=0`. Within the positive domain the constitutive functions, unique electrical closure and rectifier are locally Lipschitz, so the solution is unique.

Consequently AE4-KO IPR-only predicts zero new VRAC current, no chemical/volume response, zero net chloride concentration slope and zero SPQ uptake slope for any fixed optical calibration. Resting outflow may be nonzero but has no IPR-induced increment. Numerical errors in a cached equilibrium cannot be used as a physical swelling seed.

This theorem is frozen before a new AE4 comparison. It does not depend on extending the sampled WT no-swelling observation into an unproved continuous bound, and it does not require a conductance fit. WT/AE2 beta-containing quantitative predictions and the AE4 combined-stimulus secretion deficit remain unavailable. They are not replaced by old Task48 results and labelled new VRAC outcomes.
