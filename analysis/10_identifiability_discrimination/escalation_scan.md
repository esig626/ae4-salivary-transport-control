# Mandatory escalation scan

The scan was performed after deriving the strongest defensible result: the projected-signature rank criterion for a regular steady-state pump--leak model. Each route was tested against both the analytical structure and the reconstruction gate.

| Route | Test performed | Evidence found | Decision |
| --- | --- | --- | --- |
| A. Global identifiability / exact equivalence | Examined the affine diagnostic globally and the generic nonlinear map locally | Exact scalar equivalence lines and global injectivity of nonsingular two-output affine maps are immediate; nonlinear physiological map unavailable | No escalation. Global diagnostic facts are not physiological results |
| B. Perturbation-induced distinguishability | Derived the stacked-condition Jacobian and two-condition determinant | A perturbation is useful exactly when it rotates the sensitivity row; no closed calcium/bath model exists to nominate or validate a perturbation | No salivary theorem; retain as a standard conditional proposition |
| C. Optimal design | Compared all two-output panels by smallest singular value and condition number in the declared diagnostic | Ranking changes with arbitrary output scaling and the chosen restoring matrix | Do not optimize an uncalibrated surrogate |
| D. General pump--leak theorem | Factored the observation Jacobian into restoring, stoichiometric, activity, and measurement components | `rank[P s2,P s4]` is clean and general, but is an application of the IFT and rank criterion rather than a new pump--leak theorem | No escalation after novelty audit |
| E. Hypothesis testing | Asked whether AE2-only and AE4-only observation images are separated | Images cannot be validated without the full observation map and an error model | Do not introduce testing theory |
| F. Inference-relevant singularities | Checked the executable diagnostic Jacobian and reviewed requirements for the full model | The diagnostic restoring matrix is regular; the full `F_u` and branch are unavailable | No bifurcation/continuation work justified |

## Route B derivation

For two scalar conditions with sensitivity rows `j_1` and `j_2`, combined local identification occurs precisely when

`det([j_1; j_2]) != 0`.

If `j_2=lambda j_1`, the second condition adds replication but no structural direction. This result is useful for planning once a reproducible salivary model exists, but it is not new and does not select a perturbation by itself.

## Route C scaling check

Condition numbers depend on units and relative measurement-noise scaling. The diagnostic uses dimensionless outputs, yet its effective restoring matrix and `Q*` weights are declared rather than estimated. Reporting a “best” experiment from those numbers would elevate construction choices into evidence. The machine-readable panel table is retained only as a code/geometry audit.

## Escalation conclusion

No route supplies concrete evidence for a stronger global, perturbational, general-network, design, testing, or singularity result in the published salivary model. The reasons are substantive rather than computational: the central model is not independently executable, while the model-independent mathematics is standard.
