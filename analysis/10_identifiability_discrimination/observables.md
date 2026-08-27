# Observation panels and diagnostic geometry

## Measurement definitions

| Observable | Published/experimental status | Role in the inverse problem | Phase 10 status |
| --- | --- | --- | --- |
| `Q` fluid secretion rate | Primary published model output and experimentally accessible gland output | Scalar headline measurement | The general scalar nonidentifiability result applies; no reproduced numerical `Q` map |
| `Na_i` | Model state; measurable with ion-sensitive methods | AE4 has a direct cation signature absent from AE2 | Candidate second measurement; full-model rank unverified |
| `K_i` | Model state; measurable with ion-sensitive methods | Same direct cation distinction, with a different AE4 partition | Candidate second measurement; full-model rank unverified |
| `Cl_i` | Model state with direct relevance to secretion | Both exchangers import chloride | Candidate but may align the two effects; full-model rank unverified |
| `HCO3_i` or `pH_i` plus closure | Model acid-base state | AE2/AE4 differ by one versus two bicarbonates per cycle | Candidate; closure and buffering propagation unresolved |
| Cell volume or a lumen/osmotic readout | Model state/output | Adds a water/osmotic projection | Candidate only after the Eq. (9) orientation conflict is resolved |

## Minimum cardinality

Two independent activity parameters require at least two scalar steady-state observables locally. This is a necessary cardinality statement, not a guarantee that any arbitrary pair works. A panel is sufficient only when its `m x 2` observation Jacobian has rank two.

Biologically plausible pairs to test once the numerical model is closed are:

- `Q + Na_i` or `Q + K_i`, because the cation component is a direct stoichiometric distinction;
- `Q + Cl_i`, because it combines secretion and the intracellular chloride pool;
- `Cl_i + HCO3_i`, because the exchanger cycle ratios differ; and
- `Q + volume` only after water-flux orientation and observable definitions are fixed.

The projected-signature criterion in `theory.md` is the required test for each pair. No pair is called physiologically identifying in this audit.

## Executable stoichiometric diagnostic

The diagnostic defines normalized outputs

`Q* = 1 + (0.10,0.10,0.70,0.10)(x-1)`

and `volume*` as the equal-weight mean state ratio. `Q*` is a chloride-weighted plotting proxy, **not** a reconstructed secretion rate. The weights and effective restoring matrix are declared so that no fitted or historical information is smuggled into the calculation.

On the `41 x 41` grid `G2,G4 in [0,2]`:

- the scalar `Q*` map has an exact line through `(1,1)` on which `Q*=1`;
- all 1681 solved diagnostic states remain positive;
- all 15 two-observable panels happen to have rank two under the declared restoring matrix; and
- scaled condition numbers range from `1.62` (`Q*+K_i`) to `63.4` (`K_i+volume*`).

This broad full-rank outcome is itself a warning: an arbitrary regular restoring matrix can make almost every panel look sufficient. It cannot substitute for reconstruction of the physiological `F_u`.

## Files

- `results/10_identifiability_discrimination/panel_metrics.csv`: rank, singular values, determinant, and condition number for every one- and two-output panel.
- `results/10_identifiability_discrimination/parameter_observations.csv`: complete parameter grid and normalized/raw diagnostic outputs.
- `results/10_identifiability_discrimination/q_equivalence_set.csv`: exact scalar-equivalence line segment.
- `figures/10_identifiability_discrimination/observation_geometry.png`: scalar equivalence and the `Q*`/`Na_i` image.
- `figures/10_identifiability_discrimination/panel_rank_diagnostics.png`: smallest singular values for pairs.

## Mechanism discrimination

Competing “AE2-only” and “AE4-only” mechanism classes would generate curves or surfaces in observation space. Discrimination requires their images to be disjoint at the stated error resolution; local rank within either class is not enough. Since the full salivary observation map is not reproducible, the image intersection/separation question is not numerically answerable here. No classifier, hypothesis test, or biological mechanism-separation claim is reported.
