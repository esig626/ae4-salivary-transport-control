# Viability audit

## Classification: STOP

The required stop conditions are met.

## Evidence

| Gate | Required standard | Result |
| --- | --- | --- |
| Independent published-model reconstruction | Closed numerical steady-state model with a baseline check | Failed: critical scaling, sign, forcing, closure, and solver details remain unresolved |
| Analytical content | At least one clean result | Passed at a general level: IFT derivative, scalar equivalence manifold, projected-signature rank criterion, and stacked-condition rank test |
| Salivary minimal panel | Rank verified in the reproduced physiological observation map | Failed: only an explicitly uncalibrated diagnostic can be evaluated |
| Mechanism discrimination | Validated mechanism images shown to intersect or separate | Not justified |
| Escalation | Concrete global, perturbational, network, design, testing, or singularity evidence | None |
| Novelty | Exact main claim not already standard | Failed: the strongest results are standard IFT/rank/experimental-design consequences |

## Why the diagnostic does not rescue the project

The diagnostic is reproducible, its derivative matches finite differences to `5.43e-11`, and it displays the expected one-dimensional scalar equivalence set. It also makes every tested two-output panel full rank. Those facts validate the implementation and illustrate the mathematics, but they depend on a declared restoring matrix and proxy-output weights. Treating them as physiological evidence would be circular.

## Deliverables completed despite the stop

- independent structural reduction of AE2/AE4 contributions;
- analytical proofs and scope limits;
- candidate measurement-panel framework;
- machine-readable rank and observation geometry for a labeled diagnostic;
- automated derivative, stoichiometry, positivity, rank, and sign-invariance tests;
- mandatory escalation scan;
- targeted novelty audit; and
- checked results/decision records.

## Consequence

Do not begin a manuscript from this phase and do not claim numerical identifiability or discrimination for the 2018 salivary model. A restart is justified only if an authoritative, citable specification resolves the numerical blockers (especially AE4 scaling, current/flux conversions, calcium forcing, water orientation, CO2 sign, and luminal closure). At that point the same code interfaces and rank criteria can be applied to the reconstructed `F_u` and actual observables.
