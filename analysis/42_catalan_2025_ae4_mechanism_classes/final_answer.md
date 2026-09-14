# Task 42 final interpretation

`CATALAN 2025 STOICHIOMETRIC SOURCE CLASSES ALONE DO NOT RECOVER THE AE4 LOSS PHENOTYPE UNDER THE INHERITED TASK 40 CYCLE LAW.`

Task 42 was frozen before phenotype reveal at commit `4cb75308bd81d06f064bb80b11831985fa4bb219`. That checkpoint contains all seven predeclared mechanism classes, seven WT resting solves and 21 production attempts. There were no numerical retries, parameter sweeps, optimisation calls, genotype resting solves or phenotype informed changes.

The held out comparison is now complete. None of the Catalán 2025 source supported stoichiometric classes recovers an approximately 35% AE4 loss secretion deficit when inserted into the inherited Task 40 scalar AE4 cycle law with the rest of the model fixed.

## Held out secretion comparison

| Class | Affinity direction at Task 40 state | WT 600 s status | 5% AE4 total deficit | AE4 null total deficit | Null early deficit | Null late deficit | Null NKCC1 compensation |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| C0 | opposed | pass | 3.42% | 3.86% | 1.37% | 4.16% | +23.16% |
| C1 | supported | pass | -10.12% | -10.76% | 0.26% | -12.32% | +95.48% |
| C2 | supported | WT pH failure | unavailable | unavailable | unavailable | unavailable | unavailable |
| C3a | supported | pass | -0.76% | -0.72% | 1.05% | -0.95% | +37.59% |
| C3b | opposed | pass | 6.86% | 7.67% | 1.68% | 8.38% | +11.30% |
| C4a | supported | WT pH failure | unavailable | unavailable | unavailable | unavailable | unavailable |
| C4b | opposed | WT pH failure | unavailable | unavailable | unavailable | unavailable | unavailable |

Early is 0 to 60 s and late is 60 to 600 s, exactly as frozen in the prediction checkpoint.

The largest valid AE4 null secretion loss is only 7.67% in C3b, and that class has a forward affinity opposed to the inherited Task 40 scalar cycle direction at the reference state. The source supported classes with valid 600 s WT trajectories do worse: C1 predicts 10.76% more secretion after AE4 deletion, while C3a predicts essentially no phenotype and a slight increase in secretion.

## The Catalán 2025 NaCl inward and K plus base outward proposal

C1 is the direct 1:1:1:1 NaCl inward and K plus HCO3 outward proposal. Its transport direction is thermodynamically supported at the Task 40 reference state and its WT trajectory remains physiological for 600 s. However, AE4 deletion increases rather than decreases total secretion. NKCC1 positive chloride loading rises by 95.48% in the null relative to the class specific WT, overwhelming the lost AE4 contribution.

C2 is the corresponding carbonate class, with NaCl inward and 2K plus CO3 outward. Its forward direction is also supported, but the stimulated WT trajectory crosses the inherited intracellular pH lower limit at about 190 s. It therefore has no admissible 600 s WT reference and cannot be used for a phenotype comparison under the predeclared rules.

Thus the specific NaCl inward and K plus base outward mechanisms proposed in Catalán 2025 do not recover the salivary phenotype in this first stage whole cell test.

## Acid base balance is a primary discriminator

All three carbonate source classes, C2, C4a and C4b, fail the stimulated WT pH gate. Their WT resting states are admissible, but their production trajectories cross the inherited lower pH bound at about 190 s, 222 s and 317 s respectively.

This is not a numerical failure. It is a biological consistency failure under the inherited acid base architecture. The carbonate classes remove one total inorganic carbon unit and two alkalinity equivalents per positive AE4 cycle, so they place a stronger demand on alkalinity replacement during stimulation.

This result reinforces the earlier NBC work. Once AE4 is required to carry substantial stimulated chloride transport, bicarbonate and alkalinity supply become part of the secretion mechanism rather than a secondary pH bookkeeping issue. NBC provides the additional inward base pathway needed to support productive AE4 transport while preserving WT pH. The present mechanism panel shows that changing AE4 stoichiometry alone does not remove that constraint.

## Temporal phenotype

The experimental phenotype is not only a scalar secretion deficit. Early secretion is relatively preserved and the sustained deficit develops later.

C0 and C3b reproduce that direction qualitatively. Their null deficits grow from 1.37% to 4.16% and from 1.68% to 8.38% between the early and late windows. However, the magnitude remains far too small.

C1 and C3a fail the temporal test more strongly. Their early null deficits are small and positive, but the late window becomes a secretion advantage after AE4 deletion.

So the Task 42 panel can produce the correct direction of temporal separation in some classes, but not the required magnitude.

## Comparison with Task 40 and Task 41

Task 40, the equal Na and K routing control, gives only a 3.86% null secretion deficit and +23.16% NKCC1 compensation. Task 42 shows that replacing that source vector with the Catalán 2025 permitted stoichiometric classes does not solve the problem. In the valid source supported classes, compensation is at least as problematic: +95.48% in C1 and +37.59% in C3a.

Task 41 answers a different question. It was deliberately designed against the desired secretion range and therefore is not an independent validation result. Its added AE4 dependent stimulated CaCC recruitment gives a 30.26% null cumulative secretion loss and strongly suppresses NKCC1 compensation to +2.79%. However, its temporal pattern is not the experimental one: the percentage deficit is already large early and declines towards the 600 s value. Task 41 therefore establishes that an additional network coupling can create the required magnitude without pathological chloride depletion, but it does not identify the biological mechanism.

Taken together, Tasks 40 to 42 separate three possibilities:

1. Equal AE4 cation routing alone is insufficient.
2. The source supported Catalán 2025 stoichiometric classes alone are also insufficient under the inherited Task 40 kinetic law.
3. A substantial phenotype can be generated when AE4 loss also changes another part of the secretory network, as demonstrated constructively by Task 41.

The remaining problem is therefore not simply the Na versus K split or the formal AE4 stoichiometry. It is the way AE4 transport is coupled to whole cell regulation during sustained secretion.

## Interpretation limit

This result does not falsify Catalán et al. 2025. That work constrains ion coordination and proposes stoichiometric classes, but it does not provide a salivary whole cell kinetic law, transporter recruitment law or complete regulatory model. Task 42 deliberately kept the Task 40 scalar AE4 cycle law fixed so that only transported source coefficients changed.

A future mechanism could therefore still combine a Catalán compatible stoichiometry with class specific kinetics or regulation and produce a different whole cell phenotype. That question remains open.

What Task 42 establishes is narrower and useful: **changing the AE4 transported source vector, by itself, is not enough to recover the knockout secretion phenotype.**

The frozen predictions are in `prediction_checkpoint.md` and `results/42_catalan_2025_ae4_mechanism_classes/`. The post reveal scalar comparison is in `results/42_catalan_2025_ae4_mechanism_classes/post_reveal_comparison.csv`.
