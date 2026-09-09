# AE4 ALONE INSUFFICIENT; MISSING PATHWAY LOCALIZED

## Controlling decision

**Changing AE4 alone on the best available internally coherent seven-state
chassis does not reproduce the held-out AE4-knockout phenotype.** No tested
mechanism simultaneously satisfies the source constraints, WT/resting gates,
negligible AE2-knockout requirement, and approximately `35%` AE4-knockout
secretion decrement.

The primary ex vivo experiment used `0.3 micromolar` carbachol plus
`5 micromolar` isoproterenol for 10 minutes and found total AE4-null saliva
`35 +/- 4.7%` lower (`n=6` per genotype), giving the held-out ratio
`0.65 +/- 0.047`. The first 2-3 minutes were comparable before a sustained
deficit developed.

This is a chassis-scoped exclusion, not a claim that AE4 is biologically
unimportant. The primary evidence requires transported Na and K and shows
beta/cAMP/PKA activation of AE4. Those features are necessary parts of an
adequate AE4 mechanism, but on this chassis they do not create the required
stimulus-dependent cation-current support. The first missing model property is
localized to **Na/K-pump-led cation homeostasis and its K-channel/acid-base
support**, not to another scalar AE4 capacity, Hill factor, saturation term, or
stoichiometric prefactor.

Three distinctions control the interpretation:

1. The exact frozen-row obstruction applies only to the historical WT row. It
   does not prove that mixed-cation AE4 has no other resting root.
2. KO-blind continuation did find rest-screen-passing `2:1:3` roots for C6,
   C7, and C8-SAT. Those roots close the moving-rest loophole within the
   declared finite source-default continuation, but every tested capacity in
   their rest-compatible bands predicts `AE4-KO/WT > 1`; positive PKA
   activation makes the disagreement worse. Unsupported axes and remote
   branches are not globally excluded.
3. The rest-compatible, source-shaped robustness grid therefore attains **no
   positive secretion reduction**. The historical C1 decrement (`0.144386%`)
   and the broader closure-only envelope (about `0.1476%`) are inadmissible
   controls, not universal maxima of AE4 biology.

## Direct answers to the twelve required questions

### 1. Can the experimental AE4-knockout phenotype be reproduced by changing AE4 alone?

**No.** The non-AE4 equations and parameters were held fixed, and the
knockout observation was withheld until AE4 capacity/rest calibration was
frozen. The investigation covered the full C1-C8 family, a 768-row static
interpolation, candidate-specific resting roots, rest-compatible capacity
bands, and beta/cAMP/PKA activation folds constrained by primary assays.

The failure is stronger than a fixed-row mismatch:

- C1 exactly closes the historical WT row but predicts only
  `AE4-KO/WT = 0.99855614` at the endpoint.
- C2 and C3 have mathematical moved roots, but their intracellular chloride
  lies outside the KO-blind, fixed-before-heldout two-SEM WT resting gate; C3
  also predicts the knockout change in the wrong direction.
- Exact-volume C6 roots for all three candidate stoichiometries fail the
  resting-Cl gate. That slice is not treated as exhaustive.
- KO-blind continuation subsequently finds simultaneous chloride, pH, and
  volume-sensitivity passes for source-shaped C6, C7, and C8-SAT `2:1:3`
  candidates. Their best tested no-activation endpoint remains above one, and
  their source-bounded PKA responses move farther above one.
- With AE4 set identically to zero, the deterministic multistart screen finds
  one candidate-independent local knockout root. It closes numerically but
  predicts `Cl_i = 48.782 mM` versus the measured
  `36.50 +/- 1.60 mM` (`+7.68` reported SEM), `pH_i = 7.3785` versus
  `6.89 +/- 0.02` (`+24.4` reported SEM), and cell volume about `2.974 pL`.
  No alternative AE4 law can alter that AE4-null vector field, although global
  uniqueness of the screened equilibrium is not claimed.
- On the C1 chassis, the candidate-specific AE2-null root predicts
  `Cl_i = 49.525 mM` versus the direct AE2-knockout target
  `54.50 +/- 1.80 mM` (`-2.764` reported SEM), although its pH and secretion
  change remain small. Candidate-specific AE2-null roots were therefore also
  solved for all 15 surviving `2:1:3` capacity points: every one fails the
  direct AE2-knockout chloride target (about `-4.74` to `-4.40` reported SEM),
  despite only about one-percent dynamic secretion effects. Thus near-normal
  AE2 secretion does not validate chronic resting homeostasis.

The AE4-null result localizes a fixed-chassis acid-base/volume/cation-homeostasis
failure independently of candidate kinetics. Absolute physical time and flow
remain uncertified, so simulated outputs are reported as within-protocol
ratios; an endpoint is not relabelled as an exact physical steady state.

### 2. If yes, what is the minimal structural AE4 property required?

**Not applicable: no AE4-only mechanism succeeds.** Two AE4 properties are
independently required for biological adequacy but are not sufficient here:

- transported K, as well as Na, is necessary for the physiological direction
  of passive reversible chloride loading at the primary-paper concentrations;
  and
- beta/cAMP/PKA regulation is experimentally established, including H89
  sensitivity and loss of activation for the S173A mutant.

K transport creates a distinct K-balance demand. A positive PKA gate changes
capacity but does not rotate a source vector, move a thermodynamic zero, or
add a missing cation-current pathway.

### 3. Which AE4 mechanisms fail, and why?

| Family | Result | Decisive reason |
| --- | --- | --- |
| **C1 historical Na-only `1:1:2`** | Exact frozen-row closure; endpoint `0.99855614`, post-step integral `1.00563816` | The `0.144386%` endpoint loss is about 240-fold smaller than the nominal `35%` target, while the integral changes in the wrong direction. Direct K transport and passive Na-only thermodynamics also exclude it as physiology. |
| **C2 printed 2018 pooled Na+K** | Fixed-row AE4-vector error about `0.0328154`; a moved root fails resting Cl | Its K source conflicts with the fitted Na-only row. Allocating pooled *net* reverse flux by intracellular Na/K fractions is not a species-resolved reversible branch law. Ratios from the drifting fixed-row simulation are not predictions. |
| **C3 explicit parallel Na/K** | Fixed-row AE4-vector error about `0.0331963`; moved `Cl_i about 46.31 mM`; endpoint ratio above one | Correct branch insertion exposes the missing cation-homeostasis direction. State relaxation does not rescue resting Cl or the secretion direction. |
| **C4 cooperative occupancy** | Nonnegative fixed-row calibration selects zero AE4 capacity | A scalar Hill response cannot rotate the balance signature; the empirical occupancy law is not independently affinity-sign certified. |
| **C5 transported cation plus allosteric gating** | No fixed-row rescue and no phenotype rescue | A positive gate changes magnitude only; it cannot remove the K residual or reverse an unfavorable affinity. |
| **C6/C7 reversible `1:1:2` and `1:2:3`** | The tested KO-blind rest calibrations fail the independent resting constraints | K coupling can repair local thermodynamic direction but does not by itself repair whole-network balance. The `1:2:3` law also changes the anion source signature. |
| **C6/C7 reversible `2:1:3`** | Rest-compatible continuation bands exist; all tested endpoint and integral ratios exceed one | This is the decisive moving-rest test. The source-compatible shape can close WT rest, but knockout increases rather than decreases predicted secretion; PKA activation worsens it. Candidate-specific AE2-null roots also fail the direct resting-Cl target. |
| **C8 mass-action/saturating pairs** | Mass-action rows are C6 aliases; rest-compatible C8 `2:1:3` saturation also predicts ratios above one | A positive saturation denominator preserves equilibrium, affinity sign, and the source ray. The extra unconstrained saturation parameter does not create a new balance pathway. Its PKA calculation is a Round-3 common-gate sensitivity adapter, not a separately measured candidate. |
| **Common PKA gates and recruited-K diagnostic** | Assay folds `1.25`, `1.6`, and `3.0` do not rescue the tested representatives | Common C1 gating leaves its endpoint nearly unchanged; recruited K and source-compatible mixed-cation gates move KO/WT farther above one. Selective recruited K is additionally novel and not thermodynamically certified. |
| **Pure cation-dependent gating without transported cation** | Rejected before simulation as a complete AE4 law | Primary experiments directly measured AE4-associated Na and K transport. |

The `1.25` recombinant and approximately `3` native AE4-only folds bound the
positive AE4-specific sensitivity. The approximately `1.6` fold is retained
only as a total-WT AE2-plus-AE4 sensitivity, not as an AE4-specific constraint;
fold `1` is the no-activation control.

The source-discussed carbonate and undefined bicarbonate self-exchange modes
are **untestable**, not failed: the seven-state chassis has no carbonate state,
and the primary experiments do not identify their missing stoichiometric
coefficients.

### 4. Is the published 2018 AE4 formulation sufficient, insufficient, or impossible to assess?

**Insufficient as printed on this coherent comparison chassis.** C2 is
assessable as the exact printed pooled rate and balance insertion. At the
closed historical row it leaves a raw WT residual of about `1.14e-3` and the
AE4-vector residual above. Its reverse cation partition is not a general
reversible branch convention, and its moved root fails the independent
resting-Cl gate.

The exact publication-generating implementation remains unavailable, so this
does not claim to reproduce every numerical choice used for the 2018 figures.
It establishes that the printed formulation, inserted coherently without
hidden non-AE4 compensation, is not a successful reconstruction.

### 5. Is the historical Na-only formulation sufficient?

**No.** C1 retains the negligible endpoint AE2 *secretion* effect
(`AE2-KO/WT = 0.99905227`) and is the only one-capacity candidate collinear
with the exact historical WT source requirement. That does not validate the
AE2-null resting state: its predicted chloride is `-2.764` reported SEM from
the direct AE2-knockout target. C1's AE4 endpoint decrement is only
`0.144386%`, and its post-stimulus flow integral increases by `0.563816%`
after knockout. It is also contradicted by direct K-transport evidence and has
the wrong passive thermodynamic direction at physiological Na concentrations.

### 6. What role do Na and K play separately?

Both are transported substrates, but their energetic and network roles differ.

- At the primary-paper physiological concentrations, the forward
  `Delta G/RT` values for Na-only `1:1:2`, `1:2:3`, and `2:1:3` cycles are
  `+1.8940`, `+4.4368`, and `+1.2453`: all oppose inward chloride loading.
- The strong outward K gradient can drive a K branch forward. In the
  equal-standard-constant parallel `1:1:2` calculation, net inward chloride
  requires a K capacity fraction above `0.3895`. This is a model-derived
  threshold, not a measured selectivity.
- On the full chassis, forward K-coupled AE4 creates a separate K-balance
  load. The exact steady cation/current identity ties apical chloride delivery
  to Na/K-pump and K-efflux support. K can therefore fix local transporter
  direction yet worsen secretion unless the supporting homeostasis changes.

The measured maximal Na and K assay responses (`1.5` and `1.6`) do not
identify physiological event fractions. Using their ratio as branch weights
is an explicit minimal assumption, never a knockout-fitted choice.

### 7. Is cation cooperativity required?

**It is measured, but its mechanistic placement is not identified and it is
not sufficient for the secretion phenotype.** The primary Na and K dose
responses have Hill coefficients about `2.0` and `1.8`. These slopes can arise
from a two-cation cycle, cooperative occupancy with one transported cation, or
transport plus allosteric regulation; the source authors leave that choice
open.

A positive Hill or allosteric factor changes rate magnitude but not the source
direction or thermodynamic zero. C4, C5, and the alternative integer
stoichiometries do not identify cooperativity as the missing property. The
needed discriminator is simultaneous AE4-specific Na/K and anion flux near
reversal, not another secretion-only knockout measurement.

### 8. Is an alternative transport stoichiometry required?

**No alternative stoichiometry is established as the required phenotype
mechanism.** The primary transport study considered, but did not measure, electroneutral
Cl:cation:HCO3 ratios `1:1:2`, `1:2:3`, and `2:1:3`. All three were tested as
branch-resolved reversible candidates.

The `2:1:3` ratio is numerically necessary for rest compatibility within the
tested reversible family: it is the only ratio for which the KO-blind
robustness continuation found simultaneous resting Cl, pH, and volume-
sensitivity passes. It is neither directly measured nor sufficient for the
phenotype. C6, C7, and C8-SAT variants all predict the held-out secretion
change in the wrong direction. The exact
electroneutrality relation `HCO3 = Cl + cation` rejects charged substitutions
but does not select a unique ratio, and a Hill slope near two does not prove a
two-cation cycle.

### 9. Is cation transport itself required, or can cation-dependent gating explain the phenotype?

**Cation transport is required by direct experiment; gating alone cannot
explain the phenotype.** Reversible Na and K changes were measured directly in
AE4-expressing cells, excluding a cation-gated `Cl/HCO3` exchanger with no
transported cation as the complete mechanism.

Additional regulation can coexist with transport. The 2021 primary study
supports beta/cAMP/PKA activation, H89 sensitivity, and an S173-associated
response. But a positive allosteric or PKA multiplier preserves stoichiometry
and affinity sign. Mapping its initial-rate assay folds to a sustained,
instantaneous whole-cell input is already a favorable new assumption; under
that assumption all rest-compatible `2:1:3` variants move farther from the
held-out phenotype as the fold increases. The common gates applied to generic
C7 and C8 roots are Round-3 sensitivity adapters; they add no independent
mechanistic evidence beyond the measured positive activation.

### 10. What maximum AE4-only effect is attainable, and which pathway is missing?

The effect must be reported by domain:

The moving-root domain fixes the source-default Na:K weights `1.5:1.6`,
`K_eq = 1`, saturation `sigma = 1`, and bath bicarbonate `21 mM`; it varies
only capacity and the positive PKA fold on a finite KO-blind grid. Its maximum
is local to those declared axes, not a bound over unsupported Na:K, energetic,
saturation, or bath grids.

| Tested domain | Best endpoint effect | Meaning |
| --- | ---: | --- |
| Original exact-volume source + rest + thermodynamic comparison | **No qualifying row; no numerical maximum** | The original exact-volume table contains no qualifying row, so no bound is inferred from that slice. This is not a claim that the complete source-admissible family is empty or that mixed-cation AE4 has no resting root; the post hoc KO-blind continuation below explicitly finds such roots. |
| Rest-compatible, source-shaped `2:1:3` robustness bands with tested PKA folds | **No positive reduction; smallest source-supported-fold endpoint ratio `1.01377776`** | Even the favorable fold-`1.25` row predicts a `1.3778%` increase. Positive PKA activation worsens it. The capacity scan is KO-blind but post hoc WT-rest robustness, not a primary capacity identification or global proof. |
| Exact-closure historical C1 control | `0.144386%` reduction (`KO/WT = 0.99855614`) | Reproducible but Na-only and assay-inadmissible; its integral has no decrement. |
| Deliberately generous static Round-2 closure-only envelope | about `0.1476%` reduction (`KO/WT = 0.99852364`) | Pure Na, `K_eq = 100`, and unconstrained saturation; an unsupported diagnostic, not a universal family bound. Its integral also has no decrement. |

Thus the **maximum positive endpoint reduction over the rest-compatible,
source-shaped tested grid is `0%`**. Among the AE4-specific source-supported
folds (`1.25` and `3`), the smallest endpoint ratio is `1.01377776` and the
corresponding integral ratio is `1.02093900`, already a secretion increase.
The fold-`1` reference is slightly more favorable but remains above one and is
not evidence against measured activation. Across the five frozen capacity
points per core at the source-supported fold `3`, the production ranges are:

| Rest-compatible `2:1:3` law | KO-blind capacity interval | Endpoint KO/WT range | Post-step integral KO/WT range |
| --- | ---: | ---: | ---: |
| C6 reversible mass action | `[3.5105005e-12, 5.9165257e-12]` | `1.020180-1.026207` | `1.027228-1.033437` |
| C7 reversible affinity law | `[0.0033498713, 0.0072388250]` | `1.031710-1.054591` | `1.038475-1.061737` |
| C8 symmetric saturation | `[0.0344779857, 0.0378955063]` | `1.070946-1.079120` | `1.077435-1.085671` |

The fixed steady water map is strictly monotone in apical chloride delivery.
Within that map, a flow ratio of `0.65` requires an apical chloride-flux ratio
of `0.641489`, or a `35.851%` loss. The exact steady cation/current identity
places that missing support in the combination of Na/K-pump current and K
efflux, with acid-base/Na loading controlling the compatible state.

At the frozen historical row, replacing the Na-only AE4 source by the printed
mixed-cation source gives a pure Na-to-K correction. The Na/K-ATPase signature
has the highest one-component alignment, `|cos(theta)| = 0.98058`, but its
fixed `-3:2` Na:K ratio leaves an irreducible `19.612%` residual. NHE1 alone
cannot cancel any K residual. The smallest exact positive decomposition uses
equal code-native cycle-flux increments of the pump and an independent
Na-loading direction such as NHE1.

Round 4 therefore ends with an analytical one-component rejection and module
localization, not a phenotype-tuned simulation. Pump-only and NHE1-only scalar
corrections cannot provide exact frozen-row closure. A quantitative
two-direction test would require an independently measured flux conversion,
physical-time map, and stimulus law. The available normalized NHE fluorescence
slope provides none of these and is unchanged by genotype. The result
localizes **stimulus-dependent cation-current support, pump-led at the frozen
row and coupled to K-channel/acid-base homeostasis**; it does not assert that a
specific protein is altered in knockout or prove that a joint correction is
sufficient.

### 11. Which conclusions are exact/analytical, numerical, or assumption-dependent?

| Status | Conclusions |
| --- | --- |
| **Exact/analytical under the declared fixed-chassis equations** | Electroneutrality requires `c = a + b`; a positive gate or saturation denominator cannot change an affinity sign or equilibrium; exact preservation of the historical WT row requires the AE4 source to equal its Na-only `1:1:2` vector; with AE4 activity zero, the AE4-KO vector field is candidate-independent; the steady water map is monotone in apical chloride delivery; the steady cation/current identity places chloride support in pump plus K efflux; pump-only and NHE1-only source directions cannot each close the frozen mixed-cation residual, whereas pump plus an independent Na-loading direction spans it exactly. |
| **Numerical, reproducible, and chassis-scoped** | Capacity projections and closure residuals; C1/AE2 endpoint and integral ratios; the 768-row envelope; C2/C3 and C6 moved roots; the screened AE4-null local branch of the candidate-independent vector field and the candidate-specific AE2-null resting roots; rest-pass `2:1:3` capacity intervals; and PKA-gated trajectory ratios. Production Radau covers all 75 capacity/fold rows; independent BDF covers all 30 positive AE4-specific rows, all remain above one, and the maximum absolute BDF/Radau differences are `1.36e-9` endpoint and `3.37e-9` integral. Conservation residual is at most `4.98e-14`, and thermodynamic sign violation is zero. A separate 15-row set covers one in-band point per core across all five folds. The `0.3895` K-fraction threshold and `35.851%` chloride-flux threshold are numerical values after inserting primary concentrations and fixed water parameters. Root continuation establishes the reported branches and finite bands, not absence of every remote mathematical branch. |
| **Primary measurements, not model deductions** | Knockout saliva, resting chloride, and pH phenotypes; the negligible AE2-knockout secretion effect; direct Na and K transport; voltage independence/electroneutral behavior; EC50/Hill fits; reversal; and beta/cAMP/PKA/H89/S173 regulation. |
| **New modelling assumptions, kept explicit** | The historical seven-state comparative chassis; fixed non-AE4 parameters; pure Na and pure K branches; assay-maxima branch weights; passive `K_eq = 1`; `21 mM` bath bicarbonate; empirical occupancy/allostery and `sigma = 1` saturation closures; mapping initial-rate PKA folds to an instantaneous sustained input; using the candidate's WT root for every genotype and a 100-code-unit prestimulus segment; exact-volume calibration; the post hoc but KO-blind finite `1e-8`-to-`1e4` capacity continuation with `0.025`-log-decade local refinement and `+/-10%` volume sensitivity (a sensitivity convention, not reported experimental uncertainty); the `+/-2` reported-SEM WT/rest acceptance gates; the one-SEM held-out flow interval `[0.603, 0.697]` with both endpoint and integral required; the `+/-5%` AE2 dynamic-normality threshold used by the automated hard-gate screen (a new qualitative sensitivity, not an experimental confidence interval); endpoint/integral reporting on a time-uncertified coordinate; and pathway-signature projection as localization rather than parameter fitting. |

The frozen-row theorem is not generalized to moved resting states. Conversely,
the existence of rest-compatible roots is not promoted to mechanism success:
their held-out dynamics and the screened AE4-null branch of the
candidate-independent vector field fail. The conclusion is therefore
numerical across the declared candidate/rest bands, supported by exact pathway
constraints, and not a claim
of global exclusion over every imaginable AE4 law.

### 12. Does the result justify reopening identifiability or hypothesis testing?

**No.** A validated physiological parameter-to-observation map and a mechanism
that survives held-out phenotypes were prerequisites for reopening that
project. Task 12 instead finds a candidate-independent AE4-null resting
vector field whose screened local branch fails direct phenotypes, an
uncertified physical time/flux scale, and missing
stimulus-dependent cation/acid-base support. Physiological identifiability or
mechanism-discrimination calculations on this map would assign precision to a
falsified chassis.

The next admissible step is to reconstruct and validate the Na/K-pump-led
cation-homeostasis, K-efflux, and acid-base coupling against independent WT
and knockout Na, K, Cl, pH, volume, and time-resolved stimulated data. The
secretion knockout must remain held out. Only after that chassis closes,
reproduces the independent phenotypes, and survives audit should
identifiability, discriminating experiments, or hypothesis testing be
reopened. The conditional Task-10 algebra remains valid under its assumptions;
it is not promoted to physiological inference here.

## Evidence and reproducibility boundary

The experimental constraints are recorded in `experimental_evidence.md` from
the primary knockout study
([DOI 10.1074/jbc.M114.612895](https://doi.org/10.1074/jbc.M114.612895)), the
direct transport study
([DOI 10.1085/jgp.201611571](https://doi.org/10.1085/jgp.201611571)), and the
PKA study
([DOI 10.1152/ajpgi.00145.2021](https://doi.org/10.1152/ajpgi.00145.2021)).
Candidate definitions, numerical comparisons, the independent audit, and the
structural derivations are in `ae4_candidate_family.md`,
`model_comparison.md`, `adversarial_audit.md`, and
`structural_explanation.md`. Machine-readable values are under
`results/12_ae4_mechanism_reconstruction/`.

> **Final classification: `AE4 ALONE INSUFFICIENT; MISSING PATHWAY LOCALIZED`.**

The exclusion is exhaustive over the declared C1-C8 implementations and the
tested WT-rest-compatible capacity/PKA grids on this fixed chassis. The
localized defect is a cation-current-support module, not a uniquely identified
single protein correction. No manuscript or physiological identifiability
claim follows from this result.
