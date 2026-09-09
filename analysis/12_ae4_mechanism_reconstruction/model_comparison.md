# AE4 model comparison

## Result

**No tested AE4-only reconstruction satisfies the complete validation
contract on the fixed seven-state chassis.** The historical Na-only law closes
the historical WT row but is contradicted by direct K-transport evidence and
predicts only a `0.144386%` endpoint secretion decrement after AE4 deletion.
Source-shaped, thermodynamically valid `2:1:3` Na/K mechanisms can instead be
placed on KO-blind WT-rest-compatible branches under an explicit post-hoc
volume sensitivity. Across every frozen capacity point and the independently
supported positive PKA folds, however, they predict `AE4-KO/WT > 1`: deleting
AE4 increases, rather than decreases, secretion. Their candidate-specific
AE2-null resting roots also miss the directly measured AE2-null chloride
phenotype.

This establishes the Task-12 conclusion **`AE4 ALONE INSUFFICIENT; MISSING
PATHWAY LOCALIZED`** within the tested chassis and finite reconstruction
domains. It is not a universal exclusion of remote equilibria, untested
Na:K selectivities, or other non-AE4 chassises.

## 1. Comparison contract

### 1.1 Anti-leak calibration

The experimental AE4-knockout saliva ratio, `0.65 +/- 0.047`, was never an
optimizer input, capacity selector, equilibrium selector, or scan-bound
selector. The computation proceeded in this order:

1. freeze the non-AE4 chassis and its fingerprint;
2. use the historical WT source demand, independent AE4 transport assays, and
   source-supported resting observations to determine AE4-only quantities;
3. freeze the mechanism, capacity, resting-root branch, and stimulus fold;
4. simulate WT, AE4-null, and AE2-null trajectories;
5. join the held-out phenotype only to score predictions.

The post-hoc capacity continuation was also KO-blind. Its axes and selected
log-capacity quantiles were frozen before any genotype simulation. Machine
records expose `knockout_target_read_by_variant_selector = false` and the
calibration targets actually used by each family; they do not generically
claim that Hill or EC50 data calibrated a law that does not contain those
quantities.

### 1.2 Fixed and varied quantities

The primary comparison retains every non-AE4 equation and parameter in
`chassis.md`. Scenario changes are restricted to AE4 or AE2 genotype
multipliers and, in the explicit stimulus round, a separate PKA input. A
stored AE4 capacity is multiplied by the genotype factor, so WT preserves the
calibrated capacity and AE4 knockout gives exactly zero contribution. The
AE4-null vector field is therefore candidate-independent. Its trajectory is
identical across candidates only when the initial condition and input history
are also identical; moving-rest candidates begin from their own frozen WT
roots before genotype-specific evolution.

Candidate complexity counts distinguish fitted/free values from independent
assay constants and new assumptions. In particular, C3 has one fitted
capacity when its `1.5:1.6` Na:K weights are fixed from assay response maxima;
the C8 saturation strength remains an unconstrained modelling degree even
when fixed at one. C8 mass-action labels are exact code-level aliases of their
C6 counterparts and are retained for registry completeness, not counted as
independent confirmations.

### 1.3 Observables and gates

The principal outputs are the secretion proxy at the end of the common
protocol and its post-step integral. Both are dimensionless within-protocol
ratios because the historical model's absolute flow and time units are not
certified. The endpoint is a fixed protocol coordinate, not automatically a
certified steady state; endpoint raw-RHS norms are retained.

The cumulative 10-minute experimental phenotype cannot be mapped uniquely to
one code-time observable. Accordingly, a successful flow screen would require
both endpoint and integrated ratios to lie within the experimental interval.
Even that would only be a flow screen: success additionally requires WT rest,
AE2 knockout, genotype-specific Cl and pH, conservation, thermodynamics,
complexity review, and an external adversarial audit.

Resting WT constraints use the primary measurements `Cl_i = 50.10 +/- 1.50
mM` and `pH_i = 6.91 +/- 0.07`; the numerical rest screen uses a declared
two-reported-SEM window. A `1.3 pL` reference volume and the later `+/-10%`
volume band are modelling references/sensitivities, not measured WT
uncertainty. Direct genotype targets are kept separate:

| Genotype | Direct chloride target | Direct pH target |
| --- | ---: | ---: |
| AE4 knockout | `36.50 +/- 1.60 mM` | `6.89 +/- 0.02` |
| AE2 knockout | `54.50 +/- 1.80 mM` | `6.95 +/- 0.05` |

The AE2 study's own-control chloride and pH values are exported as additional
comparators, not substituted for the direct knockout line.

## 2. Round 1: source and lineage candidates

At the exact historical WT state the chassis requires the AE4 intracellular
source vector

`(Na, K, Cl, HCO3) = (-d, 0, +d, -2d)`, with
`d = 0.0388833386` in code-native amount/time units.

This is a chassis-compatibility constraint, not a biochemical endorsement of
Na-only AE4. Any source-supported positive K transport creates a K-balance
component that a single fitted AE4 capacity cannot cancel at the same state.

| Family | Fixed-row result | Predictive interpretation |
| --- | --- | --- |
| C1 historical Na-only | Exact WT source closure | Endpoint `AE4-KO/WT = 0.99855614`; integrated ratio `1.00563816`. It is a lineage negative control because direct K transport and passive thermodynamics reject it. |
| C2 printed pooled Na+K | AE4-vector max error about `0.0328154` | The printed mixed-cation insertion is incompatible with the Na-only-calibrated row. Its reverse pooled-flux allocation is not a species-resolved reversible branch law. |
| C3 explicit Na/K branches | AE4-vector max error about `0.0331963` | Correct branch balances expose the missing K-homeostasis direction. |
| C4 cooperative occupancy | Nonnegative capacity projection reaches zero | Hill occupancy changes magnitude, not the required balance direction; Hill slopes are not stoichiometry. |
| C5 transported cation plus gate | No fixed-row closure | A positive gate cannot rotate the transported source vector. |
| C6 alternative stoichiometries | No source-supported fixed-row closure | Mixed-cation balance rays conflict with the exact historical row. |
| C7 reversible affinity laws | No source-supported fixed-row closure | Thermodynamic sign control does not eliminate the cation-balance mismatch. |
| C8 mass action/saturation | MA rows alias C6; SAT adds one assumption | Saturation preserves the source ray and thermodynamic zero and cannot repair fixed-row closure. |

The compact primary ledger below reports the requested prediction and
residual fields. `T1` is the fixed-chassis WT AE4 balance vector, `T2` is the
primary Na:K maximal-response capacity ratio, and `T3` is the primary cation
dose-response shape. `fit/free` separates actually fitted quantities from all
counted AE4 degrees. `Delta` values are AE4-null minus WT at the protocol
endpoint. Solver status is `pass` for every row, but that means integration
completed; it does not override the WT-closure gate.

| Mechanism ID(s) | fit/free | Targets | Calibration max error | WT raw closure / gate | AE4-null endpoint ratio | AE2-null endpoint ratio | `Delta Cl_i` (mM) | `Delta pH_i` | Max conservation | Thermodynamics |
| --- | ---: | --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | --- |
| C1 | `1/1` | T1 | `0` | `1.78e-12` / pass | `0.99855614` | `0.99905227` | `-1.06145` | `+0.36679` | `4.97e-14` | empirical, not sign controlled |
| C2 | `1/1` | T1 | `0.0328154` | `1.143e-3` / fail | `1.0290560` | `0.9847905` | `-0.07751` | `+0.31901` | `4.97e-14` | empirical, not sign controlled |
| C3 | `1/1` | T1+T2 | `0.0331963` | `1.157e-3` / fail | `1.0295240` | `0.9845608` | `-0.06436` | `+0.31837` | `4.97e-14` | empirical, not sign controlled |
| C4 | `1/1` | T1+T2+T3 | `0.0777667` | `2.710e-3` / fail | `1.0000000` | `1.0016481` | `0` | `0` | `5.33e-14` | empirical, not sign controlled |
| C5 | `1/1` | T1+T2+T3 | `0.0777667` | `2.710e-3` / fail | `1.0000000` | `1.0016481` | `0` | `0` | `5.33e-14` | pass |
| C6/C7/C8-MA `112`, `123`, `213` | `1/1` | T1+T2 | `0.0777667` | `2.710e-3` / fail | `1.0000000` | `1.0016481` | `0` | `0` | `5.33e-14` | pass |
| C8-SAT `112`, `123`, `213` | `1/2` | T1+T2 | `0.0777667` | `2.710e-3` / fail | `1.0000000` | `1.0016481` | `0` | `0` | `5.33e-14` | pass |

The zero changes in C4-C8 default rows arise because the nonnegative
fixed-row projection selects zero AE4 capacity. They are not independent
evidence that those mechanisms predict no phenotype after an admissible
moving-rest calibration; the later continuation supplies that test.

The only fixed-row closure is therefore C1, which is assay-inadmissible. This
does not imply that the source-supported mixed-cation candidates lack other
closed resting states; that loophole is tested below.

Two source-discussed modes remain explicitly **untestable**, rather than
failed: chloride/carbonate exchange cannot be represented without a carbonate
state, and the chloride-free bicarbonate self-exchange observations do not
identify the missing stoichiometric coefficients. Neither was silently
assigned arbitrary parameters.

## 3. Round 2: minimal structural interpolation

Round 2 varied only declared AE4 axes: Na versus K fraction, reversible
equilibrium-bias diagnostics, saturation strength, and the explicitly labelled
AE4 bicarbonate-driving-force conventions. It did not tune a non-AE4
parameter and did not read the knockout target.

The fixed-row result is sharp. If the maximum per-species closure error must
remain below `1e-8`, the transported-K fraction must satisfy
`f_K <= 2.5718e-7`, effectively the pure-Na boundary. That contradicts direct
AE4-associated K flux. Ten grid rows close, all on that effectively pure-Na
boundary; they are therefore diagnostic, not source-admissible.

The most favorable endpoint on the deliberately generous interpolation grid
is the closure-only row
`R2_0462_C8_112_SAT_fNa1_Keq100_sat100_B21`, with
`AE4-KO/WT = 0.9985236403`, a `0.1476%` endpoint decrement. It is pure Na, uses
an unsupported `K_eq = 100` energetic bias and unconstrained saturation, and
its integrated response does not show the experimental decrement. Passive
aqueous exchange has `K_eq = 1` as the reference; decade shifts are retained
only as missing-energetic-bias localization. Similarly, injecting `40` or
`42.9 mM` bicarbonate into AE4 while the remaining acid/base baths remain at
the historical `21 mM` convention is a driving-force sensitivity, not a
coherent whole-bath success.

Thus `0.1476%` is an unsupported static diagnostic envelope, not a
source-supported maximum and not a universal AE4 bound.

## 4. Moving-rest reconstruction

### 4.1 Root method

To avoid treating the historical coordinates as compulsory, the moving-rest
screen eliminates cell height through electroneutrality and solves seven
steady residuals in six concentrations plus `log(H+)`. It uses bounded TRF
least squares, scaled raw residuals, analytical height recovery, and no
knockout target. The recorded multistart audits use fixed seed `20260827` and
find one local closed basin per screened candidate from the declared starts;
this is robustness evidence, not a proof of global uniqueness.

C2 and C3 have exact positive roots. Their approximate resting values are:

| Candidate | `Cl_i` (mM) | pH | Volume (pL) | Result |
| --- | ---: | ---: | ---: | --- |
| C2 pooled | `46.35` | `6.908` | `1.248` | WT chloride outside two-SEM gate |
| C3 parallel | `46.31` | `6.908` | `1.248` | WT chloride outside two-SEM gate; knockout flow direction wrong |

The remaining active C4-C8 families were also given KO-blind moving-root
screens. Exact `1.3 pL` normalization is a useful profile slice, not an
exclusion criterion by itself. At that slice all source-default roots close
numerically but miss the WT chloride gate; C8-MA rows are the corresponding
C6 aliases. The closest exact-volume row is C7 `2:1:3`, still at about
`-2.24` reported SEM in chloride.

### 4.2 Post-hoc finite capacity continuation

Because exact-volume normalization can hide a trade-off, a later KO-blind
robustness protocol continued each source-default capacity over factors
`1e-8` through `1e4` around its volume-normalized reference. The coarse scan
was refined at `0.025` log-decade spacing near transitions, connected pass
intervals were bracketed, and five fixed log-quantiles (`0.001`, `0.25`,
`0.5`, `0.75`, `0.999`) were selected. A pass requires the two-SEM WT Cl and
pH windows plus the new `+/-10%` volume sensitivity. This protocol was added
post hoc during adversarial review, but its axes and points were frozen before
held-out genotype simulations. It is not a measured uncertainty and not a
global search for disconnected or remote roots.

Only the three source-shaped `2:1:3` cores pass within the finite envelope:

| Core | Capacity interval | Interpretation |
| --- | ---: | --- |
| C6 reversible mass action | `[3.5105005159e-12, 5.9165256942e-12]` | Thermodynamically valid source-shaped pass band |
| C7 affinity law | `[0.0033498713, 0.0072388250]` | Thermodynamically valid source-shaped pass band |
| C8 symmetric saturation | `[0.0344779857, 0.0378955063]` | Thermodynamically valid pass band with `sigma = 1` as an unconstrained assumption |

The existence of these roots corrects the stronger fixed-row claim: the
intersection of source-supported mixed-cation structure and mathematical WT
closure is not empty after state relaxation. The roots are instead tested
predictively.

## 5. Beta/cAMP/PKA activation and held-out prediction

The static candidate set was incomplete without the primary evidence for
beta-adrenergic/cAMP/PKA stimulation. A positive common multiplier was applied
independently of calcium and outside the genotype multiplier. The tested fold
values were frozen from evidence and diagnostics:

- `1.25`: heterologous forskolin response, retained as positive
  AE4-specific evidence;
- about `3`: digitized native AE2-null/AE4-only isoproterenol response, the
  strongest direct AE4-specific constraint;
- `1.6`: total WT AE2-plus-AE4 exchange, reported only as an intermediate
  sensitivity rather than an AE4-specific fold;
- `1`: no-activation control;
- `10`: favorable inverse diagnostic outside the source-observed range.

The instantaneous step and alternative onsets are code-protocol conventions;
`120` or `130` code units are not relabelled as two or three physical minutes.
Common gating preserves stoichiometry, thermodynamic sign, and equilibrium.
The recruited-K C1 adapter is a novel flux-space localization, not a complete
thermodynamic transport law.

For the three rest-compatible capacity bands, every tested fold and capacity
predicts the wrong direction. Over the AE4-specific positive folds `1.25` and
`3`, the smallest ratio is the near-lower C6 capacity at fold `1.25`:

- endpoint `AE4-KO/WT = 1.01377776`;
- post-step integrated `AE4-KO/WT = 1.02093900`.

The separate total-WT fold-`1.6` sensitivity has minimum endpoint/integrated
ratios `1.01526971` and `1.02239373`. At fold `3`, the finite-band ranges are:

| Core | Endpoint ratio range | Integrated ratio range |
| --- | ---: | ---: |
| C6 `2:1:3` | `1.0201803` to `1.0262068` | `1.0272281` to `1.0334374` |
| C7 `2:1:3` | `1.0317096` to `1.0545913` | `1.0384749` to `1.0617369` |
| C8-SAT `2:1:3` | `1.0709456` to `1.0791195` | `1.0774353` to `1.0856713` |

Hence the maximum positive endpoint reduction over the qualified tested
rest-compatible source-shaped grid is **zero**. This statement is confined to
the finite capacity continuation, its fixed source-default Na:K weights,
`K_eq = 1`, `sigma = 1` where applicable, the historical bath convention, and
the tested positive folds. It is more defensible than using an empty strict
exact-volume set to invent a numerical maximum.

## 6. Full residual phenotype audit

### 6.1 Candidate-independent AE4-null rest

With AE4 activity exactly zero, every candidate reduces to the same non-AE4
vector field. The independently solved local AE4-null root closes to numerical
precision but gives approximately:

- `Cl_i = 48.7822 mM`, versus direct AE4-null `36.50 +/- 1.60 mM`;
- `pH_i = 7.3785`, versus direct AE4-null `6.89 +/- 0.02`;
- volume `2.974 pL` and `HCO3_i = 58.90 mM`.

Only the Cl and pH comparisons are direct genotype phenotype
falsifications. Volume and bicarbonate are model-reference diagnostics, not
measured AE4-null uncertainty. The direction is also wrong relative to the
rest-compatible C3/C6 WT roots: the chassis AE4-null root has higher, not
lower, intracellular chloride.

### 6.2 Candidate-specific AE2-null rest and dynamics

AE2 knockout is not candidate-independent because AE4 remains active. The C1
lineage AE2-null root has `Cl_i = 49.5249 mM`, or `-2.764` reported SEM from
the direct `54.50 +/- 1.80 mM` target; its pH is within the direct knockout
window. Candidate-specific AE2-null roots were also solved for every one of
the 15 surviving `2:1:3` capacity points. All close numerically and all fail
the direct chloride target:

| Core | AE2-null chloride range (mM) | Chloride residual range (reported SEM) | pH residual range |
| --- | ---: | ---: | ---: |
| C6 `2:1:3` | `46.012` to `46.289` | `-4.716` to `-4.562` | `+1.828` to `+2.412` |
| C7 `2:1:3` | `45.975` to `46.358` | `-4.736` to `-4.523` | `+1.216` to `+2.711` |
| C8-SAT `2:1:3` | `46.428` to `46.578` | `-4.484` to `-4.401` | `+2.663` to `+2.926` |

Dynamic AE2 secretion ratios remain comparatively near one, but the `5%`
normality threshold used by the automated internal screen is explicitly a new
sensitivity, not an experimental confidence interval. Near-normal secretion
cannot override the chronic chloride contradiction.

### 6.3 Conservation and thermodynamics

All numerically interpreted rows retain the explicit Na, K, Cl, and HCO3
source terms, transported charge, current balance, electroneutrality, and full
raw residual ledger. The surviving `2:1:3` rows are affinity-sign valid and
conserve to numerical tolerance. Solver completion alone is not labelled
resting closure; settledness and endpoint residuals are separate fields.

The common PKA gate does not alter the thermodynamic zero. Representative
in-band checks for each C6/C7/C8 core and all five folds were independently
cross-run with BDF and Radau; the reported bounds are below `3e-10` for
endpoint-ratio differences and below `3e-9` for integral-ratio differences.
That representative check covers one capacity per core and five folds (15
rows). A separate independent BDF run covers all 30 positive AE4-specific
fold/capacity rows. Against the frozen production Radau values, its exact
maximum endpoint and integral differences are `1.3592778e-9` and
`3.3704604e-9`; every BDF ratio is also above one. The representative reported
bounds and full-30 exact comparison are kept as distinct provenance records.

## 7. Round 3 localization and Round 4 decision

At the exact historical row, replacing the fitted Na-only AE4 source with a
published-style Na/K source creates a cation correction proportional to
`(-1,+1)`. Among existing scalar pathways, the Na/K-ATPase direction
`(-3,+2)` aligns best, with cosine `5/sqrt(26) = 0.98058`, but its best scalar
projection leaves residual fraction `1/sqrt(26) = 0.196116`. NHE1 or K efflux
alone leaves `1/sqrt(2)`. No one existing scalar pathway can exactly restore
the mixed-cation WT balance.

Exact closure requires two independent cation-homeostasis directions, for
example pump plus NHE1 (or an equivalent independent Na-loading direction),
or an unsupported direct `1:1` Na/K source. Pump plus K efflux is not the
exact positive decomposition of this correction.
The candidate-independent AE4-null resting failure separately localizes a
missing acid-base/chloride-homeostasis module, without identifying a unique
protein. The Na/K vector alignment is a resting-closure localization; it does
not prove that the pump explains the held-out 35% flow phenotype.

Round 4 was therefore **analytically screened and not simulated**. NHE1 alone
cannot cancel the K component, the pump alone has the exact `19.6116%` lower
bound above, and no independently calibrated magnitude exists for a single
component that restores every WT, AE2, Cl, pH, conservation, and
thermodynamic gate. A phenotype-tuned joint correction would require at least
two non-AE4 directions and would violate the requested one-component rule.

## 8. Evidentiary-domain summary

| Domain | Best endpoint result | Status |
| --- | ---: | --- |
| Strict exact-volume source/rest/thermodynamic table | No qualifying row | Empty slice; no numerical maximum inferred |
| Post-hoc KO-blind WT-rest robustness bands, positive AE4-specific folds | Minimum ratio `1.01377776`; maximum positive reduction `0%` | Qualified finite numerical exclusion; all rows wrong direction |
| Fold-`1.6` total-WT sensitivity | Minimum ratio `1.01526971` | Not AE4-specific evidence |
| Historical C1 exact-row control | Ratio `0.99855614`; decrement `0.144386%` | Reproducible but assay-inadmissible |
| Generous Round-2 closure-only grid | Ratio about `0.99852364`; decrement about `0.1476%` | Unsupported energetic/pure-Na diagnostic |

The tested AE4-only family therefore cannot satisfy the held-out secretion
phenotype and the independent resting phenotypes simultaneously. Reopening
identifiability or mechanism discrimination would be premature until the
non-AE4 cation and acid-base/chloride homeostasis closure is repaired from
independent data.

## 9. Reproducibility and artifacts

Authoritative machine-readable outputs are under
`results/12_ae4_mechanism_reconstruction/`:

- `model_comparison.csv/json` and `round2_envelope.csv` contain the frozen
  static comparisons;
- `resting_equilibria.csv/json`, `rest_capacity_envelope.csv/json`, and
  `thermodynamic_rest_profile.csv/json` contain the moving-root screens;
- `pka_activation_envelope.csv/json` contains the frozen stimulus predictions;
- `genotype_resting_equilibria.csv/json` contains direct genotype-target joins;
- `trajectories.csv` contains sparse solver-labelled landmarks at
  `t=0,50,100/right-limit,150,200` and both sides of any additional input
  change, while dense internal grids were used for the integrals;
- `calibration_manifest.json`, `solver_crosscheck.json`, and `summary.json`
  record provenance, solver conventions, tested-domain bounds, and decision
  logic.

The default historical regression uses BDF; the expanded moving-root and PKA
screen uses positivity-safe Radau, with the decisive rows checked against BDF.
LSODA, RK45, and DOP853 were not used for frozen claims because unconstrained
trial states can violate the chassis positivity domain. All absolute output
units, remote-root uniqueness, and physical code-time mapping remain explicit
limitations.
