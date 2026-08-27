# Adversarial audit of the AE4 mechanism reconstruction

## Verdict

**PASS for the scoped negative/localization conclusion; no AE4 candidate
passes as a successful mechanism.** The adversarial reopen was material: a
KO-blind capacity continuation found rest-pass `2:1:3` C6, C7, and saturating
C8 roots outside the earlier exact-volume slice. The earlier claim that the
complete source-admissible intersection was empty is therefore retracted.
After freezing the connected rest-pass intervals and their phenotype-blind
prediction points, however, every source-bounded positive PKA comparison has
both endpoint and stimulated-integral `AE4-KO/WT > 1`. The audited substantive
conclusion is

> **AE4 ALONE INSUFFICIENT; MISSING PATHWAY LOCALIZED**

and is valid only with the following scope:

1. the independently reimplemented, internally closed historical seven-state
   chassis is fixed;
2. the enumerated primary-source-compatible AE4 stoichiometries, empirical
   laws, reversible laws, saturation laws, candidate-specific resting roots,
   and beta/cAMP/PKA activity envelope are the tested AE4 family. The reopened
   moving-root bands fix the source-default Na:K weights `1.5:1.6`, `K_eq=1`,
   C8 saturation shape `sigma=1`, and AE4 bath HCO3 at `21 mM`; unsupported
   moving-root combinations of weight, energetic bias, saturation, or bath are
   not globally excluded; and
3. “localized” means a **cation-current/acid-base/chloride-homeostasis
   module**. At the frozen historical row the Na/K pump is the closest
   one-component cation direction, but no one-component correction is exact.
   This does not mean that a genotype-dependent pump, K channel, or NHE1
   change was experimentally established or that a joint correction was
   fitted.

The audit rejects any broader statement that `0.1476%` is a universal maximum
over every imaginable AE4 kinetic law. That number is a generous endpoint
bound only for the exact fixed-row, closure-only interpolation, and its
attaining candidate is pure-Na and energetically unsupported. It is not the
bound for the newly reopened moving-root family.

## Independence and evidence order

This audit was performed independently of the candidate ranking. The evidence
priority was:

1. conservation, charge, and dimensions;
2. primary transport and knockout experiments;
3. independently reproduced historical behavior;
4. the 2018 model text; and
5. explicitly labelled new modelling assumptions.

The held-out AE4-KO secretion ratio was unavailable to every calibration and
resting-root function. It appears only in post-freeze residual columns. In
particular, the source-derived PKA envelope `(1.25, 1.6, 3.0)` (with `1.6`
retained only as a total-WT sensitivity), the fixed Na:K maximal-response
weights `1.5:1.6`, the Round-2 structural grid, the
resting-state starts, capacity scan axes, interval refinement, and prediction
point fractions contain no `0.65` target.

The primary paper reports a `35 +/- 4.7%` reduction; the machine ledger maps
that linear complement to `KO/WT = 0.65 +/- 0.047`. The decision rule uses the
declared one-reported-SEM interval `[0.603, 0.697]` and requires **both** the
endpoint ratio and the stimulated-integral ratio to lie inside it. It does not
treat SEM as a population standard deviation or optimize an exact `0.65`.
The rest-pass `2:1:3` rows retain near-unity candidate-specific AE2
**secretion** endpoint ratios, but this does not rescue their inverse AE4
phenotype. Candidate-specific AE2-null resting roots were also solved for all
15 surviving capacities and every one fails the primary resting-Cl gate.

## Material disagreements and corrections

The audit did not merely confirm the first analysis. It blocked or reopened
the result at each of the following points:

- A static C1--C8 screen was incomplete after the 2021 primary PKA paper was
  identified. Positive beta/cAMP/PKA gates were added before an insufficiency
  conclusion could pass.
- Failure of the exact-`1.3 pL` C6 roots was initially generalized to the
  moving-root family. KO-blind continuation disproved that inference by
  finding C6/C7/C8-SAT `2:1:3` rest-pass intervals; the structural theorem was
  restricted back to the frozen historical row.
- A coarse capacity grid missed the narrow C8-SAT pass interval. The final
  scan uses `0.025`-decade resolution and Brent-refined edges.
- The initial genotype-rest export compared the C1 AE2-null root with AE4-WT
  targets and called it common. Direct AE2-null targets replace those values,
  the C1 label is candidate-specific, and all 15 surviving capacities receive
  separate AE2-null rest predictions.
- The first dual-solver description incorrectly implied all 75 rows had been
  rerun with both solvers. The certified bound applies to 15 representative
  dual-solver rows; the 30 source-fold BDF and full production Radau grids are
  separately identified.
- The audit rejects protein-specific localization. Pump/NHE/K-channel labels
  are balance directions; the evidence supports a homeostasis module, not a
  demonstrated genotype-dependent protein change.

## Hard-gate matrix

| Gate | Status | Independent finding |
| --- | --- | --- |
| E1. Held-out knockout | **PASS** | `calibrate_candidate` has no KO argument; the resting solver has no phenotype argument; PKA variant selection is source-defined; every machine row declares both calibrator and selector KO-blind. Changing the held-out value can change only post-freeze residuals. |
| E2. Fixed non-AE4 chassis | **PASS** | Every row and genotype carries fingerprint `55885266c44ecbbfe86be622c6ebd08b77151ee6ef7dfd5c8bbc8d4250b58e50`. At a fixed state, PKA changes only the AE4 source; pump, channels, water, tight junctions, AE2, NHE1, buffer, geometry, and voltages are identical. |
| E3. AE4 scale discipline | **PASS with scope** | Static capacity is fitted only to the WT balance vector. Exact-volume roots use a KO-blind WT-volume target; the later post-hoc continuation deliberately uses independent WT Cl/pH and a new `+/-10%` volume gate to select finite C6/C7/C8 capacity bands. No KO result enters either selection. PKA folds come from the independent assay envelope. Diagnostic `K_eq`, saturation, and Na/K grids are not promoted as calibrated biology. |
| E4. Common units | **PASS with limitation** | All candidates insert the same code-native amount-rate source into the same RHS and water map. Absolute flow and physical time are uncertified. No result is described as an absolute secretion rate, and code-time delays are not called minutes. |
| E5. Stoichiometric balances | **PASS** | Every active candidate satisfies its declared per-cycle Na, K, Cl, and HCO3 coefficients and zero transported charge. Parallel Na and K fluxes enter separate balances. KO annihilates every basal and recruited source exactly. |
| E6. Reversal and affinity | **PASS for labelled reversible laws** | C5-C8 reversible branches were independently placed on equilibrium and perturbed on both sides for Na and K. Positive PKA, saturation, and gating multipliers preserve zero and sign. C1-C4 and recruited-K C1 remain explicitly empirical/non-certified and cannot pass this gate as complete mechanisms. |
| E7. Primary-evidence priority | **PASS** | Direct Na and K transport, electroneutrality, reversibility, dose-response fits, the primary resting Cl/pH measurements, and PKA activation override the historical Na-only implementation and pooled 2018 prose. Candidate stoichiometries `1:1:2`, `1:2:3`, and `2:1:3` are labelled author-considered possibilities, not measured ion counts. |
| E8. Complexity and aliases | **PASS** | Fitted and independently fixed quantities are separated. C8 mass-action rows are exact aliases of their C6 rows and are not counted as independent support. Saturating C8 has one additional unconstrained shape choice. PKA fold scans remain sensitivities rather than a selected fit. |
| E9. Numerical reproduction | **PASS** | The pre-freeze 15-row representative dual-solver check is bounded by `3e-10` endpoint and `3e-9` integral. A post-freeze BDF rerun of the exact 30 source-fold rows agrees with frozen Radau within `1.360e-9` endpoint and `3.371e-9` integral; every ratio remains above one. Production Radau covers the full 75-row sensitivity grid. Conservation residuals are below `5e-14` and thermodynamic sign violations are zero. |
| E10. Convention robustness | **PASS with scope** | Endpoint and integral expose the C1 sign disagreement. Across five frozen points in each rest-pass `2:1:3` interval and the independently supported folds `1.25` and `3`, both metrics exceed one. Exact-volume, `+/-10%` volume, no-activation, and solver checks preserve the negative conclusion; alternate-start checks preserve the sign for the four original diagnostics, while the model's own AE4-null rest is independently falsified. Physical time is uncertified, so the numerical flow bound remains protocol-specific rather than a universal transient theorem. |
| E11. Stimulus/family exhaustion | **PASS for the declared source-default finite family** | The 2021 beta/cAMP/PKA evidence is included as a positive gate. C6, C7, and C8-SAT `2:1:3` rest-pass bands were all tested after KO-blind selection; every positive source-bounded gate is adverse. Unsupported moving-root weight/`K_eq`/`sigma`/bath combinations are outside this numerical exclusion. Carbonate and unresolved self-exchange remain explicitly unevaluable, not silently failed. |
| E12. Scope of structural proof | **PASS** | The exact collinearity theorem is reported only at the frozen historical row. Candidate-specific moving roots are solved separately; their existence is not represented as a contradiction of the frozen-row theorem. |

## Exact source-vector audit

Let positive `J` denote forward chloride entry. A cycle importing `a`
chlorides and exporting `b` monovalent cations and `c` bicarbonates contributes

\[
s=(-b_{Na}J,-b_KJ,+aJ,-cJ)
\]

in `(Na,K,Cl,HCO3)` order. Electroneutrality requires

\[
\boxed{c=a+b_{Na}+b_K}.
\]

This identity was checked for every active source, including PKA-gated sources
at full activation. For parallel branches the audit additionally verified

\[
s_{Na}=(-bJ_{Na},0,+aJ_{Na},-cJ_{Na}),\qquad
s_K=(0,-bJ_K,+aJ_K,-cJ_K),
\]

rather than accepting a pooled cation term.

At the fixed historical WT row the non-AE4 chassis demands

\[
s_0=J_0(-1,0,+1,-2),\qquad
J_0=0.03888333867648.
\]

Therefore exact closure by one scalar AE4 source requires collinearity with the
historical Na-only `1:1:2` ray. A nonzero net K source cannot close that same
row. At a per-species tolerance of `1e-8`, the exact bound is

\[
\frac{|J_{K,AE4}(u_0)|}{J_0}
\le 2.57179561745\times 10^{-7}.
\]

This is a bound on **net K source at that state**, not on K transport capacity.
A reversible K branch can have nonzero capacity but zero net flux at its
equilibrium. The primary direct-K evidence therefore conflicts with the
historical fixed-row source demand but does not mathematically forbid a
candidate-specific resting state.

The exact steady-state identities for a forward cycle with `c` chlorides,
`m` cations, and `b=c+m` bicarbonates are

\[
A=J_2+cJ_4,\qquad B=J_2+bJ_4,
\]

and hence

\[
J_4=\frac{B-A}{m},\qquad
J_2=\frac{bA-cB}{m}.
\]

The forward wedge is `A <= B <= (b/c)A`. Summing the cation/current balances
removes AE4 entirely and yields the common chassis identity

\[
C_{apical}=P_{NaK}+K_{out}.
\]

These identities explain why changing AE4 stoichiometry alone cannot repair a
missing cation-homeostasis relation.

## Thermodynamic audit

For a `1 Cl : 1 cation : 2 HCO3` forward cycle,

\[
Q_C=\frac{[Cl]_i[HCO_3]_e^2[C]_e}
{[Cl]_e[HCO_3]_i^2[C]_i},\qquad
A_C=\log(K_{eq}/Q_C).
\]

Every candidate called reversible was independently required to satisfy
`J_C=0` at `A_C=0` and `sign(J_C)=sign(A_C)` above and below equilibrium for
both `C=Na` and `C=K`. The analogous products were used for `1:2:3` and
`2:1:3`. Common PKA multiplication is positive and preserves these signs.

The historical and published empirical polynomials can reverse independently
of the declared passive affinity and are correctly labelled
`empirical_not_sign_controlled`. The novel PKA-recruited K branch is also
rejected as a complete mechanism because its K direction is locked to the C1
Na polynomial rather than a K affinity.

## Independent numerical reproduction

### Chassis and solver

The historical resting row was independently reproduced with maximum raw RHS
below `5.3e-12`, `(V_a,V_b)=(-50.24,-62.8)` mV, and cell volume `1.3` pL.
The AE4-KO vector field is exactly candidate-independent because genotype zero
multiplies the complete AE4 source.

| Cross-check | BDF | Radau | Absolute ratio difference |
| --- | ---: | ---: | ---: |
| C1 endpoint KO/WT | `0.998556140220` | `0.998556141804` | `1.58e-9` |
| C1 stimulated integral KO/WT | `1.005638156093` | `1.005638123330` | `3.28e-8` |
| Diagnostic Round-2 endpoint bound | `0.998523640298` | `0.998523640428` | `1.31e-10` |
| Diagnostic Round-2 integral | `1.005603798734` | `1.005603791143` | `7.59e-9` |

The last two rows are a pure-Na, `K_eq=100`, saturating diagnostic and are not
source-admissible.

### Static candidates

| Candidate group | Fixed-row status | Endpoint AE4-KO/WT | Stimulated-integral AE4-KO/WT | Audit disposition |
| --- | --- | ---: | ---: | --- |
| C1 historical Na-only | exact | `0.99855614` | `1.00563816` | Reject: direct K evidence, empirical reversal, phenotype |
| C2 published pooled | fails; AE4-vector error `0.0328154` | `1.02905601` | `1.03680886` | Reject: closure, pooled reversal convention, phenotype |
| C3 explicit Na/K | fails; AE4-vector error `0.0331963` | `1.02952399` | `1.03728261` | Reject: fixed-row closure, empirical reversal, phenotype |
| C4 and default C5-C8 | nonnegative fit sets capacity to zero | `1.0` | `1.0` | Reject: no calibrated AE4 contribution |
| Closed Round-2 pure-Na boundary | exact diagnostic only | best `0.99852364` | best `1.00560381` | Reject: direct K evidence and unsupported energetic bias |

Only C1 closes the historical row, and it predicts a `0.144386%` endpoint
loss while the stimulated integral changes in the opposite direction. The
closure-only Round-2 grid improves the endpoint loss to merely `0.147636%` and
still produces no integrated decrement.

### Candidate-specific resting roots

All roots below were solved without the knockout-flow target. The production
screen uses deterministic broad multistart, verifies no active state bound,
and enforces raw RHS `<=1e-8`. The independent primary WT gates are
`Cl_i=50.10 +/- 1.50 mM` and `pH_i=6.91 +/- 0.07` (reported SEMs); a two-SEM
limit is used.

| Core | Max raw RHS | Rest `Cl_i` (mM) | Cl z-score | Rest pH | Volume (pL) | Rest gate |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| C1 | `2.87e-13` | `50.000` | `-0.067` | `6.9100` | `1.3000` | **PASS** |
| C2 | `1.11e-13` | `46.353` | `-2.498` | `6.9078` | `1.2482` | FAIL |
| C3 | `8.56e-15` | `46.313` | `-2.525` | `6.9082` | `1.2479` | FAIL |
| C6 `1:1:2`, volume-normalized | `1.55e-13` | `43.476` | `-4.416` | `7.0073` | `1.3000` | FAIL |
| C6 `1:2:3`, volume-normalized | `8.77e-14` | `40.675` | `-6.283` | `7.0415` | `1.3000` | FAIL |
| C6 `2:1:3`, volume-normalized | `1.22e-13` | `45.739` | `-2.907` | `6.9761` | `1.3000` | FAIL |

These exact-volume roots are only one calibration slice and are **not** a
global exclusion. A broader KO-blind capacity continuation independently
reproduced the following source-default rest-pass roots:

| Core | Capacity | Max raw RHS | Rest `Cl_i` | Cl z-score | Rest pH | pH z-score | Volume |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| C6 `2:1:3` | `3.58681628117e-12` | `5.35e-13` | `47.3702` | `-1.8199` | `7.04891` | `1.9845` | `1.41137` |
| C7 `2:1:3` | `0.00605934707159` | `1.20e-13` | `47.2159` | `-1.9227` | `7.00153` | `1.3076` | `1.34874` |

Both solves succeeded from the neutral baseline, touched no bound, and read no
KO target. They sit close enough to the two-SEM/volume gates that a single
selected capacity is not adequate evidence. Brent refinement of the observed
connected pass regions gives the following finite intervals. The five
prediction points were fixed at log-interval fractions
`0.001, 0.25, 0.5, 0.75, 0.999` before any knockout run.

| Core | Finite rest-pass capacity interval | Boundary-limiting gates | Frozen prediction points |
| --- | ---: | --- | ---: |
| C6 `2:1:3` mass action | `[3.5105005e-12, 5.9165257e-12]` | pH at lower edge; Cl at upper edge | 5 |
| C7 `2:1:3` affinity | `[0.0033498713, 0.00723882497]` | pH at lower edge; Cl at upper edge | 5 |
| C8 `2:1:3` saturation | `[0.0344779857, 0.0378955063]` | pH at lower edge; Cl at upper edge | 5 |

The scan covers `1e-8` through `1e4` times each KO-blind exact-volume
capacity, with `0.025`-decade resolution around the observed transition. It
establishes these finite connected branches, not global absence of remote or
disconnected mathematical roots. The `+/-10%` volume gate is a new post-hoc,
KO-blind robustness convention; it is not an experimental error bar.

### Beta/cAMP/PKA family-exhaustion test

The 2021 primary follow-up prevents a static-family-only insufficiency claim.
It supports positive AE4 activation by beta/cAMP/PKA, with a source sensitivity
envelope of `1.25`, approximately `1.6`, and approximately `3.0`. It does not
identify whether the fold is turnover, active transporter number, or changed
stoichiometry, and its initial Cl-uptake assay does not supply a physiological
chassis time constant.

The PKA input is therefore separate from calcium, zero at rest, and applied as
a positive gate after each capacity/rest calibration is frozen. An immediate
model step and delayed steps at 20%/30% of the code-native stimulated window
are convention sensitivities only; they are not called 2-3 minutes.

Representative fold-3 calculations are shown below. The first four simulated
rows are frozen production Radau values; the reopened in-band rows are
independent BDF points inside their frozen capacity intervals.

| Candidate | Endpoint AE4-KO/WT | Stimulated-integral AE4-KO/WT | Independent gate that fails |
| --- | ---: | ---: | --- |
| C1 PKA common | `0.99855614` | `1.00563812` | Direct K assay and phenotype |
| C1 PKA K recruit | `1.03403239` | `1.04112562` | Novel direction-locked K law, thermodynamics, phenotype |
| C3 PKA common | `1.03848048` | `1.04610120` | Resting Cl, empirical reversal, phenotype |
| C6 `1:1:2` PKA common | `1.12216869` | `1.13012675` | Resting Cl and phenotype |
| C6 `1:2:3` PKA common | not simulated in exact-volume slice | not simulated | KO-blind resting Cl gate failed in that slice |
| C6 `2:1:3` PKA common, reopened G | `1.02038110` | `1.02743633` | Phenotype; independent BDF point inside frozen band |
| C7 `2:1:3` PKA common, reopened G | `1.04758775` | `1.05465508` | Phenotype; independent BDF point inside frozen band |
| C8 `2:1:3` saturating PKA common | `1.07577345` | `1.08230181` | Phenotype; closest-volume independent BDF point |

The complete high-accuracy production Radau table is authoritative; the
isolated BDF points serve only as sign-and-scale checks, while the exact
30-row BDF comparison below is the decisive independent solver audit. Over all
15 frozen interval points and the AE4-specific source folds `1.25` and `3`,
the best result is the C6 near-lower point at fold `1.25`:

| Audited source/rest/thermodynamic domain | Minimum endpoint KO/WT | Minimum integral KO/WT | Maximum positive reduction |
| --- | ---: | ---: | ---: |
| Source-default C6/C7/C8-SAT `2:1:3`; five points/core; folds `1.25,3` | `1.01377776` | `1.02093900` | `0%` |

The `1.6` total-WT exchange fold, retained only as a sensitivity because it is
not AE4-specific, has minima `1.01526971` and `1.02239373`. At fold `3`, the
full endpoint/integral ranges are C6
`1.020180--1.026207 / 1.027228--1.033437`, C7
`1.031710--1.054591 / 1.038475--1.061737`, and C8-SAT
`1.070946--1.079120 / 1.077435--1.085671`. Every value is adverse.

After the Radau artifact freeze, an independent BDF rerun used the exact same
15 capacities and folds `{1.25,3}`. Its minima are
`1.013777757954 / 1.020939000438`; the largest absolute BDF--Radau differences
over all 30 rows are `1.3593e-9` endpoint and `3.3705e-9` integral. All 30
ratios remain above one. These full-grid values are distinct from the earlier
15-row representative dual-solver bound recorded in `solver_crosscheck.json`.

Positive common activation cannot impose an arbitrary fold on the steady AE4
flux: the state relaxes until the non-AE4 balance demand is restored.
Selective K recruitment changes the C1 source relative to common gating by

\[
\Delta s=(+(F-1)J,-(F-1)J,0,0),
\]

which is precisely the adverse Na/K correction exposed by the fixed chassis.

On a fine early code-time grid, the smallest KO/WT ratios over the first three
saved post-step code units were approximately `0.96185` (C1 common), `0.96185`
(C1 K recruit), and `0.94607` (C3 common), followed by ratios above one. None
approaches `0.65`, and each fails another hard gate. At the closest-volume
fold-3 `2:1:3` points, the ratios at code coordinates `t=120/130` are
`1.03145/1.02610` (C6), `1.06035/1.05485` (C7), and
`1.08185/1.07623` (C8-SAT). These are **not minute-labelled predictions**,
because the historical time scaling is uncertified. The source paper's
initial 2--3 minute slope can constrain a fold but cannot identify a chassis
activation time constant.

The full `2:1:3` interval grid was not repeated for every delayed onset or
from genotype-specific knockout equilibria. Production initializes all
genotypes at the candidate WT root, evolves them independently for 100 code
units, and then applies the gate. This prevents a universal transient or
physical-time maximum claim. It does not reopen mechanism success: the
common tested AE4-null rest branch and every candidate-specific surviving
AE2-null rest already fail direct genotype phenotypes, while both endpoint and
integral are adverse under the declared common protocol.

The common modeled AE4-null rest closes with raw residual `5.2e-14` at
`Cl_i=48.7822 mM`, `pH_i=7.37847`, and volume `2.97418 pL`. It misses the
primary AE4-null targets by `+7.68` and `+24.4` reported SEM for Cl and pH,
respectively. Because genotype zero annihilates the entire AE4 source, its
vector field is candidate- and PKA-law independent. The reported equilibrium
is the common 20-start local branch, not a proof of unique convergence from
every initialization or absence of remote roots. The
**C1-specific** modeled AE2-null root is `Cl_i=49.5249 mM`,
`pH_i=6.9141`, and volume `1.29584 pL`.
Against the direct primary AE2-null targets `54.50 +/- 1.80 mM` and
`6.95 +/- 0.05`, these are `-2.764` and `-0.718` reported SEM: pH passes but
resting Cl fails. Even comparison with the AE2 littermate control Cl target
`53.40 +/- 1.80 mM` gives `-2.153` SEM. Thus near-unity AE2 secretion must not
be described as complete AE2 phenotype agreement. Because AE2 deletion leaves
AE4 active, this root depends on the AE4 law and capacity and is not a common
root. Independent two-start local solves at all 15 surviving capacities give
AE2-null Cl z-score ranges `-4.716..-4.562` (C6), `-4.736..-4.523` (C7), and
`-4.484..-4.401` (C8-SAT), with raw residuals below `3e-13`; all fail the
direct `54.50 +/- 1.80 mM` target. Their pH z-score ranges are
`+1.828..+2.412`, `+1.216..+2.711`, and `+2.663..+2.926`, respectively. Thus
criteria 4--5 fail independently of the inverse AE4 secretion prediction.
The negative decision does
not depend solely on mapping a particular code-time onset to the experiment;
the quantified secretion ratios nevertheless remain predictions of the
declared protocol.

## Aliases, complexity, and nonuniqueness

`C8_112_MA`, `C8_123_MA`, and `C8_213_MA` are numerically exact aliases of
`C6_112`, `C6_123`, and `C6_213`. They are retained for requested-family
traceability but collapsed for ranking and evidence counting. The default C8
saturation factor is an unconstrained modelling choice and incurs an
additional complexity unit. Non-unit `K_eq`, pure-Na capacity, and AE4-only
bicarbonate-bath overrides are diagnostic axes, not independent source facts.

The two primary-paper alternatives that lack a carbonate state or a resolved
bicarbonate self-exchange stoichiometry remain explicit inactive registry
entries. Omitting their simulation is not evidence against their chemistry;
they are not defined seven-state mechanisms.

## Localization and Round 4

Replacing the historical Na-only source by the published resting Na/K
partition produces the fixed-row correction

\[
d=(-0.03217931477,+0.03217931477,0,0).
\]

The Na/K-pump signature `(-3,+2,0,0)` has absolute cosine alignment
`5/sqrt(26)=0.98058068`, but its best one-component residual is
`1/sqrt(26)=0.19611614`. No single source-backed scalar chassis correction is
exact. The smallest exact **positive** two-direction decomposition uses equal
code-native cycle-flux increments of the Na/K pump and NHE1 (or an equivalent
Na-loading direction). Pump plus K efflux spans the cation plane only with a
signed decrease in the efflux direction, so it is not the stated positive
repair.

This is a flux-space localization only. The primary knockout study reports
unchanged NHE-dependent alkalinization and unchanged resting pH, so it does not
support a genotype-specific NHE1 correction. Pump density/whole-cell scaling
also remains unresolved. A Round-4 joint parameter fit would therefore add an
uncalibrated non-AE4 degree of freedom after seeing the phenotype and is not
auditable as held-out reconstruction. Not performing that fit is the correct
anti-leak decision.

## Candidate disposition

| Candidate/family | Complete assay | Rest/WT closure | Thermodynamics | Held-out phenotype | Final |
| --- | --- | --- | --- | --- | --- |
| Historical C1, with or without common PKA | FAIL | PASS | empirical only | FAIL | **REJECT** |
| Published pooled C2 | partial | FAIL | empirical/reversal-invalid | FAIL | **REJECT** |
| Explicit empirical C3, with or without PKA | PASS for direct Na/K | FAIL primary resting Cl | empirical only | FAIL | **REJECT** |
| Cooperative C4 | interpretive | FAIL/zero capacity | empirical only | FAIL | **REJECT** |
| Gated/thermodynamic C5 | direct-source compatible | FAIL/zero capacity | PASS | FAIL | **REJECT** |
| Reversible C6/C7/C8 cores | direct-source compatible | finite `2:1:3` moving-root bands PASS only under the new volume sensitivity | PASS | every frozen source-fold point is adverse | **REJECT as phenotype reconstruction** |
| PKA-recruited K C1 | novel diagnostic | PASS at rest | FAIL | wrong direction | **REJECT** |
| Pure-Na Round-2 closure boundary | FAIL direct K | PASS fixed row | unsupported bias | FAIL | **REJECT** |
| Carbonate/self-exchange entries | unresolved | not constructible | not testable | not testable | **UNEVALUABLE, not successful** |

## Reproducibility and release conditions

The independent audit suite in `tests/test_ae4_adversarial.py` enforces the
registry, charge/source signatures, separate branch insertion, capacity-times-
genotype semantics, exact aliases, thermodynamic zeros/signs, PKA locality and
KO annihilation, common fingerprint, anti-leak interfaces, moving roots, units,
and archive independence. It is intentionally fast; the frozen stimulus and
multistart tables provide the high-cost numerical regression surface.

The frozen artifact regression reads, rather than merely trusts, the produced
JSON. It finds `730` capacity-profile rows including `15` frozen prediction
points, `105` stimulus rows including `75` interval/fold sensitivities and
`30` source-specific positive-fold rows, and `17` genotype-rest rows including
the common tested AE4-null branch, the C1 AE2-null diagnostic, and `15`
candidate-specific AE2-null predictions. The `30` source/rest/thermodynamic
rows all have endpoint and integral above one; zero genotype-rest or complete
hard-gate rows pass.

This final audit `PASS` is issued only for the exact classification and scope
at the start of this report; it is not a `PASS` for any candidate mechanism.
The frozen tables show KO-blind capacity selection, join the one-reported-SEM
held-out compatibility rule (`0.65 +/- 0.047`) only after simulation, and show
that no source/rest/thermodynamic row passes both endpoint and integral
phenotype checks.

Two complete same-code generations produced byte-identical hashes for all 20
result-directory and six figure-directory files. A subsequent plotting-only
refactor changed only the closure-panel label placement; two renders through
the frozen-CSV helper were byte-identical with PNG SHA-256
`1fcfbca5fce3b12bea73d6891b812776e4b552ad880653fdb722db8dc9427471`.
The final manifest SHA-256 is
`6ff2aa8a53df7e4fb3cd4017bd6982d41de9d84cb4236da6c80aff6927d62486`,
and all 31 entries independently match their files.
The four frozen-artifact audit tests and the complete 126-test repository
suite pass, Python compilation succeeds, and `archive/` remains unchanged.

Any result that (i) calls code time minutes, (ii) promotes the
`0.1476%` closure-only diagnostic to a universal bound, (iii) labels an
exact-volume C6 root source-complete despite its Cl failure, (iv) counts
C6/C8 mass-action aliases independently, (v) calls the post-hoc `+/-10%`
volume convention experimental uncertainty, or (vi) attributes the localized
balance residual to a specific altered protein will revoke this audit.
