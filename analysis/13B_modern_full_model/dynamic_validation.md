# Task 13B WT dynamic validation and stimulus reconstruction

## Pre-reveal controlling disposition

The old G5 freeze is now **REJECTED_SUPERSEDED**.  Its 18 root/profile groups
passed the gates that existed when it was frozen, but two later native-source
audits exposed independent hard failures: `ABSOLUTE_GLAND_GEOMETRY_SCALE_FAIL`
and `NATIVE_CONDUCTANCE_SOURCE_MAP_FAIL`.  Its implied whole-gland cell counts
are incompatible with the Kondo et al. (2015) acinar-fraction bound, and its
common conductance multiplier does not implement the Palk/Gin native
whole-cell conductance map.  Therefore the old manifest's internal
`reveal_eligible: true` field is historical payload, not current permission.
The controlling disposition is
`results/13B_modern_full_model/g5_supersession.json`; held-out reveal remains
prohibited while the native-source reconstruction hierarchy is pending.

The historical selected resting reference was
`G2_BALANCED_APICAL_K_BIASED_G10:M00`.  The dimensional audit corrected the
historical-lineage water coefficients to `Pa=0.00432`, `Pb=0.0515`, and
`Pt=0.000260 pL s^-1 (mOsm/L)^-1` before any final G5 profile was frozen.
The two old effective-Ca/NKCC pairs remain useful failure evidence but are not
production candidates.  Neither is silently promoted to a point estimate.

The replacement native WT-rest hierarchy has completed canonical H1, both
hydraulic diagnostics, and the full three-start N/AN absolute-source grid.
Across 1,125 panels and 3,375 attempts it finds 739 numerical roots.  The first
495 H1/HW/HWQ panels have zero WT passes; HW and HWQ are diagnostic-only H>1
modes whose scale-one forms are exact H1 alias controls.  The 630-panel
absolute tier has ten preliminary WT passes, all `N_ABS_NKCC` scale 4, and
`AN_ABS_AE4_NHE` has none.  Seven-start confirmation and 33-start geometry are
both incomplete, so `eligible_dynamic_root_ids` remains empty and no native
candidate has yet reached this document's 600-second gates.

No AE4-null secretion target file, magnitude, or time-course datum was opened,
hashed, parsed, or passed to this layer.  No knockout trajectory is implemented
in `src/modern_full_model/validation.py` or
`src/modern_full_model/run_analysis.py`.

## 1. Physical protocol and distinct inputs

The source protocol and its model map are frozen separately:

| Quantity | Frozen representation | Evidence class | Limitation |
|---|---:|---|---|
| Gland protocol duration | `600 s` | primary measurement, P15 | Flow was read once per minute. |
| CCh dose | `0.3 uM` | primary measurement, P15 | Dose does not measure the intracellular Ca waveform. |
| IPR dose | `5 uM` | primary measurement, P15/P21 | Dose is mapped to a normalized beta input, not to a measured receptor occupancy. |
| Resting Ca input | `0.058 uM` | historical-lineage WT-root convention | Not a new primary measurement. |
| Stimulated Ca input | `0.55 uM` | historical-lineage amplitude sensitivity | The old code-time onset is rejected. |
| Beta input | `0 -> 1` | new normalized input convention | It is not PKA activity and has no fitted delay. |

`SecretagogueProtocol` exposes four arms: rest, CCh only, IPR only, and
CCh+IPR.  CCh changes only `Stimulus.calcium_uM`; IPR changes only
`Stimulus.beta_input`.  Thus calcium cannot accidentally stand in for PKA and
beta activation cannot directly open Ca-gated channels.

At exactly `t=0` the protocol returns the resting input so a frozen resting
root and basal regulatory state are unambiguous.  The stimulated right limit
is represented at `t=1e-6 s`; the production grid then uses physical seconds
through `t=600 s`.  No inherited `rhs_scale` or old code-time conversion is
used.

The instantaneous calcium step is a declared waveform sensitivity, not a
claim that cytosolic Ca rises instantaneously in the experiment.  No matched
P15 intracellular Ca trace was available to identify a rise, oscillation, or
decay law.  Adding a free calcium time constant would therefore add flexibility
without evidence.

## 2. WT dynamic evidence that may be used

Only these WT measurements enter the dynamic layer:

| Target | Value or constraint | Use |
|---|---|---|
| Matched Ae4-line WT resting Cl | `50.10 +/- 1.50 mM` | Upstream WT root calibration. |
| Matched Ae4-line WT resting pH | `6.91 +/- 0.07` | Upstream WT root calibration. |
| CCh-only initial Cl-uptake assay | `2.18 +/- 0.20 x10^-3 s^-1` | Arm-specific WT diagnostic. |
| IPR-only initial Cl-uptake assay | `0.40 +/- 0.07 x10^-3 s^-1` | Arm-specific WT diagnostic. |
| CCh+IPR initial Cl-uptake assay | `2.02 +/- 0.10 x10^-3 s^-1` | Combined-arm WT diagnostic. |
| Isolated-exchanger control/IPR contrast | `1.56 +/- 0.19` versus `4.33 +/- 0.91 x10^-3 s^-1` | Context-specific beta-response constraint. |
| WT whole-gland flow envelope | graphically about `9--10 uL/min`, `n=6` | One WT-only absolute scale, not a minute-wise point fit. |
| P21 native and construct perturbations | IPR/forskolin/PKAc increase AE4-associated activity; H89 prevents the IPR/forskolin increment; S173A loses, S273A retains the PKAc response | Regulatory direction and hierarchy. |

The model reports the exact concentration derivative
`(dCl_amount/dt - Cl*dV/dt)/V` and its fractional rate for each agonist arm.
It does **not** silently equate this derivative with normalized SPQ
fluorescence.  A quantitative fit to the three uptake numbers requires an
explicit SPQ observation map and assay preconditioning protocol.  The P15
combined value also shows why an additive CCh-plus-IPR target would be
incorrect.

## 3. What the 2021 evidence does and does not identify

P21 identifies the following causal perturbation hierarchy:

`beta agonist -> adenylate cyclase/cAMP perturbation -> PKA-dependent process -> S173-dependent AE4 activation`.

It does not identify direct phosphorylation of S173.  It also does not report
an acute cAMP trace, PKA-activity trace, AE4 regulatory-state trace,
phosphorylation trace, washout constant, or beta-onset AE4 trace.  The plotted
time axes are the imposed chloride-depletion/reuptake assay used to estimate an
initial transport slope; they are not an IPR-onset experiment.  The 18--20 h
PKAc coexpression contrast supplies no acute rate either.

Accordingly the regulatory nesting comparison has the following status:

| Family | Dynamic coordinates | What it can represent | Pre-reveal status |
|---|---|---|---|
| R0 | none | instantaneous/static beta-dependent capacity control | Required negative control. |
| R1 | one effective PKA/AE4 activation fraction | Smallest dynamic beta-to-AE4 model | Preferred on parsimony if dynamics must be retained; its time constant is unmeasured. |
| R2 | effective cAMP plus downstream activation | A two-stage delay | Endpoint perturbations cannot separate its two rates from R1. |
| R3 | upstream PKA activity plus an AE4 regulated/phosphorylation-like fraction | Reversible regulated AE4 fraction | S173 logic constrains response presence, not the two rates or direct phosphorylation. |

R1--R3 kinetic landmarks (`t10`, `t50`, and `t90`) are reported as
`UNIDENTIFIED_SENSITIVITY_NOT_FIT`.  They are never optimized against a
knockout secretion curve.  State-specific AE4-cycle modulation is not part of
this layer because P21 does not distinguish it from common active-capacity
modulation.

The predeclared freeze contains 20 members: R0 at each of the `1.25` and
`1.70` fully activated capacity readings, plus R1--R3 at both gains for three
unmeasured timing profiles.  `TFAST` uses upstream/downstream sensitivities of
`3/10 s`, `TREFERENCE` uses `10/30 s`, and corrected `TSLOW` uses `30/90 s`.
R1 uses only the corresponding downstream time (`10/30/90 s`), and symmetric
R3 rates preserve exact R2/R3 isomorphs.  These profiles define the current
native dynamic contract; they do not identify a unique family or kinetic
parameter.  Missing acute kinetic identification is carried as an ensemble
axis, but no ensemble is reveal-eligible until all other native WT gates pass.

## 4. Root and parameter intake rule

The dynamic suite accepts, but never selects or modifies:

1. one production WT charge-manifold root;
2. a `ModelVariant` and `WTCalibrationParameters` tuple;
3. the deterministic root-search rule and its branch identifier;
4. identical whole-cell and AE4 parameters for R0--R3; and
5. regulatory-family objects whose only difference is the nested regulatory
   subsystem.

The runner hashes the root, whole-cell parameters, AE4 parameters, regulatory
objects, input protocol, flow scale, and solver specifications.  It rejects a
family comparison if a regulatory family changes any whole-cell or AE4
transport parameter.

The old G2/G5 WT panel retained nine numerically rest-passing branches spanning
three cation routings and three apical topologies, all on the now-rejected
`G=10` conductance slice.  The predeclared historical selection rule chose
`G2_BALANCED_APICAL_K_BIASED_G10:M00`, then used converged-start membership and
dimensionless residual only as tie breakers.  It converged from all seven
calibration starts and from 30 of 33 fixed-geometry starts to one branch.  The
selection did not use any held-out result.

## 5. Numerical protocol

Radau and BDF are two different stiff algorithms exposed by SciPy.  Their
agreement is a dual-method check, not the independent implementation supplied
by Agent I.

| Run | Method | `rtol` | amount `atol` | volume `atol` | regulatory `atol` | max step |
|---|---|---:|---:|---:|---:|---:|
| production | Radau | `1e-7` | `1e-10 fmol` | `1e-12 pL` | `1e-10` | `2 s` |
| crosscheck | BDF | `1e-7` | `1e-10 fmol` | `1e-12 pL` | `1e-10` | `2 s` |
| loose | Radau | `1e-6` | `1e-8 fmol` | `1e-10 pL` | `1e-8` | `5 s` |
| tight | Radau | `1e-8` | `1e-11 fmol` | `1e-13 pL` | `1e-11` | `1 s` |

Absolute state differences remain labelled by coordinate.  The cross-family
summary uses a dimensionless coordinate-wise relative difference whose floor
is the corresponding solver's native-unit `atol/rtol` crossover.  It never
takes an absolute maximum across fmol, pL, and fractions.

Conservation residuals likewise retain separate units and tolerances for
fmol/s species or charge balances, pL/s water balance, amperes, and mM
speciation closure.  The only combined gate is the maximum of the resulting
dimensionless residual/tolerance ratios.

Nearby-state tests use `0.1%` deterministic perturbations:

- cell and lumen volume up/down with all amount charges unchanged; and
- equal intracellular and luminal Na/Cl amount pairs up/down.

Both preserve cell and lumen bulk charge exactly to numerical precision.

## 6. WT flow scaling and absolute geometry gate

The model produces single-cell local-lumen outflow in `pL/s`, whereas P15
reports whole-gland `uL/min`.  One WT-only multiplier can map the simulated
0--600 s mean to the midpoint of the graphical WT envelope.  The artifact also
reports the implied effective secretory cell count using
`1 pL/s per cell = 6e-5 uL/min per cell`.

This multiplier is not, by itself, validation: any positive trace can be
rescaled.  It is therefore labelled a Tier-2 WT observation scale.  The Kondo
et al. (2015) `66.91%` acinar fraction now supplies a binding upper-bound map:
one generous `100 mg` SMG at unit density contains at most `51,469,231`
modeled `1.30 pL` acinar cells.  Hence both the scale-independent minute-shape
test and the absolute cell-count test are binding.  A nonempty arbitrary
multiplier interval cannot waive the geometry gate.

### Frozen-reference baseline G5 result

All four reference-timing regulatory families completed 600 s with Radau and
BDF on the selected root:

| Family | Endpoint Cl (mM) | Endpoint pH | Endpoint volume (pL) | Endpoint flow (pL/s) | Cumulative flow (pL/cell) |
|---|---:|---:|---:|---:|---:|
| R0 | 20.02993 | 6.80976 | 0.96479 | 0.000635835 | 0.498934 |
| R1 | 20.03801 | 6.81026 | 0.96489 | 0.000636065 | 0.499431 |
| R2 | 20.04074 | 6.81043 | 0.96493 | 0.000636143 | 0.499590 |
| R3 | 20.04074 | 6.81043 | 0.96493 | 0.000636143 | 0.499590 |

The maximum coordinate-normalized Radau/BDF difference was
`4.24e-7`; cumulative-flow relative differences were at most `6.77e-8`.
R2 and R3 agree to numerical precision, as required by their exact default
isomorphism.  Agent I's independent R1 transcription reproduced endpoint Cl,
pH, flow, cumulative flow, and capacity to the reported numerical precision.
Thus this baseline failure is not a solver or transcription failure.

The provisional R1 minute-wise model flows fell from `0.0011026 pL/s` at
minute 1 to `0.0006361 pL/s` at minute 10, a max/min ratio of about `1.73`.
The declared WT graphical envelope permits a ratio of only `10/9 = 1.11`.
Equivalently, the intervals of multipliers that would put each minute in
`9--10 uL/min` have an empty intersection.  Scaling the 600-s mean to
`9.5 uL/min` produces about `12.58 uL/min` at minute 1 and `7.26 uL/min` at
minute 10; only one of the ten post-stimulus minute landmarks lies within the
declared envelope.  The same mean scale implies about `1.90e8` effective
secretory cells.

This is a scale-independent WT flow-shape failure for the selected root,
qualified by the fact that the source envelope is graphical rather than a
tabulated minute-wise series.  The production runner now tests multiplier-
interval feasibility before fitting a mean scale, so no normalization can
hide this failure.

### Historical targeted stimulated-NKCC repair and superseded WT ensemble

The first targeted repair retained the `0.058 -> 0.55 uM` historical calcium
anchor and nested the original stimulus-independent NKCC law (`N0`), an
algebraic CCh/Ca-dependent capacity (`N1`), and a one-state effective dynamic
capacity (`N2`).  No nominal P15-centered gain repaired the fixed flow shape
on any retained root.  Only deliberately out-of-scale, slow N2 probes passed,
which established that NKCC sustainment is decision-controlling but did not
make those probes physiologically eligible.

The second targeted generation treated calcium as an effective
channel-driving coordinate calibrated only from WT flow, not as measured free
intracellular calcium.  Each algebraic N1 law was normalized to its own
protocol endpoint.  The old joint-repair matrix contained `9 roots x 3 profiles x 20
regulatory members x 2 solvers = 1080` integrations; every integration passed
the numerical gate.

| Joint profile | Selected-root individual ratio range | One common scale across all members/solvers | Result |
|---|---:|---:|---|
| effective Ca `0.10 uM`, N1 multiplier `1.75` | `1.0111--1.0527` | feasible interval `17409.88--18375.76`; frozen scale `17682.38 uL min^-1 per pL s^-1` | pass on all 9 roots |
| effective Ca `0.15 uM`, N1 multiplier `1.75` | `1.0973--1.1589` | empty selected-root interval | fail on all 9 roots |
| effective Ca `0.12 uM`, N1 multiplier `2.0` | `1.0292--1.0827` | feasible interval `15439.59--15844.91`; frozen midpoint `15642.25 uL min^-1 per pL s^-1` | pass on all 9 roots |

For the selected root, the `0.12/M2` profile gives ensemble endpoint ranges
of `42.750--45.985 mM` Cl, `6.87368--6.89297` pH, and
`1.19542--1.23935 pL` volume.  Its onset mechanistic NKCC fractional-Cl
diagnostic is about `2.33e-3 s^-1`; this is reported with an explicit assay-
context caveat and is not equated to the P15 total SPQ slope.  Kondo et al.
(2015) support a sustained SMG calcium-response shape under CCh, but their
ratiometric signal does not identify either absolute effective-Ca value.

The shared whole-gland scales imply `294,706,403` and `260,704,116` effective
cells for the selected-root `0.10/M1.75` and `0.12/M2` profiles.  Kondo et al.
(2015) measured an SMG acinar fraction of `66.91%`.  Even the deliberately
generous transfer of one `100 mg` gland, unit tissue density, and a `1.30 pL`
modeled cell gives only `51,469,231` modeled acinar cells.  The old profiles
therefore fail the absolute-scale gate by factors of about `5.73` and `5.07`;
this is a hard rejection, not a warning.  A paired-gland interpretation cannot
rescue a protocol collected from the ipsilateral gland (Kondo et al. 2019).

Independent P15 AE2-line validation was then applied without changing the
frozen WT parameters.  On the selected `0.12/M2`/preferred-R1 trajectory,
dynamic AE2 deletion changes 10-min cumulative flow by only `+0.0604%`, with
Radau/BDF agreement.  This is qualitatively consistent with the reported
absence of a large AE2 secretion effect, not an equivalence claim: P15
reported no numerical equivalence margin.  The comparison starts from the
same frozen WT root and is not presented as an Ae2-line resting fit.

Finally, the compact H23 check perturbed cell volume, lumen volume, equal cell
Na/Cl amounts, and equal lumen Na/Cl amounts by `+/-0.1%`.  Both then-surviving
selected-root profiles passed with Radau and BDF under the same old frozen
scales; those robustness checks do not reverse supersession.

## 7. Frozen pre-reveal artifacts

The runner wrote:

| Artifact | Contents |
|---|---|
| `wt_trajectories.csv` | Tidy R0--R3 WT states, ions, pH, volume, regulation, and flow. |
| `wt_trajectory_summary.csv` | Endpoints and per-key, unit-labelled conservation maxima. |
| `wt_solver_crosschecks.csv` | Radau/BDF and loose/tight comparisons. |
| `wt_initial_state_sensitivity.csv` | Charge-preserving nearby-state results. |
| `regulatory_fit_summary.csv` | R0--R3 gains and unfit kinetic landmarks. |
| `wt_ion_uptake_diagnostics.csv` | Separate CCh, IPR, and combined onset/preconditioned diagnostics. |
| `wt_flow_scale.json` | WT-only scale, implied cell count, and geometry-gate status. |
| `wt_flow_trajectory.csv` | Minute landmarks on the scaled WT trace. |
| `wt_dynamic_gate.json` | Resting, dynamic, solver, conservation, flow, and regulatory gate ledger. |
| `dynamic_frozen_manifest.json` | Immutable root/parameter/input/tolerance freeze. |
| `dynamic_artifact_hashes.csv` | SHA-256 hashes for all preceding WT-only artifacts. |

`dynamic_artifact_hashes.csv` records hashes for the ten preceding artifacts;
those artifacts describe the superseded baseline failure.  The old ensemble
freeze consists of the following immutable historical artifacts:

- `g5_final_pre_reveal_manifest.json`, SHA-256
  `6a820b35d6f475dfaddb5f66096ef9c75661aba2a4f0ff82f2c9f496b32c30fd`;
- `g5_final_pre_reveal_gate.json`, SHA-256
  `19eb2300ca10da44e09e10138d5440787bde0c4f4a59a121e139d035231b7c26`;
- equation/source-tree SHA-256
  `0d763be527ec1d2f8dad8a5317478ba2a70326f8152989cbba2506bdd08692e1`;
  and
- `g5_final_pre_reveal_hashes.csv`, whose own SHA-256 is
  `4130de9e41490bfdc734006e5b62da93ac8db39e5d7c50f78245c36a5c9dd3ca`.

These hashes preserve the exact old payload, including its then-current
eligibility field.  They do not override the separate controlling artifact
`g5_supersession.json`, which records `REJECTED_SUPERSEDED`, both hard-failure
reason codes, and current `reveal_eligible: false`.  The old files must not be
rewritten to make their historical fields look current.

The executable tests cover input separation, 600-second timing,
coordinate-unit tolerances, charge-preserving nearby states, short Radau/BDF
agreement, WT-only flow conversion, ensemble coverage, and pre-reveal
manifest rejection.

## 8. Current gate

| Gate | Current status | Reason |
|---|---|---|
| Distinct CCh/Ca and IPR/beta inputs | pass | Executable and tested. |
| Physical seconds and 600-s protocol | pass | No inherited code-time scale. |
| R0--R3 interface | pass | All families integrate; R2/R3 exact nesting is reproduced. |
| Old G2/G5 WT resting roots | rejected, preserved | Exact numerical roots, but constructed with a non-native common conductance multiplier. |
| Old G5 Radau/BDF dynamics | historical pass | Complete old 20-member ensemble passed numerical, positivity, conservation, solver, and nearby-state checks before source supersession. |
| WT uptake comparison | partial by design | Model derivative is reported; SPQ observation map is absent. |
| Superseded baseline WT flow shape | fail, preserved | No single scale fits the `0.55 uM`-Ca baseline into the graphical envelope. |
| Old joint WT flow shape | historical within-profile pass | Two profiles passed the old scale-free shape gate; that cannot compensate for their two later hard failures. |
| Absolute gland scale | **fail** | Selected old profiles require 294,706,403 or 260,704,116 cells versus the generous 51,469,231-cell one-SMG bound. |
| Native conductance source map | **fail** | Old `G=10` scaling is not the literal 31.4 nS CaCC plus 14 nS total-K map with the Palk Ca gate. |
| Independent AE2 validation | qualitative pass | No conspicuous modeled secretion effect; source has no equivalence margin. |
| Cell+lumen nearby-state robustness | pass | Both selected-root profiles pass all eight `0.1%` perturbations with both solvers. |
| Regulatory kinetic identification | unidentified, non-blocking ensemble axis | P21 contains no acute regulatory onset trace; all declared sensitivities remain frozen. |
| Old effective-Ca/NKCC pair identification | rejected non-unique profiles | `0.10/M1.75` and `0.12/M2` passed old shape gates but neither is measured directly and both belong to superseded roots. |
| Strict-native H1 resting gate | **fail, complete** | 45 panels/135 attempts/35 numerical roots; zero pass WT Cl/K gates. |
| HW hydraulic diagnostic | **fail, diagnostic-only** | 225 panels/675 attempts/175 numerical roots; zero WT-passing roots. H>1 cannot rescue production. |
| HWQ hydraulic+outflow diagnostic | **fail, diagnostic-only** | 225 panels/675 attempts/175 numerical roots; zero WT-passing roots. At every tested H its displayed intracellular ranges match HW, so coupled outflow scaling does not repair Cl/K. |
| Combined native rest checkpoint | **fail, complete through hydraulics** | H1+HW+HWQ: 495 panels/1,485 attempts/385 numerical roots/zero WT-passing roots or panels. Scale-one HW/HWQ aliases canonicalize exactly to H1. |
| Absolute-source three-start grid | **complete, preliminary** | 630 panels/1,890 attempts/354 numerical roots. Ten WT passes occur only for `N_ABS_NKCC` scale 4: all five topologies at pump 1/`AE4NA05` and pump 2/`AE4NA20`. `AN_ABS_AE4_NHE` has zero passes. |
| Combined native-source hierarchy | **complete through preliminary screen** | 1,125 panels/3,375 attempts/739 numerical roots/10 preliminary WT passes. Every pass is pre-seven-start and pre-33-start. |
| Confirmation and root geometry | **pending; engineering-gated** | A code audit paused promotion and hardened panel-atomic reruns, both-quartile start coverage, fail-closed ID/eligibility/upstream checks, scale-one aliases, and geometry intake. No scientific parameter or WT gate was retuned. |
| Eligible for native WT dynamics | **no** | Seven-start confirmation and 33-start geometry are incomplete; no root ID is dynamics-eligible. |
| Eligible for held-out reveal | **no** | `g5_supersession.json` is controlling until a new native-source pre-reveal freeze passes every WT gate. |

The separate beta-associated volume-sensitive apical anion-exit requirement
also remains unresolved.  Catalan et al. (2015) license the topology but not
an absolute whole-cell conductance.  The diagnostic law
`gV*beta*max(Vi/Vrest-1,0)` stays closed because the tested IPR-only trajectory
shrinks from `1.3000` to `1.273730 pL` and never swells; it is therefore a
deadlocked, production-ineligible rescue rather than a fitted pathway.

The decisive missing dynamic measurement remains a simultaneous acute WT
beta/cAMP-or-PKA/AE4 activity time course under the matched CCh+IPR protocol.
It would identify a preferred kinetic member; until then, the corrected
`10/30/90 s` timing panel remains an explicit ensemble axis.  A simultaneous
absolute calcium and NKCC-activity time course under the matched `0.3 uM` CCh
+ `5 uM` IPR protocol would additionally constrain the source-class hierarchy
without reviving the two rejected old joint profiles.

## Final native-source dynamic result

The pending confirmation/geometry checkpoint described above is now complete.
All ten scale-4 N-class roots retained a single full-rank branch under 33
starts and were independently reproduced. The frozen native dynamic matrix
then completed 800 production and 320 nearby-state cases. Every numerical,
conservation, solver/tolerance, sustainment, co-stimulation, regulatory, and
nearby-state gate passed.

No root passed the absolute one-SMG WT gate. The required shared observation
scale ranged from `6362.3123` to `8210.7187`, while the generous one-SMG
ceiling is `3088.15386`. The final gate is therefore
`SCIENTIFIC_GATE_FAIL`, `phenotype_reveal_permitted=false`. A complete
high-capacity diagnostic stress tier produced zero new WT-rest passers.
