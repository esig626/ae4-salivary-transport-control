# Round-6 dynamic reconstruction

## Controlling result

There is no admissible Task-13 WT trajectory from which to form a predictive
AE4-null secretion ratio. All five conditional state-cycle survivors fail the
unchanged-chassis WT calibration/physiology gate, and all five nested Round-5
cation-topology attempts also fail before dynamics. The correct Round-6 result
is therefore **blocked/N/A**, not a manufactured endpoint ratio.

This is a model-gate result, not a statement that AE4 biology is unimportant.
It says only that no tested state-resolved AE4 mechanism has yet been embedded
in an independently calibrated, WT-valid whole-cell chassis.

The strict knockout secretion magnitude and time course were never passed to
an optimizer, root finder, solver setup, parameter bound, model selection, or
trajectory denominator. A candidate-independent AE4-null code-coordinate
trajectory was integrated as a numerical diagnostic; it is explicitly not a
held-out prediction.

## 1. Prerequisite gate ledger

| Gate | Required before a held-out dynamic claim | Frozen result | Consequence |
| --- | --- | --- | --- |
| D0: target firewall | No `E15-02` or `E15-03` value used in calibration | pass | strict targets remain unfit |
| D1: transporter compatibility | Conditional transporter-assay survivor | five rows | necessary, not sufficient |
| D2: interior WT capacity | Exact WT closure with an admissible calibrated capacity | `0/5` Stage-A rows; `0/5` Round-5 nested rows | fail |
| D3: WT physiology | Valid WT state including resting Cl/pH gates | `0` Gate-2 passes | fail |
| D4: independently calibrated extension | Frozen pump/K or acid-base extension calibrated without knockout secretion | none | fail |
| D5: source-constrained dynamics | Physical time map and independently constrained Ca/PKA inputs | absent | fail |
| D6: delayed-state mechanism | Explicit kinetic state capable of a sourced delay | AE4 occupancies are QSS; no PKA time constant | fail for an AE4/PKA delay claim |
| D7: dual-solver trajectory | Positive, conserving WT/KO/AE2 runs with a valid denominator | KO-only Radau/BDF diagnostic passes; predictive pair N/A | N/A |
| D8: model freeze before reveal | Freeze an eligible model before joining strict targets | no eligible model exists | N/A |
| D9: held-out scoring | Endpoint, integral, onset, and late-shape evaluation | no eligible row | N/A |

Failure at D2 is already decisive. Later failures are reported because they
identify what would remain even if static closure were repaired.

## 2. Stimulus and time audit

The experiment and implementation do not yet share a calibrated time axis or
input map.

| Item | Primary experiment | Current implementation | Status |
| --- | --- | --- | --- |
| Stimulation | `0.3 uM` CCh plus `5 uM` IPR for 10 min | Ca changes `0.05 -> 0.55` at `t_code=100`; PKA changes `0 -> 1` at the same code coordinate for SR6 | historical comparison input, not a dose-response map |
| Physical time | flow recorded each minute | `rhs_scale=10000`; code-time conversion unresolved | cannot map 2--3 min to code time |
| PKA onset | no activation time constant reported | immediate sustained step | sensitivity only |
| AE4 state dynamics | molecular relaxation not measured | stationary CTMC occupancy is re-solved on every whole-cell RHS call | state competition without carrier memory |
| Absolute flow | about `9--10 uL/min` WT in the matched gland protocol | chassis-native water/flow unit | absolute comparison uncertified |

The immediate PKA step cannot be relabelled as the experimental delayed
divergence. Fitting a new delay to `E15-03` would violate the holdout.

## 3. Candidate delay mechanisms

| Proposed origin of divergence | What the current equations contain | Round-6 verdict |
| --- | --- | --- |
| delayed PKA activation of AE4 | edge-specific PKA regulation, but an immediate input and QSS occupancy | not testable without independently measured PKA/AE4 kinetics |
| intracellular Cl depletion | dynamic `Cl_i` exists | code-coordinate sensitivity is possible, but an invalid WT chassis and absent physical-time map prevent phenotype scoring |
| changing Na/K gradients | dynamic `Na_i` and `K_i` exist | same limitation; native stimulated trajectories were not measured |
| pump/K-channel recruitment | Round 5 redistributes fixed total capacities | no recruitment state or time constant; fractions are unmeasured SMG sensitivities |
| acid-base adaptation | algebraic proton/CO2 reduction and dynamic bicarbonate | no conserved carbon pool or independently calibrated buffer/NHE dynamics |
| lumen/osmotic feedback | lumen Na/K, height, and historical water closure | retained in the numerical diagnostic, but flow and time scales are uncertified |

Thus a code trace can show relaxation of the historical equations but cannot
decide which physical process caused a 2--3 minute delay.

## 4. Round-5 intake

Round 5 tested the smallest conditionally transporter-compatible SR2 `1:2:3`
row with nested cation-current topologies. Pump/K fractions are declared
sensitivities, not measured submandibular fractions. WT pH remained a gate and
was not included in the WT-Cl capacity objective.

| Topology | Attempted apical pump/K fraction | Height (um) | Max raw closure | pH | Valid WT capacity/Gate 2 |
| --- | ---: | ---: | ---: | ---: | --- |
| M0 control | `0/0` | `199.999995` | `3.05048e-5` | `7.52055` | no/no |
| pump only | `0.2/0` | `199.999989` | `3.86550e-5` | `7.52127` | no/no |
| K only | `0/0.2` | `199.999987` | `3.14718e-5` | `7.52022` | no/no |
| coupled | `0.2/0.2` | `199.999995` | `4.37303e-5` | `7.52293` | no/no |
| coupled | `0.4/0.4` | `199.999984` | `5.40240e-5` | `7.53393` | no/no |

All attempts hit the derived-height boundary and miss the `1e-8` closure
requirement. Round-5 null roots and target-manifold evaluations remain useful
ionic sensitivities, but none supplies a WT denominator. The structural
acid-base audit has rank three in its declared full flux space, rank one in
the released Cl/pH observation, and flux nullity two. Consequently the
released pair cannot calibrate a conserved-carbon/NHE/buffer module either.

These results are machine-readable in
`results/13_state_resolved_ae4/chassis_reconstruction_scan.csv`.

## 5. Executable code-coordinate diagnostic

The only integration performed for this round deliberately avoids a failed WT
denominator:

- model: `FixedChassis(ZeroAE4())`;
- scenario: candidate-independent AE4 knockout;
- initial state: the frozen common bounded null root;
- input: historical Ca step at `t_code=100`, with PKA fixed at zero;
- grid: 452 points from `stage_a_time_grid()`;
- solvers: Radau and BDF, each with `rtol=1e-8`, `atol=1e-10`.

| Diagnostic | Radau | BDF | Interpretation |
| --- | ---: | ---: | --- |
| endpoint `q_total` | `0.001358119887856` | `0.001358119887848` | chassis-native flow, not a secretion prediction |
| integral over `t_code>100` | `0.1369138922591` | `0.1369138922229` | flow x code time, not a 10-min total |
| endpoint max raw RHS | `1.57e-12` | `9.27e-12` | settled numerically at the code endpoint |
| max conservation check | `4.26e-14` | `4.26e-14` | numerical conservation |

The maximum state-coordinate difference over the shared grid is `2.36e-6`;
because the seven coordinates have mixed meanings, that number is a numerical
crosscheck rather than a physical norm. The maximum `q_total` difference is
`1.28e-9`. All states remain positive.

The diagnostic establishes solver applicability to this KO-only IVP. It does
**not** establish:

- an admissible KO/WT endpoint ratio;
- a 10-minute cumulative ratio;
- first-2--3-minute comparability;
- a sustained physical late deficit;
- a valid absolute WT flow; or
- a mechanism for the observed delay.

Landmarks from both solvers are in
`results/13_state_resolved_ae4/dynamic_diagnostic_trajectory.csv`. The gate and
metric disposition is in
`results/13_state_resolved_ae4/dynamic_trajectory_status.csv`. The executable
crosscheck is `tests/test_state_resolved_ae4_dynamic.py`.

## 6. Specific information needed to reopen dynamics

The next admissible dynamic reconstruction needs one matched experiment, not
a free delay fit. In WT mouse submandibular acini under the exact `0.3 uM CCh
+ 5 uM IPR` protocol, measure a calibrated cAMP/PKA activity time course and
AE4-associated uptake together with time-resolved intracellular Na, K, Cl,
pH, and volume, plus luminal or collected-effluent Na and K. In parallel,
estimate apical and basolateral pump and Ca-activated K currents in the same
preparation and retain minute-resolved gland flow for scale mapping. To fully
separate cation-current feedback from acid-base compensation, co-measure total
inorganic carbon and an NHE- or buffer-flux readout. Without that carbon/flux
measurement, the cation/current protocol is the first decisive branch test,
not a complete mechanism identifier. Freeze all WT-derived fractions,
capacities, time constants, and code-to-physical-unit maps before running the
AE4-null flow prediction.

Until those static and dynamic calibrations exist, the rigorous disposition
is **phenotype unresolved: WT denominator and stimulus-dependent
cation/acid-base time scale specifically unidentified**.
