# Task 48 parameter admissibility advisory

Advisory only. The orchestrator must accept all scientific decisions. This worker ran no model inference, sweeps, genotype resting solves or Task 47 reproduction.

## Source scope

The common Task 48 starting SHA is `5c945c4a7857657b5f7c76d5bf9e3d3abc7d5721`. Task 43 was read at `547f113d9123ab1976774639a69483faf40cef34`; Task 44 at `e5fa9840bf96be3147c6118daee4942468c78f8b`. Those branches are intentionally unmerged. Their original provenance has not been promoted to measured salivary parameters. The complete 137 record source inventory and the condensed machine readable classification accompany this note.

## Decisive audit result

**The inherited model already has explicit AE4 activation driven by beta adrenergic input.** Task 47 `common.py` instantiates `R1EffectiveActivation`; `camp_pka.py` implements `dx/dt=(beta-x)/tau` and `G=basal+increment*coupling*x`; the model passes that gain to `equal_routing_adapter`, which passes it to the QSS AE4 evaluator. That evaluator multiplies carrier amount by G. Deleting AE4 removes its transport contribution. No duplicate PKA state is justified.

The inherited values are basal=1, increment=0.25, coupling=1 and tau=30 s. They are effective modelling assumptions. The 2021 primary study gives the approximate 25% forskolin increase in heterologous CHO cells, not a native SMG maximal gain. It supports beta/cAMP/PKA activation, H89 inhibition, a response to constitutively active PKA and loss of that response in S173A; it does not directly demonstrate phosphorylation of S173. Reidentifying the existing activation magnitude using the native isolated exchanger data is admissible if their observation map is valid. Adding activation as a newly missing mechanism is factually incorrect. [Primary 2021 study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8887885/).

## Exact screening set and restrictions

The following 33 shared quantities are provenance admissible for a local full-vector sensitivity screen. This list is not a claim that all are identifiable, nor permission for high dimensional unconstrained optimisation. Each must be removed or combined if the mapped-data Jacobian has zero or collinear sensitivity. No finite physiological bounds are supplied by Task 43 for these effective settings. Positive domains, fractions between 0 and 1 and logarithmic parameterisation are mathematical domains, not confidence intervals.

| Parameter | Inherited value | Class | Initial priority |
| --- | ---: | --- | --- |
| `acid_base.cell_buffer_pka` | 7.0 | unconstrained_effective_parameter | ancillary_shared |
| `geometry.cell_buffer_total_fmol` | 30.0 | unconstrained_effective_parameter | ancillary_shared |
| `geometry.fixed_cell_anion_equivalents_fmol` | 78.06176220594485 | model_determined | ancillary_shared |
| `geometry.cell_impermeant_osmoles_fmol` | 124.69936937124379 | model_determined | ancillary_shared |
| `homeostasis.nhe1_cha_carrier_amount_fmol` | 2.339370005697548e-05 | model_determined | primary_transport |
| `homeostasis.ae2_capacity_fmol_s` | 0.005 | unconstrained_effective_parameter | primary_transport |
| `homeostasis.thermodynamic_saturation_log_width` | 2.0 | unconstrained_effective_parameter | ancillary_shared |
| `homeostasis.co2_basolateral_permeability_fmol_s_mM` | 0.02 | unconstrained_effective_parameter | ancillary_shared |
| `homeostasis.co2_apical_permeability_fmol_s_mM` | 0.01 | unconstrained_effective_parameter | ancillary_shared |
| `membranes.nak_capacity_fmol_s` | 0.08 | unconstrained_effective_parameter | ancillary_shared |
| `membranes.nak_na_half_mM` | 10.0 | unconstrained_effective_parameter | ancillary_shared |
| `membranes.nak_k_half_mM` | 1.5 | unconstrained_effective_parameter | ancillary_shared |
| `membranes.apical_pump_fraction` | 0.075075 | unconstrained_effective_parameter | ancillary_shared |
| `membranes.apical_k_fraction` | 0.3 | unconstrained_effective_parameter | ancillary_shared |
| `membranes.g_k_total_S` | 1.4e-08 | unconstrained_effective_parameter | ancillary_shared |
| `membranes.g_cl_apical_S` | 3.14e-08 | unconstrained_effective_parameter | ancillary_shared |
| `membranes.calcium_half_uM` | 0.26 | unconstrained_effective_parameter | ancillary_shared |
| `membranes.calcium_hill` | 1.46 | unconstrained_effective_parameter | ancillary_shared |
| `membranes.g_para_na_S` | 3e-10 | unconstrained_effective_parameter | ancillary_shared |
| `membranes.g_para_k_S` | 5e-11 | unconstrained_effective_parameter | ancillary_shared |
| `membranes.g_para_cl_S` | 2.5e-10 | unconstrained_effective_parameter | ancillary_shared |
| `membranes.g_para_hco3_S` | 5e-11 | unconstrained_effective_parameter | ancillary_shared |
| `water.lumen_dead_volume_pL` | 0.1 | unconstrained_effective_parameter | ancillary_shared |
| `water.outflow_rate_s` | 0.2 | unconstrained_effective_parameter | ancillary_shared |
| `ae4.carrier_amount_fmol` | 0.1501264265070657 | model_determined | primary_transport |
| `ae4.na_loaded_attempt_rate_s` | 0.1 | model_determined | primary_transport |
| `ae4.k_loaded_attempt_rate_s` | 1.9 | model_determined | primary_transport |
| `nbc.capacity_fmol_s` | 0.11570913197464398 | model_determined | primary_transport |
| `nbc.log_width` | 2.0 | unconstrained_effective_parameter | ancillary_shared |
| `nbc.nhe1_stimulated_multiplier` | 2.3 | unconstrained_effective_parameter | primary_transport |
| `nkcc.alpha_eff_fmol_s` | 0.017334746894096052 | model_determined | primary_transport |
| `regulation.fully_activated_increment` | 0.25 | unconstrained_effective_parameter | primary_transport |
| `nkcc_activation.fully_activated_multiplier` | 1.75 | model_determined | primary_transport |

The primary transport subset is only a priority ordering. Pump and channel settings must not be selected because changing them moves saliva towards 35%. Their admissibility follows from their provenance; infer them only if joint observations independently constrain them. Failure in an arbitrarily reduced numerical subset does not prove that all shared admissible parameterisations fail.

`regulation.tau_ae4_s` is a conditional additional candidate only after isolated native exchanger exposure times and transient information are verified. Without such data, preserve the inherited 30 s as an openly unidentifiable assumption and do not infer delay from saliva. A single timed endpoint generally identifies a gain/time combination, not both.

Published Cha microscopic coefficients, published Palk coefficients and hydraulic conversion values remain fixed here: no source uncertainty ranges have been established. Bath composition, temperature, agonist schedule, inhibitor action and deletion are source prescribed inputs. The historical initial state coordinates are replaced by solved genotype resting states, not estimated separately. Structural zero background channels and the inactive AE4 cooperative gate cannot be turned on under ordinary parameter reidentification. Calcium step values and recruitment thresholds are retained as declared inputs unless independently measured protocol data replace them.

## Identifiability deductions before numerical fitting

1. In the QSS AE4 law, multiply every microscopic attempt rate by a common positive factor lambda and divide carrier amount by lambda. QSS occupancy is unchanged; each turnover scales by lambda; every whole cell AE4 flux is unchanged. The data cannot identify all four separately. Fix the common chloride attempt rate to 1 per second as a gauge; the estimable parameters are an effective carrier/rate scale and two rate ratios. The molecular carrier count remains unidentified.
2. AE4 gain is `carrier*(basal+increment*coupling*x)`. Coupling and increment are exactly multiplicatively redundant. Basal and carrier are also redundant after rescaling the increment. Fix basal=coupling=1 and beta input=1 when IPR is on. Fit an effective native gain only. No independent receptor occupancy is measured.
3. At fully recruited constant stimulation, NKCC effective amount and NKCC stimulated multiplier enter only as their product. NHE carrier amount and NHE stimulated multiplier have the same degeneracy. Task 44 found identical stationary sensitivity columns. Rest and isolated protocols may break these degeneracies; a stimulated stationary dataset alone cannot.
4. The AE4 activation time has exactly zero equilibrium sensitivity. It requires timed data and cannot be identified from resting means or constant stimulus equilibria.
5. Isolated initial slope ratios can remove common fluorescence scale, but they also remove absolute flux scale. Absolute transporter carrier capacities cannot be extracted from such ratios without an independently supported fluorescence/volume map.
6. A common acinus to gland multiplier cancels in genotype saliva or bicarbonate ratios; that multiplier is not identifiable and should not be added to the model.
7. The fixed 0.5/0.5 realised AE4 Na/K routing is a structural assumption, not an experimentally measured fraction. Reidentifying microscopic Na and K attempt rate ratios changes scalar concentration dependence but does not identify this imposed source split.
8. Native isolated control minus AE4 KO uptake is not automatically pure AE4 activity: genotype resting and depleted states can differ. The tempting calculation `(4.33-1.08)/(1.56-0.95)` is a descriptive difference ratio, not a source measured AE4 activation factor.

## Protocol representability limits that parameters cannot fix

The selected Palk NKCC law is `alpha*activity*(A1-A2*Na_i*K_i*Cl_i^2)/(A3+A4*Na_i*K_i*Cl_i^2)`. It has no extracellular concentration arguments. It therefore cannot literally represent chloride bath depletion/readdition or sodium withdrawal with the same source interpretation. Changing alpha cannot fix missing bath dependence. The inherited AE4 symmetric rates also require strictly positive external ions and bicarbonate, so exact ion-free protocols need a source justified boundary treatment rather than an arbitrary fitted epsilon. These are issues for the orchestrator protocol and feasibility decision; this advisory neither accepts nor implements an equation correction.

The Task 46 and 47 results establish no numerical NKCC cap. Task 47 rejected a limiting inherited K recycling block at its frozen point and exposed uncertain capacities; it did not prove no admissible parameterisation can match the full joint data. Those results were read, not reproduced.

## Required subsequent decision

The orchestrator must publish the exact executable observation vector and its local Jacobian before deciding a numerical fit subspace, unless an exact structural contradiction already rules out the full shared parameter class. If essential assays cannot be represented, no fully joint admissible fit should be claimed. A restricted available-data fit is conditional evidence and must be labelled accordingly. There is presently no newly missing AE4 activation mechanism, no verified finite physiological prior on its effective gain or timing, and no source basis for genotype-specific capacity multipliers.

## Exact contrast obstruction reported during independent review

The feasibility worker identifies a stronger shortcut than numerical sensitivity screening: in AE4 KO the only beta driven transport contribution vanishes. With identical CCh input, rest and observation map, CCh only and CCh plus IPR therefore give identical chemical trajectories for every shared admissible parameterisation. Altering the existing AE4 gain or its timing cannot affect this KO contrast. Its parameter Jacobian is identically zero.

The reported KO uptake means are 2.30 ± 0.10 and 0.90 ± 0.09. Under independent Gaussian observation errors, forcing a common predicted mean yields a minimum contribution of 108.287292818 to chi squared, with a contrast of 10.406118047 standard errors. The best common mean would be 1.526519337 in the reported uptake units. These are algebraic data calculations, not a model fit. The protocol/observation workers must verify that all other inputs and the fluorescence normalisation coincide before the orchestrator accepts this argument. If confirmed, a high dimensional fit cannot remove this exact discrepancy and need not be run to establish rejection.


Orchestrator integration: reviewed and accepted for the 48A evidence freeze. This document supplies evidence and restrictions, not a completed fit.
