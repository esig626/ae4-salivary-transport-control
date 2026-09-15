# Independent Task 48 feasibility advisory

Frozen input: `5c945c4a7857657b5f7c76d5bf9e3d3abc7d5721`.

This advisory contains equation inspection, saved output inspection and a few direct evaluations at a valid state. It contains no authoritative inference, parameter fitting, resting solve, trajectory, scientific correction, branch modification, commit or push. The accompanying `check_structure.py` reproduces `structural_checks.json` using the frozen Task 47 parameter loader. A saved WT state is solely a fixture for a pointwise algebraic identity, never a surrogate established knockout rest.

## Principal finding

The inherited architecture predicts exactly the same AE4 knockout core dynamics under CCh alone and CCh plus IPR whenever the initial core state and the remaining protocol history are the same. This identity holds for every shared admissible parameter vector, including any AE4 specific activation gain or time constant. The two Table 1 knockout uptake means strongly disagree with that identity. Genotype specific rest cannot remove an invariance within the same genotype.

The protocol worker independently checked the primary Table 1 caption and the Results description: these measurements use physiological bicarbonate solution, CCh at 0.3 micromolar and IPR at 5 micromolar, with no documented additional bath or inhibitor difference between the CCh and combined comparisons. Initial rates are fitted over indicated regions of separate traces; there is no licence to introduce protocol specific fitted windows or calibration multipliers.

## Exact dependency proof

Write the twelve conserved cell and lumen coordinates as `x` and the effective AE4 regulatory fraction as `r`. For the inherited model:

1. `validation.SecretagogueProtocol.calcium_input` supplies the same calcium under `CCH_ONLY` and `CCH_IPR`.
2. `SecretagogueProtocol.beta_input` differs between the arms.
3. `camp_pka.R1EffectiveActivation` changes `r` and the AE4 capacity multiplier only.
4. `nkcc_stimulation.CompositeAe4NkccRegulation` passes calcium, rather than beta or the AE4 regulatory state, to `N1AlgebraicNkcc1`.
5. `nbc_minimal.MinimalNbcModel` obtains NBC recruitment and the NHE stimulation multiplier from calcium only.
6. `model.ModernFullModel.evaluate` passes `genotype.ae4_expression` to the AE4 contribution; a zero value removes all AE4 sources exactly.

Consequently, for AE4 knockout,

`dx/dt = f(x, calcium(t), common bath and inhibitors; theta)`

contains neither beta nor `r`. On the physical domain where the ODE solution is unique, a common genotype resting state and common CCh/bath/inhibitor history give identical `x(t)` under both arms. This implies identical chloride concentration, volume, pH, secretion and fluorescence under any common causal observation rule. The proof does not equate fluorescence slope to transporter flux.

The direct check used calcium 0.25 micromolar, AE4 regulatory fraction zero for CCh and one for combined stimulation, and exact AE4 deletion. The maximum difference of all twelve core RHS entries was **zero**, bit for bit. NHE and NKCC were also bit identical. This check verifies the implementation dependency; the algebra establishes its parameter independence.

## Quantitative consequence of the identity

The primary data worker verified these source values in Table 1:

| AE4 knockout protocol | Mean initial uptake | Reported SE |
| --- | ---: | ---: |
| CCh, n=6 | 2.30 | 0.10 |
| CCh plus IPR, n=9 | 0.90 | 0.09 |

Units are `10^-3 s^-1` for the SPQ observation, not fmol/s.

If a common observation scale is not otherwise known, allow the common predicted mean to take any value `m`. Even this relaxation leaves

`min_m [(m-2.30)^2/0.10^2 + (m-0.90)^2/0.09^2] = 108.2872928177`.

The minimising common mean is `1.5265193370`; its two standardised residuals are approximately `-7.7348` and `6.9613`. This is an analytic lower bound under independent Gaussian sampling errors for the reported means, not a model fit and not an invented hard tolerance. The difference is `10.4061` combined standard errors under independence.

Unknown covariance cannot erase the contrast: for any correlation in `[-1,1]`, the standard error of the difference is at most `0.10+0.09=0.19`, so the difference is at least **7.3684** standard errors. This uses the reported SEs as mean uncertainty, not biological standard deviation. A formal small sample p value would need the underlying sampling assumptions and is not asserted.

This bound remains valid if the full joint objective contains additional observations with nonnegative squared residual contributions. It also persists when shared parameters, genotype resting states or an AE4 only PKA activation law are reidentified. It depends on a common observation operator for the two documented comparable protocols. Allowing independent arbitrary stimulus observation scales would remove the test by adding unsupported freedom.

## AE4 PKA activation is already present

`analysis/47_dynamic_potassium_recycling_reconstruction/common.py:load_model` reconstructs `R1EffectiveActivation` within `CompositeAe4NkccRegulation`. The frozen gain increment is 0.25, the basal multiplier is one, and the activation time is 30 seconds. `camp_pka.py` explicitly labels acute kinetics as unmeasured and gain choices as sensitivity assumptions.

Therefore Task 48 must not claim that beta/cAMP/PKA activation is absent, add a duplicate state, or misread the existing class name `R1` as proof that Task 48's R1 reconstruction has already been performed. Reidentifying the existing gain from isolated exchanger data is an inference operation; it does not repair the knockout invariance. Isolated endpoint measurements do not uniquely identify the activation time constant.

## Additional protocol representability limits

| Issue | Inherited code | Consequence |
| --- | --- | --- |
| Exact zero Na or total carbon bath | `parameters.FullModelParameters.__post_init__` rejects each | Na removal and bicarbonate free assays cannot be passed through unchanged |
| Positive activities assumed throughout | `transporters.AE4Environment`; `membranes.evaluate_homeostasis`; `nbc_minimal` electrical closure | Genotype/inhibitor masks do not universally bypass invalid zero activities |
| Carbonic anhydrase inhibition | `acid_base.speciate` is instantaneous equilibrium and has no CA kinetics | Ethoxyzolamide has no represented target; changing a transporter cannot impersonate its action |
| Finite bath NKCC specificity | `nkcc1_palk2010.palk_cycle_flux_fmol_s` accepts only intracellular Na, K and Cl | Active Palk law does not respond directly to the source's 4 mM depletion versus high chloride readdition |
| Ohmic electrodiffusion at zero bath ions | `membranes.nernst_voltage_V` and paracellular reversal terms | Simply accepting zero concentration would encounter infinite Nernst potentials; a finite boundary flux extension is not merely parameter validation |

The NKCC audit held the intracellular state fixed and changed extracellular chloride from the inherited 126.1616 mM bath to the source's 4 mM depletion concentration. NKCC remained exactly `0.1283531881` cycles fmol/s. This is a dependency test, not a physiological assay simulation. The active reduced law supplies no unique bath varying extension. Its unchanged inward transport rate across the large bath change is an audit limitation even though no zero chloride concentration is used in the actual source assay.

The script also demonstrates a zero chloride validator rejection, but the source protocol uses 4 mM chloride and that rejection is not claimed as an exact assay blocker. Source high chloride differs by solution type: 128.3 mM in bicarbonate solution and 153.3 mM without bicarbonate. External charge and osmolarity bookkeeping also require the documented replacement ions; changing chloride alone is not a complete bath reconstruction.

An inherited protocol cannot be called exact by replacing missing ions with a small arbitrary floor, silently suppressing other transporters, or treating ethoxyzolamide as deletion of bicarbonate transport. Analytic zero limits, if finite and unique, are numerical extensions; a new extracellular substrate response or finite carbonic anhydrase kinetics changes equations and must be declared. None is authorised here by this worker.

## AE4 constitutive Na response has the opposite fixed state direction

The active `AE4Parameters` has `cooperative_gate=None`. Let `a` and `c` denote the chloride inward and outward edge rates; let `D` be the sum of Na and K outward loaded rates and `U` the corresponding inward sum. The inherited QSS total chloride flux is

`J4 = C (a D - c U)/(a+c+D+U)`.

It has derivatives

`dJ4/dD = C (a+c)(a+U)/(a+c+D+U)^2 > 0`,

`dJ4/dU = -C (a+c)(c+D)/(a+c+D+U)^2 < 0`.

At fixed other concentrations the Na outward rate is proportional to `Na_o^(-1/2)` and the inward rate to `Na_o^(1/2)`. Thus the ungated inherited law gives `dJ4/dNa_o < 0`: removing extracellular Na increases the predicted chloride uptake. Its Na zero limit for total chloride is finite, `C*a`, although individual rate expressions diverge and the API rejects zero.

Direct evaluations at the saved valid fixture give J4 approximately 0.00493393, 0.01632837, 0.03501276 and 0.20259652 fmol/s as external Na decreases from 145 to 49 to 1 to `10^-6` mM. No parameter was changed. The native AE2 knockout assay instead reports lower uptake after Na removal. This is a constitutive direction conflict at fixed intracellular state. It is not, on its own, a proof about the complete assay because Na removal may alter intracellular pH, ions and volume before measurement. It should not be promoted into a universal dynamical impossibility statement without reconstructing that history.

## Resting states and AE2 control

At resting calcium, the inherited NBC recruitment is identically zero. The cell alkalinity balance then imposes at stationarity

`NHE = AE2 + 2 AE4`.

For AE4 knockout this becomes `NHE = AE2`. In the saved WT rest, the fluxes are NHE `0.007880876`, AE2 `-0.001986984` and AE4 `0.004933930` fmol/s. The inherited AE2 capacity is only 0.005 fmol/s, so unchanged knockout rest cannot retain that WT NHE flux. Its pH, Na or other state variables must change. This is a useful explanation for a resting pH residual, not a proof against a new shared admissible AE2 capacity.

The saved WT rest has chloride `60.5723 mM`, already about `6.98` reported SE above the measured `50.10 mM` mean. Its pH `6.91098` agrees with the WT central pH. The historical physiological gate is therefore not equivalent to the new data likelihood.

Use the ten coordinate charge manifold in `task30_nhe1_repair.state_from_coordinates` and ten independent conserved rows for new genotype resting solves. Keep the geometry, fixed charge and transporter capacities common. Do not run the old WT calibration objective separately on each genotype or change fixed cell charge to make each desired chloride concentration a root. Its broad numerical coordinate bounds are not measured parameter uncertainty.

AE2's large isolated exchanger contribution and small whole gland knockout secretion effect are not an algebraic contradiction by themselves: physiological compensation differs from the isolated transporter assay. Those constraints must be compared at their actual protocols. Likewise, unchanged NHE mediated pH slope is not an assertion that a saved physiological NHE flux must be exactly equal across genotypes.

## Timing and admissible decision

The inherited calcium input is a step; AE4 has a nominal 30 second regulatory time, while ionic and water modes can be slower. A lower resting chloride state does not by itself prove an immediate large saliva deficit, since fluorescence normalisation, ion amounts, electrical driving and volume all matter. The required initial two to three minute similarity and later sustained deficit need explicit common minute bins or a source justified qualitative statement; no arbitrary percentage band is justified.

The exact knockout beta invariance is the central structural rejection. Changing only AE4 gain or kinetics cannot repair it, even before considering secretion. A further beta effect outside AE4 is needed to represent the comparison, but the observation alone does not identify which pathway or how it acts. The primary worker found speculative discussion of distinct calcium and cAMP effects on AE4/AE2, not an independently measured unique beta inhibitory route that fixes this knockout contrast. Without uniquely identifying independent evidence, selecting a beta target would be another mechanism search. The orchestrator should therefore distinguish an analytical failure of the inherited family from unresolved assay mappings and from fixed parameter residuals, and publish the negative result without claiming that all physiological models have been excluded.

Primary source for the independently verified data and protocol: [Peña-Münzenmayer et al. 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4409235/), Table 1 and Results describing Figure 1.

Task 43 and Task 44 directories are absent from the frozen checkout. This worker inspected their retained Task 46 replay and Task 47 provenance/handoff summaries and did not fetch or modify any refs. Any full provenance audit requiring their original reports remains with the orchestrator/provenance worker.
