# 49G measured resting-state balance and source geometry

## Decision

The complete fixed source basis spans the conditional measured-state residual. It does **not** identify a unique shared correction. At the declared AE4-KO state, an exact alkalinity separating functional proves that the required source is outside the native, affinity-oriented addition cone. Increasing basal inward NBC fails the same sign test. No model correction is accepted, so 49H and 49I are inapplicable.

This is a conditional equation-level diagnosis, not a proof that every hidden state compatible with the measured Cl/pH is impossible. The distinction is essential: the experiment does not supply a complete conserved state. No secretion observation, genotype-specific multiplier, pathway subset or parameter fit is used.

## States and measurements

The analysis began after the remotely verified 49E commit `d6e0404a9b87279b507b998edea01ed0820e2b84`. It uses the four frozen Task48 resting cohorts, with no strain-specific model parameters:

| Cohort | Cl mean ± SE, mM | pH mean ± SE | Conditional cell volume, pL |
|---|---:|---:|---:|
| AE4 WT | 50.10 ± 1.50 | 6.91 ± 0.07 | 1.175491 |
| AE4 KO | 36.50 ± 1.60 | 6.89 ± 0.02 | 1.041972 |
| AE2 control | 53.40 ± 1.80 | 6.87 ± 0.01 | 1.200314 |
| AE2 KO | 54.50 ± 1.80 | 6.95 ± 0.05 | 1.241651 |

The predeclared algebraic lift retains the corresponding cached genotype's intracellular Na concentration and complete lumen. Both control cohorts use the same WT reference and laws. Cell CO2 follows the necessary stationary identity `C=H`, where C is total CO2 influx and H is NHE turnover. Carbonate equilibrium supplies TIC at measured pH; cell water stationarity and charge/osmotic accounting supply volume and K. Both compartment charges, measured Cl/pH, `TIC_dot-TA_dot=0`, and cell water balance hold numerically. Lumen water also remains balanced because its saved state is retained. These are conditional states, **not solved resting equilibria**.

All retained and derived values and cached-state hashes are in `output/resting_geometry.json`. Positivity is checked without clipping. Nothing is taken from a WT root and labelled an AE4-null equilibrium.

For measured chloride c and pH h, define carbonate fractions α₀, α₁, α₂, a=α₁+2α₂, buffer deprotonation b and w=[OH]−[H]. The observation/charge manifold is

`nCl_i=c V_i`,

`nTA_i=a nTIC_i+B_fixed b+w V_i`,

`nNa_i+nK_i=(c+w)V_i+a nTIC_i+B_fixed b+X_fixed`.

These leave three cell and five lumen degrees of freedom. The two imposed stationary identities eliminate two of them; six remain fixed to cached nuisance values in this calculation. They are not experimental measurements. The two distinct observed control means are not asserted to be two exact equilibria of the same WT equations.

## Full conserved residual

Positive entries below are inherited net sources at the conditional states, in fmol/s. The required correction is the **negative** of each entry. These residuals are not inferred individual transporter rates.

| Cohort | Na_i | K_i | Cl_i | TIC_i | TA_i |
|---|---:|---:|---:|---:|---:|
| AE4 WT | 0.100564 | 0.129989 | 0.270357 | −0.039804 | −0.039804 |
| AE4 KO | 0.458414 | 0.545235 | 0.994255 | 0.009395 | 0.009395 |
| AE2 control | 0.062123 | 0.077822 | 0.134758 | 0.005186 | 0.005186 |
| AE2 KO | 0.021766 | 0.035409 | 0.127586 | −0.070410 | −0.070410 |

| Cohort | Na_l | K_l | Cl_l | TIC_l | TA_l |
|---|---:|---:|---:|---:|---:|
| AE4 WT | −0.006264 | −0.015685 | −0.022993 | 0.001448 | 0.001044 |
| AE4 KO | −0.021040 | −0.052684 | −0.077231 | 0.000859 | 0.003507 |
| AE2 control | −0.003624 | −0.009075 | −0.013304 | 0.000396 | 0.000604 |
| AE2 KO | −0.002382 | −0.005965 | −0.008745 | 0.001368 | 0.000397 |

Cell/lumen volume residuals are at most 2.0×10⁻¹³ pL/s. The charge-rate and current audit rows are at most 4.7×10⁻¹⁶ fmol-charge/s. Their near-zero values do not remove the nonzero chemical balances. `output/resting_residuals.csv` contains all 64 rows, required corrections and measurement-error propagation.

## One fixed source basis and its projection

The predeclared 21 signatures per cohort include NKCC1, NHE1, AE2, AE4, NBC, both pump and K locations, separate TMEM16A and VRAC chloride columns, both CO2 interfaces, four paracellular ion paths, three water paths and advective lumen outflow. Deleted transporters have zero genotype columns. Dormant NBC/VRAC chemistry is retained in the span but its native availability is recorded as zero. Water laws, pump partitions and other capacities remain fixed.

The 64×84 matrix stacks the four cohorts' 12 conserved rows, two charge rows and two membrane-current rows. Charge/current rows are algebraic audits, not additional independent observations. The block construction permits arbitrary flux coefficients only as a relaxed source-space calculation; it does not fit independent genotype multipliers or claim those coefficients satisfy a shared constitutive law.

For the SVD metric, water rows are multiplied by the bath osmolarity (320.884441 mOsm/L); other rows have unit scale. Nonzero columns are then normalised. The saved matrix has **rank 48, nullity 36**, tolerance 4.04×10⁻¹⁴. The smallest retained singular value is 0.352928; discarded values are below 2.76×10⁻¹⁶. The target norm is 1.290687 in this declared metric; its orthogonal residual norm is **1.36×10⁻¹⁵**. Native source decomposition reproduces the inherited RHS within 1.12×10⁻¹⁶.

There is therefore **no orthogonal missing conserved-coordinate direction** in this full basis. This is expected because `F=Sj` and `−F=S(−j)`. Cancelling every native flux or taking the minimum-norm SVD coefficient representative is not an admissible physiological correction. That representative minimises the norm in column-normalised coefficient coordinates, not the unweighted raw flux-coefficient norm. The full nullspace and projection are saved in `output/source_svd_projection.json`; no constrained least-squares solve is needed.

Two exact dependencies explain why source coefficients cannot identify a transporter:

`s_AE4 = 2 s_AE2 − ½ s_NKCC1`,

`s_NBC = 2 s_NHE1 + 2 s_CO2,b + ⅓ s_pump,b + ⅔ s_K,b`.

The second includes the same electrical-current signature. TMEM16A and VRAC also have identical chloride chemistry. These are identities in the full basis, not selected subsets or proposed replacements.

## Exact constitutive obstruction

For all genotypes the unchanged equations give

`TA_dot = H + 2B − E − 2A`,

`TIC_dot = C + 2B − E − 2A`,

where B is inward NBC, E is inward-Cl AE2 and A is inward-Cl AE4. Consequently any rest requires `C=H`. The inherited resting NBC gate is exactly zero, making the additional conditions

`H=E+2A` in WT, `H=E` in AE4 KO, and `H=2A` in AE2 KO.

The AE4-KO lift evaluates to

`H = +0.0087632634`, `E = −0.0006315750`, `TA_dot = +0.0093948384 fmol/s`.

NHE adds alkalinity, and reverse AE2 also adds it. Every other natively available signature has zero direct cell-TA source in this knockout. If ℓ selects its cell-TA row and C_native is the same source matrix oriented by the actual native flux signs, then

`ℓᵀ C_native ≥ 0`, but `ℓᵀ(−F) = −0.0093948384 < 0`.

This is an explicit Farkas certificate: the target is outside that addition cone, with at least 0.0093948384 separation in the scaled metric. The distance is a lower bound, not a computed nearest-cone projection. Current closure and shared-law constraints can only shrink the cone, so they cannot evade this certificate at the same states and native directions.

Allowing reductions in existing rates does not rescue a complete positive-capacity reweighting at this fixed chemical state: total TA balance would read `a_H H + a_E |E|=0`, with nonnegative capacities. It can hold only if both nonzero rates are suppressed to zero. That degenerate state is not a reconstruction of the inherited nonzero laws. This claim does not cover changes that move the unmeasured state or reverse affinities.

The same obstruction has an explicit scalar equation. Set

`D=P_CO2,b+P_CO2,a`,

`Cbar=(P_CO2,b CO2_b+P_CO2,a CO2_l)/D`,

`Ccrit=(Cl_i HCO3_o/Cl_o) 10^(pKa1−pH_i)`.

For the native AE4 knockout, carbon balance fixes `CO2_i=Cbar−H/D` and AE2 balance demands

`H = G2 tanh(log((Cbar−H/D)/Ccrit)/w_AE2)`.

Necessary conditions include `0<H<min(G2, D(Cbar−Ccrit))`. Here `Cbar=1.186858 mM`, `Ccrit=1.153470 mM`, `G2=0.005 fmol/s`, and `D(Cbar−Ccrit)=0.001001641 fmol/s`, whereas H is 0.008763263 fmol/s. The frozen NHE turnover exceeds even AE2's absolute capacity, independently of the retained lumen; at this actual lift AE2 also has the wrong direction. Increasing AE2's capacity at the same state cannot reverse its sign. The result does not select an NHE replacement or any excluded mechanism.

## Why resting NBC is not the correction

Only after the negative cell-TA direction was identified was the NBC affinity audited. Its dimensionless inward affinity at the unchanged AE4-KO voltage is **+2.7330**. A nonnegative basal capacity would therefore initially supply inward NBC and add alkalinity, opposing the required correction. No basal-NBC trajectory or altered electrical closure was run.

The conditional NBC cycle required to cancel just the TA row would be `(E+2A−H)/2`: +0.0199020 in AE4 WT, **−0.0046974 in AE4 KO**, −0.0025931 in AE2 control and +0.0352052 fmol/s in AE2 KO. The unchanged native voltage gives positive inward NBC affinity in all four cohorts. These signed demands are a diagnostic of incompatibility, not a four-genotype fit or an instruction to force outward NBC. Even a hypothetical TA-only repair would leave the other residual coordinates to satisfy.

Task36's zero-rest NBC choice remains a nesting assumption, not experimental evidence of absence. However, this calculation does not support simply switching it on. NBC adds no new independent source direction and its local inward sign fails the knockout requirement.

## Uncertainty and identifiability limit

Local centred differences propagate the eight measured Cl/pH SEs without rebuilding the source matrix or fitting anything. For the AE4-KO TA residual, derivatives are 6.7400×10⁻⁵ fmol/s per mM and −0.0561551 fmol/s per pH unit. The independent-error working SE is **0.00112827 fmol/s**. With unknown Cl–pH correlation, the linearised SE is bounded above by **0.00123094 fmol/s**. These are conditional measurement-error summaries, not p-values, hard physiological bounds, or proof of a globally persistent sign.

The corresponding TA working SEs for AE4 WT, AE2 control and AE2 KO are 0.0679021, 0.0102152 and 0.0457186 fmol/s. Their central signs must not be described as certain. The full propagated covariance and the six retained-nuisance derivatives per cohort are in `output/resting_uncertainty.json`. A derivative's units are its raw residual-row unit divided by its coordinate unit; covariance entries carry the product of the two row units. Roundoff-level derivatives/SEs of imposed water, charge and current identities are numerical noise, not biological uncertainty. No covariance or prior is invented for unmeasured Na or lumen chemistry. For example, at AE4-KO the TA derivatives with respect to retained intracellular Na and luminal TIC are −0.0003028 and −0.0009153 fmol/s per mM. Those observation-invisible coordinates can affect the obstruction; a local derivative is not a licence for extrapolation or adjustment.

There is no uniquely identified shared correction after considering source nullspace, flux signs and unmeasured states. The precise unresolved requirement for the next task is the **established-genotype carbon/alkalinity compatibility equation above**, together with the full conserved residual. Any proposed correction must supply the missing negative KO TA/TIC rates with charge-consistent companion sources at independently justified resting states. It must not be selected by a coefficient norm, basal-NBC availability or the secretion phenotype. The present data do not authorise naming one transporter or editing one law. The requested stronger outcome of a uniquely implementable correction has not been achieved.

## Execution record

One predeclared source basis and one set of central states were used. A first execution reached export after the SVD and failed because a diagnostic dictionary key was written `ae2` instead of `AE2`. The identical calculation was replayed once to recover unsaved results; both SVD invocations are disclosed in `execution_provenance`. No basis, parameter, state-selection rule or scientific gate changed. This is a deviation from a literal single numerical invocation, not a second basis or a search.

The accepted execution has four central and 64 local derivative RHS evaluations; the failed export attempt had 20 evaluations. There were zero trajectories, zero constrained least-squares/parameter fits, zero pathway subsets and zero model edits. Previously published work was not recomputed. Source/hash verification and advisory checks accompany the checkpoint.
