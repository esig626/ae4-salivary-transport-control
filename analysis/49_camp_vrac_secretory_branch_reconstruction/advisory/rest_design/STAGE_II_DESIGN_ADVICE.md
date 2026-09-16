# Advisory Stage-II equation and identifiability audit

Status: source-code/algebra advisory only. No model RHS evaluation, measured-state residual, rank/SVD, constrained solve or trajectory has been run. No canonical files were modified. The orchestrator must accept or reject this design and execute Stage II only after the verified remote 49E checkpoint.

Inspected: `AGENTS.md`; Task 49 entry, exclusion ledger, source seed, waterfall and 49a addendum; Task 48 protocol/observation/report/data and saved genotype rests; active `protocol_layer.py`, `model.py`, `membranes.py`, `water.py`, `acid_base.py`, `ae4_equal_cation_routing.py`, `nbc_minimal.py`; Task 36 design.

## Critical distinctions

1. The measurements are four pairs of intracellular Cl and pH, not four complete 12-coordinate cell/lumen states. A numerical missing-source vector is conditional on an explicit embedding of all unmeasured coordinates.
2. The complete stoichiometric source basis necessarily represents the inherited vector field: `F(x)=S(x)j(x)`. Hence the unconstrained claim that the inherited residual is in the span of all inherited pathways is tautological. The substantive questions are the thermodynamically admissible cone, genotype sharing of a constitutive law, and uniqueness modulo null directions and unmeasured-state uncertainty.
3. A minimum-norm or NNLS coefficient vector is one representative, not evidence of unique transporter identification. Analytic null identities below exist before any numerical calculation.
4. Multiplicative flux-law sensitivities and stoichiometric signatures differ. NBC has a nonzero physical signature but an exactly zero resting capacity derivative because its existing gate is zero. Making it nonzero is a model edit, not a numerical basis convention.

## Coordinates and full source dictionary

Use the model's exact 12-coordinate order:

`(Na_i,K_i,Cl_i,TIC_i,TA_i,V_i, Na_l,K_l,Cl_l,TIC_l,TA_l,V_l)`.

Let `e_X` denote the unit vector of that coordinate. Positive directions and signatures are:

| Column | Positive direction | Full signature |
|---|---|---|
| NKCC1 | bath to cell cycle | `e_Na_i+e_K_i+2e_Cl_i` |
| NHE1 | Na influx, H efflux | `e_Na_i+e_TA_i` |
| AE2 | Cl influx, bicarbonate efflux | `e_Cl_i-e_TIC_i-e_TA_i` |
| AE4 | Cl influx, fixed equal Na/K export | `-0.5e_Na_i-0.5e_K_i+e_Cl_i-2e_TIC_i-2e_TA_i` |
| NBC | inward 1Na:2HCO3 cycle | `e_Na_i+2e_TIC_i+2e_TA_i` |
| pump apical | cell Na to lumen | `-3e_Na_i+2e_K_i+3e_Na_l-2e_K_l` |
| pump basolateral | cell Na to bath | `-3e_Na_i+2e_K_i` |
| K apical | cell to lumen | `-e_K_i+e_K_l` |
| K basolateral | cell to bath | `-e_K_i` |
| TMEM16A Cl | cell to lumen | `-e_Cl_i+e_Cl_l` |
| Stage-I VRAC Cl, if represented | cell to lumen | `-e_Cl_i+e_Cl_l` (identical chemical signature; resting activity governed by the frozen Stage-I gate) |
| CO2 basolateral | bath to cell | `e_TIC_i` |
| CO2 apical | lumen to cell | `e_TIC_i-e_TIC_l` |
| paracellular Na | lumen to bath | `-e_Na_l` |
| paracellular K | lumen to bath | `-e_K_l` |
| paracellular Cl | lumen to bath | `-e_Cl_l` |
| paracellular bicarbonate | lumen to bath | `-e_TIC_l-e_TA_l` |
| water basolateral | bath to cell | `e_V_i` |
| water apical | cell to lumen | `-e_V_i+e_V_l` |
| water paracellular | bath to lumen | `e_V_l` |
| lumen outflow | lumen to sink | `-C_Na_l e_Na_l-C_K_l e_K_l-C_Cl_l e_Cl_l-C_TIC_l e_TIC_l-C_TA_l e_TA_l-e_V_l` |

Solute coefficients are fmol/s; water and outflow coefficients are pL/s. The outflow signature is state dependent and includes advective solute export. Do not add a bare water-only outflow column and omit its solutes. Finite buffer sites and fixed anions have no source columns: they are constant amounts in this model. Do not add the AE4 regulatory state to the conserved source analysis; it is exactly zero at rest and is not a chemical conservation coordinate.

Keep the two apical chloride columns visibly duplicated if both are represented; dropping one silently hides a structural inability of resting chemical data to distinguish conductances. KO masks zero the deleted transporter's *constitutive flux*, not its stoichiometric chemistry.

## Exact charge/current dependencies

The bulk charge functionals are

`q_i = Na_i+K_i-Cl_i-TA_i-X_fixed`,

`q_l = Na_l+K_l-Cl_l-TA_l`.

Therefore charge-rate rows are exact linear combinations of amount-rate rows. They are useful audit rows, not extra independent observations or rank. Current closure imposes

`I_a-I_para=0`, `I_b+I_para=0`.

In fmol-charge/s units (multiply by `F*1e-15` for A), the relevant current coefficients are:

| Flux | `I_a` | `I_b` | `I_para` |
|---|---:|---:|---:|
| NBC inward | 0 | +1 | 0 |
| apical pump | +1 | 0 | 0 |
| basolateral pump | 0 | +1 | 0 |
| apical K outward | +1 | 0 | 0 |
| basolateral K outward | 0 | +1 | 0 |
| apical Cl outward | -1 | 0 | 0 |
| paracellular Na/K outward | 0 | 0 | +1 |
| paracellular Cl/HCO3 outward | 0 | 0 | -1 |

All other columns have zero current. Both outward pump current and inward NBC current are positive; inward NBC imports one net negative charge. If flux coefficients are manipulated independently, impose both current equalities. Prefer calculating the response of the exact electrical closure when discussing realizable conductance/capacity changes. A single-channel source added at fixed voltage generally violates that closure.

## Algebraic degeneracies that should be disclosed

The active equal-routing AE4 signature obeys exactly

`s_AE4 = 2 s_AE2 - (1/2) s_NKCC1`.

The NBC signature obeys exactly

`s_NBC = 2 s_NHE1 + 2 s_CO2,b + (1/3) s_pump,b + (2/3) s_K,b`.

The second identity also preserves the outward basolateral current. It uses nonnegative coefficients in the listed positive directions. Thus NBC does not add an independent conserved-coordinate direction to the complete inherited source set. Neither identity means the corresponding shared nonlinear laws can realise those flux changes at measured states.

Further bookkeeping degeneracies are

`s_K,a = s_K,b - s_para,K`,

`s_pump,a = s_pump,b - 3 s_para,Na + 2 s_para,K`,

`s_water,a = -s_water,b + s_water,para`,

and `s_TMEM16A = s_VRAC` for the Cl-only branch.

These are analytic identities, not the result of a subset search. Keeping every column in a single SVD makes the same issue visible numerically.

## The exact resting acid/base obstruction to inspect

In all genotypes the active equations imply

`dTA_i = J_NHE1 + 2J_NBC - J_AE2 - 2J_AE4`,

`dTIC_i = J_CO2,b + J_CO2,a + 2J_NBC - J_AE2 - 2J_AE4`.

Thus every stationary solution, irrespective of a basal NBC extension, must satisfy

`J_CO2,b + J_CO2,a = J_NHE1`.

For AE4 KO under the inherited resting zero gate these reduce further to

`J_NHE1 = J_AE2 = J_CO2,b + J_CO2,a`.

For AE2 KO under that same gate, `J_NHE1=2J_AE4`.

The exact gate is `u_sec=clip((Ca-0.058)/(0.25-0.058),0,1)`, followed by the literal `activation == 0.0` return in `MinimalNbcModel.evaluate`. Task 36 chose it for nesting, not from evidence of zero basal transport.

Do not interpret this known gate as the selected correction. At a fixed measured-state embedding, a positive TA residual cannot be repaired by *inward* NBC, which adds `+2` TA. An outward NBC coefficient also needs the correct affinity sign and current closure. The model's CO2 paths alter TIC but never TA; they cannot directly cancel a pure TA residual. Changes in pump, K and water affect acid/base balance only through state/current changes, not through a direct TA source.

## Measured-state embeddings must be honest

The four cohorts retain their own Cl/pH means and SE. AE4 WT and AE2 control are distinct observed controls; the model has no source-supported strain multiplier. Both controls therefore use the same WT constitutive law while retaining distinct observation uncertainty. Treating the two observed means as exact distinct equilibria of an identical genotype is itself a conditional diagnostic, not evidence of two biological roots.

At fixed measured `c=[Cl]_i` and `h=pH_i`, write carbonate fractions `alpha_0,alpha_1,alpha_2`, `a=alpha_1+2alpha_2`, deprotonated buffer fraction `b=1/(1+10^(pKa_buffer-h))`, and `w=[OH]-[H]`. The exact observation/charge constraints are

`nCl_i=c V_i`,

`nTA_i=a nTIC_i+B_fixed b+w V_i`,

`nNa_i+nK_i=(c+w)V_i+a nTIC_i+B_fixed b+X_fixed`.

They leave three independent cell coordinates, for example `(nNa_i,nTIC_i,V_i)`, plus five independent lumen coordinates. Accordingly each measured Cl/pH pair leaves **eight** core state degrees of freedom even after both compartment charge constraints. A full missing-source vector is not a directly observed quantity.

Recommended single conditional embedding for a physically consistent diagnostic, if the orchestrator chooses to avoid an unreported latent-state fit:

1. Take the saved Task 48 genotype root's Na concentration and its complete lumen state as explicit nuisance reference values. AE2 control uses the WT saved root. Do not call these experimental measurements.
2. Evaluate inherited NHE1 at that Na and measured pH. The necessary stationary difference balance fixes CO2 concentration algebraically:

   `CO2_i = (P_b CO2_b + P_a CO2_l - J_NHE1)/(P_b+P_a)`.

   Then `TIC_i=CO2_i/alpha_0(h)`. Reject the embedding if this is nonpositive; do not clip.
3. Cell water stationarity requires `O_i=(L_b O_b+L_a O_l)/(L_b+L_a)`. Combining this with charge and the fixed osmole accounting gives

   `V_i=(X_fixed+B_fixed b+M_other+B_fixed)/(O_i-2c-(a+1)TIC_i-w)`.

   Then `K_i=c+a TIC_i+w+(B_fixed b+X_fixed)/V_i-Na_i`, and all amounts follow. Require positive volume and concentrations. No optimisation is needed.
4. This construction enforces only charge, measured Cl/pH, the carbon-minus-alkalinity identity, and cell water stationarity. It does **not** enforce all solute/lumen balances and must not be described as a measured equilibrium or a solved genotype rest. Evaluate every remaining residual once and publish all retained and derived nuisance values.

This is one transparent embedding, not proof that hidden state choices are correct. More generally use the exact eight-dimensional observation-compatible manifold above to explain the non-identifiability; do not reinterpret a conditional residual as a unique empirical source vector. A simpler saved-root replacement embedding is mathematically allowed if labelled, but its large water/CO2 residuals may mostly report arbitrary hidden-state retention.

## Single fixed matrix and cone protocol

Predeclare the full dictionary above, signs, units, genotype masks, row scales, and the one embedding before evaluating data. Include all four cohorts in one saved matrix construction. Verify that the scalar source decomposition reproduces each inherited full core RHS; this is a bookkeeping identity check, not a new trajectory.

Publish the complete raw signature matrix, the oriented cone convention, and any fixed shared-law tying matrix in the same artefact. The stoichiometric matrix alone is a generous feasibility test. A shared-law coefficient has to act across cohorts through the actual state dependence and genotype mask. It must not become four independently fitted pathway multipliers. Pump partition and K-channel capacity partition are fixed inherited laws: treating apical and basolateral pump multipliers independently would be an extra model change.

Use full SVD with a declared numerical tolerance and disclose nullspace directions, orthogonal projection, and singular values. Scale species and water rows explicitly; raw fmol/s and pL/s cannot form a meaningful unqualified Euclidean distance. Positive row scaling preserves rank but changes the least-squares projection, so metric-dependent vector magnitudes are not invariants.

Cone membership must distinguish: (i) addition of forward flux, (ii) reduction of an existing flux, and (iii) reversal of a reversible law. Pump turnover is irreversible; channel and exchanger signs are constrained by driving force at the declared state. A physically represented dormant NBC signature is not an available positive-capacity resting flux under the inherited zero gate.

If needed, use one continuous constrained least-squares solve with all predeclared columns and equality constraints simultaneously. Do not subsequently drop columns, refit subsets, or select a mechanism by the largest coefficient. Unrestricted capacity scaling has a trivial all-capacities-zero stationary solution; do not present that as physiological closure. Direct additive cancellation `delta j=-j` is the same pitfall. Preserve source-fixed laws and independent nonzero requirements explicitly, or state that the relaxation cannot identify a correction.

## Uncertainty and the correction gate

Resting Cl/pH SE are observation uncertainty, not hard state bounds or SD. Carry first-order sensitivities of residuals to the two measurements, with a clearly labelled independent-error working covariance. The Cl–pH covariance is unreported; for any residual functional the unknown-correlation one-SE envelope is bounded by `|dR/dCl| SE_Cl + |dR/dpH| SE_pH`. This is a covariance bound on linearised uncertainty, not a deterministic physiological interval.

Also expose sensitivities to the eight hidden coordinates (or the smaller subset retained in the declared embedding). If an asserted source direction can be reproduced by an observation-invisible state perturbation, the measurements have not identified it. No additional numerical mechanism search resolves that structural issue.

A valid uniquely identified correction requires that the desired joint source direction survives the source nullspace, shared-law constraints, measurement error and latent-state ambiguity, with a thermodynamically consistent sign. If that condition fails, the defensible output is the precise balance identities and nullspace certificate above plus the computed conditional residual/projection; no basal NBC or other model edit is justified. This is an equation-level non-identifiability result, not a menu of candidate transporters.

