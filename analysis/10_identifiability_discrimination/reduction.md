# Independent steady-state reduction and reconstruction gate

## Outcome

The published transporter stoichiometry and the form of the steady-state inverse problem can be reduced independently. The numerical 2018 salivary model cannot be reconstructed uniquely from the published specification, so the required full-model baseline reproduction did **not** pass. This document derives the recoverable structure and declares exactly what the executable diagnostic does and does not represent.

Primary source: Vera-Sigüenza et al. (2018), *A Mathematical Model Supports a Key Role for Ae4 (Slc4a9) in Salivary Gland Secretion*, [doi:10.1007/s11538-017-0370-6](https://doi.org/10.1007/s11538-017-0370-6), including its appendices.

## Recoverable published balance structure

Let the intracellular balance coordinates relevant to AE2/AE4 be ordered as

`u_c = (Na_i, K_i, Cl_i, HCO3_i)`.

Using positive signs for intracellular accumulation, the published cycle signatures are

`s2 = (0, 0, +1, -1)^T`

and

`s4 = (-rho_Na, -rho_K, +1, -2)^T`,

where `rho_Na = Na_i/(Na_i+K_i)`, `rho_K = K_i/(Na_i+K_i)`, and `rho_Na+rho_K=1`. Thus AE2 imports one chloride while exporting one bicarbonate; AE4 imports one chloride while exporting two bicarbonates and one monovalent cation. These signs agree with the amount balances in published Eqs. (1)--(4).

The published activity factors, before their unknown whole-cell scaling is resolved, are

`phi2(u) = Cl_e/(Cl_e+K_Cl) HCO3_i/(HCO3_i+K_B) - Cl_i/(Cl_i+K_Cl) HCO3_e/(HCO3_e+K_B)`

and

`phi4(u) = k_+ Cl_e HCO3_i^2 (Na_i+K_i) - k_- Cl_i HCO3_e^2 (Na_e+K_e)`.

All remaining transport, buffering, electrical, luminal, and water terms may be collected in `f0(u)`. The smallest independently defensible steady-state form is therefore

`F(u; G2, G4) = f0(u) + s2 G2 phi2(u) + s4 G4 phi4(u) = 0`.

This is a structural reduction, not a closed numerical model: `f0` contains the published Na/K pump, NKCC1, NHE1, channels, tight-junction currents, buffer/CO2 terms, lumen balances, membrane constraints, and water/volume equations.

## Why the published numerical model is not closed

| Blocking item | Consequence for a steady-state solve | Why it cannot be silently chosen |
| --- | --- | --- |
| AE4 density is reported in `amol/um^3`, while Eq. (31) and amount balances require a compatible whole-cell rate | `G4 phi4` cannot be placed on the same scale as AE2 and other fluxes | Any density/area/volume multiplier changes both the baseline and identifiability Jacobian |
| Currents `I` and molar fluxes `J` are interchanged without complete Faraday/area conversion | Ion and voltage residuals do not have a unique dimensional implementation | A guessed conversion changes `F_u` and therefore every projected signature |
| The numerical calcium input used for stimulation is not supplied | CaCC and CaKC conductances at the target secretory state are unknown | A resting-state value is not a substitute for the stimulated secretory steady state |
| Eq. (9) prints both `q_a-q_b` and `q_b-q_a` conventions, with prose/formula orientation conflicts | Cell volume and the water-flow observable cannot both be assigned uniquely | The physical orientation, not merely an equation-row sign, determines `Q` |
| Net CO2 flux is printed with a sign inconsistent with its influx description | Acid-base steady state may reverse | This affects both AE activity factors through bicarbonate |
| Luminal acid-base closure and DAE solver details are incomplete | The luminal state and algebraic potentials are not closed reproducibly | Different closures can give different regular branches |

These are the same critical provenance gaps identified in Phase 00; Phase 10 found no published information that resolves them.

## Executable baseline-anchored diagnostic

To verify analytical derivatives and illustrate observation geometry without inventing a full-model reproduction, the code uses dimensionless concentration ratios `x = u_c/u_c0`, activity multipliers `theta=(G2,G4)`, and

`F_d(x;theta) = A(x-1) - S(theta-1) = 0`.

Here:

- `u_c0=(25,120,50,12.1) mM` is the published resting concentration vector;
- the columns of `S` are the normalized published signatures above, multiplied by declared diagnostic levers `(0.08,0.12)`;
- `A` is a declared symmetric positive-definite effective restoring matrix; and
- `theta0=(1,1)` is an anchor, not an estimate of transporter density.

The exact solution and sensitivity are

`x*(theta)=1+A^{-1}S(theta-1)` and `D_theta x*=A^{-1}S`.

The matrix and all choices are versioned in `src/identifiability_discrimination/reduced_model.py`. The diagnostic is intentionally simple enough that its global geometry is known exactly.

## Numerical checks

| Check | Result | Interpretation |
| --- | --- | --- |
| Reduced residual at the published concentration anchor | maximum absolute residual `0` | Construction/implementation check only |
| Full published baseline reproduced | `false` | Reconstruction gate failed |
| Implicit derivative vs central difference (`h=1e-6`) | maximum absolute error `5.43e-11` | Derivative implementation verified |
| Positive restoring Jacobian | eigenvalues `0.867, 0.946, 1.052, 1.134` | Diagnostic branch is regular |
| Grid `G2,G4 in [0,2]` | minimum state ratio `0.853` over 1681 points | Diagnostic plotting domain remains positive |

As an additional transcription audit, direct substitution of the printed resting concentrations and table constants gives an AE2 dimensionless bracket of `-2.6954e-3` and an AE4 kinetic bracket of `4.1646e4` in its printed compound units. Multiplying by the printed `G2` and `G4` yields `-4.8707e-5` and `2.7487e4`, respectively. These numbers are intentionally not compared or inserted into one residual: the stated units do not put them on a common whole-cell flux scale. The check makes the scaling obstruction concrete.

Machine-readable records: `baseline_check.json`, `derivative_verification.json`, and `domain_summary.json` under `results/10_identifiability_discrimination/`.

## Sign-convention result

If an ambiguity merely changes a complete residual equation by a sign, `F` is replaced by `D F` with diagonal `D_ii in {-1,+1}`. Then

`-(D F_u)^{-1} D F_theta = -F_u^{-1}F_theta`.

The implicit sensitivity is invariant; the code verifies this to machine precision. This does **not** resolve the water-balance ambiguity, because swapping the physical orientation of `q_a` and `q_b` changes the content of the equation and the definition of `Q`, not just a complete row sign.

## Reduction verdict

The general rank/equivalence theory in `theory.md` is valid. Numerical salivary identifiability, minimal physiological panels, and mechanism-image separation remain unvalidated because the full observation map was not reproduced.
