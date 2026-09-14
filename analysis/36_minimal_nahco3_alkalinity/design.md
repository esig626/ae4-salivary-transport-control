# Task 36 design: stimulus-recruited electrogenic 1 Na : 2 HCO3 pathway

## Status

This document supersedes the earlier 1 Na : 1 HCO3 electroneutral draft.

The corrected structural surrogate is an electrogenic NBCe-like basolateral pathway with

`Na_o + 2 HCO3_o <-> Na_i + 2 HCO3_i`.

It is a mechanism test, not an isoform assignment and not a measured transporter-abundance claim.

## Why this exists

The accepted Task 31 R09 WT fluxes are

- `J_NKCC1,Cl = 0.2562444047732706 fmol/s`;
- `J_NHE1 = 0.007890580489555009 fmol/s`;
- `J_AE4,Cl = 0.004937067441546469 fmol/s`;
- `J_AE2,Cl = -0.001983554393538192 fmol/s`.

The old architecture has intracellular alkalinity balance

`dTA_i/dt = J_NHE1 - J_AE2 - 2 J_AE4 + ...`

so large productive AE4 chloride loading is impossible unless another mechanism supplies alkalinity without exporting chloride.

The corrected NBC source adds exactly that missing degree of freedom.

## Conserved-coordinate source signature

For positive inward NBC cycle flux `J_B`, one cycle adds

- `+1 Na_i`;
- `0 K_i`;
- `0 Cl_i`;
- `+2 TIC_i`;
- `+2 TA_i`.

Thus the cell source vector is

`(+J_B, 0, 0, +2 J_B, +2 J_B)`.

The net charge added to the cell is

`(+1) - 2 = -1`

charge-equivalent per inward cycle. Therefore the transporter is electrogenic and must enter the basolateral current closure.

With `V_b = phi_i - phi_e`, the positive-inward dimensionless affinity is

`A_B = log((Na_e HCO3_e^2)/(Na_i HCO3_i^2)) + V_b / V_T`,

where `V_T = RT/F`.

The bounded reversible flux is

`J_B = u_sec G_B tanh(A_B / w)`

with the existing generic width `w = 2`.

Positive inward NBC flux carries conventional current

`I_B = F * J_B`

in the cell-to-bath direction after conversion from fmol/s to mol/s.

## Secretory recruitment and exact REST nesting

The NBC implemented here is the stimulus-recruited increment, not a claim that native basal NBC transport is absent.

Use the parameter-free protocol coordinate

`u_sec = clip((Ca - 0.058)/(0.25 - 0.058), 0, 1)`.

Therefore

- at Task 31 REST (`Ca = 0.058 uM`), `u_sec = 0` exactly and the accepted Task 31 resting model is unchanged;
- at the standard central stimulated condition (`Ca = 0.25 uM`), `u_sec = 1`.

The same coordinate gives the already-discussed source-fixed NHE1 stimulation

`G_NHE(t) = 1 + (2.3 - 1) u_sec`.

No NHE1 kinetic constant, affinity, stoichiometry, proton modifier, state or time constant is changed.

## Corrected whole-cell intracellular balances

Let

- `J_N` be the NKCC1 cycle flux, so NKCC1 chloride loading is `2 J_N`;
- `J_H` be NHE1 Na influx / H extrusion;
- `J_2` be AE2 Cl influx / HCO3 efflux;
- `J_4` be AE4 Cl influx;
- `r` be the Na fraction of the AE4-exported monovalent cation;
- `P_a`, `P_b` be apical and basolateral Na/K pump cycle fluxes;
- `J_Ka`, `J_Kb` be K-channel effluxes;
- `J_Cla` be apical CaCC Cl efflux;
- `J_B` be the new NBC cycle flux.

Then

`dnNa_i/dt = J_N + J_H + J_B - r J_4 - 3(P_a + P_b)`

`dnK_i/dt = J_N - (1-r)J_4 + 2(P_a + P_b) - J_Ka - J_Kb`

`dnCl_i/dt = 2 J_N + J_2 + J_4 - J_Cla`

`dnTIC_i/dt = J_CO2,b + J_CO2,a + 2 J_B - J_2 - 2 J_4`

`dnTA_i/dt = J_H + 2 J_B - J_2 - 2 J_4`

`dV_i/dt = q_b - q_a`.

The existing lumen equations are unchanged except indirectly through the new voltages and fluxes.

## Current closure

The apical QSS equation is unchanged:

`I_apical - I_para = 0`.

The basolateral equation becomes

`I_basolateral,old + I_B + I_para = 0`.

Because `I_B` depends on `V_b` through the NBC affinity, the first equation is eliminated analytically and the second is solved as one monotone scalar voltage equation. This is a voltage root, not a parameter search.

At `u_sec = 0`, the implementation delegates literally to the existing membrane closure so REST is an exact nesting limit.

## Charge and carbon checks

The transporter source signatures are

- NKCC1: `+1 +1 -2 = 0` charge;
- NHE1: `+1 Na` plus `+1 TA`, hence `+1 -1 = 0`;
- AE2: `+1 Cl`, `-1 TA`, hence `-1 -(-1) = 0`;
- AE4: `-1 cation +1 Cl -2 HCO3`, hence zero net charge;
- NBC: `+1 Na +2 HCO3`, hence `-1` charge-equivalent into the cell and the matching electrogenic current.

At sustained intracellular acid-base balance,

`0 = J_H + 2 J_B - J_2 - 2 J_4`

and TIC balance gives

`0 = J_CO2 + 2 J_B - J_2 - 2 J_4`.

Subtracting yields

`J_CO2 = J_H`.

Thus the existing CO2/NHE acid-base cycle remains intact; NBC supplies the extra bicarbonate equivalents consumed by AE4.

## Algebraic WT stimulated feasibility scale

Use the experimentally motivated stimulated positive chloride partition only as a WT architecture constraint:

`J_NKCC1,Cl / (J_NKCC1,Cl + J_AE4,Cl) = 0.70`

with AE2 negligible in the positive WT loading pool.

On the Task 31 R09 chloride-loading scale,

`J_NKCC1,Cl = 0.2562444047732706 fmol/s`,

so the required AE4 chloride loading is

`J_AE4,Cl = (3/7) J_NKCC1,Cl = 0.10981903061711597 fmol/s`.

With source-fixed stimulated NHE1

`J_H = 2.3 * 0.007890580489555009 = 0.01814833512597652 fmol/s`,

the 1:2 NBC alkalinity balance requires the NBC **cycle** flux

`J_B = J_AE4 - J_H/2 = 0.1007448630541277 fmol/s`.

This is the important correction relative to the discarded 1:1 draft: each NBC cycle supplies two bicarbonates, so the required cycle rate is about half the required bicarbonate-equivalent influx.

The electrogenic current associated with that cycle rate is about 9.72 pA.

Solving the full two-membrane current closure at the Task 31 R09 reference chemical state and `Ca=0.25 uM` gives the fixed reference capacity

`G_B = 0.11570913197464398 fmol/s`.

This is a derived structural feasibility scale, not a measured molecular capacity and not a knockout-phenotype fit.

## Held-out AE4 expression consequence

If NKCC1 remains the 70% positive chloride-loading arm and no compensatory refit is allowed, then

`R(x) = (C + x A)/(C + A) = 0.70 + 0.30 x`,

where `x` is AE4 expression relative to WT.

Therefore the algebraic non-compensated predictions are

- `x = 1.00`: `R = 1.000`;
- `x = 0.05`: `R = 0.715`;
- `x = 0.00`: `R = 0.700`.

These ratios are not fitted to the AE4 knockout secretion phenotype. They are a direct consequence of a 70/30 WT stimulated chloride-loading architecture.

## What remains to be tested numerically

The stoichiometric, charge, current, carbon and REST-nesting architecture is now explicit. The remaining question is only whether the full nonlinear ODE trajectory, with the existing NKCC1, AE4, AE2, pump, channel, lumen and water feedbacks, remains on a physiological branch and retains substantial AE4-dependent secretion.

Do not search parameter space before that direct test.
