# Task 40 source contract: equal Na/K AE4 cation routing

## Objective

Task 39 showed that replacing the generic NKCC1 concentration-response law with the source-fixed Palk/Benjamin law does not solve AE4-loss compensation. With the existing donor-concentration AE4 routing, AE4 loss increases integrated NKCC1 loading by about 29.65% at AE4=0.05 and 33.10% at AE4=0.00, leaving only small secretion deficits.

Task 40 tests one mechanism change only: keep the existing net AE4 cycle law and all anion/base stoichiometry unchanged, but route the single monovalent-cation component equally through Na and K rather than weighting it by donor-side Na/K concentrations.

This is a structural hypothesis test, not a fitted routing fraction.

## Frozen parent model

Start from the completed Task 39 model merged to `main` at commit `afd101448763439f369f2682c467f069ea18a442`.

Keep fixed:

- Palk/Benjamin NKCC1 concentration-response law and `alpha_eff`;
- inherited NKCC1 stimulus activity multiplier;
- Cha NHE1 kinetics, carrier amount, and stimulus multiplier;
- electrogenic 1 Na : 2 HCO3 NBC, capacity, recruitment and current closure;
- total AE4 net cycle law, total capacity, 1 Cl : 1 monovalent cation : 2 HCO3 stoichiometry, and beta/cAMP/PKA regulation;
- AE2;
- Na/K pump;
- K and Cl channels;
- CO2 and acid-base closure;
- paracellular transport;
- water laws;
- bath and geometry;
- standard CCh/IPR/Ca/beta stimulus;
- all numerical tolerances unless a purely numerical solver retry is explicitly allowed.

The only scientific change is the Na/K partition of the one AE4 monovalent-cation component.

## Fixed equal-cation routing law

Let `J4` denote the existing routed AE4 net chloride source in the cell, with positive `J4` meaning chloride loading into the cell.

Preserve `J4` exactly as calculated by the inherited AE4 evaluator. Do not change its affinity, attempt rates, occupancy, capacity, expression scaling, or regulation.

For every forward or reverse cycle use the ensemble-average cation partition

- Na branch fraction = `0.5`
- K branch fraction = `0.5`.

Thus the intracellular conserved sources are

`S_Na_AE4 = -0.5 * J4`

`S_K_AE4  = -0.5 * J4`

`S_Cl_AE4 = +1.0 * J4`

`S_TIC_AE4 = -2.0 * J4`

`S_TA_AE4  = -2.0 * J4`.

For negative `J4`, all signs reverse automatically. Do not branch on donor concentrations.

The transported charge-equivalent source must remain exactly zero:

`S_Na + S_K - S_Cl - S_TA = 0`.

The total cation source must remain exactly the inherited one-cation-per-cycle value:

`S_Na + S_K = -J4`.

No new AE4 state, regulator, capacity, or free routing parameter is introduced.

## Algebraic preflight that must be reproduced

At the accepted Task 31 R09 reference state:

- `Na_i = 11.636125748680639 mM`
- `K_i = 116.76320932871538 mM`
- resting routed AE4 chloride source `J4 = 0.004937067441546469 fmol/s`.

The current donor-weighted Na fraction is

`r_rest = Na_i / (Na_i + K_i) = 0.09062450161183984`.

Changing only the cation partition to 50:50 therefore changes the old reference-state AE4 cation sources by

`delta S_Na = -(0.5 - r_rest) * J4 = -0.0020211144444590447 fmol/s`

`delta S_K  = +(0.5 - r_rest) * J4 = +0.0020211144444590447 fmol/s`.

Direct AE4 Cl, TIC and TA sources are unchanged, and the net charge perturbation is zero.

Therefore exact Task 31 REST nesting is not expected under this mechanism change. Task 40 must verify this algebraic residual shift before solving a new WT rest. It must not falsely require bitwise REST nesting.

## Resting-state rule

Because the mechanism changes Na/K source partition even at rest, Task 40 may perform exactly one intended WT stationary solve using the accepted Task 31 R09 state as the sole initial guess.

- No genotype-specific resting solve.
- No parameter calibration.
- No alternate routing fraction.
- No multistart.
- No continuation family.
- No search over roots.

One purely numerical retry of the same frozen scientific problem is allowed only if the first stationary method fails for a numerical reason. If no admissible WT rest is obtained, stop.

The accepted Task 40 WT rest, if obtained, is frozen before any stimulation or genotype evaluation.

## Dynamic validation

Only after the new WT rest passes physiology and conservation gates, run exactly three intended 600 s trajectories from that same frozen Task 40 WT resting state:

1. WT, AE4 expression = 1.0
2. AE4 expression = 0.05
3. AE4 expression = 0.00

All other scientific inputs are identical.

Do not solve genotype-specific rests.

## Mechanistic questions to report, not fit

The task must answer:

1. Does AE4 loss now leave intracellular Na higher than WT during sustained stimulation, particularly at 600 s?
2. Does the large Task 39 NKCC1 compensation (`+29.65% / +33.10%`) shrink, disappear, or reverse?
3. Does intracellular chloride remain lower after AE4 loss rather than being restored by NKCC1?
4. Does AE2 remain negligible as a positive chloride-loading route?
5. What happens to NBC, NHE1, pump, CaCC and paracellular chloride return as secondary consequences?
6. Does the secretion deficit persist or still collapse with time?
7. What cumulative secretion ratios/deficits are predicted at AE4=0.05 and AE4=0.00?

The approximately 35% experimental null deficit is held-out context only. Do not tune the routing fraction, AE4 capacity, NKCC1, NBC, NHE1, pump, channels, or stimulus to approach it.

## Absolute anti-search rule

No optimisation, parameter sweep, alternate Na/K split, one-at-a-time routing study, grid, random search, multistart, Bayesian/global/evolutionary search, Shapley analysis, second mechanism, or phenotype fit is permitted.

The 50:50 split is the only routing hypothesis evaluated in Task 40.
