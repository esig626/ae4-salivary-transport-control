# Task 39 source contract: published two-state NKCC1 law

This file records the source-fixed NKCC1 equation to be implemented and tested in Task 39. It is not a fitting specification.

## Primary model source

Vera-Sigüenza E, Catalán MA, Peña-Münzenmayer G, Melvin JE, Sneyd J. *A Mathematical Model Supports a Key Role for Ae4 (Slc4a9) in Salivary Gland Secretion.* Bulletin of Mathematical Biology 80, 255-282 (2018). DOI: 10.1007/s11538-017-0370-6.

Appendix 1 states that the salivary model uses the Palk et al. (2010) two-state reduction of the Benjamin and Johnson (1997) NKCC1 model. Equation 17 is

`J_NKCC1 = alpha_NKCC1 * (a1 - a2 * Na_i * K_i * Cl_i^2) / (a3 + a4 * Na_i * K_i * Cl_i^2)`.

Positive `J_NKCC1` is one inward NKCC1 transport cycle, with the intracellular amount source signature

`(+1 Na, +1 K, +2 Cl, 0 TIC, 0 TA)`.

The printed parameter values are:

- `a1 = 157.5 s^-1`
- `a2 = 2.0096e7` in the paper's printed inverse-fourth-concentration convention
- `a3 = 1.0306 s^-1`
- `a4 = 1.3852e6` in the paper's printed inverse-fourth-concentration convention

The archived historical MATLAB implementation stores cell concentrations in mM and explicitly multiplies the two fourth-order coefficients by `1e-12`. Therefore the executable coefficients for the current mM concentration convention are exactly

- `a2_mM = 2.0096e-5 mM^-4 s^-1`
- `a4_mM = 1.3852e-6 mM^-4 s^-1`.

Task 39 must test the equivalence between evaluating the source coefficients on concentrations converted from mM to M and evaluating these converted coefficients directly on mM concentrations. Do not silently choose a different unit convention.

## Effective whole-cell prefactor

The paper's `alpha_NKCC1` is a model-determined transporter-density scale in a different historical geometry/amount convention. Do not transplant its absolute density directly into the modern conserved-amount model.

Instead, determine one whole-cell effective prefactor algebraically from the already accepted Task 31 R09 WT resting state. This is a nesting calibration only, not a phenotype fit and not an optimiser.

Frozen Task 31 R09 values:

- `Na_i = 11.636125748680639 mM`
- `K_i = 116.76320932871538 mM`
- `Cl_i = 60.30496692587399 mM`
- accepted NKCC1 cycle influx `J_NKCC1_rest = 0.1281222023866353 fmol/s`

For the published shape law, these give

- `X_rest = Na_i * K_i * Cl_i^2 = 4941065.334966328 mM^4`
- source-law shape factor `F_rest = 7.391062769440939`
- effective whole-cell prefactor `alpha_eff = J_NKCC1_rest / F_rest = 0.017334746894096052 fmol/s`

Task 39 must recompute these values independently and reject the implementation if they do not agree to numerical precision.

At REST the existing NKCC1 stimulus multiplier is exactly 1.0, so using this `alpha_eff` must reproduce the accepted Task 31 NKCC1 source exactly and therefore nest the accepted resting RHS without a new root solve.

## Current-architecture isolation rule

Task 39 replaces only the concentration-response law of NKCC1.

Keep the already-declared stimulus-dependent NKCC1 activity multiplier exactly unchanged. In particular, do not remove, refit, or replace the inherited N1 algebraic stimulation pathway. The point of Task 39 is to ask whether the published Palk/Benjamin concentration dependence prevents the large compensatory increase seen with the generic thermodynamic `tanh` law, while holding the rest of the current architecture fixed.

The old generic NKCC1 law must remain available as a comparator in source code but must not be active in the Task 39 production model.

## Source-level verification requirements

Before any full dynamics, test all of the following:

1. Direct evaluation of Eq. 17 agrees with an independent literal implementation.
2. M-space and mM-space evaluations agree after the exact `1e-12` fourth-order conversion.
3. The source signature is exactly `1 Na : 1 K : 2 Cl` inward per positive cycle.
4. The law reverses at `Na_i * K_i * Cl_i^2 = a1/a2_mM`; do not clip negative flux.
5. A multiplicative NKCC1 activity multiplier scales the whole cycle flux without moving that reversal.
6. The algebraically derived `alpha_eff` reproduces `0.1281222023866353 fmol/s` at the accepted Task 31 R09 resting state.
7. With the Task 39 model at REST, the complete accepted Task 31 R09 RHS, voltages, and water outflow remain unchanged within the existing production tolerances, without solving another root.

## Held-out biological context

The 2018 paper reports that AE4 knockout produces only a small NKCC1 increase in the model, and Peña-Münzenmayer et al. (2015) reported no detectable difference in NKCC1 activity between control and AE knockouts. These observations are validation context only. They must not be used to tune `alpha_eff`, the source coefficients, the stimulus multiplier, AE4, NBC, NHE1, or any other parameter.

Task 38 is the frozen generic-law comparator:

- AE4=0.05 cumulative secretion ratio to WT: `0.984315354219`
- AE4=0.00 cumulative secretion ratio to WT: `0.981604521729`
- integrated/declared NKCC1 loading increase after AE4 loss: approximately `+24.13%` at 5% expression and `+26.89%` at null.

Task 39 must report whether the source-fixed Palk law reduces that compensatory NKCC1 response. It must not tune the model to force such a reduction.
