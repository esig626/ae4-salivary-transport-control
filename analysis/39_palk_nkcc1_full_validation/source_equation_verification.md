# Task 39 source and nesting verification

The final source gate and exact REST gate both passed before integration.
The only scientific replacement is the NKCC1 concentration-response core.

## Primary equation and concentration convention

[Vera-Sigüenza et al. (2018), Appendix 1, Eq. 17](https://pmc.ncbi.nlm.nih.gov/articles/PMC5792321/)
identifies the Palk/Gin two-state reduction of Benjamin–Johnson NKCC1 and gives

\[
X=[Na]_i[K]_i[Cl]_i^2,\qquad
F(X)=\frac{a_1-a_2X}{a_3+a_4X}.
\]

Its parameter table prints 157.5, 2.0096e7, 1.0306, and 1.3852e6, respectively.
The printed concentration-unit labels require the historical implementation
convention to resolve their numerical interpretation. This is documented,
not silently corrected or treated as a free parameter choice.

The immutable archived [Parameters.m](https://github.com/esig626/ae4-salivary-transport-control/blob/b9d4e769747e636ad79d44a23c64c729f31b506d/archive/legacy-2017/Ae4_Dynamics_Project/Original/Parameters.m)
labels a2/a4 as inverse fourth powers of M and defines concentrations in mM.
[Saliva_Ae4.m](https://github.com/esig626/ae4-salivary-transport-control/blob/b9d4e769747e636ad79d44a23c64c729f31b506d/archive/legacy-2017/Ae4_Dynamics_Project/Original/Saliva_Ae4.m)
sets `s=1e-3` and inserts `(Na*s)*(K*s)*(Cl*s)^2` in both numerator and denominator.
[Par.m](https://github.com/esig626/ae4-salivary-transport-control/blob/b9d4e769747e636ad79d44a23c64c729f31b506d/archive/legacy-2017/Ae4_Dynamics_Project/Par.m)
stores the corresponding converted a2/a4 values directly.

Thus `X_M = (1e-3)^4 X_mM = 1e-12 X_mM`. The production coefficients are exactly

| Coefficient | Task 39 value |
| --- | ---: |
| a1 | 157.5 s^-1 |
| a2_mM | 2.0096e-5 mM^-4 s^-1 |
| a3 | 1.0306 s^-1 |
| a4_mM | 1.3852e-6 mM^-4 s^-1 |

The historical code's a1 is 157.55; Task 39 follows the required published
157.5 exactly. No alternate coefficient set was evaluated.

At accepted REST, evaluating source coefficients with M concentrations gives
`F = 7.391062769440935`; direct mM evaluation gives `7.391062769440939`.
The relative difference is below the declared `2e-15` precision tolerance.
Exact decimal arithmetic also verifies both fourth-order coefficient conversions.

## Reversal, source signature, and scale

The denominator is positive for positive concentrations. Reversal is at
`X_rev = a1/a2_mM = 7837380.573248408 mM^4`. The tests construct this product
algebraically, verify zero within roundoff, and verify both signs around it.
Negative flux is preserved through the actual homeostasis source assembly.
This is the fixed source-law reversal; no new bath-dependent coefficient or
voltage term is introduced.

A signed inward cycle supplies `(+1 Na,+1 K,+2 Cl,0 TIC,0 TA)`. Its charge is
`1+1-2=0`. Tests independently remove the unchanged NHE1/AE2 contributions
from homeostasis sources and recover this signature in both directions.

The shape is dimensionless in the contract's rate convention. The effective
whole-cell cycle prefactor has units fmol/s. It is obtained by one division:

| Algebraic quantity | Recomputed value |
| --- | ---: |
| X_rest (mM^4) | 4941065.334966328 |
| F_rest | 7.391062769440939 |
| J_rest (fmol cycles/s) | 0.1281222023866353 |
| alpha_eff = J_rest/F_rest (fmol/s) | 0.017334746894096052 |

A separate 60-digit Decimal evaluation agrees within floating-point precision.
The accepted Na/K/Cl tuple is recovered exactly from the saved conserved amounts.
The full freeze payload and SHA-256 are in `results/39_palk_nkcc1_full_validation/alpha_eff_freeze.json`.
It was written before any trajectory and never refitted.

## Implementation and gate evidence

`Nkcc1Kinetics` selects `palk_benjamin_2010_eq17` explicitly. The default
`generic_tanh` branch remains available. The selector is forwarded through
the existing core, NKCC1 wrapper, and NBC wrapper. Inherited parameter
dataclasses and their hashes are preserved. The NKCC1 affinity diagnostic
under Palk is `log(X_rev/X)`, consistent with the selected core's reversal.

The original N1 expression multiplies the entire cycle flux: 1.0 at REST and
1.75 for the prescribed stimulus. Both multipliers, the full regulatory
object, stimulus, AE4, NBC, and whole-cell parameter hashes are unchanged.
No new dynamical coordinate exists.

All seven focused source tests pass. One permitted literal test correction was
used: an initial bitwise coefficient-conversion assertion encountered the
one-ULP difference between `1.3852e-6` and `1.3852e6*1e-12`. The corrected test
checks exact decimal conversion and floating-point agreement at `2e-15`.
The scientific coefficients and implementation were unchanged. Both initial
and final test logs are retained; the initial coding-test failure token is
superseded by the final passing source gate.

The saved Task 31 R09 state was evaluated directly under Task 31, frozen Task
38, and Task 39. Task 31 and Task 38 REST RHSs agree bitwise. Task 39 also
agrees bitwise on all 13 rows. Amount, volume, voltage and outflow differences
are exactly zero; regulatory RHS is unchanged. The inherited amount `1e-12`,
volume `1e-14`, voltage `1e-12`, outflow `1e-12`, and conservation tolerances
all pass. A fixed-state onset evaluation also passes bookkeeping and checks
that the NKCC1 source is the only permitted RHS difference.

There were no stationary/resting solves. The frozen acid-base speciation and
NBC algebraic current-closure routines retain their inherited internal
`brentq` evaluations; these were neither replaced nor used to find a new
resting state. Full evidence is in `source_verification.json` and
`rest_nesting.json`.
