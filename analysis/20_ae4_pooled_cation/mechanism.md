# Task 20 pooled-cation mechanism

`src/modern_full_model/pooled_ae4.py` exposes the selectable
`POOLED_CATION_112_NO_SLIP` mechanism through the existing production evaluator
interface. `SR2_SHARED_112_QSS` remains selectable with its original parameters.
The production model, chemistry, membrane equations, water equations and
regulation files are unchanged.

The assumed positive cycle is

\[
\mathrm{Cl_o+C_i+2HCO3_i\rightleftharpoons Cl_i+C_o+2HCO3_o},
\qquad C_s=\mathrm{Na}_s+\mathrm K_s.
\]

Its affinity and rates are

\[
A=\log\frac{\mathrm{Cl_o}}{\mathrm{Cl_i}}+
\log\frac{C_i}{C_o}+2\log\frac{\mathrm{HCO3_i}}{\mathrm{HCO3_o}},
\quad k_+=k_0e^{A/2},\quad k_-=k_0e^{-A/2},
\quad J=E g(k_+-k_-)=2E g k_0\sinh(A/2).
\]

Here `E` is carrier amount in fmol, `g` is the inherited active-capacity
regulation gain, and `k0` is in s^-1. The primary gauge fixes `k0=1/s` for
both legacy provenance families. The inherited carrier is the reference for
capacity log-fold accounting. Neither gauge nor carrier amount is a measured
cross-architecture kinetic equivalence. Infeasibility from the sign of `A`
does not depend on this gauge.

There is exactly one scalar cycle. The donor Na fraction is `Na_i/(Na_i+K_i)`
for positive `J` and `Na_o/(Na_o+K_o)` for negative `J`. Its complement is
the K fraction. At `J=0` both cation sources vanish; the chosen partition has
no effect. No interpolation or extra slip degree of freedom is used. The
source vector is continuous at zero, although its cation derivatives can have
different one-sided limits.

Intracellular Na, K, Cl, HCO3, TIC and alkalinity sources are
`(-p J, -(1-p) J, J, -2J, -2J, -2J)`. Na and K have the same sign or vanish,
and their absolute source sum is `abs(J)`. The charge identity is exactly
`-p-(1-p)-1+2=0`. Numerical charge residuals are retained rather than silently
discarded. `J*A>=0` follows because `sinh(A/2)` has the sign of `A` for all
nonnegative carrier amounts and gains. This is the declared pooled affinity
and entropy accounting, not identification of microscopic cation cycles.

Single-cation absence is supported explicitly: the other cation takes the
whole donor partition. Each pooled side must have positive total cation
activity; chloride and bicarbonate activities must be positive. Negative or
nonfinite activities/capacities are rejected.

The 14 assay tests include all 12 required checks, exact rational charge
accounting, independently computed mass-action ratios over broad random
activity panels, pure and near-pure cations, regulation scaling and exact
legacy production-output reproduction at all ten inherited WT states.
Assay inputs are independent of secretion phenotypes.
