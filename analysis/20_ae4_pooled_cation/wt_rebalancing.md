# Fixed-WT feasibility result

**All 20 requested pooled WT parameterizations are structurally infeasible.**
All ten conserved WT states are retained byte-for-byte in both P0 and P10
decision payloads. None supports either requested positive AE4 chloride flux
under the prescribed affinity and nonnegative capacity constraints.

| Quantity across ten roots | Minimum | Maximum |
|---|---:|---:|
| Pooled affinity `A` | -2.638392096 | -1.650095042 |
| Inherited intracellular HCO3, mM | 4.429886222 | 7.188876785 |
| HCO3 at pooled reversal, with other activities fixed, mM | 16.31639310 | 16.67804707 |
| Cl at pooled reversal, with other activities fixed, mM | 3.646512178 | 10.18204976 |

These reversal values are algebraic diagnostics only. No conserved state,
concentration, acid-base parameter or bath value was replaced by them.

For every root, the independent 50-digit mass-action calculation gives

\[
q=\frac{\mathrm{Cl_o}C_i\mathrm{HCO3_i}^2}
{\mathrm{Cl_i}C_o\mathrm{HCO3_o}^2}<1,
\qquad A=\log q<0.
\]

Local detailed balance therefore requires `k_forward<k_reverse` and
`J=E*g*(k_forward-k_reverse)<=0` for every `E,g>=0`. Both prescribed targets
are strictly positive: P0 equals the inherited legacy AE4 chloride source;
P10 equals `0.10*Lplus`. Consequently `J_target*A<0` at every requested
candidate. This contradicts the required `J*A>=0` before any balance equation
is solved. This proof covers every positive rate gauge with the prescribed
forward/reverse ratio, not just the chosen symmetric gauge.

The inherited positive-loading pool consists of AE4 and NKCC1; AE2 is a
negative chloride counterflux. For P10, the required remaining-loader factor
is recorded as `(Lplus-J_target)/J_NKCC1_Cl`. P0 retains that factor at unity
within floating-point rounding. These are target constraints, not accepted
capacity changes.

Task 18's fixed-state flux-coordinate reduction and full-production closure
checks are reused for the feasible path. The canonical objective remains
minimax absolute capacity log-fold, followed by squared log-fold. Here its
feasible set is empty, so both optimum values and accepted capacity payloads
are null. Running a local optimizer cannot resolve an analytically inconsistent
sign constraint. No optimizer failure is called infeasibility.

Changing K conductance, CaCC, pumps, NHE1, CO2 exchange or paracellular
capacities cannot change `A` at a fixed conserved state with unchanged bath,
chemistry and geometry. Algebraic voltage reclosure cannot help because the
pooled cycle is electroneutral. No finite amount of allowed capacity
rebalancing suffices. A formal signed carrier that would match the requested
flux is negative in every case; it is recorded only as a contradiction
diagnostic and is never passed to the production model.

Each payload records inherited resting closure, its exact state digest,
reference-only pooled sources/RHS, the sign certificate, and null accepted
parameters/closure. Reference pooled evaluations are deliberately not marked
as accepted WT roots. Conservation accounting still holds pointwise, while
the unrebalanced reference state is not a pooled steady state.

The frozen manifest retains every decision, including all infeasible ones.
It contains zero feasible models. New WT dynamics and genotype continuation
are consequently inapplicable under the mandated fixed-state protocol. The
post-freeze report must retain explicit unavailable rows rather than simulate
from a nonsteady reference state, fit replacement states or compare unpaired
genotypes.
