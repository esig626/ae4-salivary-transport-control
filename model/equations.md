# Baseline equations

## Canonicalization status

This is a published-equation skeleton, not an implementation specification. It records the 2018 paper’s balances and constraints and marks conflicts that must be resolved from published evidence before Phase 01 code is written.

Let `e`, `i`, and `l` denote interstitium, cytoplasm, and lumen. Intracellular balances are amount balances because cell volume `ω_i` varies.

## Intracellular amount balances

From §2.3, Eqs. (1)–(6):

\[
\frac{d([\mathrm{Cl}^-]_i\omega_i)}{dt}
=2J_{Nkcc1}+J_{Ae4}+J_{Ae2}-\frac{I_{CaCC}}{Fz_{Cl}}.
\]

\[
\frac{d([\mathrm{Na}^+]_i\omega_i)}{dt}
=J_{Nkcc1}+J_{Nhe1}
-\frac{[\mathrm{Na}^+]_i}{[\mathrm{Na}^+]_i+[\mathrm{K}^+]_i}J_{Ae4}
-3J_{NaK}.
\]

\[
\frac{d([\mathrm{K}^+]_i\omega_i)}{dt}
=J_{Nkcc1}+2J_{NaK}-\frac{I_{CaKC}}{Fz_K}
-\frac{[\mathrm{K}^+]_i}{[\mathrm{Na}^+]_i+[\mathrm{K}^+]_i}J_{Ae4}.
\]

\[
\frac{d([\mathrm{HCO}_3^-]_i\omega_i)}{dt}
=J_{Buffer}-2J_{Ae4}-J_{Ae2}.
\]

\[
\frac{d([\mathrm{H}^+]_i\omega_i)}{dt}=J_{Buffer}-J_{Nhe1}.
\]

\[
\frac{d([\mathrm{CO}_2]_i\omega_i)}{dt}=J_{CO2}-J_{Buffer}.
\]

These equations require a single dimensional convention for every `J` and `I/(Fz)` term. The paper states in §2.3 that balance fluxes must have units of moles/time, but Appendix 7 describes `J_Ae4` as concentration/time.

## Membrane potentials

From §2.4, Eqs. (7)–(8):

\[
C_m\frac{dV_b}{dt}=-I_{CaKC}-I_{NaK}+I_K^t+I_{Na}^t,
\]

\[
C_m\frac{dV_a}{dt}=-I_{CaCC}-I_K^t-I_{Na}^t.
\]

The simulations replace the derivatives with zero (§2.9), giving two algebraic QSS constraints. The convention is that `V_a` is measured lumen-to-cytoplasm and `V_b` interstitium-to-cytoplasm.

## Water and volume

The constitutive water relations in §2.5 are:

\[
q_a=P_a\left(\sum[c]_l+\Psi_l-\sum[c]_i-\frac{x_i}{\omega_i}\right),
\]

\[
q_b=P_b\left(\sum[c]_i+\frac{x_i}{\omega_i}-\sum[c]_e\right),
\]

\[
q_t=P_t\left(\sum[c]_l+\Psi_l-\sum[c]_e\right),
\qquad q_{tot}=q_a+q_t.
\]

The published cell-volume balance is internally inconsistent:

- §2.5 Eq. (9) prints `dω_i/dt = q_a - q_b`.
- §2.7 prints `dω_i/dt = q_b - q_a`.
- The formulas imply that `q_a` uses the apical/lumen gradient and `q_b` the basolateral/interstitial gradient, while the sentence immediately after Eq. (9) swaps the membrane labels.

Phase 10 resolves the physical equation as

\[
\frac{d\omega_i}{dt}=q_b-q_a,
\]

because `q_a` is positive from cell to lumen and `q_b` is positive from bath
to cell. The same sign appears in the collected system and in the upstream
Palk et al. (2010) model. At steady state, either printed sign imposes the
same condition `q_a=q_b`, so the exact Phase 10 steady identities are
independent of the correction.

The ionic sums printed after Eq. (11) include K+, Na+, H+, Cl−, and HCO3− in lumen; those species plus CO2 in cytoplasm and interstitium. `x_i/ω_i` accounts for impermeant intracellular negative charge and `Ψ_l` for untracked neutral luminal osmolytes/proteins.

## Luminal balances

From §2.6, Eqs. (14)–(16), with fixed `ω_l`:

\[
\omega_l\frac{d[\mathrm{Na}^+]_l}{dt}
=\frac{I_{Na}^t}{Fz_{Na}}-q_{tot}[\mathrm{Na}^+]_l,
\]

\[
\omega_l\frac{d[\mathrm{K}^+]_l}{dt}
=\frac{I_K^t}{Fz_K}-q_{tot}[\mathrm{K}^+]_l,
\]

\[
\omega_l\frac{d[\mathrm{Cl}^-]_l}{dt}
=\frac{I_{CaCC}}{Fz_{Cl}}-q_{tot}[\mathrm{Cl}^-]_l.
\]

The final subscript `l` in the Cl− loss term is supplied by the §2.7 collected system and the surrounding prose; Eq. (16) itself omits it typographically. This correction must be recorded explicitly if adopted in Phase 01.

## Electroneutrality and tight-junction constraints

Appendix 9, Eqs. (32)–(35), imposes electroneutrality in all three compartments. For the cell:

\[
[\mathrm{K}^+]_i+[\mathrm{Na}^+]_i+[\mathrm{H}^+]_i
-[\mathrm{Cl}^-]_i-[\mathrm{HCO}_3^-]_i-\frac{x_i}{\omega_i}=0.
\]

Appendix 10 uses linear current-voltage relations for tight-junction Na+ and K+, their Nernst potentials, and

\[
V_t=V_a-V_b.
\]

## Constitutive equations

The full published constitutive forms are located rather than duplicated here:

| Term | Published location |
|---|---|
| `J_Nkcc1` | Appendix 1 Eq. (17) |
| `J_NaK` | Appendix 2 Eq. (18) |
| CaKC open probability, flux/current, Nernst potential | Appendix 3 Eqs. (19)–(21) |
| CaCC open probability, flux/current, Nernst potential | Appendix 4 Eqs. (22)–(24) |
| `J_Ae2` | Appendix 5 Eq. (25) |
| `J_Nhe1` | Appendix 6 Eq. (26) |
| Ae4 state system and simplified flux | Appendix 7 Eqs. (27)–(31), Fig. 10 |
| CO2 equilibrium, membrane diffusion, and `J_Buffer` | Appendix 8 |
| Tight-junction Na+ and K+ relations | Appendix 10 Eqs. (36)–(37) and following Nernst equations |

See `fluxes.md` for direction, provenance, and unit warnings.

## Numerical form reported by the paper

With QSS potentials, the paper reports a DAE with ten differential equations and two algebraic constraints, solved with MATLAB `ode15s`. It does not report solver tolerances, consistent-initialization procedure, Jacobian handling, or a numerical table for the prescribed Ca2+ drive. These are Phase 01 blockers, not values to infer here.
