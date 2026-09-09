# Published-source steady-state reduction and reconstruction gate

## Status and scope

This note is an independent derivation from the published model:

> Vera-Sigüenza E, Catalán MA, Peña-Münzenmayer G, Melvin JE, and Sneyd J (2018), “A Mathematical Model Supports a Key Role for Ae4 (Slc4a9) in Salivary Gland Secretion,” *Bulletin of Mathematical Biology* 80(2):255–282. [DOI 10.1007/s11538-017-0370-6](https://doi.org/10.1007/s11538-017-0370-6); [open author manuscript, PMCID PMC5792321](https://pmc.ncbi.nlm.nih.gov/articles/PMC5792321/).

No unpublished manuscript, historical code, or archived figure was used. Equation numbers below refer to the published paper.

**This is not an executable reconstruction of the 2018 model.** The published balances support exact symbolic steady-state reductions, but the published record does not specify a dimensionally closed numerical model. In particular, the water, carbon-dioxide, current/flux, and density scalings contain contradictions that prevent an honest baseline regression test. The purpose of this note is therefore to:

1. record the exact consequences of the published balances;
2. resolve the cell-volume sign by conservation;
3. state every additional closure needed for a small steady-state system;
4. expose `G_Ae2` and `G_Ae4` without concealing unit assumptions; and
5. identify the gate that must be passed before any numerical result is called a reproduction.

The following labels distinguish levels of support:

- **Published**: transcribed directly from the 2018 article.
- **Exact reduction**: algebraic consequence of the published balances, assuming their flux terms are commensurate.
- **Conservation resolution**: one of two inconsistent printed signs is selected from mass conservation.
- **Closure assumption**: not completely specified by the paper and must be declared before use.
- **Unit-dependent**: algebraically valid only after an explicit conversion of every current, concentration, volume, density, and permeability to one unit system.

## Notation and positive directions

Write intracellular concentrations as

\[
n=[\mathrm{Na}^+]_i,\qquad
k=[\mathrm{K}^+]_i,\qquad
c=[\mathrm{Cl}^-]_i,
\]

\[
b=[\mathrm{HCO}_3^-]_i,\qquad
h=[\mathrm{H}^+]_i,\qquad
d=[\mathrm{CO}_2]_i.
\]

Subscripts `e` and `l` denote the fixed interstitium and lumen, respectively. The variable cell volume is `omega_i`; the lumen volume is fixed. The paper prescribes intracellular calcium rather than giving it a state equation.

For the steady-state algebra, use positive particle fluxes

\[
N=J_{\mathrm{Nkcc1}},\quad
P=J_{\mathrm{NaK}},\quad
A_2=J_{\mathrm{Ae2}},\quad
A_4=J_{\mathrm{Ae4}},
\]

\[
H=J_{\mathrm{Nhe1}},\quad
B=J_{\mathrm{Buffer}},\quad
C=\frac{I_{\mathrm{CaCC}}}{Fz_{\mathrm{Cl}}},\quad
K_c=\frac{I_{\mathrm{CaKC}}}{Fz_K}.
\]

Here `C > 0` means chloride leaves the cell for the lumen, and `K_c > 0` means potassium leaves the cell for the interstitium. Since `z_Cl=-1`, an outward chloride particle flux has a negative electrical current. This distinction matters because Appendix 4 calls `G P(V-E)/F` a `J`, although it is `I/F`, not the positive particle flux `I/(F z_Cl)` used in the amount balances.

Let

\[
f_n=\frac{n}{n+k},\qquad f_k=\frac{k}{n+k}=1-f_n.
\]

## Published transporter functions with the two activities explicit

The paper's Appendix 5 gives

\[
A_2=G_{\mathrm{Ae2}}\,\phi_2(c,b),
\]

where

\[
\phi_2(c,b)=
\frac{c_e}{c_e+K_{\mathrm{Cl}}}
\frac{b}{b+K_B}
-
\frac{c}{c+K_{\mathrm{Cl}}}
\frac{b_e}{b_e+K_B}.
\tag{R1}
\]

Appendix 7 gives the simplified Ae4 expression

\[
A_4=G_{\mathrm{Ae4}}\,\phi_4(n,k,c,b),
\]

\[
\phi_4=
k_+c_e b^2(n+k)
-k_-c b_e^2(n_e+k_e).
\tag{R2}
\]

Equation (R2) is used here because it is the explicit activity form printed as Eq. (31), but it is not an exact algebraic consequence of the preceding displayed Markov flux as written. The preceding expression has a concentration-dependent denominator; it cannot reduce to a constant-`G_Ae4` bilinear difference unless that denominator is absorbed into a state-dependent effective `G_Ae4` or an additional approximation is made. In addition, Eq. (29) prints the gain term `k_3 cl_i B_i` in the exterior-facing state equation, whereas the cycle diagram and Eq. (31) require an exterior complex and exterior chloride. Therefore this note treats Eq. (31) as a published phenomenological constitutive choice, not a verified reduction of Eqs. (27)--(30).

Thus the inferential parameter vector can be written

\[
\theta=(G_2,G_4)=(G_{\mathrm{Ae2}},G_{\mathrm{Ae4}})
\]

without changing the published transporter kinetics.

The other constitutive functions required by the reduction are

\[
N(n,k,c)=\alpha_{\mathrm{Nkcc1}}
\frac{a_1-a_2nkc^2}{a_3+a_4nkc^2},
\tag{R3}
\]

\[
P(n)=\alpha_{\mathrm{NaK}}r
\frac{k_e^2n^3}{k_e^2+\alpha_1n^3},
\tag{R4}
\]

and

\[
H(n,h)=G_{\mathrm{Nhe1}}
\left[
\left(\frac{h}{h+K_H}\right)^2
\frac{n_e}{n_e+K_{\mathrm{Na}}}
-
\frac{n}{n+K_{\mathrm{Na}}}
\left(\frac{h_e}{h_e+K_H}\right)^2
\right].
\tag{R5}
\]

The buffer expression is

\[
B(d,h,b)=k_1d-k_{-1}hb.
\tag{R6}
\]

Equations (R1)--(R6) are published functional forms, not yet a common-unit implementation. The numerical tables describe `G_Ae2` as fmol/s, `G_Ae4`, `alpha_Nkcc1`, and `alpha_NaK` as amol/micrometre cubed, and the main balances as total moles/time. No membrane-area or whole-cell conversion is supplied.

## Published steady balances

Setting the derivatives in Eqs. (1)--(6) to zero gives

\[
0=2N+A_4+A_2-C,
\tag{S1}
\]

\[
0=N+H-f_nA_4-3P,
\tag{S2}
\]

\[
0=N+2P-K_c-f_kA_4,
\tag{S3}
\]

\[
0=B-2A_4-A_2,
\tag{S4}
\]

\[
0=B-H,
\tag{S5}
\]

\[
0=J_{\mathrm{CO2}}-B.
\tag{S6}
\]

For nonzero total lumen flow `Q=q_tot`, Eqs. (14)--(16), with the missing lumen subscript in Eq. (16) restored from section 2.7, give

\[
\frac{I_{\mathrm{Na}}^t}{F}=Qn_l,\qquad
\frac{I_K^t}{F}=Qk_l,\qquad
C=Qc_l.
\tag{S7}
\]

The QSS membrane equations in section 2.9 are

\[
-I_{\mathrm{CaKC}}-I_{\mathrm{NaK}}+I_K^t+I_{\mathrm{Na}}^t=0,
\tag{S8}
\]

\[
-I_{\mathrm{CaCC}}-I_K^t-I_{\mathrm{Na}}^t=0.
\tag{S9}
\]

The conversion `I_NaK=F P` is consistent with the pump's net movement of one positive charge out of the cell per cycle, but it is not written explicitly in the paper. It is therefore a closure assumption whenever (S8) is combined with (R4).

## Exact symbolic reductions

This section uses only the steady balances and stoichiometry. It does not require a choice of water sign, carbon-dioxide sign, or numerical unit conversion beyond the assumption that terms added in one balance are commensurate.

### Acid/base and chloride identities

Equations (S4)--(S6) imply

\[
B=H=J_{\mathrm{CO2}}=A_2+2A_4.
\tag{E1}
\]

Equation (S1) gives

\[
C=2N+A_2+A_4.
\tag{E2}
\]

Solving the two-by-two stoichiometric system (E1)--(E2) gives

\[
\boxed{A_4=B-C+2N},
\tag{E3}
\]

\[
\boxed{A_2=2C-4N-B}.
\tag{E4}
\]

Consequently, wherever both kinetic factors are nonzero,

\[
G_4=\frac{B-C+2N}{\phi_4},\qquad
G_2=\frac{2C-4N-B}{\phi_2}.
\tag{E5}
\]

This identity is the clearest link from the published reduction to the identifiability problem. A steady measurement of flow and luminal chloride determines `C=Q c_l`, but it does not by itself determine `B` and `N`. Therefore it cannot separate `G_2` and `G_4` without additional measurements, trusted nuisance kinetics, or perturbations.

### Charge-balance identity and redundant potassium balance

Adding (S2) and (S3), then using (E1)--(E2), yields

\[
\boxed{C=P+K_c}.
\tag{E6}
\]

The same identity follows from the two QSS current balances. Once (S1), (S2), (S4), (S5), and (E6) hold, the potassium balance (S3) is the negative of the sodium residual and is redundant. This is a conservation redundancy, not a dropped physiological mechanism.

### Luminal electroneutrality consequence

Equations (S7) and (S9) give

\[
Q(n_l+k_l)=Qc_l.
\]

For `Q != 0`,

\[
\boxed{c_l=n_l+k_l}.
\tag{E7}
\]

The paper also imposes lumen electroneutrality,

\[
n_l+k_l+h_l-c_l-b_l=0.
\]

Together with (E7), this requires

\[
h_l=b_l.
\tag{E8}
\]

The rounded Table 1 values satisfy this approximately: `pH_l=6.81` gives `h_l=1.55e-4 mM`, while the table lists `b_l=1.5e-4 mM`. If lumen acid/base values are not fixed in this way, the published ten-state system lacks the additional lumen acid/base balances needed for exact closure.

## Conservation resolution of the volume sign

The paper defines

\[
q_a=P_a(\Omega_l-\Omega_i),
\qquad
q_b=P_b(\Omega_i-\Omega_e),
\tag{W1}
\]

where `Omega_l` includes `Psi_l` and `Omega_i` includes `x_i/omega_i`.

When positive:

- `q_a` moves water from cell to lumen and therefore removes cell volume;
- `q_b` moves water from interstitium to cell and therefore adds cell volume.

Mass conservation therefore requires

\[
\boxed{\frac{d\omega_i}{dt}=q_b-q_a}.
\tag{W2}
\]

This is the sign printed in the collected system in section 2.7, not the opposite sign printed as Eq. (9) in section 2.5. It also gives restoring osmotic feedback: a decrease in volume raises `x_i/omega_i`, increasing basolateral influx and decreasing apical efflux. The Eq. (9) sign would reverse this feedback and violate the stated flux directions.

The upstream Palk et al. model from which the water relations were inherited
prints the same conservation equation; see Palk et al. (2010),
[DOI 10.1016/j.jtbi.2010.06.027](https://doi.org/10.1016/j.jtbi.2010.06.027),
[open full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC2954280/).

At steady state either printed sign reduces to `q_a=q_b`, so the exact identities (E1)--(E8) do not depend on this correction. Transient volume dynamics do depend on it.

## Conservation resolution of the carbon-dioxide sign

Appendix 8 first defines two influxes into the cell:

\[
J_{\mathrm{CO2},a}=P_{\mathrm{CO2}}(d_l-d),
\qquad
J_{\mathrm{CO2},b}=P_{\mathrm{CO2}}(d_e-d).
\]

Their conservation-consistent sum is

\[
J_{\mathrm{CO2}}^{\mathrm{cons}}
=P_{\mathrm{CO2}}(d_l+d_e-2d).
\tag{C1}
\]

The paper instead prints

\[
J_{\mathrm{CO2}}^{\mathrm{pub}}
=P_{\mathrm{CO2}}(2d-d_l-d_e)
=-J_{\mathrm{CO2}}^{\mathrm{cons}},
\tag{C2}
\]

while Eq. (6) adds `J_CO2` to the cell. Both conservation of the two individually printed membrane influxes and the sign in Eq. (6) select (C1). Thus (C1) is the defensible primary convention, while (C2) should be retained only as a printed-sign sensitivity branch. Unlike the volume case, the Table 1 baseline does not support the conservation resolution numerically: `d=6.6 mM` exceeds the listed `d_e=1.9 mM` and the equilibrium lumen value, so passive diffusion is outward although the forward buffer term must be replenished at steady state. This is a baseline/acid-closure inconsistency, not a reason to reverse Fickian transport.

Results relying only on (E1)--(E8) are independent of the choice. A numerical `F(u;theta)` should use (C1) as primary and carry (C2) far enough to report whether a conclusion depends on the paper's printed aggregate sign.

Appendix 8 also prints an intracellular chemical-equilibrium relation

\[
d=\frac{k_{-1}}{k_1}bh.
\]

Imposing this relation simultaneously with the dynamic mass-action expression (R6) makes `B=0` identically. At a nontrivial secretory steady state, (E1) would then force `H=A_2+2A_4=0`. The paper does not explain whether the equilibrium relation is used only to initialize or close non-dynamic compartments, or is enforced on the intracellular state. The conditional reduction below treats `(d,h,b)` through the dynamic mass-action balance and does **not** impose the intracellular equilibrium relation as an additional algebraic constraint.

## A smallest conditional steady-state system

This section gives a useful analytical form of `F(u;theta)=0`. It is **conditional**, not a claim that the published numerical model has been reconstructed.

### Additional assumptions

The reduction below assumes:

1. all transporter and channel terms have first been converted to total amount/time;
2. membrane potentials are QSS variables as in section 2.9;
3. the interstitium is fixed;
4. `h_l`, `b_l`, `d_l`, and `Psi_l` are fixed, with `h_l=b_l`;
5. `Q>0`;
6. impermeant intracellular charge `x_i` is fixed in moles and has the unit valence used in Eq. (32);
7. `I_NaK=F P` and the Appendix 10 `g` values are the Table 2 conductances after unit conversion;
8. (C1) is used as the conservation convention, with (C2) reserved for a printed-sign sensitivity check; and
9. the intracellular mass-action buffer is dynamic, so the separate intracellular chemical-equilibrium relation is not imposed simultaneously.

Use the six unknowns

\[
u=(n,c,b,h,n_l,k_l).
\tag{M1}
\]

The following quantities can then be eliminated algebraically.

### Buffer and osmotic eliminations

From `B=H`,

\[
d=\frac{H(n,h)+k_{-1}hb}{k_1}.
\tag{M2}
\]

Cell electroneutrality gives

\[
\frac{x_i}{\omega_i}=n+k+h-c-b.
\tag{M3}
\]

The intracellular osmolarity printed after Eq. (11) therefore simplifies to

\[
\Omega_i=n+k+h+c+b+d+\frac{x_i}{\omega_i}
=2(n+k+h)+d.
\tag{M4}
\]

Using (E7)--(E8),

\[
c_l=n_l+k_l,
\qquad
\Omega_l=2(n_l+k_l+h_l)+\Psi_l.
\tag{M5}
\]

In accordance with the ionic sums printed after Eq. (11), the fixed bath osmolarity used here is

\[
\Omega_e=n_e+k_e+h_e+c_e+b_e+d_e,
\]

whereas the paper's lumen sum does not include `d_l`.

The steady volume condition `q_a=q_b` gives

\[
\Omega_i^*=
\frac{P_a\Omega_l+P_b\Omega_e}{P_a+P_b}.
\tag{M6}
\]

Thus

\[
k=\frac{\Omega_i^*-d}{2}-n-h,
\tag{M7}
\]

and, after solving,

\[
\omega_i=\frac{x_i}{n+k+h-c-b}.
\tag{M8}
\]

Positivity of `k`, `omega_i`, all concentrations, and the denominator in (M8) is part of the admissible-state definition.

If the tabulated transporter “densities” must be multiplied by cell volume or membrane area to obtain total fluxes, the corresponding factor must be evaluated using (M8) and inserted into (R2)--(R4). This does not increase the number of reduced unknowns, but it changes the numerical residuals and their parameter derivatives. The paper does not state which scaling is intended.

The steady total flow becomes

\[
Q=\left(
\frac{P_aP_b}{P_a+P_b}+P_t
\right)(\Omega_l-\Omega_e).
\tag{M9}
\]

Equation (M9) is unit-dependent: if osmolarities are stored in mM while permeabilities use mol/L, the `1e-3` conversion must be explicit.

### Electrical eliminations

For fixed calcium, define the published open probabilities

\[
p_{\mathrm{Cl}}=
\left(\frac{[\mathrm{Ca}]_i}{[\mathrm{Ca}]_i+K_{\mathrm{CaCC}}}\right)^{\eta_1},
\qquad
p_K=
\left(\frac{[\mathrm{Ca}]_i}{[\mathrm{Ca}]_i+K_{\mathrm{CaKC}}}\right)^{\eta_2}.
\]

Define Nernst potentials

\[
E_{\mathrm{Cl}}=\frac{RT}{z_{\mathrm{Cl}}F}\log\frac{c_l}{c},
\qquad
E_{K,i}=\frac{RT}{F}\log\frac{k_e}{k},
\]

\[
E_{\mathrm{Na}}^t=\frac{RT}{F}\log\frac{n_l}{n_e},
\qquad
E_K^t=\frac{RT}{F}\log\frac{k_l}{k_e}.
\]

With

\[
C=2N+A_2+A_4,
\tag{M10}
\]

the channel and tight-junction relations give

\[
V_a=E_{\mathrm{Cl}}-\frac{FC}{G_{\mathrm{CaCC}}p_{\mathrm{Cl}}},
\tag{M11}
\]

\[
V_b=E_{K,i}+\frac{F(C-P)}{G_{\mathrm{CaKC}}p_K},
\tag{M12}
\]

\[
V_t=E_{\mathrm{Na}}^t+\frac{FQn_l}{g_{\mathrm{Na}}^t}
=E_K^t+\frac{FQk_l}{g_K^t}.
\tag{M13}
\]

### Six residual equations

The conditional six-dimensional steady problem is

\[
F(u;G_2,G_4)=0,
\]

with residuals

\[
F_1=H-A_2-2A_4,
\tag{M14}
\]

\[
F_2=J_{\mathrm{CO2}}(d)-H,
\tag{M15}
\]

\[
F_3=N+H-3P-f_nA_4,
\tag{M16}
\]

\[
F_4=C-Qc_l,
\tag{M17}
\]

\[
F_5=
E_{\mathrm{Na}}^t+\frac{FQn_l}{g_{\mathrm{Na}}^t}
-E_K^t-\frac{FQk_l}{g_K^t},
\tag{M18}
\]

\[
F_6=
E_{\mathrm{Cl}}-\frac{FC}{G_{\mathrm{CaCC}}p_{\mathrm{Cl}}}
-E_{K,i}-\frac{F(C-P)}{G_{\mathrm{CaKC}}p_K}
-E_{\mathrm{Na}}^t-\frac{FQn_l}{g_{\mathrm{Na}}^t}.
\tag{M19}
\]

Equations (M2), (M5), (M7), and (M9)--(M13) recover `d`, `c_l`, `k`, `Q`, and the membrane potentials. Equation (M8) recovers cell volume. The potassium balance is redundant by (E6).

This is the smallest useful reduction found without solving a model-specific nonlinear constitutive equation in closed form. It is suitable for symbolic rank and equivalence arguments. It is not yet suitable for numerical baseline claims because assumptions 1, 7, and 8 are unresolved by the published source.

## Published baseline targets

Table 1 and section 3.1 report the following nominal resting values. Rows marked “model-derived” in the paper are calibration outputs, not independent validation data.

| Quantity | Published nominal value |
| --- | ---: |
| `[Ca2+]_i` | 58 nM |
| `n`, `k`, `c`, `b` | 25, 120, 50, 12.1 mM |
| `pH_i`, `d` | 6.91, 6.6 mM |
| `x_i/omega_i0` | 82.8 mM |
| `V_a`, `V_b` | -50.24, -62.8 mV |
| `n_l`, `k_l`, `c_l` | 118.7, 5.6, 124.3 mM |
| `b_l`, `pH_l`, `Psi_l` | `1.5e-4 mM`, 6.81, 48.8 mM |
| `n_e`, `k_e`, `c_e`, `b_e` | 140.2, 5.3, 102.6, 42.9 mM |
| `pH_e`, `d_e` | 7.4, 1.9 mM |
| `omega_i0`, `omega_l/omega_i0` | 1.3 pL, 0.02 |
| Resting total flow | `0.24e-8 microlitre/min` |
| Maximum stimulated total flow | `3.2e-8 microlitre/min` |

The paper reports approximately 27.3% cell shrinkage at maximal stimulation, approximately 24% lower simulated flow for Ae4 knockout, and no significant simulated flow change for Ae2 knockout.

The prescribed calcium drive is only a figure, not an equation or data table. The text gives a resting value of 58 nM and stimulation from minute 6 to minute 12. The plotted curve is approximately 0.05 micromolar at rest, spikes to 0.30 micromolar, falls to 0.25 and then 0.20 micromolar, and returns to rest. Therefore no exact stimulated steady state or transient can be reconstructed from the article alone.

Relevant activity and conductance values are:

| Parameter | Published value and unit |
| --- | ---: |
| `G_Ae2` | 0.01807 fmol/s |
| `G_Ae4` | 0.66 amol/micrometre cubed |
| `k_+`, `k_-` | `1.92e-2`, `1.3e-5 mM^-4 s^-1` |
| `alpha_Nkcc1`, `alpha_NaK` | 2.15, 4.84 amol/micrometre cubed |
| `G_Nhe1` | 0.0305 fmol/s |
| `G_CaCC`, `G_CaKC` | 71.3, 30.4 nS |
| `g_Na^t`, `g_K^t` | 12.46, 0.9 nS |
| `P_a`, `P_b`, `P_t` | `4.32e-12`, `5.15e-11`, `2.6e-13 L^2 mol^-1 s^-1` |

## Numerical baseline inconsistencies

These checks use only the published equations and displayed values. They show why an exact baseline reproduction is presently impossible; they are not a proposed recalibration.

### 1. The reported osmolalities do not satisfy steady water balance

Section 2.5 reports

\[
\Omega_e=292.6,\qquad
\Omega_i=296.6,\qquad
\Omega_l=297.4\quad\text{mM}.
\]

Steady volume requires

\[
P_a(\Omega_l-\Omega_i)=P_b(\Omega_i-\Omega_e).
\]

With the published permeabilities, the right-hand side is about 59.6 times the left-hand side:

\[
\frac{P_b(4.0)}{P_a(0.8)}\simeq 59.6.
\]

This discrepancy is far larger than rounding of Table 1.

### 2. The water equations do not give the reported resting flow

Converting mM to mol/L, the displayed baseline gives

\[
q_a=4.32\times10^{-12}(0.8\times10^{-3})
=3.456\times10^{-15}\ \mathrm{L/s},
\]

\[
q_t=2.6\times10^{-13}(4.8\times10^{-3})
=1.248\times10^{-15}\ \mathrm{L/s}.
\]

Therefore

\[
q_{\mathrm{tot}}=q_a+q_t
=2.82\times10^{-7}\ \mathrm{microlitre/min},
\]

whereas section 3.1 reports

\[
2.4\times10^{-9}\ \mathrm{microlitre/min}.
\]

The discrepancy is approximately 118-fold.

### 3. The displayed bath values do not reproduce the reported bath osmolarity

Using `pH_e=7.4`, so that `h_e=10^(3-7.4) mM`, the printed Table 1 bath entries sum to approximately 292.90 mM, not 292.6 mM. This smaller discrepancy may reflect rounding or an unreported unrounded value, but it does not repair checks 1 or 2.

### 4. Acid/base equilibrium does not give the reported intracellular carbon dioxide

Appendix 8 states

\[
d=\frac{k_{-1}}{k_1}bh.
\]

Using `pH_i=6.91`, hence `h=1.2303e-4 mM`, together with `b=12.1 mM`, `k_1=11`, and `k_-1=2.6e4`, gives

\[
d\simeq 3.52\ \mathrm{mM},
\]

not the Table 1 value 6.6 mM. The printed reverse-rate unit is also `s^-1` even though the reverse term multiplies two concentrations.

Moreover, Table 8 gives `P_CO2` units of `s^-1`, so both the diffusion and buffer expressions are concentration/time. Equation (6) is an amount balance. A factor of cell volume, and possibly membrane geometry, is required before these terms can be added to the transporter amount fluxes.

### 5. Tight-junction current and reported flow use incompatible normalizations

The Table 1 potentials give `V_t=V_a-V_b=12.56 mV`. With the Table 2 conductances and published Nernst relations, the Na and K tight-junction currents correspond to a combined particle flux of about 2.30 fmol/s. The lumen balances would then imply a resting flow of roughly

\[
1.1\times10^{-6}\ \mathrm{microlitre/min},
\]

more than 400 times the reported resting flow. This calculation is conditional on interpreting the printed nS values as whole-cell conductances, exactly the interpretation the article appears to invite. It demonstrates that a missing geometry or normalization cannot be ignored.

### 6. The displayed baseline does not satisfy the apical QSS equation

Using the reported `58 nM` calcium value, the CaCC gate, conductance, concentrations, and `V_a=-50.24 mV` give an apical chloride current of approximately `-154.1 pA`. The Table 1 potentials and Appendix 10 tight-junction relations give

\[
I_{\mathrm{Na}}^t+I_K^t\simeq 221.9\ \mathrm{pA}.
\]

The apical QSS residual from (S9) is therefore about `-67.8 pA`, not zero. With the other displayed values fixed, calcium near `79.5 nM`, rather than `58 nM`, would close this one equation. This mismatch cannot be attributed to the unknown NaK current because that current does not appear in the apical QSS balance.

### 7. Transporter density and concentration scaling is unspecified

The intracellular equations explicitly require moles/time. Nevertheless:

- `G_Ae2` and `G_Nhe1` are total fmol/s;
- `G_Ae4`, `alpha_Nkcc1`, and `alpha_NaK` are densities in amol/micrometre cubed;
- the prose after Eq. (31) instead calls `G_Ae4` fmol and calls `J_Ae4` concentration/time;
- no membrane area or density-to-whole-cell conversion is given; and
- literal insertion of concentrations displayed in mM into the fourth-order Nkcc1 and Ae4 polynomials produces implausible scales or directions.

One part of the NKCC1 issue is source-resolvable. The 2013
[corrigendum to Palk et al.](https://www.sciencedirect.com/science/article/pii/S0022519312005668)
gives `a_1=157.55 s^-2`, `a_2=2.0096e7 M^-4 s^-2`,
`a_3=1.0306 s^-1`, and `a_4=1.3852e6 M^-4 s^-1`. For mM states,
the corresponding numerical `a_2` and `a_4` values are `2.0096e-5` and
`1.3852e-6`. At the Table 1 state the corrected bracket is approximately
`0.598 s^-1`; literal use of the 2018 table's large numbers with mM states
gives approximately `-14.51 s^-1` and reverses direction. This correction is
implemented as a regression diagnostic.

The correction does not supply the missing NKCC1 density-to-whole-cell
conversion. No equivalent explicit corrigendum was found for the NaK values,
and the Ae2 `K_B=10^4 mM` value remains source-ambiguous because it makes the
printed resting Ae2 bracket negative. The original implementation may also
have used implicit volume conversions, but these choices are not published.
Inferring them by forcing the baseline would be a new calibration, not a
reproduction.

### 8. The Ae4 Markov and simplified flux displays are not algebraically equivalent as written

As noted after (R2), the Markov steady flux preceding Eq. (31) contains a concentration-dependent denominator, while Eq. (31) uses a constant `G_Ae4` multiplying a bilinear forward-minus-reverse expression. The state equation printed as Eq. (29) also uses intracellular chloride and the intracellular complex in a gain term for the exterior-facing state. These are source-level blockers to reconstructing the four-state Ae4 mechanism. Using Eq. (31) is a declared phenomenological choice, not a verified elimination.

### 9. The baseline is partly the fitting target

The paper marks channel and tight-junction conductances, transporter densities/activities, all Ae4 parameters, `x_i`, `Psi_l`, and several concentrations as determined from the model. It does not publish the fitting objective, weights, bounds, parameter covariance, or unrounded fitted vector. Matching the displayed state after refitting these quantities would therefore not independently verify the equations.

## Reconstruction decision

The following conclusions are defensible now:

1. The steady stoichiometric identities (E1)--(E8) are exact consequences of the published balances.
2. Conservation resolves the volume equation as `d omega_i/dt = q_b-q_a`.
3. The potassium balance is redundant in the steady QSS system.
4. Subject to stated closure and unit assumptions, the model reduces to the six residuals (M14)--(M19) with `G_2` and `G_4` explicit.
5. The identities (E3)--(E5) expose which additional flux information would separate Ae2 and Ae4 at steady state.

The following conclusion is not defensible:

> The published 2018 baseline has been numerically reconstructed.

An executable reconstruction requires, at minimum:

1. a complete amount/current unit ledger, including density-to-cell scaling
   and source-backed or explicitly inferred NaK/Ae2 concentration conventions;
2. implementation of the conservation-resolved carbon-dioxide sign, an acid/base closure, and a printed-sign sensitivity check;
3. a corrected water normalization consistent with the baseline;
4. exact fixed lumen acid/base values and a numerical calcium input;
5. an unrounded, algebraically consistent initial state; and
6. a parameter vector that distinguishes published calibration targets from independently fixed inputs.

Until those items are obtained from a published correction, upstream source audit, or author clarification, any numerical implementation must be labelled an **independently calibrated surrogate based on the published topology**, not a reproduction. Symbolic identifiability results should state explicitly whether they use only (E1)--(E8) or also the conditional numerical closure (M14)--(M19).
