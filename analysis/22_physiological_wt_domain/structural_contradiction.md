# Task 22 necessary-condition certificate

**Failure class: `PROVEN_STRUCTURAL_CONTRADICTION` within the predeclared
Task 22 domain and unchanged production chassis.** This is a coupled
physiological/thermodynamic conflict in the carbon and alkalinity balances,
not an inference from local optimizer failure and not a general rejection of
pooled/no-slip AE4. All ten immutable root backgrounds give the same bounds.

## Production identities

Let $T_i,T_l$ be total inorganic carbon (TIC), $A_i,A_l$ total alkalinity,
$C_i,C_l,C_b$ dissolved CO2, $Q\geq0$ the production lumen outflow, $H$ the
inward NHE1 cycle flux, and $J_a,J_b$ the lumen-to-cell and bath-to-cell CO2
fluxes. The $R$ symbols below are **amount derivatives** in fmol/s, not
concentration derivatives. These identities follow by subtracting the
corresponding rows of `ModernFullModel.evaluate`:

$$R_{T_i}-R_{A_i}=J_b+J_a-H,$$
$$J_a=Q(A_l-T_l)+R_{A_l}-R_{T_l}.$$

The cell membrane has no TIC or alkalinity source. AE2 and pooled AE4
contribute equally to cell TIC and alkalinity, so they cancel in the first
identity. The paracellular bicarbonate source contributes equally to lumen
TIC and alkalinity, so it cancels in the second. Finite carbon speciation
does not create carbon or alkalinity. These identities were also checked
directly against production evaluations away from steady state.

## The luminal domain prevents a positive apical CO2 supply at rest

For carbon fractions $\alpha_0$ (CO2) and $\alpha_2$ (CO3),

$$A_l-T_l=T_l(\alpha_2-\alpha_0)+B_l^-+[OH^-]_l-[H^+]_l+\epsilon_{\rm spec}.$$

As pH increases, $\alpha_0$ decreases and $\alpha_2$ increases. The frozen
carbon pKa values are 6.1 and 10.3. At the maximum allowed lumen pH, 8.0,
$\alpha_2-\alpha_0<0$. With $T_l\geq1$ mM, lumen buffer amount
$10^{-9}$ fmol, volume at least 0.02 pL, and the inherited alkalinity
speciation residual bound $|\epsilon_{\rm spec}|\leq10^{-9}$ mM,

$$A_l-T_l\leq
 [\alpha_2(8)-\alpha_0(8)](1\ {\rm mM})
 +\frac{10^{-9}}{0.02}\ {\rm mM}
 +1000(10^{-6}-10^{-8})\ {\rm mM}+10^{-9}\ {\rm mM}
 <-0.00645\ {\rm mM}.$$

The unrounded upper bound is approximately -0.00645626719353 mM. This
bound covers the entire luminal pH/TIC/volume box, including the most
favorable buffer and speciation error; it is not a grid-search observation.

Since $Q\geq0$, the second production identity gives
$J_a\leq\varepsilon_{A_l}+\varepsilon_{T_l}$ for any accepted resting state.
This remains true at zero outflow, including lumen volumes below the
0.10 pL outflow dead volume.

## The cell balances force CO2 below bath

The fixed NHE1 law is

$$H=V_{\rm NHE1}\tanh\!\left(
\frac{\log(Na_b/Na_i)+\log(10)(pH_b-pH_i)}{2}\right).$$

The most permissive bounds are $Na_i\leq60$ mM and $pH_i\leq7.05$;
the bath has $Na_b=145$ mM and $pH_b=7.4$. The capacity lower bound
is 0.01 times the inherited NHE1 reference. For every root,
$H\geq5.11938969945\times10^{-5}$ fmol/s.

The inherited accepted residual allowances are:

| Production row | Absolute allowance (fmol/s) |
|---|---:|
| Cell TIC | $2\times10^{-9}$ |
| Cell alkalinity | $10^{-9}$ |
| Lumen TIC | $2\times10^{-9}$ |
| Lumen alkalinity | $10^{-9}$ |

The TIC rows use the inherited 0.02 fmol/s scale times the 1e-7 scaled
root tolerance; the two alkalinity rows use the 1e-9 omitted-row gate.
Using both production identities with all four errors in their most
permissive directions,

$$J_b\geq H_{\min}-(2+1+2+1)10^{-9}
 >5.118\times10^{-5}\ {\rm fmol/s}>0.$$

Production defines $J_b=P_b(C_b-C_i)$ with strictly positive bounded
basolateral CO2 permeability. Therefore $C_i<C_b$. This conclusion
survives the inherited numerical tolerances by a margin of thousands of
times the total four-row residual allowance.

## Positive pooled AE4 affinity is then impossible

The shared carbon equilibrium gives exactly

$$\frac{[HCO_3^-]_i}{[HCO_3^-]_b}
=\frac{C_i}{C_b}10^{pH_i-pH_b}<10^{7.05-7.4}.$$

The unchanged pooled law's dimensionless productive chloride-loading
affinity is

$$\mathcal A=log\frac{Cl_b}{Cl_i}
+\log\frac{Na_i+K_i}{Na_b+K_b}
+2\log\frac{[HCO_3^-]_i}{[HCO_3^-]_b}.$$

Using $Cl_i\geq47.10$ mM, $Na_i+K_i\leq60+200=260$ mM,
the frozen $Cl_b=126.16159280366088$ mM and $Na_b+K_b=150$ mM,

$$\mathcal A<
\log\frac{126.16159280366088}{47.10}
+\log\frac{260}{150}+2\log(10)(7.05-7.4)
=-0.0764726613570\ldots<-0.0764.$$

Task 22 requires $\mathcal A>0$ and positive AE4 chloride loading. The
strictly negative upper bound contradicts that requirement. The extrema
used above need not be jointly realizable: allowing them simultaneously
makes an upper bound more permissive, which strengthens the exclusion.

## Scope and independent checks

`src/modern_full_model/task22_certificate.py` evaluates these monotonic
analytic bounds with 60-digit decimal arithmetic and reports conservative
outward-rounded margins. `structural_certificate.json` records every
constant, residual allowance, bound and root background. The accompanying
tests verify the production carbon/alkalinity identities, carbon-fraction
monotonic bound, numeric margins and fail-closed behavior when the luminal
pH assumption is changed in an isolated test.

The proof does not require the capacity **upper** limits, the osmotic
screens or the anti-cancellation screen. Expanding only upper capacity
limits would not remove this particular contradiction. It does rely on
the fixed source topology/chemistry and the declared physiology, including
the very small inherited lumen buffer. It does not establish that pooled
AE4 cannot work in a different physiological domain or a different chassis.
No model term, bound, buffer or tolerance was changed to test alternatives.

The 27 bounded WT restorations cover every inherited root, all seven
Task 21 resting witnesses, and ten additional WT-only starts. Their failure
is preserved separately as numerical search evidence. The analytic proof,
not those optimizer outcomes, establishes that the accepted WT ensemble is
empty. Thus no Task 22 WT trajectory can be initialized from an admissible
resting state, and genotype evaluation remains prohibited.
