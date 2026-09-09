# Fixed non-AE4 comparison chassis

## Decision and scope

The primary comparison uses an independently reimplemented version of the
Task-11 **seven-state historical branch** as its fixed non-AE4 chassis. This is
the smallest available branch that satisfies all of the following:

1. its WT resting row is internally closed (`R11-03`);
2. it retains dynamic luminal Na and K, cell volume/height, and intracellular
   Na, K, Cl, and HCO3, so candidate cation stoichiometries reach every relevant
   balance;
3. AE2 knockout remains negligible under its own protocol (`R11-04`); and
4. AE4 can be replaced at one explicit four-species insertion boundary.

The implementation is
`src/ae4_mechanism_reconstruction/chassis.py`. It does not read, import, or
execute a historical file. The two membrane potentials are rederived as the
solution of a two-by-two linear QSS system rather than copied from the archived
closed-form expressions.

This is a **comparative chassis**, not a recovered 2018 physiological model.
Its absolute flow and time units remain undocumented, its non-AE4 equations
differ from the publication in several places, and its fitted resting row is
not independent biological validation. Dimensionless candidate-to-WT and
knockout-to-WT ratios under one fixed protocol are the defensible outputs.

## Provenance labels

- **Historical implementation evidence:** equation or number documented by
  Task 11 from the internally closed seven-state branch.
- **Published/public correction:** independently supported by the 2018 paper,
  Palk et al. (2010), or the Palk corrigendum.
- **New modelling decision:** a Task-12 interface, numerical safeguard, or
  comparison convention not asserted as historical biology.

Every retained numeric chassis value is historical implementation evidence.
The conservation volume sign and corrected NKCC concentration convention also
have public support. The swappable AE4 vector and immutable fingerprint are new
modelling/software decisions.

## State, sign, and dimensional conventions

The seven states are

\[
y=(N_l,K_l,H,N_i,K_i,C_i,B_i),
\]

where `N`, `K`, `C`, and `B` denote Na, K, Cl, and HCO3 concentrations,
respectively; `l` and `i` denote lumen and cell; and `H` is cell height. All
concentrations are mM, height is micrometres, calcium is numerically in
micromolar, and potentials are mV. Cell volume is recovered as

\[
V_i=\delta H\times10^{12}\quad\text{pL}.
\]

All transporter, channel, tight-junction, and amount-balance terms share the
historical internal amount/area/time scale. Water terms share its height/time
scale. The global factor `10000` defines the historical integration coordinate
but has no established physical-time conversion.

For AE4, positive cycle flux means basolateral chloride entry. Candidate source
terms are signed changes of intracellular amounts:

\[
s_4=(s_N,s_K,s_C,s_B).
\]

Electroneutrality requires

\[
s_N+s_K-s_C-s_B=0.
\tag{1}
\]

The chassis rejects a candidate violating (1) before simulation.

## Algebraic reductions

### Cellular and luminal electroneutrality

The eliminated intracellular proton concentration is

\[
h_i=B_i+C_i+X_i/H-N_i-K_i.
\tag{2}
\]

The lumen chloride concentration is

\[
C_l=N_l+K_l.
\tag{3}
\]

Thus the cellular and luminal electroneutrality residuals are zero to floating
point precision. Equation (2) is historical implementation evidence and is an
algebraic consequence of the published charge constraint when `X_i/H` is the
fixed impermeant-charge concentration. Equation (3) neglects the nearly
cancelling luminal H/HCO3 pair and is the historical reduction.

### CO2 quasi-steady closure

The historical branch eliminates CO2 using

\[
D_i=
\frac{p_{CO2}(D_l+D_e)-g_B k_- B_i h_i}{a_D}.
\tag{4}
\]

Here `a_D` is the fitted denominator stored by the historical parameter pack.
This is historical implementation evidence, not conservation- or
thermodynamically validated chemistry. It replaces the paper's dynamic CO2
balance and retains the paper's disputed aggregate sign.

## Fixed non-AE4 fluxes

### Na/K ATPase

The source convention uses molar concentrations inside the kinetic factor:

\[
J_{NaK}=A_{NaK}r
\frac{(sK_e)^2(sN_i)^3}{(sK_e)^2+\alpha_1(sN_i)^3},
\qquad s=10^{-3}.
\tag{5}
\]

One cycle removes three intracellular Na, adds two intracellular K, and carries
one net positive charge outward.

### Channel gates, Nernst terms, and QSS potentials

The historical gates are

\[
p_C=\frac{1}{1+(K_C/[Ca])^{\eta_C}},\qquad
p_K=\frac{1}{1+(K_K/[Ca])^{\eta_K}}.
\tag{6}
\]

Define

\[
v_C=RT/F\log(C_l/C_i),\quad
v_K=RT/F\log(K_e/K_i),
\]

\[
v_N^t=RT/F\log(N_l/N_e),\quad
v_K^t=RT/F\log(K_l/K_e).
\tag{7}
\]

With `g_C=G_C p_C`, `g_K=G_K p_K`, and `g_t=g_N^t+g_K^t`, the independently
derived QSS system is

\[
\begin{pmatrix}
g_C+g_t&-g_t\\
-g_t&g_K+g_t
\end{pmatrix}
\binom{V_a}{V_b}
=
\binom{-g_Cv_C+g_N^tv_N^t+g_K^tv_K^t}
{-FJ_{NaK}+g_Kv_K-g_N^tv_N^t-g_K^tv_K^t}.
\tag{8}
\]

The signed current-equivalent fluxes are

\[
j_C=g_C(V_a+v_C)/F,\quad
j_K=g_K(V_b-v_K)/F,
\]

\[
j_N^t=g_N^t(V_a-V_b-v_N^t)/F,\quad
j_K^t=g_K^t(V_a-V_b-v_K^t)/F.
\tag{9}
\]

Outward chloride has `j_C<0`; outward potassium has `j_K>0`. The two regression
identities are

\[
j_C+j_N^t+j_K^t=0,
\qquad
J_{NaK}+j_K-j_N^t-j_K^t=0.
\tag{10}
\]

### Water and volume

\[
q_a=b_a\{2(N_l+K_l-N_i-K_i-h_i)-D_i+\Psi_l\},
\tag{11}
\]

\[
q_b=b_b\{2(N_i+K_i+h_i)+D_i-I_e\},
\tag{12}
\]

\[
q_t=b_t\{2(N_l+K_l)+\Psi_l-I_e\},
\qquad Q=q_a+q_t,
\tag{13}
\]

\[
\dot H=q_b-q_a.
\tag{14}
\]

The sign in (14) is fixed by conservation and agrees with Palk et al. (2010),
the collected 2018 system, and every historical branch. The numeric water
coefficients are historical fitted values and are not a common-unit conversion
of the published permeability table.

At steady volume, equations (11)-(13) imply the exact chassis reduction

\[
Q=\left(\frac{b_ab_b}{b_a+b_b}+b_t\right)
\{2C_l+\Psi_l-I_e\}.
\tag{15}
\]

### NKCC1

\[
J_{NKCC}=A_{NKCC}
\frac{a_1-a_2N_iK_iC_i^2}{a_3+a_4N_iK_iC_i^2}.
\tag{16}
\]

The values of `a2` and `a4` are the public-corrigendum M-to-mM conversion and
historical implementation convention. A cycle adds 1 Na, 1 K, and 2 Cl.

### AE2

\[
J_2=G_2\left[
\frac{C_e}{C_e+K_C}\frac{B_i}{B_i+K_B}
-\frac{C_i}{C_i+K_C}\frac{B_e}{B_e+K_B}
\right].
\tag{17}
\]

The retained historical values are `B_e=21 mM` and `K_B=10^-4 mM`. The first
is different from both the historical water bath and the published bath; the
second is likely a missing-exponent-sign issue in the paper but remains
historical implementation evidence here. AE2 knockout changes only the
scenario multiplier on (17), never a chassis number.

### NHE1 and buffer

\[
J_{NHE}=G_{NHE}\left[
\frac{N_e}{N_e+K_N}\frac{h_i}{K_H+h_i}
-\frac{N_i}{N_i+K_N}\frac{h_e}{K_H+h_e}
\right],
\tag{18}
\]

\[
J_B=g_B(k_+D_i-k_-B_ih_i).
\tag{19}
\]

The first-power proton dependence in (18) and fitted scaling in (19) are
historical chassis choices; the 2018 paper prints a squared NHE proton factor.

## Differential balances and AE4 insertion

The lumen and volume equations are

\[
\dot N_l=j_N^t-QN_l,\qquad
\dot K_l=j_K^t-QK_l,\qquad
\dot H=q_b-q_a.
\tag{20}
\]

Before product-rule dilution, the four cell amount-rate numerators are

\[
R_N=J_{NKCC}-3J_{NaK}+J_{NHE}+s_N,
\tag{21}
\]

\[
R_K=J_{NKCC}+2J_{NaK}-j_K+s_K,
\tag{22}
\]

\[
R_C=2J_{NKCC}+J_2+j_C+s_C,
\tag{23}
\]

\[
R_B=J_B-J_2+s_B.
\tag{24}
\]

The concentration equations are

\[
\dot X_i=(R_X-\dot H X_i)/H,
\quad X\in\{N,K,C,B\}.
\tag{25}
\]

Equations (21)-(24) are the **only direct AE4 entry points**. With (1) and the
QSS current identities, differentiating cellular electroneutrality gives the
eliminated proton balance exactly:

\[
\frac{d(Hh_i)}{dt}=J_B-J_{NHE}.
\tag{26}
\]

AE4 has no direct membrane-current term because every admitted cycle is
electroneutral. Its indirect paths are:

1. altered Na/K/Cl/HCO3 change `h_i` through (2);
2. `h_i` and HCO3 alter CO2, NHE1, AE2, and buffer fluxes;
3. Na/K/Cl alter Nernst potentials, QSS voltages, channel and tight-junction
   fluxes;
4. ionic and acid/base changes alter water flux, cell height, lumen ions, and
   secretion; and
5. the changed volume dilutes every cell concentration through (25).

## Complete retained parameter ledger

All values in this table are fixed across AE4 candidates. “Historical” means
the exact number is retained as historical implementation evidence; it does
not mean a validated biological constant.

| Code name | Value | Role / unit convention | Provenance |
|---|---:|---|---|
| `buffer_scale` | 0.199652558943136 | fitted buffer multiplier | Historical |
| `co2_transport_scale` | 1970 | CO2 aggregate scale | Historical |
| `buffer_forward` | 0.132 | forward buffer rate | Historical scaled value |
| `buffer_reverse` | 312 | reverse buffer rate | Historical scaled value |
| `impermeant_charge_amount` | 2439.50353087137 | `X_i`, mM·micrometre in height form | Historical calibration |
| `co2_elimination_denominator` | 3939.9145458622197 | denominator in (4) | Historical calibration |
| `co2_external_lumen_sum` | 13.2 | `CO2_l+CO2_e`, mM | Historical bath |
| `nak_capacity` | 0.00138020583534274 | Na/K-pump scale | Historical calibration |
| `nak_rate` | 1,305,000 | Na/K kinetic rate | Historical/source family |
| `k_e` | 5.3 | external K, mM | Historical; matches paper |
| `nak_alpha` | 0.641 | pump saturation coefficient | Historical/source family |
| `mm_to_molar` | 0.001 | explicit M conversion | Historical convention |
| `thermal_voltage_mv` | 26.7137302360965 | `RT/F`, mV | Historical; physical constants match paper |
| `faraday` | 96485.3365 | C/mol | Historical; matches paper |
| `cacc_half_activation` | 100.26 | historical CaCC gate constant | Historical, differs from paper |
| `kcc_half_activation` | 100.3 | historical K-channel gate constant | Historical, differs from paper |
| `cacc_hill` | 1.49 | CaCC Hill exponent | Historical, differs from paper |
| `kcc_hill` | 1.7 | K-channel Hill exponent | Historical, differs from paper |
| `tight_na_conductance` | 357.863312079523 | pre-scaled tight Na conductance | Historical calibration |
| `tight_k_conductance` | 25.893178360928 | pre-scaled tight K conductance | Historical calibration |
| `cacc_conductance` | 20,468,119.8347071 | pre-scaled CaCC conductance | Historical calibration |
| `kcc_conductance` | 87,495,178.2340436 | pre-scaled K conductance | Historical calibration |
| `nkcc_capacity` | 0.0063812 | NKCC scale | Historical calibration |
| `nkcc_forward` | 157.55 | `a1` | Historical; public source-corrected |
| `nkcc_reverse` | 2.0096e-5 | `a2` for mM states | Historical; public corrigendum |
| `nkcc_denominator_forward` | 1.0306 | `a3` | Historical; public source family |
| `nkcc_denominator_reverse` | 1.3852e-6 | `a4` for mM states | Historical; public corrigendum |
| `nhe1_capacity` | 0.674737308094733 | NHE1 scale | Historical calibration |
| `ae2_capacity` | 0.399084872489615 | AE2 scale | Historical calibration |
| `water_apical` | 5.79346121973713e-4 | `b_a` | Historical calibration |
| `water_basolateral` | 5.45083716910752e-5 | `b_b` | Historical calibration |
| `water_paracellular` | 7.31969562708724e-6 | `b_t` | Historical calibration |
| `lumen_impermeant_osmolyte` | 48.8001357236744 | `Psi_l`, mM | Historical calibration |
| `external_osmotic_sum` | 288.1 | folded external osmotic term, mM | Historical bath convention |
| `area_to_litre_per_height` | 4.52961672473868e-14 | geometry conversion | Historical geometry |
| `cl_e` | 102.6 | external Cl, mM | Historical; matches paper |
| `na_e` | 140.2 | external Na, mM | Historical; matches paper |
| `hco3_e_transport` | 21 | external exchanger HCO3, mM | Historical hard-coded value |
| `ae2_cl_half` | 5.6 | AE2 Cl half-saturation, mM | Historical/source family |
| `ae2_hco3_half` | 1e-4 | AE2 HCO3 half-saturation, mM | Historical convention |
| `nhe1_na_half` | 15 | NHE1 Na half-saturation, mM | Historical/source family |
| `nhe1_h_half` | 4.5e-4 | NHE1 H half-saturation, mM | Historical/source family |
| `h_e` | 3.8904514499428e-5 | external H, mM | Historical pH 7.41 |
| `rhs_scale` | 10,000 | hidden integration-coordinate scale | Historical; physical time unresolved |

The SHA-256 fingerprint of this complete non-AE4 vector is

`55885266c44ecbbfe86be622c6ebd08b77151ee6ef7dfd5c8bbc8d4250b58e50`.

The fingerprint is stored in every simulation result. A scenario can multiply
AE4 or AE2 activity for genotype definition, but cannot override a chassis
parameter.

The comparison harness also exposes an independent dimensionless
beta-adrenergic/cAMP/PKA input, `u_PKA` in `[0,1]`. Its default is exactly zero
and only explicitly stimulus-aware AE4 candidates read it; calcium is never
used as its proxy. The PKA protocol is an experimental input convention, not a
non-AE4 parameter and therefore not part of the chassis fingerprint. No
physical-time conversion is inferred from the historical integration
coordinate.

## AE4-specific reference law kept outside the chassis

The default negative-control object `HistoricalSodiumOnlyAE4` has
`capacity=18.9502607221547`, `k_forward=0.0192`, and
`k_reverse=0.00159305881398359`. These are explicitly AE4-specific and are not
part of the non-AE4 fingerprint. It returns

\[
s_4=(-J_4,0,+J_4,-2J_4),
\]

where

\[
J_4=G_4[k_+C_eB_i^2N_i-k_-C_i(21)^2N_e].
\tag{27}
\]

All other candidates implement the same callable contract but may return a
different electroneutral four-species vector.

## WT closure demand and the published-partition incompatibility

Setting AE4 to exactly zero at the historical WT row, while retaining every
non-AE4 value, gives the required AE4 amount-source vector

\[
s_4^{required}=
(-0.0388833386254,\;-1.67\times10^{-16},\;
+0.0388833386765,\;-0.0777666773911).
\tag{28}
\]

The small charge residual, `8.92e-11`, comes from the rounded stored baseline;
the full historical law closes the raw seven-state RHS to
`4.94e-12`. This demand uses **only WT closure** and no knockout observation.

At the same cycle rate, the printed 2018 intracellular mole-fraction partition
uses

\[
\alpha_N=25/145,\qquad \alpha_K=120/145,
\]

and therefore inserts

\[
s_4^{published}=(-\alpha_NJ_4,-\alpha_KJ_4,+J_4,-2J_4).
\]

Relative to (28), the amount-balance residual is

\[
(+0.03217931472,-0.03217931477,0,3.81\times10^{-11}).
\tag{29}
\]

Thus the published cation allocation cannot simply be dropped into the fitted
historical resting row. It leaves Cl/HCO3 closure unchanged but breaks Na/K
closure by equal-and-opposite terms. Candidate comparison must recalibrate only
AE4-specific capacity under a common WT objective and report any irreducible
four-row closure error; it must not retune Na/K pump, NKCC, channels, or baths.

## AE4-specific versus chassis differences from the 2018 publication

| Difference | Classification | Treatment in Task 12 |
|---|---|---|
| Na-only historical AE4 versus pooled/partitioned Na+K | AE4-specific | Swapped candidate law and source vector |
| Historical folded reverse coefficient and 21 mM reverse HCO3 | AE4-specific kinetics, with bath convention shared by AE2 | Candidate parameters vary only when source-justified; chassis bath stays fixed in the primary comparison |
| Missing four-state Markov denominator/derivation | AE4-specific | Candidate family reconstructs empirical and reversible alternatives explicitly |
| AE4 capacity/geometry scale | AE4-specific | WT-only activity normalization; KO held out |
| HCO3_i 10 rather than 12.1 mM | Chassis | Fixed across candidates |
| External HCO3/CO2 split and acid/base QSS | Chassis | Fixed across candidates; limitation reported |
| Historical channel gates/conductances | Chassis | Fixed across candidates |
| First-power NHE1 proton term | Chassis | Fixed across candidates |
| Fitted water coefficients and undocumented absolute flow scale | Chassis | Fixed; only dimensionless flow ratios interpreted |
| Seven-state reduction and CO2 QSS | Chassis | Fixed across candidates |
| 0.05-to-0.55 calcium step and hidden time multiplier | Chassis/protocol | Fixed historical comparison protocol; not called the primary experimental or 2018 protocol |

## Numerical protocol and verified regressions

The historical comparison protocol starts from

`(118.7, 5.6, 28.7, 25, 120, 50, 10)`,

uses calcium `0.05` through coordinate time 100 and `0.55` thereafter, splits
the solve at the discontinuity, and retains relative and absolute tolerances
`1e-6`. BDF is the exact Task-11 lineage solver substitute. Radau is permitted
for high-throughput candidate screening only after endpoint cross-checks; for
historical AE4 knockout its endpoint flow differs from BDF by less than
`1e-8` relative.

Verified independent chassis results are:

| Check | Result |
|---|---:|
| WT maximum raw baseline RHS | `4.94e-12` |
| WT resting `q_total` | `0.000531415476760683` chassis units |
| WT resting `(V_a,V_b)` | `(-50.24,-62.8)` mV |
| Historical AE2-KO endpoint flow / WT | `0.999052` |
| Historical AE4-KO endpoint flow / WT | `0.998556` |
| Cell/lumen electroneutrality at WT | below `1.5e-14` |
| QSS current residuals at WT | below `4.1e-17` |
| Proton-balance identity residual at WT | below `2.8e-17` |
| Luminal-Cl redundancy residual at WT | below `1.1e-17` |

Tests are in `tests/test_ae4_chassis.py`. The test harness also rejects an
electrogenic AE4 source, freezes the non-AE4 fingerprint, verifies the demand
vector and published-partition mismatch, and reproduces the historical
knockout ratios without fitting to either knockout.

## Interpretation boundary

This chassis answers a controlled counterfactual:

> With every non-AE4 historical equation and number fixed, can a
> source-justified AE4 law, calibrated only to WT closure and independent
> transport evidence, predict the held-out knockout phenotype?

It cannot establish absolute physiological flow, physical response time, or
identity with the publication-generating model. If all defensible AE4 source
vectors fail, the exact balance identities and fixed fingerprint make that a
localization result: the missing effect lies in an indirect non-AE4 pathway or
in a chassis convention, rather than in an untested retuning of the same
comparison system.
