# Resolution of the Task 10 reconstruction blockers

## Decision rule and recomputation basis

The Task 10 diagnostics substituted the published Table 1 row into the
published equations. This audit separately evaluates the archived outer model
at its own initial row:

`(Na_l,K_l,h_i,Na_i,K_i,Cl_i,HCO3_i) = (118.7,5.6,28.7,25,120,50,10)`.

All historical numbers below were recomputed from the immutable MATLAB
expressions and the values in `Par.mat`, before the final global multiplier of
the right-hand side. They show whether the historical implementation is
internally closed. They do not turn historical implementation evidence into a
published correction. The full machine-readable calibration, residual, and
parameter-concordance record is in
`results/11_forensic_reconstruction/reproduction_summary.json` and
`reproduction_residuals.csv`.

The required classification phrases are used exactly:

- `resolved by historical implementation convention`;
- `publication typo resolved by code and conservation`;
- `implementation-specific calibration not stated clearly in paper`;
- `historical code also inconsistent`; and
- `still unresolved`.

Where a blocker has separable parts, the table gives a primary overall status
and the detailed section identifies any component that is resolved.

## Blocker disposition

| Task 10 blocker | Published-equation result | Historical recomputation / evidence | Classification |
| --- | --- | --- | --- |
| Opposite signs for `d omega_i/dt` | Eq. (9) prints `q_a-q_b`; the collected system and conservation require `q_b-q_a`. | Every archived dynamics file uses `q_b-q_a`. The resting zero set is sign-invariant. | **publication typo resolved by code and conservation** |
| Approximately `59.6 x` steady water-balance mismatch | Published values give `q_b/q_a=59.606`. | Historical osmotic brackets are `Delta_a=0.7997659713`, `Delta_b=8.5003697524`; fitted `b1/b2=10.62857143=Delta_b/Delta_a`, giving `q_a=q_b=0.000463341314` and ratio `0.9999999999998`. | **implementation-specific calibration not stated clearly in paper** |
| Published resting-flow mismatch (about `117.6 x`) | Literal published permeabilities give `2.82e-7 microlitre/min`, not `2.4e-9`. | Historical `q_total=0.000531415477` has no declared physical unit. The natural area conversion gives `1.4443e-9 microlitre/min`, only 60.18% of the published resting flow, and the code never applies that conversion. | **still unresolved** |
| Printed/recomputed bath osmolarity | The Table 1 entries sum to `292.9000398` mM, not the stated 292.6 mM. | Historical water code instead uses `Ie=288.1` mM, omitting bath H⁺/CO₂ and using HCO₃ₑ=40 mM; it does not recover either published number. | **historical code also inconsistent** |
| Tight-junction-current-implied flow | Published tight currents imply `1.11008e-6 microlitre/min`, `462.53 x` the reported resting flow. | Historical `GtNa` and `GtK` are solved so `(JtNa+JtK)=q_total(Na_l+K_l)` in internal units; height division nearly recovers Table 2 conductances, but no physical current/flow normalization is stated. | **implementation-specific calibration not stated clearly in paper**; absolute flow **still unresolved** |
| Approximately `67.8 pA` apical electrical residual | Published Table 1 voltages, conductances, gate and Ca=58 nM give about `-67.839 pA`. | Historical calibration gives `V_a=-50.24`, `V_b=-62.8` mV and an apical molar-equivalent residual `-4.03e-17` in its own units. It uses different gates, calcium, and transformed conductances. | **implementation-specific calibration not stated clearly in paper** |
| Basolateral current/NaK closure | Published current-to-cycle conversion and fitted density are not explicit. | Historical QSS uses `I_NaK=F J_NaK`; at rest `J_K=0.0453088161`, `J_NaK=0.0207461277`, and tight flux sum `0.0660549438`, leaving `+4.03e-17`. | **resolved by historical implementation convention** internally; published scale is **still unresolved** |
| Intracellular acid/base mismatch (`CO2_eq=3.51857` versus 6.6 mM) | The published equilibrium relation at Table 1 does not give its CO2 row. | The historical row uses HCO3=10 mM and algebraic `J_CO2=J_Buffer` to force CO2=`6.60012370`; its own chemical-equilibrium value is only `2.90790800` mM, a factor `2.26972` lower. | **historical code also inconsistent** |
| CO2 membrane-flux sign | The printed aggregate is the negative of the sum of the two individually defined inward fluxes. | The code adopts the printed aggregate sign and coefficient `1.99997`; it therefore confirms an implementation branch but does not resolve conservation. | **historical code also inconsistent** |
| Buffer/CO2 rate and amount units | Published `kminus1` and `P_CO2` units do not close the amount balances. | Code multiplies both buffer rates by `0.012`, then by fitted `gb=0.19965256`, and uses `P_CO2=1970` rather than `1.97e-13`; no documented common conversion exists. | **implementation-specific calibration not stated clearly in paper**; dimensions **still unresolved** |
| AE4 fourth-order M-versus-mM ambiguity (`10^12`) | The paper labels rates mM-based but gives inconsistent concentration/time, density, and amount prose. | Historical code unequivocally uses mM states. It converts published volumetric activity by height: `18.95026072/28.7=0.66028783 amol/um^3`. However it implements a Na-only law and an effective `kminus=0.00159306`, not the published Na+K law and `1.3e-5`. | Concentration/geometry component **resolved by historical implementation convention**; published AE4 model overall **still unresolved** |
| AE4 monovalent-cation allocation | Published Eqs. (2)–(3) split Ae4 between Na and K. | All archived variants put the full Ae4 cycle in Na and none in K. At rest the omitted K share is `0.0321793147` historical flux units. | **historical code also inconsistent** with the published model |
| NKCC1 coefficient units/direction | Literal 2018-table use gives about `-14.5077 s^-1`; the Palk corrigendum gives about `+0.598095 s^-1`. | `Par.mat` contains the corrected mM-state numbers `2.0096e-5`, `1.3852e-6`, reproducing `+0.598095 s^-1`. | **publication typo resolved by code and conservation** (also independently resolved by the Palk corrigendum) |
| NKCC1 density-to-cell scale | Paper gives 2.15 amol/um3 without an executable geometry rule. | Code uses `aNkcc1=0.0063812`; combining it with height and a `10^4` factor gives 2.2234, not exactly 2.15, and the factor is not documented. | **implementation-specific calibration not stated clearly in paper** |
| NaK concentration and density scaling | Paper labels `r` and `alpha1` mM-based and gives density 4.84 amol/um3. | Code explicitly converts mM states to M with `s=0.001` and uses fitted `aNaK=0.00138020584`; no source-backed mapping of the fitted scale is present. | M basis **resolved by historical implementation convention**; whole-cell scale **still unresolved** |
| NaK reaction stoichiometry display | Appendix 2's displayed reaction produces only one intracellular K, contrary to the stated 3 Na:2 K pump cycle. | Every archived K balance uses `+2 J_NaK`; the Na/K and current balances close only with the net +1 outward charge of a 3:2 cycle. | **publication typo resolved by code and conservation** |
| Ae2 `K_B` and direction | Printed `K_B=10^4` mM gives a small reverse resting factor (`-0.00269545` in Task 10). | Both MATLAB parameter layers use `1e-4` mM and generate forward Ae2 flux; however, the later historical TeX repeats `10^4`, and the upstream source was not recovered here. | Historical branch **resolved by historical implementation convention**; publication error **still unresolved** |
| Ae4 Markov-to-Eq. (31) reduction | The printed rational Markov steady flux has a state-dependent denominator; Eq. (31) replaces it by a constant activity, and Eq. (29) has an interior/exterior index conflict. | No archived MATLAB file implements the four Markov states or supplies the missing approximation; all use only a different Na-only Eq. (31)-type polynomial. | **still unresolved** |
| Current versus particle-flux sign | Paper alternates `I` and `J` and sometimes divides by `F` in an expression already called a flux. | Code consistently stores channel/tight terms as signed `I/F`, adds negative Cl `I/F` to the cell Cl row, uses `-I/F` for luminal Cl, and converts NaK cycle flux to current with `F`. | **resolved by historical implementation convention** for signs; absolute units **still unresolved** |
| Missing membrane area / volume convention | Published model gives volume but not explicit areas. | Historical geometry gives apical area `36.09256 um2`, basolateral `54.49977 um2`, mean `45.29617 um2`, height 28.7 um, and exactly 1.3 pL. Ae2, Nhe1, Ae4 and tight-junction published values are recovered to less than 0.3% by area/height transformations. | **resolved by historical implementation convention** for those families; other families remain **still unresolved** |
| Luminal Cl and acid/base closure | Published ten-state model has no luminal acid/base ODE and Eq. (16) drops the Cl subscript. | Code enforces `Cl_l=Na_l+K_l`, making the Cl row redundant under apical QSS; it omits explicit luminal acid/base variables. | Cl subscript/reduction **publication typo resolved by code and conservation**; general lumen acid/base dynamics **still unresolved** |
| Intracellular H and CO2 closure | Paper reports dynamic Eqs. (5)–(6). | H elimination follows from electroneutrality plus QSS charge balance; CO2 is additionally forced to steady balance. | H row **resolved by historical implementation convention**; CO2 transient **historical code also inconsistent** |
| Missing numerical calcium input | Paper gives a qualitative curve, rest 58 nM, stimulation during minutes 6–12. | Historical code supplies a different single step: 50 nM to 550 nM at `t>100`. | **still unresolved** for the published input |
| Solver and time scale | Paper names `ode15s` but not tolerance/initialization. | Code uses `ode15s`, `RelTol=AbsTol=1e-6`, but multiplies the full RHS by `10^4`; the five-state variant uses `10^3`. | Solver tolerance **resolved by historical implementation convention**; physical time **still unresolved** |
| Knockout implementation | Paper says set `G_Ae2` or `G_Ae4` to zero. | `Salivary_ex.m` explicitly sets `g2=0`; the function API permits `g4=0`, but no archived driver runs all WT/Ae2-KO/Ae4-KO cases. | Definition **resolved by historical implementation convention**; complete reproduction **still unresolved** |
| Baseline is partly a fitting target | The paper states that conductances/densities were chosen to require a physiological steady state but gives no objective, weights, bounds, or unrounded vector. | `Original/Saliva_Ae4.m` explicitly solves conductances/activities from the target row; all 45 `Par.mat` fields derive from that calibration to at most `2.51e-15` relative discrepancy, and its 12 residuals close to `1.39e-17`. | **implementation-specific calibration not stated clearly in paper**; independent validation **still unresolved** |

## Numerical details

### Water: why `59.6 x` disappears

At the historical resting row, electroneutrality and the algebraic CO2 solve
give

\[
[H^+]_i=1.2302687701\times10^{-4}\ \mathrm{mM},\qquad
[CO_2]_i=6.6001236986\ \mathrm{mM}.
\]

The archived reduced water brackets are

\[
\Delta_a=0.7997659713,\quad
\Delta_b=8.5003697524,\quad
\Delta_t=9.3001357237.
\]

With `b1=5.7934612197e-4`, `b2=5.4508371691e-5`, and
`b3=7.3196956271e-6`,

\[
q_a=4.6334131397\times10^{-4},\quad
q_b=4.6334131397\times10^{-4},\quad
q_t=6.8074162788\times10^{-5}.
\]

The volume residual is `-8.51e-17` and `q_b/q_a=0.9999999999998`.
This closure is deliberate: `Original/Parameters.m` constructs the three `b`
values with independent factors 74.4, 7, and 0.94. It is not a common-unit
conversion of the published permeabilities. In fact, the historical ratio
`b1/b2=10.6286` has the opposite hierarchy from published
`P_a/P_b=0.08388`. Therefore the previous mismatch was produced by correctly
using the published table; it disappears only after replacing that table with
an implementation-specific calibration.

The historical total is `q_total=5.3141547676e-4` in an unlabeled code unit.
If it is interpreted as a height velocity and multiplied by the recovered mean
area (`delta=4.5296167247e-14 L/um`), the result is
`1.44426506e-9 microlitre/min`, versus the paper's
`2.4e-9 microlitre/min`. Since the scripts only plot `q_total/q_total(1)` and
never declare this conversion, absolute flow remains uncertified.

### Bath sum and tight-junction flow normalization

Direct summation of the published Table 1 bath species, including the proton
concentration implied by pH 7.4, gives `292.900039811` mM rather than the
stated 292.6 mM. This `0.30004` mM discrepancy is small enough to be rounding
or an unreported fitted value, but the historical code does not identify which:
its water expression uses `Ie=140.2+5.3+102.6+40=288.1` mM and deliberately
omits bath CO2 and H. The historical value therefore cannot be used as a
correction to the published total.

The printed tight-junction currents provide an independent absolute-flow
failure. Their sum at Table 1 is `221.8881092 pA`, or
`2.299708093e-15 mol/s` after division by Faraday's constant. Dividing by
`Na_l+K_l=124.3 mM` gives

\[
q=1.850127187\times10^{-14}\ \mathrm{L/s}
=1.110076312\times10^{-6}\ \mathrm{microlitre/min},
\]

which is 462.532 times the reported resting flow. The historical calibration
forces `JtNa+JtK=q_total(Na_l+K_l)` exactly in its internal units. It does so
with `GtNa=357.8633121` and `GtK=25.89317836`; division by height produces
12.4691 and 0.9022, close to the printed 12.46 and 0.9. This is strong evidence
of a hidden normalization, but without a declared current, area, and time
conversion it does not certify the absolute published flow.

### Electrical closure: why `67.8 pA` disappears

At historical Ca=`0.05 microM`, the archived gates and calibrated
conductances give

| Historical term | Value in signed `I/F` or cycle-flux units |
| --- | ---: |
| CaCC | `-0.0660549437613` |
| Tight Na | `+0.0630790170915` |
| Tight K | `+0.00297592666986` |
| CaKC | `+0.0453088160576` |
| NaK | `+0.0207461277037` |

Thus

\[
-J_{CaCC}-J_{Na}^t-J_K^t=-4.03\times10^{-17},
\]

and

\[
-J_{CaKC}-J_{NaK}+J_{Na}^t+J_K^t=+4.03\times10^{-17}.
\]

The analytic voltage solve returns the published voltages to machine
precision. But `Original/Saliva_Ae4.m` computes the conductances from those
same target voltages and flows. The code also changes both channel gates and
uses enormous pre-scaled conductances. Hence this is exact historical internal
closure, not evidence that the published Table 2 values satisfy the published
QSS equation. Task 10's `-67.839 pA` result remains the correct audit of the
printed table.

### Acid/base closure: numerical match without chemical equilibrium

The historical CO2 algebra is the steady equality between its chosen printed-
sign aggregate and the scaled buffer. It yields

\[
[CO_2]_i=6.6001236986\ \mathrm{mM},\quad
J_{Buffer}=0.09730516118,\quad
J_{CO2}=0.09730516119.
\]

The equilibrium relation printed in Appendix 8 would instead give, using the
historical HCO3=10 mM and pH=6.91,

\[
[CO_2]_{i,eq}=\frac{312}{0.132}(10)(1.23026877\times10^{-4})
=2.907908002\ \mathrm{mM}.
\]

The code value is 2.269715 times larger. The historical implementation thus
chooses a dynamic-buffer steady flux and does not impose the paper's chemical
equilibrium. It also changes the published HCO3 row from 12.1 to 10 mM and the
impermeant concentration from 82.8 to 85.000123 mM. This explains how the
historical row was made internally steady, but it does not resolve the
published acid/base contradiction.

### AE4: unit clues do not rescue model identity

Two useful historical conventions are recoverable:

1. concentrations are inserted as mM into `kplus=1.92e-2`; and
2. the volumetric Ae4 activity is converted to an areal coefficient by cell
   height: `18.95026072/28.7=0.66028783 amol/um3`, matching the published 0.66.

Those clues remove the `10^12` M/mM branch and identify a plausible geometry
conversion for the historical law. They do not recover the paper's law. The
archived turnover is

\[
\phi_4^{hist}=0.002051862990,
\qquad J_4^{hist}=18.95026072\phi_4^{hist}=0.03888333862,
\]

using Na alone, external bicarbonate 21 mM, and effective
`kminus=0.001593058814`. The paper uses Na+K, external bicarbonate 42.9 mM,
and `kminus=1.3e-5`. Moreover, the historical cation balances allocate all of
`J4` to Na. Restoring the published fractions at the same point allocates
`0.00670402390` to Na and `0.03217931472` to K, producing equal-and-opposite
residuals relative to the calibrated historical balance. This is a substantive
model difference, not a unit conversion.

The archive also does not repair the paper's internal Markov derivation. The
rational steady Markov expression contains a state-dependent denominator,
whereas Eq. (31) uses a constant activity, and the gain term printed in Eq.
(29) refers to an interior complex in an exterior-facing equation. No archived
file carries the four Markov states or the missing approximation. The archived
polynomial is therefore evidence for the expression actually used in this
historical branch, not an algebraic proof of published Eq. (31).

### NKCC and NaK

For Nkcc1, the code and the upstream Palk corrigendum agree. At
`Na_i=25`, `K_i=120`, `Cl_i=50` mM,

\[
\frac{157.55-(2.0096\times10^{-5})Na_iK_iCl_i^2}
{1.0306+(1.3852\times10^{-6})Na_iK_iCl_i^2}
=0.5980945042\ \mathrm{s}^{-1}.
\]

Literal use of the large printed 2018 numbers with mM states gives
`-14.50765232 s^-1`. This blocker is genuinely resolved at the coefficient
level. The fitted density/whole-cell scale is not.

For NaK, the archived code uses `s=0.001` to convert the mM states before the
kinetic expression, giving a turnover factor `15.03118388` at rest. Literal
use of the same numerical constants with mM states would give
`5.702796786e7`, a factor `3.793977e6` larger. The historical M convention is
clear, but unlike Nkcc1 it is not backed by a published corrigendum, and the
mapping from `aNaK` to the printed density remains family-specific.

For Ae2, the code's `K_B=1e-4` mM is eight orders of magnitude below the
printed `10^4` mM. It makes the historical Ae2 factor forward and, with the
fitted activity, gives `J_Ae2=0.0195384838` in internal units. In contrast,
Task 10's literal printed evaluation gives factor `-0.0026954453`. Two
independent MATLAB parameter layers support `1e-4`, so a missing minus sign is
plausible. The archived `Ae4_Basis.tex` parameter table nevertheless repeats
`10^4`, and the upstream Falkenberg convention is not established by these
historical files. The code branch is clear; the publication typo is not
independently proven.

### Current/flux convention and remaining dimensional gap

The code resolves the sign bookkeeping as follows:

- K, Na, and Cl channel/tight expressions are signed `I/F` quantities;
- because `z_Cl=-1`, outward Cl has negative `I/F` and is added to the
  intracellular Cl RHS, which equals `-I/(F z_Cl)`;
- luminal Cl receives the opposite sign; and
- one NaK cycle carries one net positive charge outward, so the voltage solve
  uses `F J_NaK`.

This convention is conservation-consistent and explains why no second Faraday
division appears in the historical mass balances. It does not establish a
single physical unit for `gcl`, `gk`, `gtna`, `gtk`, transporter activities,
water coefficients, and the global `10^4` multiplier. Several transformed
values closely match the published table, but with different family-specific
powers of ten. The sign blocker is resolved; the absolute dimensional ledger
is not.

Finally, the original calibration is not an independent baseline test.
`Original/Saliva_Ae4.m` computes the conductances and transporter activities
from the target row itself. The forensic reproduction recovers a maximum
12-row calibration residual of `1.38778e-17`, and every one of the 45 fields
stored in `Par.mat` follows from that calibration or an exact transform to
within `2.51e-15` relative discrepancy. This precisely explains the small
historical residual and confirms lineage, while also confirming the Task 10
warning that a fitted resting row cannot validate its own fitted parameters in
the absence of the objective and unrounded source vector.

## Consequence for the previous `STOP` gate

The historical audit changes the interpretation of several Task 10 failures:
the 59.6 water ratio and 67.8 pA current residual are not failures of the
archived resting row; they arise because the published tables omit the
historical fitting/scaling conventions. The archive also supplies coherent
channel-current signs, a geometry, a QSS voltage elimination, solver
tolerances, and the corrected Nkcc1 concentration basis.

It does not remove the executable-model gate. The archive's central Ae4
mechanism is different, the acid/base equilibrium remains violated, the
published calcium curve is absent, absolute flow and time units are not
recoverable, and several parameter families require undocumented and mutually
different scale factors. Consequently, these two reports do not establish a
trustworthy published-2018 physiological parameter-to-observation map. Any
later decision to reopen identifiability must depend on independent output and
phenotype reproduction of a validated model, not on the small residual of this
historically calibrated but scientifically different system.
