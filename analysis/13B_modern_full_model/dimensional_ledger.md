# Task 13B dimensional and physical-time ledger

## Controlling result

The modern full-model scaffold is dimensionally closed on a declared physical
coordinate system:

- finite-compartment amounts are stored in `fmol`;
- cell and local-lumen volumes are stored in `pL`;
- concentrations are in `mM`;
- conventional currents are in `A` and potentials in `V`;
- all amount and volume rates are per physical second;
- regulatory fractions are dimensionless and their derivatives must be in
  `s^-1`;
- osmolarities are numerically in `mOsm L^-1` and water flows in `pL s^-1`.

This closes the *dimensions*. It does **not** certify the present default
numbers as physiological kinetics. Most transporter capacities and
local-outflow parameters remain executable placeholders.  The native-source
panel now fixes the separately sourced whole-cell channel maxima (`31.4 nS`
CaCC and `14 nS` total K), while their cross-gland transfer and the effective
calcium input remain explicit assumptions. The three water coefficients use the exact Palk/2018
published-model conversion, but their transfer to the reconstructed native
unit has not been independently remeasured. A trajectory produced from these
defaults is therefore on a seconds-labelled axis by construction, but it is
not yet a physiologically calibrated time course.

The inherited `rhs_scale=10000` is explicitly rejected as a physical-time
conversion. It has no established mapping to seconds or minutes and is absent
from `src/modern_full_model/`. No historical code-coordinate stimulus time is
relabelled as experimental time.

No AE4-null secretion magnitude or timing datum was inspected or used in this
audit. The dimensional choices, tests, defects, and acceptance rules are all
independent of that held-out validation.

The exhaustive machine-readable inventory is
`results/13B_modern_full_model/dimensional_ledger.csv`. Its 331 rows enumerate
every current core and regulatory state, parameter, acid-base observable,
transporter flux, current, water term, RHS row, conversion, scaling decision,
and relevant diagnostic.

## 1. Exact base conversions

### Amount and concentration

The coordinate identity is exact in SI prefixes:

\[
(1\ \mathrm{mM})(1\ \mathrm{pL})
=10^{-3}\frac{\mathrm{mol}}{\mathrm L}\,10^{-12}\mathrm L
=10^{-15}\mathrm{mol}
=1\ \mathrm{fmol}.
\]

Thus, for every explicitly represented solute,

\[
c\ [\mathrm{mM}]=\frac{n\ [\mathrm{fmol}]}{V\ [\mathrm{pL}]},
\qquad
n\ [\mathrm{fmol}]=c\ [\mathrm{mM}]V\ [\mathrm{pL}].
\]

There is no hidden factor of `1000`, `1e-12`, or `1e15` in this conversion.
Because the dynamic coordinates are amounts, concentration dilution follows
from `c=n/V`; it must not be added as a second source term.

This amount convention also resolves the historical fourth-order-rate hazard.
A concentration monomial written for molar activities cannot be reused with
mM-valued variables: a fourth-order term changes by `1000^4=10^12`. The
modern reversible coarse laws and AE4 affinities use ratios of like-unit
activities, while every absolute rate is declared separately in `fmol/s` or
`s^-1`.

### Current and Faraday conversion

For conventional current positive from the first named compartment to the
second and ionic valence `z`, the signed particle flux is

\[
J_{1\rightarrow2}[\mathrm{fmol/s}]
=\frac{I_{1\rightarrow2}[\mathrm A]}{zF[\mathrm{C/mol}]}
10^{15}.
\]

At `F=96485.33212 C/mol`:

- `+1 pA` carried by a monovalent cation is
  `+0.0103642696566 fmol/s` in the named direction;
- the same positive conventional current on chloride is
  `-0.0103642696566 fmol/s` chloride in the named direction;
- `1 fmol/s` of a Na/K-pump cycle produces
  `9.648533212e-11 A`, because `3 Na` out and `2 K` in is one net positive
  charge equivalent outward per cycle.

The code uses volts with siemens, so `S*V=A`; it does not mix historical mV
with SI conductance. Both membrane current identities and their chemical
source equivalents are tested to roundoff.

## 2. State and derivative units

| State rows | Stored unit | RHS unit | Interpretation |
|---|---|---|---|
| `na_i_fmol`, `k_i_fmol`, `cl_i_fmol` | fmol particles | fmol particles/s | finite cell ionic amounts |
| `tic_i_fmol` | fmol carbon | fmol carbon/s | cell total inorganic carbon |
| `alk_i_fmol` | fmol charge equivalents | fmol equivalents/s | cell total alkalinity |
| `volume_i_pL` | pL | pL/s | cell volume |
| `na_l_fmol`, `k_l_fmol`, `cl_l_fmol` | fmol particles | fmol particles/s | finite local-lumen ionic amounts |
| `tic_l_fmol` | fmol carbon | fmol carbon/s | local-lumen total inorganic carbon |
| `alk_l_fmol` | fmol charge equivalents | fmol equivalents/s | local-lumen total alkalinity |
| `volume_l_pL` | pL | pL/s | local-lumen volume |
| regulatory suffix | dimensionless fractions or normalized activities | s^-1 | family-specific cAMP/PKA/AE4 regulatory states |

`pH`, free H, HCO3, CO2, CO3, and effective-buffer occupancy are algebraic
observables. pH, pKa, and pKw are dimensionless logarithmic p-scale values;
they are not physical units named “pH”. The explicit H and OH concentrations
are converted from molar p-scale definitions to mM with the required factor
`1000`.

Total alkalinity is stored in equivalents, not carbon moles:

\[
TA=HCO_3+2CO_3+B^-+OH-H.
\]

Consequently the modeled mobile charge equivalent is exactly
`Na + K - Cl - TA`. This is why adding one HCO3 changes both TIC and TA by one,
whereas neutral CO2 changes TIC alone.

## 3. Parameter audit

All fields returned by `parameter_records(FullModelParameters())` have an
explicit unit and provenance tag. The machine-readable ledger lists every
field individually. Their dimensional and numerical status is summarized
below.

| Family | Declared units | Dimensional result | Numerical/evidence result |
|---|---|---|---|
| physical constants | K, J mol^-1 K^-1, C mol^-1 | pass | `R` and `F` certified; temperature remains protocol-specific |
| acid-base constants | dimensionless p-scale | pass after annotation correction | effective pKa/pKw and buffer pKa require temperature/preparation audit |
| bath composition | mM, dimensionless p-scale, mOsm L^-1 | pass | present values are representative or reference-derived, not frozen protocol data |
| buffer/fixed pools | fmol sites, fmol equivalents, fmol osmoles | pass as coordinates | values and overlap among physical pools are not identified |
| NKCC1/NHE1/AE2 capacities | fmol cycles/s | pass | placeholders; no physiological time claim |
| CO2 exchange | fmol s^-1 mM^-1 | pass | whole-cell area is absorbed; values uncalibrated |
| pump capacity | fmol cycles/s | pass | irreversible ATP-replete open-energy closure; not a reversible thermodynamic pump law |
| half-saturation constants | mM or uM | pass | the Palk Ca gate uses `0.26 uM`; other transporter values remain historical/test starting values |
| localization/gating factors | dimensionless fractions/exponents | pass | pump `0.056878/0.075075/0.094586` is an area-density-derived transfer; K `0.20/0.30/0.40` remains a modeling bracket |
| conductances | S | pass | native-source literals are `31.4 nS` CaCC and `14 nS` total K; they are absolute whole-cell maxima and cannot receive another area or common-scale factor |
| hydraulic coefficients | pL s^-1 (mOsm L^-1)^-1 | pass after annotation and source-unit corrections | exact Palk/2018 published-model conversion; physiological transfer remains uncertified |
| local lumen/outflow | pL and s^-1 | pass | local compliance closure, not a measured gland secretion law |
| initial conditions | mM, pL, dimensionless p-scale | pass | positive executable reference only, not an accepted WT root |
| AE4 carrier amount | fmol carrier | pass | no source-certified absolute amount |
| AE4 attempt rates | s^-1 per carrier | pass | no source-certified absolute rates; QSS time scale remains unlicensed |
| AE4 EC50/Hill summaries | mM and dimensionless | pass | only a declared observation-map sensitivity in the original assay context |

The parameter registry originally labelled pKa/pKw/pH as `pH` units and
hydraulic coefficients as `pL s^-1 mOsm^-1`. Those annotations have been
corrected by the whole-cell owner to `dimensionless (p-scale)` and
`pL s^-1 (mOsm L^-1)^-1`, respectively. Field names were retained for API
stability.

### Dynamic beta/cAMP/PKA regulation

The mandatory `camp_pka.py` module is dimensionally closed across its nested
R0--R3 families:

| Family | Dynamic states | Parameters with time dimension | RHS unit |
|---|---|---|---|
| R0 | none | none | none; instantaneous static control |
| R1 | effective AE4 activation fraction | `tau_activation_s` in s | s^-1 |
| R2 | effective cAMP fraction; AE4 activation fraction | `tau_camp_s`, `tau_activation_s` in s | s^-1 for each row |
| R3 | effective PKA activity fraction; AE4 regulated fraction | `tau_pka_s` in s; forward/reverse rates in s^-1 | s^-1 for each row |

`beta_input`, every regulatory state, activity, construct-response scale, and
capacity multiplier are dimensionless. In particular,
`camp_effective_fraction` is a normalized pathway coordinate, not a cAMP
concentration, and `pka_effective_activity_fraction` is not a biochemical
activity unit. R3's `forward_regulation_rate_s` field name is retained for API
stability, but its numerical unit is inverse seconds, not seconds.

The equations are

\[
\dot x=(\beta-x)/\tau_x,
\]

for R1;

\[
\dot c=(\beta-c)/\tau_c,\qquad
\dot x=(c-x)/\tau_x,
\]

for R2; and

\[
\dot p=(\beta-p)/\tau_p,\qquad
\dot q=k_f p(1-q)-k_r(1-p)q,
\]

for R3. At every zero boundary the corresponding derivative is nonnegative,
and at every unit boundary it is nonpositive, so the unit hypercube is an
analytic invariant set. At constant input, R1 equilibrates at `x=beta`, R2 at
`c=x=beta`, and R3 at

\[
p=\beta,\qquad q=\frac{k_f p}{k_f p+k_r(1-p)}.
\]

Implicit Radau/BDF evaluations can query an unaccepted internal stage just
outside the analytic invariant set. The module therefore declares a
dimensionless `INTEGRATION_TRIAL_TOLERANCE=1e-4` and clamps only regulatory
*state* evaluations within that numerical band. Larger excursions fail;
`beta_input`, parameters, and analytic equilibria retain exact `[0,1]`
validation. This is a solver-stage guard, not biological variability and not
an absolute ODE-solver tolerance.

When `k_f=k_r=1/tau_activation`, the R3 downstream row reduces exactly to
`(p-q)/tau_activation`; the default R2 and R3 families are therefore
dynamically isomorphic after renaming their states. That is a structural
non-discrimination result, not kinetic identification.

The common observation map is dimensionless:

\[
M_{AE4}=M_0+\Delta M\,s_{coupling}\,s_{construct}\,a.
\]

It changes AE4 capacity without changing the thermodynamic reversal. The
default `Delta M=0.25` and the alternative visual reading `0.70` are preserved
as separate 2021 sensitivities, not averaged or treated as exact targets.
The qualitative `S173A=0` and `S273A=1` response scales encode dependence of
the detected increment on S173 while explicitly not claiming direct S173
phosphorylation. A transition-specific R4 coupling is rejected because the
data do not identify which AE4 edge would be modulated.

Most importantly, the 2021 evidence establishes beta/adenylate-cyclase/cAMP,
PKA dependence, and construct logic, but reports no acute cAMP, PKA-activity,
AE4-state, phosphorylation, activation, washout, or deactivation trace that
identifies these kinetic constants.  The current contract uses R1 downstream
times `10/30/90 s`, R2 upstream/downstream pairs
`(3,10)/(10,30)/(30,90) s`, and symmetric R3 rates that reproduce the R2
rows exactly.  These are dimensionally valid sensitivity assumptions only. R0 is a
useful negative control but cannot by itself satisfy the mandatory dynamic
regulation requirement. No regulatory time or gain was selected from the
AE4-null secretion holdout.

### Calibration and dynamic-validation numerics

The WT calibration uses ten independent core balance rows on exact cell and
lumen charge manifolds. Each raw row is divided by a same-unit scale: amount
rows by an explicit fmol/s scale and volume rows by an explicit pL/s scale.
Dynamic regulatory rows initially entered that objective as unscaled `s^-1`
values beside dimensionless core residuals. The retained correction requires
one explicit positive `s^-1` conditioning scale per regulatory row, making
every least-squares residual dimensionless. Those conditioning scales are
numerical choices, not measurements of regulatory kinetics.

The root-Jacobian SVD initially used columns in mixed mM, pL, p-scale, fmol,
and fmol/s coordinates. Its singular spectrum and relative numerical-rank
threshold were therefore unit-label dependent. The retained code multiplies
each column by its declared bound span after the residual rows have been made
dimensionless. Reported numerical rank is now invariant to relabelling mM as
M or pL as L. The bundled diagnostic `RestingRoot.max_abs_raw_rhs` remains an
explicitly debug-only heterogeneous scalar. Production WT-calibration
artifacts do not emit that scalar: they report separate maxima for amount
rows in fmol/s, volume rows in pL/s, regulatory rows in s^-1, and dependent
alkalinity rows in fmol-equivalents/s. Acceptance uses the dimensionless
scaled residuals, never a cross-unit maximum.

Dynamic solver absolute tolerances are now a vector with separate native-unit
values for amount states (`fmol`), volume states (`pL`), and regulatory
fractions. Conservation diagnostics retain their individual units and are
gated only after division by a per-key same-unit tolerance. The solver
comparison reports per-coordinate absolute differences plus dimensionless
coordinate-normalized differences; it no longer forms one absolute maximum
over fmol, pL, and regulatory fractions.

Radau and BDF remain two algorithms in the same SciPy implementation. Their
agreement is a dual-method check, not the independent reimplementation
required of Agent I.

## 4. Flux ledger

### Electroneutral homeostasis

`J_NKCC1`, `J_NHE1`, and `J_AE2` are in fmol **cycle events** per second.
Their exact cell-source signatures are

\[
\begin{array}{c|ccccc}
&Na&K&Cl&TIC&TA\\ \hline
NKCC1&1&1&2&0&0\\
NHE1&1&0&0&0&1\\
AE2&0&0&1&-1&-1
\end{array}
\]

so each column has zero charge source under `Na+K-Cl-TA`. Their log activity
ratios are dimensionless and their `tanh` arguments divide by a dimensionless
log-affinity width. The source dimensions and reversal signs are exact; the
bounded constitutive shapes and capacity values are modeling decisions.

### AE4

AE4 attempt rates are in `s^-1` per carrier, carrier amount is in fmol, and
the QSS branch event rates are therefore in fmol cycles/s. The retained
working branches return the exact conserved-coordinate sources

\[
J_{Na}(-1,0,+1,-2,-2),\qquad
J_K(0,-1,+1,-2,-2).
\]

Carrier occupancies and regulatory gains are dimensionless. The summed edge
rates are in `s^-1`, and the reported carrier relaxation time is in seconds.
The algebra and units are complete, but no independent carrier relaxation
measurement currently licenses the QSS approximation or any numerical
whole-cell delay. Test values such as `1` or `10 s^-1` are fixtures only.

The quantity `entropy_production_over_r_fmol_s` has units fmol/s because it is
cycle flux times dimensionless affinity. Converting it to entropy production
requires fmol-to-mol conversion and multiplication by `R`.

### Pump, channels, and paracellular transport

Pump outputs are fmol cycles/s. Channel and paracellular currents are first
formed in amperes from absolute conductance and voltage, then converted once
through `I/(zF)`. No second area or volume factor is permitted after that
conversion. Chemical source vectors account for the 3Na/2K pump stoichiometry
and anion signs exactly.

The native channel map uses absolute whole-cell maxima, converted once as
`31.4 nS = 31.4e-9 S` for CaCC and `14 nS = 14e-9 S` for total Ca-activated
K.  The K value is partitioned between apical and basolateral membranes; it is
not duplicated.  Both K background conductances are zero in the production
source map.  The Palk-lineage activation coordinate is explicitly
dimensionful:

\[
p(C_{eff})=\frac{C_{eff}^{1.46}}
{(0.26\ \mathrm{uM})^{1.46}+C_{eff}^{1.46}},
\qquad C_{eff}\ [\mathrm{uM}].
\]

Thus `0.058 uM` is a source-linked resting transfer, while
`0.10/0.12/0.15 uM` are assumed spatially averaged WT stimulus sensitivities.
They are neither dimensionless activation fractions nor matched-protocol
free-Ca measurements.  None was selected from an AE4-null output.  Romanenko
et al.'s approximately `0.126 uM` mouse-SMG CaCC sensitivity stays a
voltage/protocol-dependent diagnostic unless an explicit observation map is
provided.

The current closure is algebraic and neglects capacitive charge. This is a
fast-electrical modeling assumption, not a units problem. If membrane
capacitance is later retained, its state equation must use
`C_m [F] dV/dt [V/s] = -I [A]` with a separately certified whole-cell
capacitance.

The present Na/K-pump law is ATP-replete, irreversible, and independent of
ATP/ADP/Pi and voltage. It may be used as a declared open-energy homeostatic
closure, but it cannot support a claim of pump thermodynamic reversal or zero
cycle affinity at global equilibrium.

### CO2 and luminal outflow

Neutral CO2 exchange multiplies a whole-cell coefficient
`fmol s^-1 mM^-1` by a concentration difference, giving fmol carbon/s. It
changes TIC but not TA or current.

For every luminal amount,

\[
J_{out,X}[\mathrm{fmol/s}]
=q_{out}[\mathrm{pL/s}]c_{l,X}[\mathrm{mM}],
\]

again using `1 mM*pL=1 fmol`. TIC and TA outflow units are carbon fmol/s and
equivalent fmol/s, respectively.

## 5. Water audit

The sign convention is

- `q_b`: bath to cell;
- `q_a`: cell to lumen;
- `q_t`: bath to lumen;
- `q_out`: lumen to ductal sink.

The water coefficients multiply osmolarity differences in `mOsm L^-1` and
produce `pL/s`. Therefore

\[
\dot V_i=q_b-q_a,\qquad
\dot V_l=q_a+q_t-q_{out}
\]

and internal `q_a` cancels from the total finite-volume balance exactly.

### Published coefficient conversion

The Palk/2018 lineage reports

- `P_a = 4.32e-12 L^2 mol^-1 s^-1`;
- `P_b = 5.15e-11 L^2 mol^-1 s^-1`;
- `P_t = 2.60e-13 L^2 mol^-1 s^-1`.

For a concentration difference whose numerical value is supplied in mM,

\[
P\left[\frac{L^2}{mol\,s}\right]
\Delta c[\mathrm{mM}]
\frac{10^{-3}\ mol/L}{1\ \mathrm{mM}}
\frac{10^{12}\ \mathrm{pL}}{1\ L}
=P\,10^9\Delta c\ [\mathrm{pL/s}].
\]

The exact executable coefficients are therefore `0.00432`, `0.0515`, and
`0.000260 pL s^-1 (mOsm L^-1)^-1`, respectively. The final unit uses the
model's ideal-particle convention in which one solute particle at `1 mM` is
numerically `1 mOsm L^-1`. The historical equation is already `q=P*delta_osm`;
no additional membrane-area factor is licensed. These absolute coefficients
therefore carry the historical geometry/area convention.

The earlier executable placeholders (`2e-5`, `2e-5`, `1e-5`) were smaller by
factors of `216`, `2575`, and `26`. Because those factors differ, correcting
them changes `P_b/P_a` and `P_t/P_a`, not merely a common water time scale.
Every WT resting root and dynamic trajectory produced with the placeholders
is invalidated and must be recalibrated. `outflow_rate_s` and
`lumen_dead_volume_pL` are separate new-model closures; this source-unit
conversion provides no reason to change them. Published-model provenance
certifies the conversion, not its independent physiological transfer to the
new reconstruction.

The audit initially found a material osmotic omission/ambiguity. The finite
effective cell-buffer pool contains a fixed number of buffer particles, so it
contributes one osmole per buffer site regardless of protonation. The first
implementation added `cell_impermeant_osmoles_fmol` but did not add
`cell_buffer_total_fmol`, with no provenance statement that one included the
other. At the executable reference (`30 fmol` buffer, `1 pL` cell), this was a
`30 mOsm L^-1` distinction.

The retained correction defines `cell_impermeant_osmoles_fmol` as *other*
impermeant particles and adds the finite buffer explicitly once through
`water.cell_impermeant_osmoles_fmol`. The other-impermeant default was reduced
from `86.1615928` to `56.1615928 fmol`, preserving the reference total of
`86.1615928 fmol` instead of disguising the correction as a changed reference
osmolarity. Regression tests verify the exact sum and show that adding
`10 fmol` buffer at `1 pL` increases osmolarity by `10 mOsm L^-1` when the
charge-derived reference quantity is updated consistently. The dimensional
defect is therefore corrected; both component values remain uncalibrated.

Free H and OH are omitted from osmolarity as a declared negligible correction.
Nonideal osmotic coefficients and reflection coefficients are also not
represented; any later addition must preserve the stated units.

## 6. Area, volume, and tissue scaling

The modern implementation deliberately uses absolute whole-cell lumps:

- pump and transporter capacities in fmol/s per modeled cell;
- conductances in S per modeled cell;
- CO2 exchange coefficients in fmol s^-1 mM^-1 per modeled cell;
- hydraulic coefficients in pL s^-1 (mOsm L^-1)^-1 per modeled cell/local
  unit; the current values inherit an absolute Palk/2018 geometry convention.

Therefore no membrane area factor appears in an RHS evaluation. This is
dimensionally valid only if fitted/imported quantities are already absolute.
If a source reports a density, the missing conversion must be explicit:

\[
G=g_A A_m,\qquad J=j_A A_m.
\]

For conductance normalized by capacitance, a whole-cell capacitance or a
specific-capacitance-plus-area map is required.  The native channel values are
already whole-cell conductances.  The apical K fraction partitions its
absolute total.  The pump fraction is unusual: it is derived from rat-parotid
luminal/basolateral area densities under an explicitly assumed uniform pump
surface density, giving `0.075075` and conservative `0.056878/0.094586`
sensitivities.  These fractions route pump capacity; the source areas must not
also multiply that capacity in the RHS.

Cell and local-lumen volumes are dynamic states, so amount-to-concentration
scaling is explicit. The present local `q_out` is not whole-gland saliva.
Scaling pL/s to gland uL/min additionally requires a source-supported number
of contributing cells/secretory units and any acinus/duct collection map:

\[
Q_{gland}[\mu\mathrm L/\mathrm{min}]
=Q_{local}[\mathrm{pL/s}]N\frac{60}{10^6}.
\]

Kondo et al. (2015) provide an independent upper-bound construction.  Their
mouse-SMG acinar fraction is `66.91%`; using a deliberately generous one-gland
mass of `100 mg`, `1 mg approximately 1 uL` tissue density, and the model's
`1.30 pL` cell gives at most `51,469,231` modeled acinar cells.  At the WT
`9--10 uL/min` envelope, a viable model must therefore supply at least
`0.00291436--0.00323818 pL/s` per cell, or
`0.00224182--0.00249091 s^-1` specific flow.  This bound does not identify a
unique cell count or acinus/duct collection efficiency, but it makes the old
G5 counts (`294,706,403` and `260,704,116`) an absolute geometry failure.
Kondo et al. (2019) describe ipsilateral-gland collection, so silently doubling
the production gland count is not licensed; paired glands remain a sensitivity
only.  No part of this bound uses the knockout secretion holdout.

## 7. Physical-time decision

Task 13's inherited chassis used a numerical `rhs_scale=10000`, a calcium
transition at code time `100`, and flux/water units whose absolute whole-cell
mapping was not documented. The following inferences are forbidden:

- `t_code` is seconds;
- `t_code` is minutes;
- dividing or multiplying `t_code` by `rhs_scale` yields physical time;
- the old step at `100` corresponds to any experimental early/late landmark;
- agreement of historical code-coordinate trajectories establishes a
  physical AE4, PKA, ionic, or water time scale.

The modern model instead defines `t_s` in seconds and requires every rate to
carry a per-second unit. This is the correct forward construction, but it
licenses physiological timing only after the following independent WT-facing
quantities are frozen:

1. AE4 attempt rates or an independent bound justifying carrier QSS;
2. the retained cAMP/PKA/AE4 activation and deactivation rates;
3. NKCC1, NHE1, AE2, pump, channel, and CO2 whole-cell rate scales;
4. transfer/validation of the published hydraulic coefficients and calibration
   of the local-outflow coefficients;
5. calcium and beta-input waveforms on the same seconds axis;
6. cell/local-lumen geometry and, for absolute gland flow, tissue scaling.

The 100-second factor used only to nondimensionalize resting-root residuals is
a numerical row scale. It is not a measured relaxation time. Likewise,
`RestingRoot.max_abs_raw_rhs` is a maximum across heterogeneous fmol/s, pL/s,
and regulatory/s rows; it is debug-only and cannot be interpreted as a
physical norm. Root decisions must use explicitly scaled, per-row residuals.

## 8. Uncertified quantities that block acceptance

The following are dimensionally declared but numerically uncertified:

- absolute AE4 carrier abundance and all AE4 attempt rates;
- whether AE4 carrier QSS is fast on the independently calibrated WT time
  scale;
- native mapping of heterologous AE4 EC50/Hill summaries;
- all dynamic cAMP/PKA activation/deactivation rates and the physiological
  common-capacity gain; the module is dimensionally audited, but the 2021
  source does not identify acute kinetics and gives gain sensitivities rather
  than one exact universal multiplier;
- NKCC1, NHE1, AE2, Na/K-pump, paracellular, and CO2 absolute whole-cell
  capacities, plus cross-gland/protocol transfer of the source-mapped channel
  maxima and Palk gate;
- direct mouse-SMG validation of the area-derived pump-fraction transfer and
  quantitative identification of the assumed apical-K bracket;
- cell buffer amount/pKa, TIC, fixed charge, and impermeant-particle mapping;
- the numerical sizes and biological overlap of the now-explicit finite-buffer
  and other-impermeant osmotic pools;
- cell and local-lumen volumes and any membrane-area conversion;
- physiological transfer of the published Palk/2018 hydraulic coefficients,
  plus the lumen dead volume and local-outflow rate;
- calcium and beta-input waveforms in physical seconds;
- exact single-cell/local-unit to whole-gland scaling and collection efficiency
  within the independent Kondo-derived upper bound;
- a fully independent dynamic implementation beyond the two SciPy methods.

Radau and BDF agreement is a useful solver-method crosscheck, but two methods
from one library are not the Task 13B independent reimplementation required
for final acceptance.

## 9. Executable audit and inherited failures

`tests/test_modern_dimensions.py` now certifies:

- the exact `1 mM*pL=1 fmol` identity in both SI arithmetic and the state API;
- the signed Faraday conversion and chloride sign;
- the numerical `1 pA` monovalent conversion;
- the exact `P_source*10^9` conversion of all three Palk/2018 water
  coefficients and their published-model provenance;
- water-outflow `s^-1*pL=pL/s` and `mM*pL/s=fmol/s` identities;
- R1--R3 fraction-per-second RHS arithmetic, zero RHS at each exact
  constant-input equilibrium, nonnegative capacity output, and the explicit
  `UNMEASURED_ACUTE_KINETICS` status;
- absence of an executable `rhs_scale` name or attribute anywhere in the
  modern package (documentation may mention its explicit rejection);
- an explicit unit for every registered whole-cell parameter.

Together with the existing modern tests, the focused suite passes:

```text
PYTHONPATH=src python -m unittest \
  tests.test_modern_dimensions \
  tests.test_modern_camp_pka \
  tests.test_modern_ae4 \
  tests.test_modern_full_model_whole_cell \
  tests.test_modern_dynamic_validation \
  tests.test_modern_independent \
  tests.test_modern_adversarial

92 tests passed
```

The frozen Task 13 state-resolved suite also passes independently of the new
model:

```text
PYTHONPATH=src python -m unittest discover \
  -s tests -p 'test_state_resolved_ae4_*.py'

72 tests passed
```

An earlier audit-stage repository-wide discovery run executed 216
pre-regulatory tests and reported four failures and six errors, all caused by
absent Task 11 forensic artifacts or absent Task 12 frozen-result files. These
are pre-existing materialization failures, not Task 13 or modern dimensional
failures; this agent did not fabricate or edit those historical artifacts.

## 10. Audit disposition

The modern architecture passes the basic dimensional ledger for amounts,
concentrations, acid-base equivalents, electrical current, Faraday conversion,
water-flow units, explicit seconds, and volume dilution. The inherited
physical-time ambiguity is not carried forward.

It does **not** yet pass a physiological time-scale, absolute-flow, or WT
dynamic gate.  The geometry audit now supplies a hard upper bound and rejects
the old G5 scale, but does not identify a unique gland map. The corrected buffer-osmole bookkeeping must eventually
receive calibrated component values. The dynamic regulatory module now passes
the dimensional and positivity audit, but every retained per-second value
must still receive independent WT/regulatory/transporter support or an
explicitly justified sensitivity bound before a stimulated trajectory can be
called physiologically timed.
