# Task 13B acid-base and carbon reconstruction

## Decision

Use total inorganic carbon (TIC) and total alkalinity (TA) as the two dynamic
intracellular acid-base coordinates, with a finite non-carbonate buffer pool.
Use the same TIC/TA coordinates in the finite lumen. Recover H, pH, CO2,
HCO3, CO3, and buffer occupancy algebraically at every RHS evaluation.

This formulation was selected because it makes carbon and proton-equivalent
bookkeeping exact while avoiding a numerically ill-scaled free-H amount state.
It is not an infinite-buffer shortcut: buffer amount is finite, and its
concentration changes with cell volume.

## Species and conserved coordinates

For each finite compartment,

\[
C_T=[CO_2]+[HCO_3^-]+[CO_3^{2-}],
\]

\[
TA=[HCO_3^-]+2[CO_3^{2-}]+[B^-]+[OH^-]-[H^+].
\]

The finite monoprotic buffer satisfies

\[
B_T=[BH]+[B^-],\qquad
\frac{[B^-]}{B_T}=\frac{1}{1+10^{pK_B-pH}}.
\]

`tic_i_fmol`, `alk_i_fmol`, `tic_l_fmol`, and `alk_l_fmol` are amount
states. Buffer totals are fixed amounts, not fixed concentrations.

## Algebraic reconstruction

For \(h=10^{-pH}\) in molar units, \(K_1=10^{-pK_1}\), and
\(K_2=10^{-pK_2}\), the carbon fractions are

\[
\alpha_0=\frac{h^2}{h^2+K_1h+K_1K_2},\quad
\alpha_1=\frac{K_1h}{h^2+K_1h+K_1K_2},\quad
\alpha_2=\frac{K_1K_2}{h^2+K_1h+K_1K_2}.
\]

Then

\[
[CO_2]=\alpha_0 C_T,\quad [HCO_3^-]=\alpha_1 C_T,\quad
[CO_3^{2-}]=\alpha_2 C_T.
\]

`acid_base.speciate` solves the monotone TA equation for pH on a declared
numerical bracket, then reports explicit species and carbon/alkalinity
residuals. A requested TIC/TA pair outside that bracket raises an error rather
than clipping pH to a convenient bound.

The effective pKa values and one-site buffer pKa are not asserted as uniquely
identified physiology. They are registered as published-lineage starting
values or new modeling decisions and must be audited against temperature,
preparation, and WT buffer-capacity evidence.

## Why TA is the charge coordinate

The explicit acid-base charge is

\[
[H^+]-[HCO_3^-]-2[CO_3^{2-}]-[B^-]-[OH^-]=-TA.
\]

Therefore the complete modeled bulk-charge identities are

\[
Q_i/F=n_{Na,i}+n_{K,i}-n_{Cl,i}-n_{TA,i}-n_X,
\]

\[
Q_l/F=n_{Na,l}+n_{K,l}-n_{Cl,l}-n_{TA,l}.
\]

This identity is exact; it does not discard free H because H is already in
TA with the correct sign.

## External flux map

For positive addition of explicit species to a compartment,

| Added species | TIC source | TA source |
|---|---:|---:|
| HCO3 | +1 | +1 |
| CO3 | +1 | +2 |
| CO2 | +1 | 0 |
| H | 0 | -1 |
| OH | 0 | +1 |

This gives the following whole-cell couplings.

* NHE1 imports Na and extrudes H: intracellular TA increases by one per cycle,
  while TIC is unchanged. Its combined source charge is zero.
* AE2 imports Cl and exports HCO3: intracellular TIC and TA each decrease by
  one per cycle. Its combined source charge is zero.
* Each retained AE4 1:1:2 branch imports Cl while exporting two HCO3 and one
  intracellular cation: TIC and TA each decrease by two per positive cycle,
  and the complete source charge is zero.
* Neutral CO2 exchange changes TIC but not TA.
* Paracellular HCO3 changes luminal TIC and TA together.
* Luminal outflow removes TIC and TA at `q_out * concentration`.

## Closed-reaction conservation proof

For species rows `(CO2,HCO3,H,BH,B-)`, use the internal reaction columns

\[
CO_2\rightarrow HCO_3^-+H^+,\qquad BH\rightarrow B^-+H^+.
\]

The source matrix is

\[
S=\begin{pmatrix}
-1&0\\
+1&0\\
+1&+1\\
0&-1\\
0&+1
\end{pmatrix}.
\]

The carbon row `(1,1,0,0,0)`, buffer-site row `(0,0,0,1,1)`, and charge row
`(0,-1,+1,0,-1)` are exact left-null vectors of `S`. The executable test uses
the same matrix and requires all six residual entries to be exactly zero.

Carbonate dissociation and water ionization add no new TIC or TA source; they
are already included in the algebraic TA definition.

## Open bath and two finite compartments

The bath composition and pH are fixed open-reservoir inputs. Basolateral CO2,
AE2, AE4, NHE1, and NKCC1 exchanges may therefore change finite-cell amounts
without adding a bath state. Apical CO2 exchange has equal and opposite cell
and lumen TIC sources. The whole-model diagnostic explicitly subtracts bath
exchange and luminal outflow before asserting carbon accounting.

This distinction prevents a false claim of closed-system carbon conservation
when the experimental preparation is bath-coupled.

## Relation to the inherited acid-base block

The new block differs structurally from a single HCO3 concentration plus an
effective buffer flux:

1. carbon cannot appear or disappear inside the finite compartments;
2. buffer sites are finite and conserved;
3. NHE, AE2, AE4, CO2 exchange, paracellular HCO3, and outflow have separate
   TIC/TA signatures;
4. cell-volume dilution is automatic in amount coordinates;
5. pH is coupled to carbon, buffer amount, and volume rather than assigned by
   one residual source law;
6. charge follows directly from TA.

No claim is made that the historical published model was wrong. This is a
modern conservation-explicit reconstruction needed to test WT pH closure.

## Fast-equilibrium assumption and escalation rule

The current implementation assumes carbonate and effective-buffer equilibrium
is fast relative to the modeled secretion/regulatory dynamics. This is a
`NEW_MODELING_DECISION`, not a measured carbonic-anhydrase time constant.

Retain it only if WT perturbation data and solver checks show no required
finite relaxation. If independent WT acid-load or CO2-step data show a
resolved lag, add exactly one independent carbon species coordinate (for
example CO2 amount) while retaining TIC and TA conservation. Do not choose a
hydration delay from the AE4-null secretion curve.

## Calibration implications

WT pH alone cannot separately identify TIC, finite buffer amount, buffer pKa,
NHE1 capacity, AE2 capacity, AE4 capacity, and CO2 permeability. Staged fitting
must therefore use, where available:

* bath CO2/HCO3 and pH;
* WT intracellular pH and resting chloride;
* independent buffer-capacity/acid-load data;
* independent NHE/AE2 perturbations;
* geometry and cell-volume data;
* transporter assays for AE4 rather than whole-cell phenotype compensation.

The AE4-null secretion magnitude and timing are prohibited for all such
choices. Released knockout resting pH/Cl may be used only under the staged
localization rule and must not be silently promoted to WT calibration data.

## Tests and present conclusion

Focused tests certify:

* TIC/TA-to-pH round-trip accuracy;
* positive explicit CO2, HCO3, CO3, BH, and B- at the reference state;
* carbon, buffer-site, and charge left-null identities;
* correct external H/HCO3/CO2/CO3 source mapping;
* integrated cell/lumen/open-bath carbon accounting;
* exact bulk-charge closure using TA;
* finite-buffer dilution through amount/volume coordinates.

These analytical and executable results establish a conservation-correct
acid-base architecture. They do not yet establish a physiological WT root or
identify the buffer and CO2 parameters.
