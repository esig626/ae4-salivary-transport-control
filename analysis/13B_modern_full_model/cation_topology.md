# Task 13B modern cation topology

## Evidence-constrained decision

The modern architecture permits Na/K-ATPase and Ca-activated K pathways at
both apical and basolateral membranes. This is required by the post-historical
topology evidence assigned in Task 13B, including Almássy et al. (2018), DOI
`10.1007/s00424-018-2109-0`. The evidence ledger must retain the exact gland,
preparation, assay, and localization context.

The architecture does **not** treat an apical fraction as a direct SMG
measurement.  For the pump, uniform surface density plus the Poulsen--Bundgaard
rat-parotid membrane areas gives the declared transfer
`f_P=0.075075`, with conservative area-error combinations
`0.056878/0.094586`.  For K, `f_K=0.20/0.30/0.40` remains an explicit modeling
bracket rather than a measured interval.  Total pump capacity and the literal
`14 nS` whole-cell Ca-activated K maximum are partitioned without changing
their totals. No genotype-specific pump or K-channel recruitment is included.

## Capacity partition

Let \(f_P,f_K\in[0,1]\). Then

\[
P_a=f_P P_T,\qquad P_b=(1-f_P)P_T,
\]

\[
G_{K,a}=f_K G_{K,T}p_K(C_{eff}),\qquad
G_{K,b}=(1-f_K)G_{K,T}p_K(C_{eff}).
\]

Both K background conductances are zero in the native-source production map;
adding a background would be a new mechanism, not part of the `14 nS` source
literal.  The production calcium gate is the Palk-lineage law

\[
p_K(C_{eff})=\frac{C_{eff}^{1.46}}
{(0.26\ \mathrm{uM})^{1.46}+C_{eff}^{1.46}}.
\]

`C_eff` is dimensionful and must be supplied in `uM`.  The `0.058 uM`
rest value is a source-linked rat-lineage transfer; `0.10/0.12/0.15 uM` are
assumed spatially averaged WT sensitivity inputs, not direct free-Ca
measurements under the matched protocol and not genotype-fit coordinates.

The native CaCC maximum is separately fixed at `31.4 nS` and uses the same
declared Palk reference gate.  At `C_eff=0.058/0.10/0.12/0.15 uM`, the gate is
`0.10062056/0.19860330/0.24437206/0.30936970`, giving CaCC conductances
`3.1595/6.2361/7.6733/9.7142 nS`.  These are already whole-cell
conductances: neither a common multiplier nor a second membrane-area factor is
permitted.

The exact identities `f + (1-f) - 1 = 0` are reported at every membrane
evaluation. Capacity is conserved by partition. Realized pump turnover need
not be identical after redistribution because the apical pump sees luminal K
whereas the basolateral pump sees bath K. That distinction is physiological
and is not a failure of capacity conservation.

The effective pump turnover starting law is

\[
J_P=P\frac{Na_i^3}{Na_i^3+K_{Na}^3}
       \frac{K_o^2}{K_o^2+K_K^2}.
\]

This is a published-lineage saturation form with unfrozen whole-cell
parameters. It is irreversible in the present ATP-replete closure; a
thermodynamically reversible pump must replace it if ATP/ADP/Pi or pump
reversal becomes a relevant gate.

## Pump and K source signatures

For rows `(Na_i,K_i,Na_l,K_l)`, one apical pump cycle and one apical K efflux
have source columns

\[
S_{a,cation}=\begin{pmatrix}
-3&0\\
+2&-1\\
+3&0\\
-2&+1
\end{pmatrix}.
\]

This matrix has rank two. Therefore an apical pump split and an apical K split
add two independent cation/lumen source directions. They add no direct Cl,
TIC, TA, pH, carbon, or water source row. Effects on those variables must be
indirect through voltage, osmolarity, volume, and the other transporters.

This structural result is important for generation decisions: redistributing
pump/K capacity cannot by itself repair an acid-base residual whose coupling
through the rest of the model is absent.

## Local electrochemical driving forces

For a membrane from cell to side `s`, the Nernst voltage is

\[
E_{x,s}=\frac{RT}{z_xF}\log\frac{[x]_s}{[x]_i}.
\]

Conventional outward current is

\[
I_{x,s}=G_{x,s}(V_s-E_{x,s}),
\]

and physical amount flux is `I/(zF)`. The apical K path uses luminal K and
`V_a`; the basolateral K path uses bath K and `V_b`. The apical pump likewise
uses luminal K. This is the required local topology; evaluating every path
with bath K or a shared voltage would silently erase the apical distinction.

The apical Cl current uses the same charge-correct law with `z=-1` and luminal
chloride. It is Ca gated separately from the K partition.

## Paracellular paths

The finite lumen exchanges Na, K, Cl, and HCO3 with the open bath. For each
species,

\[
I_{x,t}=G_{x,t}(V_t-E_{x,t}),\qquad
V_t=V_b-V_a,
\]

with positive current defined lumen to bath. The luminal amount source is
`-I/(zF)`. Including Cl and HCO3 as well as the historical-lineage Na/K paths
gives charge-complete paracellular bookkeeping. Their conductances remain
modeling decisions until the evidence/dimensional audit supplies defensible
values.

## Coupled current closure

Define

\[
I_a=I_{K,a}+I_{Cl,a}+FJ_{P,a},
\]

\[
I_b=I_{K,b}+FJ_{P,b},
\]

and \(I_t=\sum_x I_{x,t}\). The quasi-electrical closure is

\[
I_a-I_t=0,\qquad I_b+I_t=0.
\]

All conductance currents are affine in \(V_a,V_b\), so the implementation
solves one two-by-two linear system exactly at each RHS call. It reports:

* `v_apical_V`, `v_basolateral_V`, and `v_transepithelial_V`;
* each channel, pump, and paracellular current in A;
* both current residuals;
* cell and lumen molar-source-versus-current residuals;
* capacity-partition residuals.

No code unit is relabelled as volts or seconds: conductance is in S, current in
A, Faraday conversion in C/mol, sources in fmol/s, and time in seconds.

## Exact nesting and generation logic

The following limits are exact.

1. `apical_pump_fraction = 0` makes apical pump capacity and source terms zero.
2. `apical_k_fraction = 0` and the source-fixed zero background make the apical K
   current zero.
3. `(f_P,f_K)=(0,0)` is the basolateral-only historical cation-topology limit
   inside the new conservation-explicit architecture.
4. Nonzero fractions preserve total assigned capacity/conductance.

The default fractions are only executable nonzero examples. They must not be
accepted because they yield a convenient root or genotype ratio. The retained
generation must either bound them from independent topology/current evidence
or show that conclusions are robust across their full defensible range.

## Relation to G2 and earlier Task 13 split tests

Task 13 showed that a pump/K split on the unchanged seven-state chassis did not
produce a valid WT calibration. Task 13B therefore retains the scientifically
supported topology but places it inside a larger conservation-correct system:

* intracellular Na and K are conserved amount states;
* luminal Na and K are finite amount states;
* both membrane voltages are reclosed;
* luminal chloride and carbon/alkalinity are explicit;
* water and outflow alter finite volumes;
* acid-base chemistry no longer resides in one HCO3 residual row;
* dynamic beta/cAMP/PKA regulation couples independently to AE4.

This is not an undocumented patch to the Task 13 split chassis. Nor does it
imply that the earlier negative result was caused solely by cation topology.

## Focused tests

Current tests verify:

* the apical pump/K source matrix has rank two;
* zero apical fractions remove both apical pathways exactly;
* total fraction identities are exact;
* local luminal K enters apical pump turnover;
* both membrane-current equations close to numerical roundoff;
* channel/pump sources agree with Faraday-converted currents;
* cell and lumen charge sources close;
* the new cation directions have no direct acid-base row;
* short Radau and BDF whole-cell integrations remain positive and agree.

## Outstanding measurements and gates

Before the topology is accepted in a WT full model, the lead reconstruction
must resolve or bound:

1. whether the localization/current evidence transfers to the modeled gland,
   preparation, and stimulus;
2. transfer robustness of the area-derived pump fractions
   `0.056878/0.075075/0.094586` from rat parotid to mouse SMG;
3. whether functional pump surface density is sufficiently uniform for that
   area-derived transfer;
4. robustness across the assumed `0.20/0.30/0.40` apical-K partition while
   holding the source-mapped `14 nS` total fixed;
5. local luminal K during the relevant WT protocol;
6. transfer of the dimensionful Palk gate and whether apical/basolateral K
   channels share that gate;
7. paracellular Na/K/Cl/HCO3 conductances;
8. resting and stimulated apical/basolateral voltage constraints;
9. robustness of the WT root over the allowed fraction rectangle;
10. whether the simpler basolateral-only nested model fails independent WT
    gates before extra topology is retained.

None of these quantities may be chosen from the AE4-null saliva magnitude or
time course. The present result is a source-informed, current-closed topology
architecture, not a validated quantitative localization estimate or a
WT-passing production root.
