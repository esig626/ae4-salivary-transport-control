# Stage-B missing-balance localization

## Gate, evidence release, and result

Stage B opened only after `stage_a_results.md` and
`stage_a_reveal_marker.json` were frozen. Five transporter-level conditional
survivors produced zero valid WT capacity calibrations and zero WT Gate-2
passes on the unchanged seven-state chassis. None used AE4-null ionic or
secretion evidence. This is a fixed-chassis failure, not evidence that AE4 is
biologically dispensable.

Only two Stage-B localization measurements are used here:

- AE4-null resting `Cl_i = 36.50 +/- 1.60 mM` (`E15-06`); and
- AE4-null resting `pH_i = 6.89 +/- 0.02` (`E15-07_KO`).

There is no measured AE4-null acinar volume, native `Na_i` or `K_i`, luminal
Na/K composition, or membrane-partitioned current in the 2015 experiment.
Null uptake kinetics, initial Cl exit, NKCC/NHE readouts, secretion magnitude,
and the secretion time course are not localization targets in this analysis.

The controlling conclusion is:

> **The released Cl/pH pair does not determine one missing source vector. It
> defines a five-dimensional target-state manifold whose required source
> varies in all five free directions. Consequently neither a unique one-module
> correction nor a unique module pair can be named from Stage-B ionic evidence.**

The independently supported apical Na/K-pump plus apical K-conductance
topology is the smallest omitted topology that should be tested first. It adds
two exact lumen/current directions that the old intracellular projection could
not see. It is a candidate chassis extension, not a localized or calibrated
repair: the released measurements do not observe those two directions or
identify their submandibular fractions.

## 1. Candidate-independent AE4-null field and root

Every state-resolved adapter has an exact genotype-zero short circuit. With a
common non-AE4 parameter fingerprint,

\[
F^{KO}_{SRj}(y)=F_0(y)
\]

for every state graph, barrier gauge, cation branch, and PKA placement. Thus
state resolution cannot change the knockout vector field on this fixed
chassis. It can change a WT denominator only after a valid WT state exists;
Stage A found none.

The frozen 20-start bounded solve gives

\[
y^0_{KO}=(118.5652767,5.55934618,65.66070347,
21.13228225,123.70751315,48.78221503,58.90444836)
\]

in chassis order `(Na_l,K_l,height,Na_i,K_i,Cl_i,HCO3_i)`. Its raw residual is
`5.77e-14`, pH is `7.37847391`, and model-reference volume is `2.97417821 pL`.
The Cl and pH errors are `+7.676` and `+24.424` reported SEM from the released
targets. Volume is a model diagnostic, not a knockout target.

An implementation independent of both the production state-cycle code and
the Task-12 root solver re-solved the documented null equations from perturbed
positive coordinates. Its maximum state difference from the frozen CSV root
was `3.34e-9`, with residual `1.58e-12`; direct production evaluation of the
saved row retained residual `5.77e-14`.

## 2. The exact five-dimensional target manifold

Let `X=2439.50353087137` be the fixed impermeant-charge amount and let `H`
denote cell height. The reduced chassis eliminates free intracellular proton
through

\[
h_i=B_i+C_i+X/H-N_i-K_i.
\tag{1}
\]

At the released pH,

\[
h_i^*=1000\,10^{-6.89}=1.2882495517\times10^{-4}\ {\rm mM}.
\]

Therefore every state compatible with the two released means is

\[
\begin{aligned}
C_i&=36.50,\\
B_i&=h_i^*-36.50-X/H+N_i+K_i,
\end{aligned}
\tag{2}
\]

with five free coordinates

\[
(N_l,K_l,H,N_i,K_i).
\]

In `u=1/H` coordinates, Eq. (2) is the affine hyperplane

\[
-N_i-K_i+Xu+B_i=h_i^*-36.50=-36.4998711750.
\tag{3}
\]

Positivity restricts its domain but does not supply the five missing
measurements. The same dimension holds at every one-reported-SEM Cl/pH corner.
No historical or WT coordinate is silently substituted as an observed
knockout value.

## 3. Required correction is a family, not a vector

For a chosen member `y` of Eq. (2), define `d(y)` as the added amount/current/
water source required to make the AE4-null field stationary. The expanded
coordinate order is

\[
\begin{split}
(&Na_i,K_i,Cl_i,HCO_{3i},H_i,CO_{2i};\\
 &Na_l,K_l,Cl_l,HCO_{3l};q_a,q_b,q_t;V_i),
\end{split}
\tag{4}
\]

where `q_a,q_b,q_t` are positive-charge equivalents in the directions cell to
lumen, cell to bath, and lumen to bath. They are not water fluxes. `V_i` is a
separate cell-water/volume row.

The old reduction has no dynamic `CO2_i` or luminal `HCO3` state, so those two
entries are undetermined. Existing membrane potentials already re-solve the
two QSS current closures; the net additional current demand is zero, meaning
an electrogenic addition must recruit a compensating current and re-solve the
voltage closure. Chemical amount rates and water rates have uncertified,
different units and are never combined in one cosine.

Two transparent slices demonstrate why plugging in free coordinates is not
innocent. The `model_anchored` slice holds every unmeasured coordinate at the
common null root and changes only Cl/pH-consistent `HCO3`. The
`wt_reference_cations` slice uses the same null lumen coordinates but the
historical reference `H=28.7`, `Na_i=25`, and `K_i=120`. Neither slice is an
estimate of the knockout state.

| Conditional slice | Derived `HCO3_i` (mM) | Required cell source `(Na,K,Cl,HCO3,H)` | Required lumen source `(Na,K,Cl)` | Required `dH/dt` |
| --- | ---: | --- | --- | ---: |
| model anchored | `71.1867504` | `(-0.15772665,-0.08097604,-0.16687851,+0.42974694,+0.50157111)` | `(+0.00658162,+0.00047633,+0.00705795)` | `-4.30711e-8` |
| WT-reference cations | `23.5000058` | `(-0.11490889,-0.07882192,-0.12406083,+0.04707089,+0.11674087)` | `(-0.01734824,-0.00069445,-0.01804269)` | `-2.03197e-4` |

The two required cell vectors are materially different, and even the lumen
source changes sign. Differentiating the determined entries of `d(y)` with
respect to the five free coordinates at the model-anchored point gives five
nonzero singular values

\[
(0.40223543,0.01887636,0.00934814,0.00163360,0.000449923),
\]

hence numerical rank `5`. That rank persists at all four one-SEM target
corners. The ambiguity is not one nuisance scale along a fixed direction; the
source itself rotates across the observationally equivalent target family.

## 4. Exact module signatures and topology rank

The independently audited per-event signatures are stored verbatim in
`module_projection.csv`. Their exact rational structural ranks are:

| Projection | Exact rank | Interpretation |
| --- | ---: | --- |
| all 14 expanded signatures | `12` | the declared source set has two exact dependencies/aliases |
| six intracellular chemical rows | `6` | the candidate set spans that full displayed chemical space only with several directions |
| 12 determined reduced rows | `11` | omits the absent `CO2_i` and `HCO3_l` demands |

Most importantly, apical and basolateral pump cycles are identical after
projection onto intracellular chemistry, as are apical and basolateral K
channels. Redistributing fixed total capacity adds

\[
\delta\sigma_P=\sigma_{P,a}-\sigma_{P,b}
=(0_i;+3,-2,0,0;+1,-1,0;0),
\tag{5}
\]

\[
\delta\sigma_K=\sigma_{K,a}-\sigma_{K,b}
=(0_i;0,+1,0,0;+1,-1,0;0).
\tag{6}
\]

Equations (5)-(6) have exact rank `2` in the expanded space and rank `0` in
the six intracellular chemical rows. They alter lumen Na/K and apical versus
basolateral current without directly adding a new intracellular stoichiometry.
This proves both why Task 12 could not test the membrane split and why released
intracellular Cl/pH cannot determine its two fractions. Only a full re-solved
chassis or measurements in the extended coordinates can test it.

## 5. One- and two-module source projections

One/two-module projections were performed only within the five cell amount
rows, which share the chassis-native chemical flux unit. The candidate set was
`NKCC1`, `AE2`, `NHE1`, total Na/K pump, total K channel, apical Cl channel,
and carbon hydration. A positive coefficient means the displayed
physiological event direction. No secretion quantity enters the calculation.

| Slice and sign rule | Best one direction | Residual fraction | Best two directions | Coefficients | Residual fraction |
| --- | --- | ---: | --- | --- | ---: |
| model anchored, signed | carbon hydration | `0.353326` | NKCC1 + carbon hydration | `(-0.095410,+0.465659)` | `0.120909` |
| model anchored, nonnegative | carbon hydration | `0.353326` | apical Cl + carbon hydration | `(+0.166879,+0.465659)` | `0.261994` |
| WT-reference cations, signed | NKCC1 | `0.597997` | NKCC1 + NHE1 | `(-0.059278,-0.086186)` | `0.297902` |
| WT-reference cations, nonnegative | apical Cl | `0.834351` | apical Cl + carbon hydration | `(+0.124061,+0.081906)` | `0.656699` |

No one- or two-module cell signature is exact on either slice, and the best
support changes when unmeasured coordinates change. Signed improvements often
require reversal or a capacity decrease, not addition of the displayed
physiological module. These projections localize a coupled cation/anion/
acid-base deficit in the fixed equations; they do not identify a protein pair.

The older frozen-row cation result remains exact but answers a different
question. Replacing the historical WT Na-only AE4 event by its fixed Na/K
partition requires

\[
d_{WT}=(-0.03217931477,+0.03217931477,0,0).
\]

The pump direction `(-3,+2)` has alignment `5/sqrt(26)=0.980581` and irreducible
one-direction residual `1/sqrt(26)=0.196116`; pump plus an independent Na-load
direction such as NHE1 gives an exact positive algebraic decomposition. That
is a WT fixed-row mixed-cation closure theorem, not the Stage-B knockout
source. It must not be relabelled as evidence for genotype-specific pump or
NHE retuning.

## 6. Local Jacobian/source-response projection

To distinguish a source correction from an equilibrium displacement, the
AE4-null reduced-coordinate Jacobian was differentiated at the common bounded
root. Each existing chassis parameter was perturbed logarithmically, membrane
potentials were re-solved inside the RHS, and the implicit response

\[
\frac{dy^*}{d\log\theta}=-F_y^{-1}F_{\log\theta}
\tag{7}
\]

was projected to `(Cl_i,pH_i)`. The state Jacobian condition number is about
`7.03e4`, so this is a local diagnostic at a biologically failed root, not a
finite parameter reconstruction.

The `2 x 11` existing-parameter response matrix has rank `2` and nullity `9`.
Therefore many noncollinear pairs span the two observed discrepancies exactly
by dimension alone. Examples using Cl-SEM/pH-SEM scaling are:

| Pair | Exact linear coefficients in `d log(parameter)` | Implied factors | Status |
| --- | --- | --- | --- |
| total pump + total K conductance | `(-13.200,-18.338)` | `(1.85e-6,1.09e-8)` | exact rank fit; far outside local validity |
| total pump + NHE1 | `(-1.775,-1.696)` | `(0.169,0.183)` | exact rank fit; large and not an independently supported genotype retuning |
| total K conductance + NHE1 | `(+2.849,-1.960)` | `(17.28,0.141)` | exact rank fit; large |
| total K conductance + AE2 | `(+2.930,+1.833)` | `(18.73,6.25)` | exact positive rank fit; large and contradicts the intended minimal topology test |

No exact two-output fit above is a nonlinear re-solved repair. Even the best
one-direction ranking is coordinate-scale dependent: unscaled `(mM,pH)`
selects pump with residual `0.0258`, whereas reported-SEM scaling selects an
NHE1 decrease with residual `0.2998`. Requiring an explicit scale in code and
reporting rank separately prevents either cosine from becoming biological
evidence.

## 7. Localization decision and Round-5 handoff

The Stage-B result is a proof of nonidentifiability at the released evidence
level:

1. the exact knockout field is common to every state-resolved AE4 candidate;
2. released Cl/pH define a five-dimensional state family, not a complete
   coordinate;
3. the associated required source family has local rank five;
4. no one/two cell-source signature is exact across the two declared slices;
5. many two-parameter equilibrium responses fit two outputs algebraically,
   with unsupported finite changes; and
6. the source-supported pump/K membrane split is invisible in intracellular
   source projection and uncalibrated in matched submandibular tissue.

Round 5 should therefore test nested, independently sourced modules in this
order without using null secretion:

1. redistribute fixed total pump and K capacity into apical/basolateral copies
   and complete lumen/current closure;
2. if that does not jointly repair independently constrained WT and released
   null Cl/pH, replace the historical acid-base reduction with conserved
   carbon plus genotype-invariant NHE1 regulation; and
3. test their smallest coupled combination before adding total capacity or
   genotype-specific gains.

The narrowly localized missing measurement is the **stimulus-dependent,
membrane-partitioned cation-current support tuple**

\[
(P_a,P_b,J_{K,a},J_{K,b})
\]

in the matched gland, together with the native Na/K/volume coordinates needed
to map it to the chemical balance. An executable discriminator is: in WT and
`Ae4-/-` mouse **submandibular** acini under the exact `0.3 uM CCh + 5 uM IPR`
protocol, measure time-resolved intracellular Na and K and cell volume, plus
luminal or collected-effluent Na and K so that the free `Na_l,K_l`
coordinates are observed rather than fixed. Use local apical versus basal
pharmacological/electrophysiological dissection to estimate ouabain-sensitive
pump current and Ca-activated K current at each membrane. Retain Cl and pH
imaging in the same preparation, and ideally co-measure total inorganic carbon
plus an NHE- or buffer-flux readout. Without that carbon/acid-base arm, the
experiment is a decisive first test of cation topology, not a complete
mechanism identifier.
The test decides M1 if the measured extended fluxes determine nonzero
`delta sigma_P` and `delta sigma_K` support sufficient to predict the released
Cl/pH displacement under a frozen carbon model; failure with a normal cation
tuple sends the reconstruction specifically to conserved carbon/acid-base
flux rather than to another unconstrained AE4 law.

All exact signatures, conditional source vectors, local response columns,
scalings, coefficients, ranks, and limitations are recorded in
`results/13_state_resolved_ae4/module_projection.csv` and implemented in
`src/state_resolved_ae4/localization.py`.
