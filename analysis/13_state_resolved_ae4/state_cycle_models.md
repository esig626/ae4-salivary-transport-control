# State-resolved AE4 candidate-cycle families

## Scope, evidence boundary, and controlling result

This note derives transporter-level AE4 cycles only. It does not use, inspect,
or tune against any AE4-null whole-gland outcome. No parameter is fitted here.
The state graphs are candidates for subsequent calibration against
transporter-level experiments.

The main structural conclusion is:

> **The smallest source-supported skeleton is SR2: Na- and K-resolved
> 1:1:2 alternating-access branches that share the empty and Cl-loaded
> conformational states of one transporter.** The branches require distinct
> cation-specific kinetic parameters. SR1 is the corresponding pure-cation
> submodel. SR5 is the smallest extension that can mechanistically generate
> the reported cooperative cation response while transporting only one cation
> per completed cycle.

The fixed-species SR3 `Cl + Na in / HCO3 + K out` cycle is electroneutral and
thermodynamically feasible, but it is not an evidence-established replacement
for the 2016 Na branch. As a sole cycle it disagrees with the directly observed
Na direction in 2016 and cannot account for K-only exchange. It is therefore
retained only as a 2025-motivated add-on hypothesis.

### Evidence labels

| Label | Meaning in this note |
| --- | --- |
| **EXP16** | Direct ion-sensitive-dye or electrophysiological result in Peña-Münzenmayer et al. 2016, [DOI 10.1085/jgp.201611571](https://doi.org/10.1085/jgp.201611571) |
| **THERM16** | Thermodynamic calculation or mechanistic interpretation in the 2016 paper, not a measured stoichiometry |
| **REG21** | Direct regulation experiment in Peña-Münzenmayer et al. 2021, [DOI 10.1152/ajpgi.00145.2021](https://doi.org/10.1152/ajpgi.00145.2021) |
| **MUT25** | Direct functional mutagenesis or membrane-expression control in Catalán et al. 2025, [DOI 10.1152/ajpcell.00346.2024](https://doi.org/10.1152/ajpcell.00346.2024) |
| **MD25** | Molecular-dynamics or homology-model inference in the 2025 paper |
| **HYP25** | A state sequence or stoichiometry explicitly proposed, but not established, in the 2025 discussion |
| **EXACT** | Conservation, charge, cycle-affinity, or Markov-chain consequence |
| **NEW** | A modeling choice introduced here to make a proposed mechanism executable |

The prior constraints and coarse-law limits are taken from
`analysis/12_ae4_mechanism_reconstruction/experimental_evidence.md`,
`analysis/12_ae4_mechanism_reconstruction/ae4_candidate_family.md`, and
`analysis/12_ae4_mechanism_reconstruction/structural_explanation.md`. Task 10
and Task 11 establish that the four-state display in the 2018 paper does not
reduce algebraically to its printed scalar law as written. None of the graphs
below imports that unresolved reduction.

## Source constraints before choosing a graph

1. **Na is transported.** In the 2016 low-external-Cl experiment, AE4-linked
   HCO3 entry was accompanied by an increase in intracellular Na. Na therefore
   moved with HCO3 and opposite Cl under that protocol (**EXP16**). A purely
   Na-gated anion exchanger is not a complete mechanism.
2. **K is transported.** In Na-free cells, restoring external Cl produced
   AE4-linked Cl entry and K loss (**EXP16**). A purely K-gated mechanism is
   likewise incomplete.
3. **The macroscopic process is electroneutral and reversible.** AE4 did not
   create a detectable current; rates were insensitive to voltage over the
   tested range, and imposed-gradient responses reversed (**EXP16**). These
   observations constrain net charge and require reverse edges, but do not
   select one stoichiometry.
4. **Na and K cannot be represented by one undifferentiated ligand.** The 2025
   T756A-T448I double mutant lost detectable Na-supported exchange while
   retaining substantial K-supported exchange (**MUT25**). A single pooled
   `[Na]+[K]` occupancy with one common rate constant cannot express this
   hierarchy.
5. **One putative shared coordination region is supported, not a complete
   atomic cycle.** Functional mutations implicate S446, T448, D709, and T756.
   MD places cation/HCO3 coordination near the TM3-TM10 interface and suggests
   participation by G449, P450, T713, and K879 (**MUT25**, **MD25**). The MD
   model is not an AE4 structure, did not simulate Cl in the proposed step, and
   resolved one HCO3 at the modeled site, not the second site required by a
   1:1:2 bicarbonate cycle.
6. **The cation response is cooperative but its molecular origin is unknown.**
   The reported pure-cation fits were Na `EC50=49 mM`, `nH=2.0`, `Rmin=0.3`,
   `Rmax=1.5`, and K `EC50=62 mM`, `nH=1.8`, `Rmin=0.4`, `Rmax=1.6` in the
   paper's alkalinization-rate convention (**EXP16**). No fit uncertainties
   were reported. A Hill coefficient is not a transported-ion count.
7. **PKA increases macroscopic AE4 activity, but the altered microscopic step
   is unknown.** Isoproterenol/forskolin/PKAc increased activity; H89 prevented
   activation; S173A retained basal function but lost the PKAc response,
   whereas S273A retained it (**REG21**). Direct phosphorylation of S173, the
   activation time constants, and whether PKA changes turnover, active surface
   number, or stoichiometry were not measured.

## Common notation and thermodynamic contract

Let `o` denote extracellular/interstitial solution and `i` cytosol. Define

\[
B=\mathrm{HCO_3^-},\qquad Q=\mathrm{CO_3^{2-}},\qquad
a_X=[X]/c^\circ,
\]

where every activity is dimensionless and uses the same standard
concentration. `E_o` and `E_i` are outward- and inward-facing empty
transporters. A state such as `E_i C_s B_2` is an effective loaded state; it
does not assert simultaneous elementary binding or identify the order of the
three ligands. This lump is deliberate because neither the second bicarbonate
site nor a binding order has been measured.

Every displayed edge is reversible. For an edge `x <-> y`, local detailed
balance requires

\[
\log\frac{q_{xy}}{q_{yx}}
=-\frac{G_y^\circ-G_x^\circ}{RT}
+\sum_r \nu_r\log a_r .
\tag{1}
\]

Around a passive transport cycle, internal state energies cancel. Hence

\[
\log\frac{\prod_j u_j}{\prod_j w_j}=\mathcal A,
\tag{2}
\]

where `u_j` and `w_j` are forward and reverse edge rates and `A` is the
dimensionless chemical affinity. For a unicyclic continuous-time Markov model,
the exact stationary current has the form

\[
J=T_{\rm tot}\frac{\prod_j u_j-\prod_j w_j}
{\mathcal D(u,w)},\qquad \mathcal D(u,w)>0,
\tag{3}
\]

with `D` the positive matrix-tree/King-Altman denominator. Therefore

\[
\operatorname{sign}J=\operatorname{sign}\mathcal A,
\qquad J=0\Longleftrightarrow \mathcal A=0
\tag{4}
\]

for each isolated unicycle. An equilibrium constant different from the one
fixed by the identical aqueous standard states would amount to an unmodeled
energy source. Conformational asymmetry may alter individual rate constants,
but it may not change their required product ratio.

## SR1. Minimal pure-cation alternating-access 1:1:2 cycle

### Exact graph

For each `s in {Na,K}`, define a six-state reversible ring:

\[
\begin{aligned}
E_o
&\rightleftarrows E_oCl
\rightleftarrows E_iCl
\rightleftarrows E_i\\
&\rightleftarrows E_iC_sB_2
\rightleftarrows E_oC_sB_2
\rightleftarrows E_o .
\end{aligned}
\tag{SR1-G}
\]

The ligand-dependent edges are specified by

| Edge in the positive Cl-loading direction | Forward dependence | Reverse dependence |
| --- | --- | --- |
| `E_o <-> E_oCl` | `k^+_Cl,o a_Cl,o` | `k^-_Cl,o` |
| `E_oCl <-> E_iCl` | `k^+_Cl,flip` | `k^-_Cl,flip` |
| `E_iCl <-> E_i` | `k^-_Cl,i` | `k^+_Cl,i a_Cl,i` |
| `E_i <-> E_i C_s B_2` | `k^+_s,i a_s,i a_B,i^2` | `k^-_s,i` |
| `E_i C_s B_2 <-> E_o C_s B_2` | `k^+_s,flip` | `k^-_s,flip` |
| `E_o C_s B_2 <-> E_o` | `k^-_s,o` | `k^+_s,o a_s,o a_B,o^2` |

The composite `C_s B_2` association edge is **NEW**. Expanding it into an
ordered series is not justified until a second bicarbonate site or binding
order is measured.

### Stoichiometry, source, affinity, and reversal

One positive event is

\[
Cl_o+C_{s,i}+2B_i
\rightleftharpoons
Cl_i+C_{s,o}+2B_o .
\tag{SR1-R}
\]

In intracellular source order `(Na,K,Cl,B,Q)`,

\[
\nu_{Na}=(-1,0,+1,-2,0),\qquad
\nu_K=(0,-1,+1,-2,0).
\tag{SR1-S}
\]

Each vector has zero transported charge. The branch affinity is

\[
\mathcal A_s^{112}=
\log\frac{a_{Cl,o}a_{s,i}a_{B,i}^2}
{a_{Cl,i}a_{s,o}a_{B,o}^2}.
\tag{SR1-A}
\]

`A_s=0` is the exact reversal condition. At the physiological concentration
tuple used in the 2016 thermodynamic calculation, the Na branch has
`A_Na=-1.894` and runs opposite physiological Cl loading, whereas the K branch
has `A_K=+4.101` and can load Cl. This is a conditional ideal-activity
calculation (**THERM16**, **EXACT**), not a measurement of the branch mixture.

### Cation structure and mutant implications

- `k_Na,i/o`, `k_Na,flip` and the corresponding K parameters are distinct.
  Equalities between them are testable restrictions, not defaults.
- T756A-T448I is represented by a near-background effective Na branch while
  retaining a nonzero K branch. A common capacity multiplier cannot represent
  this limit (**MUT25**).
- T448I/T448G reduce the effective branch conductance by about one half in the
  Na assay; D709A, S446A, and T756A produce smaller reductions, while D709N,
  S447A, T754A, I758G, and I758R do not show the same loss. These observations
  constrain mutant/WT *lumped conductance ratios*, not a unique binding edge.
- The K-supported double-mutant result can arise from a retained K binding
  rate, a retained K-loaded flip rate, or both. The experiment does not choose
  among those placements.

### Identifiable combinations and Task 12 limit

Pure-cation steady rates can at most constrain an effective capacity
`V_s=T_tot/tau_s`, an apparent cation affinity/slope, and a cycle equilibrium
product. They do not identify the six forward rates, six reverse rates,
transporter abundance, or the binding order separately.

If the state occupancies are quasi-steady and `D` in (3) is approximated by a
constant, SR1 reduces to a branch-resolved Task 12 `C6_112` mass-action law.
A symmetric transition-state parametrization gives the `C7_112` affinity law;
retaining a positive occupancy denominator gives a `C8_112_SAT`-type law.
Thus state resolution does not create a new steady source vector by itself.
It adds relaxation modes, saturation, and an internally constrained
denominator.

**Disposition:** source-compatible as a pure-Na or pure-K assay submodel, but
not yet a mixed-cation model of one AE4 molecule.

## SR2. Na/K branches sharing conformational states

### Exact graph

SR2 joins the two SR1 rings at the empty and Cl-loaded states:

\[
E_o\rightleftarrows E_oCl\rightleftarrows E_iCl\rightleftarrows E_i,
\tag{SR2-common}
\]

and closes the common path through either one of two cation-resolved branches:

\[
\begin{aligned}
E_i&\rightleftarrows E_iNaB_2
\rightleftarrows E_oNaB_2\rightleftarrows E_o,\\
E_i&\rightleftarrows E_iKB_2
\rightleftarrows E_oKB_2\rightleftarrows E_o.
\end{aligned}
\tag{SR2-branches}
\]

The edge dependences are those in the SR1 table. There is one conserved pool

\[
T_{tot}=E_o+E_oCl+E_iCl+E_i+
E_iNaB_2+E_oNaB_2+E_iKB_2+E_oKB_2.
\tag{SR2-pool}
\]

This shared-pool constraint is the minimal mathematical expression of one
protein with a shared TM3-TM10 coordination region (**MUT25**, **MD25**) and is
what distinguishes SR2 from adding two independent scalar fluxes.

### Stoichiometry and thermodynamics

The two fundamental transport cycles have the SR1 source vectors and
affinities `A_Na^112` and `A_K^112`. Local detailed balance must hold around
both. At global equilibrium both affinities vanish. Away from equilibrium the
total entropy-production rate obeys

\[
\dot S_{prod}/R=J_{Na}\mathcal A_{Na}^{112}
+J_K\mathcal A_K^{112}\ge0 .
\tag{SR2-EP}
\]

The graph also permits a fundamental Na/K exchange circulation obtained by a
forward Na cycle plus a reverse K cycle. Its net source is

\[
(-1,+1,0,0,0)
\tag{SR2-slip}
\]

and its affinity is `A_Na-A_K`. Thus zero net Cl current
`J_Na+J_K=0` need not mean zero AE4-mediated cation exchange. This is an exact
state-network distinction hidden by one scalar `J_AE4`.

### Cation-specific structure and mutant implications

- Na and K compete for the same finite pool, but their association,
  dissociation, and loaded-flip rates remain separate.
- The T756A-T448I hierarchy is encoded by a large decrease in the Na branch
  conductance and a smaller decrease in the K branch. The preserved T713-K
  interaction and disrupted T713-Na interaction in the double-mutant MD are a
  plausible placement for that asymmetry (**MD25**), not a direct rate
  measurement.
- D709 and T756 show cation-dependent functional effects. A model with only a
  common conformational multiplier fails this mutant test even if it has
  separate Na and K binding constants.
- A molecular mutation must change one or more microscopic rates while
  preserving the local-detailed-balance product for the passive transport
  cycle. Arbitrarily changing one rate and thereby shifting the equilibrium
  is inadmissible unless a new energy source is supplied.

### Alternative 2016 stoichiometries within the shared graph

The same topology can be instantiated with effective loaded states
`E C_s^m B^b` for the three 2016 candidate counts

\[
(a,m,b)\in\{(1,1,2),(1,2,3),(2,1,3)\},\qquad b=a+m.
\tag{SR2-stoich}
\]

Their branch reaction and affinity are

\[
aCl_o+mC_{s,i}+bB_i
\rightleftharpoons
aCl_i+mC_{s,o}+bB_o,
\]

\[
\mathcal A_s^{amb}=
\log\frac{a_{Cl,o}^{a}a_{s,i}^{m}a_{B,i}^{b}}
{a_{Cl,i}^{a}a_{s,o}^{m}a_{B,o}^{b}}.
\tag{SR2-general-A}
\]

All are electroneutral by `b=a+m`. The 1:2:3 member requires either two
transported cations or repeated use of a site; the 2:1:3 member requires two
Cl events/sites. Neither microscopic arrangement is established by 2025.
These are discrete stoichiometric hypotheses from **THERM16**, not fitted
continuous exponents.

### Identifiability and reductions

Pure-Na data and pure-K data make SR1 and SR2 exactly observationally
equivalent: the unused branch has zero ligand-dependent entry rate. Mixed Na/K
titrations are required to identify the shared-pool denominator and
competition.

In the rare-occupancy or effectively constant-free-pool limit,

\[
J_s\simeq V_s(P_s^+-P_s^-),
\]

so SR2 reduces to Task 12 C3/C6 branch-resolved scalar laws. It reduces to the
published pooled C2 balance only under additional restrictive conditions:
identical branch kinetics and a reverse cation allocation consistent with the
donor-side composition. The printed intracellular mole-fraction split is not
the general reduction on reversal.

**Disposition:** minimal supported mixed-cation skeleton. It earns additional
states over SR1 because it expresses finite-pool competition and the
Na-specific/K-specific mutant hierarchy in one transporter.

## SR3. Asymmetric sequential Na/K cycle

### Exact fixed-species graph

The literal 2025 proposal can be formalized as

\[
\begin{aligned}
E_o
&\rightleftarrows E_oNaCl
\rightleftarrows E_iNaCl
\rightleftarrows E_i\\
&\rightleftarrows E_iKB
\rightleftarrows E_oKB
\rightleftarrows E_o .
\end{aligned}
\tag{SR3-G}
\]

Composite `Na+Cl` and `K+B` binding is the smallest graph that does not invent
an unmeasured binding order. The 2025 paper proposes these orientation-specific
roles but did not simulate Cl binding and did not observe a completed cycle
(**HYP25**, **NEW** formalization).

### Stoichiometry, charge, affinity, and feasibility

One positive event is

\[
Cl_o+Na_o+K_i+B_i
\rightleftharpoons
Cl_i+Na_i+K_o+B_o.
\tag{SR3-R}
\]

Its intracellular source vector is

\[
\nu_{SR3}=(+1,-1,+1,-1,0)
\tag{SR3-S}
\]

in `(Na,K,Cl,B,Q)` order. The charge sum is
`(+1)-(+1)+(-1)-(-1)=0`: `1:1:1:1` is exactly
electroneutral and contains no membrane-potential term.

With positive affinity defined to favor the displayed Cl-loading direction,

\[
\boxed{\mathcal A_{SR3}=
\log\frac{a_{Cl,o}a_{Na,o}a_{K,i}a_{B,i}}
{a_{Cl,i}a_{Na,i}a_{K,o}a_{B,o}}},
\qquad \Delta G_{forward}=-RT\mathcal A_{SR3}.
\tag{SR3-A}
\]

The reversal condition is `A_SR3=0`, equivalently

\[
[Cl]_i^*=[Cl]_o
\frac{[Na]_o}{[Na]_i}
\frac{[K]_i}{[K]_o}
\frac{[B]_i}{[B]_o}
\tag{SR3-reversal}
\]

under ideal activities. Using the matched concentration tuple from the 2016
thermodynamic calculation, supplemented by its inferred
`K_i=139.5 mM`, `K_o=3.4 mM`, gives

\[
\mathcal A_{SR3}=6.6434,
\qquad \Delta G_{forward}\simeq-17.13\;\mathrm{kJ\,mol^{-1}}
\quad (310.15\;\mathrm K).
\tag{SR3-number}
\]

Thus the proposed cycle is thermodynamically feasible and strongly favors Cl
loading under that conditional tuple. This calculation is **EXACT given the
ideal-activity inputs**, but it is not evidence that AE4 uses the cycle.

### Explicit conflict with direct observations

SR3 must not be promoted as a replacement for the 2016 mechanism:

1. **Na-direction conflict.** In SR3, positive Cl entry carries Na inward; the
   reverse event carries Cl and Na outward together. In the 2016 low-external-
   Cl protocol, Cl left and HCO3 entered while intracellular Na increased.
   The direct observation therefore places Na with HCO3 and opposite Cl
   (**EXP16**), whereas fixed SR3 places Na with Cl. The sign convention in
   (SR3-A) makes this disagreement explicit.
2. **K-only conflict.** Fixed SR3 requires extracellular Na on one half-cycle.
   Both 2016 and 2025 observe exchange when K is the main extracellular cation,
   and the T756A-T448I mutant retains substantial K-supported exchange
   (**EXP16**, **MUT25**). Fixed SR3 cannot be the sole route in those baths.
3. **Generalizing is not free.** Allowing K to substitute at the outward Na
   limb or Na at the inward K limb creates additional cycles, including
   same-cation catalytic anion exchange and Na/K exchange. Those routes need
   their own local-detailed-balance constraints and return the model toward an
   SR2-like branched network.

SR3 may therefore be tested only as an additional mixed-cation subcycle next
to a source-supported SR2 route. The added SR3 current has a source vector not
collinear with any Task 12 C1-C8 `Cl in / cation+B out` branch. It imports Na,
exports K, and has `B/Cl=1`; no kinetic or saturation limit can reduce that
source to a Task 12 scalar law without deleting one of its defining species.

### Identifiability

No available pure-cation alkalinization curve identifies the SR3 current.
Required evidence is simultaneous, AE4-specific Na, K, Cl, and inorganic-carbon
flux in a mixed Na-out/K-in condition. Until such a measurement exists, all
SR3 microscopic rates and even its contribution relative to SR2 are
unidentified.

**Disposition:** thermodynamically admissible 2025 hypothesis/add-on; rejected
as a source-complete stand-alone family.

## SR4. Explicit carbonate variants

Carbonate is not a relabeling of bicarbonate. Define the acid-base reaction

\[
B^-\rightleftharpoons H^++Q^{2-},\qquad
K_{a2}=\frac{a_Ha_Q}{a_B}.
\tag{SR4-chem}
\]

At minimum the cell model must carry `Q_i` explicitly, or impose the stated
fast-equilibrium relation while retaining its proton dependence. A dynamic
version needs a reversible reaction rate `R_2` contributing
`(-R_2,+R_2,+R_2)` to `(B,H,Q)`. The extracellular carbonate activity must be
specified consistently as well.

### SR4a: one-cation carbonate exchange

Use the SR1 graph with `E_i C_s B_2` replaced by `E_i C_s Q` and the
corresponding outward state. One event is

\[
Cl_o+C_{s,i}+Q_i
\rightleftharpoons
Cl_i+C_{s,o}+Q_o.
\tag{SR4a-R}
\]

The Na and K source vectors are

\[
(-1,0,+1,0,-1),\qquad (0,-1,+1,0,-1),
\]

and are exactly electroneutral. The affinity is

\[
\mathcal A_{s,Q}=
\log\frac{a_{Cl,o}a_{s,i}a_{Q,i}}
{a_{Cl,i}a_{s,o}a_{Q,o}}.
\tag{SR4a-A}
\]

Under fast carbonate equilibrium,

\[
\frac{a_{Q,i}}{a_{Q,o}}
=\frac{a_{B,i}a_{H,o}}{a_{B,o}a_{H,i}},
\]

so the reduced affinity contains one bicarbonate ratio and one proton ratio,
not the squared bicarbonate ratio in SR1. SR4a therefore does not reduce to a
Task 12 bicarbonate law unless pH and carbonate are artificially clamped and
their dependence is hidden in a fitted constant.

### SR4b: sequential Na/K carbonate exchange

The carbonate counterpart of SR3 is

\[
\begin{aligned}
E_o
&\rightleftarrows E_oNaCl
\rightleftarrows E_iNaCl
\rightleftarrows E_i\\
&\rightleftarrows E_iK_2Q
\rightleftarrows E_oK_2Q
\rightleftarrows E_o,
\end{aligned}
\tag{SR4b-G}
\]

with reaction, source, and affinity

\[
Cl_o+Na_o+2K_i+Q_i
\rightleftharpoons
Cl_i+Na_i+2K_o+Q_o,
\tag{SR4b-R}
\]

\[
\nu_{SR4b}=(+1,-2,+1,0,-1),
\tag{SR4b-S}
\]

\[
\mathcal A_{SR4b}=
\log\frac{a_{Cl,o}a_{Na,o}a_{K,i}^2a_{Q,i}}
{a_{Cl,i}a_{Na,i}a_{K,o}^2a_{Q,o}}.
\tag{SR4b-A}
\]

The source charge is `+1-2-1+2=0`. This is the executable interpretation of
the 2025 proposed `1:1:1:2 (Cl:Na:CO3:K)` count (**HYP25**, **NEW**).
It inherits both stand-alone conflicts listed for SR3.

### Evidence and mutant implications

- The 2025 MD coordinated one HCO3 in the modeled pocket. The authors infer
  either a second bicarbonate site for 1:1:2 or one carbonate for 1:1:1;
  neither transported species/count was directly distinguished (**MD25**,
  **HYP25**).
- The proposed dependence on D709 protonation is a hypothesis. D709A reduced
  activity whereas D709N did not, but those mutations do not report the
  transported carbon species (**MUT25**).
- BCECF alkalinization alone is particularly unable to discriminate HCO3
  entry from carbonate entry followed by proton chemistry.

**Disposition:** retain SR4a as the minimal carbonate alternative and SR4b as
the sequential carbonate hypothesis. Neither may be implemented on the
seven-state bicarbonate-only chassis.

## SR5. Transported cation plus cation-dependent conformational catalysis

SR5 tests the 2016 interpretation in which one cation is transported while an
additional cation occupancy produces cooperativity. The extra cation must not
be silently counted as transported.

### Exact minimal catalytic graph

Start from SR2. For a selected rate-limiting loaded flip
`L_i,s <-> L_o,s`, where `L_r,s=E_r C_s B_2`, add a cation-bound catalytic
path:

\[
\begin{aligned}
L_{i,s}+R_{t,r}
&\rightleftarrows L_{i,s}R_t
\rightleftarrows L_{o,s}R_t
\rightleftarrows L_{o,s}+R_{t,r},\\
L_{i,s}&\rightleftarrows L_{o,s}.
\end{aligned}
\tag{SR5-G}
\]

`R_t` is Na or K at an allosteric site accessible to the *same reservoir* `r`
before and after the conformational transition, for example a cytosolic
regulatory site. It is a catalyst and cancels from the net chemical reaction.
The reservoir/site location is **NEW and unmeasured**; it must be fixed as a
discrete model choice, not fitted continuously.

Local detailed balance around the catalytic subcycle requires

\[
\prod_{catalytic\ loop}\frac{u_j}{w_j}=1,
\tag{SR5-LDB}
\]

because `R_t` appears on both sides. The catalytic route may change the cycle
conductance but cannot move the SR2 reversal condition. If `R_t=C_s`, low
cation availability enters once in transported-cation loading and once in the
catalytic path, allowing a near-quadratic low-concentration response with only
one net transported cation.

An equivalent explicit loading-catalysis version is

\[
E_r+R_{s,r}\rightleftarrows E_rR_s,
\qquad
E_rR_s+C_{s,r}+2B_r
\rightleftarrows E_rC_sB_2+R_{s,r}.
\tag{SR5-loading}
\]

It has the same net stoichiometry and provides a discrete alternative if the
data favor altered loading rather than altered flipping. Current experiments
do not distinguish the two placements.

### Stoichiometry, cation specificity, and reductions

SR5 retains the SR2 1:1:2 source vectors and affinities. Na and K may differ in

- transported-site affinity;
- allosteric-site affinity;
- catalytic enhancement of the loaded flip; and
- mutant effects on any of those three combinations.

The T756A-T448I result requires the effective Na catalytic/transport path to
collapse much more than the K path. It does not identify a second site; the
TM3-TM10 evidence could instead affect the transported site while a separate
allosteric site remains unknown.

Fast equilibrium of the catalytic states produces a positive cation-dependent
multiplier on a thermodynamic SR2 core. This is the state-resolved limit of
Task 12 C5. It is not Task 12 C4 unless the chemical cation term is incorrectly
replaced by a Hill occupancy; SR5 retains direct cation transport and the
correct reversal.

### Identifiability

The 2016 dose-response curves identify only apparent `(EC50_s,nH_s,Rmin_s,
Rmax_s)`. On one concentration sweep, the following remain observationally
equivalent:

1. two transported cations in a 1:2:3 cycle;
2. one transported cation plus SR5 catalytic occupancy;
3. a mixture of transporter subpopulations with different affinities; and
4. state-dependent saturation in SR2 over a restricted concentration range.

The smallest discriminating measurement is simultaneous cation and anion
particle flux near reversal. A one-cation-plus-catalyst model retains the
1:1:2 net count; a two-transported-cation model has 1:2:3.

**Disposition:** smallest mechanistic extension if the Hill slopes are treated
as a structural constraint rather than an empirical interpolation.

## SR6. PKA-regulated state transitions

### Exact minimal phosphorylation graph

Choose a surviving SR2 or SR5 core with states `X_j`. Duplicate it into
unphosphorylated and regulated layers:

\[
X_0^U\rightleftarrows X_1^U\rightleftarrows\cdots
\rightleftarrows X_{N-1}^U\rightleftarrows X_0^U,
\]

\[
X_0^P\rightleftarrows X_1^P\rightleftarrows\cdots
\rightleftarrows X_{N-1}^P\rightleftarrows X_0^P,
\tag{SR6-layers}
\]

with the minimal regulatory connection at a cytosol-accessible state

\[
X_0^U
\xrightleftharpoons[k_{deph}]
{k_{ph}u_{PKA}}
X_0^P.
\tag{SR6-reg}
\]

The selected state and regulated transport edge are discrete hypotheses. To
model a PKA-induced barrier change without inventing active ion pumping, alter
both directions of one transition by the same factor:

\[
u_j^P=\gamma_Pu_j^U,
\qquad w_j^P=\gamma_Pw_j^U,
\qquad \gamma_P>0.
\tag{SR6-barrier}
\]

The forward/reverse ratio and the chemical affinity are then unchanged. A
model that changes this ratio must include the chemical affinity of ATP
hydrolysis and demonstrate energy coupling; the 2021 data do not support that
extra claim.

### Source constraints and mutant limits

- `u_PKA=0` or H89 exposure returns the unregulated layer in the minimal
  representation (**REG21**). H89 is nonselective, so it is a pathway control,
  not proof of the exact kinase edge.
- S173A is represented by `k_ph=0` or `gamma_P=1` while retaining basal
  transport. S273A retains the regulated branch (**REG21**).
- No 2021 measurement assigns regulation to the Cl-loaded flip, cation-loaded
  flip, binding, trafficking, or active transporter number. Those placements
  must remain separate nested candidates.
- No activation/deactivation time constants were measured in a protocol that
  identifies their mapping to a whole-cell stimulation time axis.

### Identifiable combinations and scalar-gate limit

Steady pretreatment experiments constrain only the regulated effective-rate
ratio

\[
F_{PKA}=\frac{J_{eff}(u_{PKA}=1)}{J_{eff}(u_{PKA}=0)}
\tag{SR6-fold}
\]

under each assay. They do not separately identify `k_ph`, `k_deph`, the
phosphorylated fraction, or `gamma_P`. The native AE4-only, total native
exchange, and recombinant folds are different experimental contexts and may
not be pooled into one microscopic constant.

If phosphorylation equilibrates rapidly and changes all relevant cycle
barriers uniformly, SR6 reduces to the Task 12 positive common PKA multiplier.
If it alters one state transition, the full response generally changes the
state occupancy and relaxation spectrum and does **not** reduce to a scalar
gate. State-transition placement can be identified only from transporter-level
activation/relaxation time courses or a perturbation that isolates the
regulated conformation.

**Disposition:** source-supported regulation layer, but its microscopic edge
and kinetics are unresolved. The minimal thermodynamically conservative test
is matched forward/reverse barrier scaling.

## Cross-family parameter identifiability

### Defensible fitted/lumped quantities

The available transporter data can support, at most, the following quantities
without further time-resolved experiments:

| Lumped quantity | What constrains it | What remains confounded |
| --- | --- | --- |
| `V_Na`, `V_K` effective pure-cation capacities | Initial alkalinization/ion-flux rates | Transporter abundance versus turnover; individual flips and binding rates |
| `Kapp_Na`, `Kapp_K`, apparent slopes | 2016 dose responses | Binding order, transported count, allostery, and state saturation |
| `Kcycle` or the product of forward/reverse rates | Reversal plus passive thermodynamics | Distribution over microscopic equilibrium constants |
| `rho_s^mut=V_s^mut/V_s^WT` | 2025 mutant/WT rates in the same cation bath | Binding versus translocation defect |
| Shared-pool competition coefficient | A mixed Na/K titration, not currently available | SR1 population mixture versus SR2 shared states |
| `F_PKA` effective fold | 2021 paired assay context | Phosphorylated fraction, regulated edge, trafficking, and microscopic efficacy |
| Relaxation eigenvalues | Would require resolved transporter time courses | Individual rates remain nonunique without multiple protocols |

Microscopic rate tables should therefore not be reported as estimated values.
Use constrained priors or symbolic parameters, enforce their product ratios,
and record all results in terms of the lumps above.

### Exact observational equivalences

1. **SR1 versus SR2 in a pure-cation bath:** exact equivalence after deleting
   the unavailable branch.
2. **SR2 versus Task 12 C3/C6 at low occupancy:** same branch source and
   forward-minus-reverse numerator; SR2's shared denominator becomes
   effectively constant.
3. **SR2 zero net anion flux versus no transport:** not equivalent. Opposing Na
   and K branch currents can leave a pure Na/K exchange current.
4. **SR5 one transported cation plus catalyst versus a two-cation cycle:**
   indistinguishable from a single alkalinization dose-response curve, but not
   from simultaneous particle stoichiometry.
5. **SR4 bicarbonate versus carbonate in one pH-clamped protocol:** rates can
   be reparameterized to agree locally; varying pH or measuring carbonate/
   proton balance breaks the equivalence.
6. **SR6 fast regulation versus a scalar activity multiplier:** equivalent
   only in the rapid/uniform-barrier limit. State-specific modulation changes
   relaxation modes.
7. **Different binding orders within a composite loaded state:** steady
   initial-rate data identify only King-Altman combinations. Binding order is
   not recoverable from the current experiments.

## Evidence disagreements that must remain explicit

| Issue | Direct evidence | Structural proposal/inference | Modeling consequence |
| --- | --- | --- | --- |
| Direction of Na relative to Cl | 2016: low external Cl caused Cl exit/HCO3 entry with intracellular Na increase; Na moves opposite Cl (**EXP16**) | 2025 discussion proposes extracellular Na stabilizes and enters with Cl (**HYP25**) | Fixed SR3 cannot replace the 2016 Na branch; test only as an added mixed-cation subcycle |
| K-only operation | 2016 and 2025: K as main external cation supports exchange; direct 2016 PBFI supports K transport (**EXP16**, **MUT25**) | Fixed SR3 assigns Na specifically to the outward limb | A stand-alone fixed SR3 graph is source-incomplete |
| Bicarbonate count/species | Voltage independence requires electroneutrality; no direct stoichiometric count (**EXP16**) | MD resolves one HCO3; authors propose a second HCO3 site or one carbonate (**MD25**, **HYP25**) | Retain SR1/SR2 and SR4 as discrete alternatives; do not hide carbonate in `B` |
| Na/K selectivity | 2016 WT pure-cation maxima and apparent affinities are similar (**EXP16**) | 2025 double mutant separates Na-supported from K-supported activity (**MUT25**) | WT similarity does not justify pooled kinetics; mutant-specific branch rates are required |
| Cooperative cation response | Hill slopes about 2 and 1.8 (**EXP16**) | More than one transported cation or an extra regulatory cation were both retained by the authors | Keep two-cation SR2 stoichiometry and one-cation SR5 allostery as observationally equivalent candidates |
| PKA mechanism | PKA-pathway activation and S173 dependence (**REG21**) | Direct S173 phosphorylation and regulated transport step were not measured | SR6 barrier, availability, and trafficking placements remain separate nested models |

## Candidate ranking before calibration

| Rank/status | Family | Reason |
| --- | --- | --- |
| **Minimal supported skeleton** | **SR2-112** | Smallest one-transporter graph with reversible alternating access, direct Na and K branches, shared states, and cation-specific mutant behavior |
| **Pure-assay submodel** | **SR1-112** | Necessary building block, but lacks a uniquely defined mixed-cation completion |
| **Minimal cooperative extension** | **SR5 on SR2-112** | Preserves one transported cation and thermodynamic reversal while allowing a near-quadratic response |
| **Discrete stoichiometric alternatives** | **SR2-123, SR2-213** | Electroneutral 2016 possibilities; microscopic extra-ion sites/order are unestablished |
| **Regulatory layer** | **SR6 on a surviving core** | PKA activation is experimentally required for stimulated comparisons; microscopic placement is unresolved |
| **Chemistry alternative** | **SR4a** | Smallest explicit carbonate cycle; requires carbonate/proton chemistry absent from the historical chassis |
| **Add-on hypothesis only** | **SR3, SR4b** | Electrically and thermodynamically admissible but conflict with direct direction/K-only evidence as stand-alone mechanisms |

## Required discrimination experiments

The most informative single protocol is a reversal-clamp experiment in
AE4-expressing cells with independently controlled intracellular and
extracellular Na, K, Cl, HCO3/CO3, and pH, while measuring AE4-specific Na,
K, and Cl flux simultaneously. Include WT and T756A-T448I, and repeat with and
without PKA activation.

This one design separates:

- SR2 same-cation branches from SR3 cross-cation cycling by flux direction;
- one versus two transported cations by particle stoichiometry;
- HCO3 from carbonate by pH/carbonate dependence;
- independent branches from shared-state competition in mixed Na/K; and
- scalar PKA gating from a changed relaxation spectrum.

Until such data exist, the defensible implementation order is SR2-112, SR5 as
a nested cooperative extension, SR6 as a nested regulation layer, and SR4 as a
separate chemistry branch. SR3 should never be used as the only AE4 route.

## Structural choices that remain unresolved

1. Net count: SR2-112 versus the discrete SR2-123 or SR2-213 alternatives.
2. Origin of cooperativity: a second transported cation versus an SR5
   catalytic occupancy versus ordinary state saturation.
3. Microscopic order: the ligand-binding sequence, the second bicarbonate or
   chloride site, and which loaded transition contains the TM3-TM10
   coordination region.
4. Carbon species: bicarbonate versus explicit carbonate/proton chemistry.
5. Mixed-cation routing: whether SR3 contributes any parallel current at all;
   direct evidence excludes only its use as the sole replacement route.
6. PKA placement: altered turnover barrier, active transporter availability,
   or trafficking, plus the phosphorylation/dephosphorylation time scales.
