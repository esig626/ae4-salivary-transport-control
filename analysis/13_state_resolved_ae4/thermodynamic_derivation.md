# Thermodynamic and mathematical audit of state-resolved AE4 cycles

## Scope, evidence boundary, and result in one paragraph

This is the independent Agent-C derivation for Task 13.  It uses conservation,
electroneutrality, local detailed balance, the transporter-level experiments of
Peña-Münzenmayer et al. (2016), the PKA experiments of Peña-Münzenmayer et
al. (2021), the structural/mutagenesis study of Catalán et al. (2025), and the
published/historical model audits from Tasks 10--12.  No AE4-null whole-gland
secretion result is used in an equation, parameter choice, model ranking, or
rejection below.

The principal audit result is deliberately two-part.  The 2025-proposed
`1:1:1:1` cycle is **exactly electroneutral and strongly favorable** for
chloride loading under the salivary concentration tuple, provided its stated
directions are respected: extracellular Cl and Na move inward while
intracellular HCO3 and K move outward.  Its dimensionless forward affinity is
about `+6.643`, or `Delta G = -17.12 kJ/mol` at 310 K.  It is nevertheless
**not admissible as the sole, fixed AE4 transport mode** for the salivary model:
the direct 2016 SBFI experiment found that Na moves with HCO3 and opposite Cl,
whereas reversing the proposed `1:1:1:1` cycle in the low-external-Cl assay
would move Na outward with Cl.  The cycle may be retained only as a clearly
labelled secondary or condition-dependent hypothesis, with a second mode and a
specific experiment capable of establishing the switch.  This is an
experimental-direction rejection, not a thermodynamic one.

## 1. Conventions

Let `o` denote the extracellular/interstitial side and `i` the cytoplasmic
side.  Write normalized activities as

\[
 a_{X,x}=\frac{\gamma_{X,x}[X]_x}{c^\circ},
 \qquad c^\circ=1\ {\rm M},
 \qquad x\in\{o,i\}.
\]

The numerical checks below use ideal activity ratios, so the common standard
concentration cancels.  Let

\[
 B=\mathrm{HCO_3^-},\qquad R=\mathrm{CO_3^{2-}}.
\]

A positive cycle current is defined as the chloride-loading direction.  The
intracellular source order is

\[
 s=(s_{Na},s_K,s_{Cl},s_B,s_R).
\]

For a bicarbonate-only cycle, omit `s_R`.  The charge added to the cell per
event, in units of the elementary charge, is

\[
 q_s=s_{Na}+s_K-s_{Cl}-s_B-2s_R. \tag{1}
\]

Thus `q_s=0` is the exact electroneutrality test.  A colon-separated list of
species is not enough to apply (1); the direction of every species is part of
the stoichiometry.

The electrochemical potential is

\[
 \widetilde\mu_{X,x}=\mu_X^\circ+RT\log a_{X,x}+z_XF\psi_x.
\]

For a forward event with free-energy change `Delta G`, define the dimensionless
affinity

\[
 \mathcal A=-\frac{\Delta G}{RT}.
\]

If `q_s` positive charges move from `o` to `i`, the electrical contribution is

\[
 \mathcal A_{elec}=-q_s\frac{F(\psi_i-\psi_o)}{RT}. \tag{2}
\]

It vanishes for every accepted AE4 cycle below.  Net-cycle electroneutrality
does not imply that every microscopic transition is voltage insensitive: a
charged loaded complex can have voltage-dependent partial transitions even
when their electrical work cancels around the full cycle.

## 2. Local detailed balance and exact Markov constraints

### 2.1 Edge and cycle parameterization

For a reversible edge `r -> s`, let `u_rs` and `w_sr` be its forward and
reverse rates.  A local-detailed-balance parameterization is

\[
 \frac{u_{rs}}{w_{sr}}=\exp(\ell_{rs}),
 \qquad
 u_{rs}=\gamma_{rs}e^{\ell_{rs}/2},
 \qquad
 w_{sr}=\gamma_{rs}e^{-\ell_{rs}/2},
 \qquad \gamma_{rs}>0. \tag{3}
\]

Binding activities are included in `ell_rs` or, equivalently, as mass-action
factors in the corresponding rate.  For every oriented graph cycle `c`, local
detailed balance requires

\[
 \boxed{
 \log\frac{\prod_{e\in c}u_e}{\prod_{e\in c}w_e}
 =\sum_{e\in c}\ell_e
 =\mathcal A_c
 }. \tag{4}
\]

Equation (4) is the Wegscheider identity.  A passive transporter that returns
to the same protein state cannot acquire a free equilibrium bias from its
binding energies: those state energies telescope around the cycle.  With the
same aqueous standard states on both sides and no coupled chemical reaction,
the cycle standard factor is one.  Any fitted `K_eq != 1` must therefore be
identified as a nonideal standard-activity correction or a specified coupled
reaction; it cannot be used as an unexplained phenotype-adjustment parameter.

For an `n`-state unicycle with clockwise rates `u_j` and reverse rates `w_j`,
the stationary cycle current has the form

\[
 J=\frac{\prod_j u_j-\prod_j w_j}{D(u,w)},
 \qquad D(u,w)>0. \tag{5}
\]

Consequently,

\[
 \operatorname{sign}J=\operatorname{sign}\mathcal A,
 \qquad J=0\Longleftrightarrow\mathcal A=0. \tag{6}
\]

A unicyclic candidate that runs opposite its chemical affinity is rejected.
For a multicyclic graph, individual branch currents may run uphill by coupling
to another downhill cycle; the correct condition is total entropy production,

\[
 \frac{\dot S_{prod}}{R}=\sum_c J_c\mathcal A_c\ge 0, \tag{7}
\]

not a separate `J_c A_c >= 0` test imposed on an arbitrary cycle basis.

### 2.2 Carrier and mass invariants

Every state-resolved implementation must conserve transporter amount,

\[
 E_T=\sum_r E_r,\qquad \frac{dE_T}{dt}=0,\qquad E_r\ge0. \tag{8}
\]

At quasi-steady transporter occupancy, edge currents reduce exactly to the
completed cycle currents and their stoichiometric source vectors.  If the
transporter states are instead integrated dynamically, transient binding
stores ions on the transporter.  Exact conservation then requires edge-level
solution fluxes and the time derivative of bound substrate.  Integrating
binding states while also inserting only the completed-cycle source into the
cell balances silently loses this storage term.  That approximation is
permitted only after showing that transporter occupancy is fast and negligible
on the whole-cell amount scale.

## 3. Exact source vectors and charge checks

The accepted forward directions and intracellular sources are as follows.

| Family or mode | Forward chemical event | Intracellular source per event `(Na,K,Cl,B,R)` | `q_s` |
| --- | --- | --- | ---: |
| SR1 Na, `1:1:2` | `Cl_o + Na_i + 2B_i -> Cl_i + Na_o + 2B_o` | `(-1,0,+1,-2,0)` | 0 |
| SR1 K, `1:1:2` | `Cl_o + K_i + 2B_i -> Cl_i + K_o + 2B_o` | `(0,-1,+1,-2,0)` | 0 |
| SR2 shared Na/K | Sum of the two preceding branch currents | `(-J_Na,-J_K,J_Na+J_K,-2(J_Na+J_K),0)` | 0 |
| SR3 sequential HCO3, `1:1:1:1` | `Cl_o + Na_o + B_i + K_i -> Cl_i + Na_i + B_o + K_o` | `(+1,-1,+1,-1,0)` | 0 |
| SR4 one-cation carbonate, `1:1:1` | `Cl_o + C_i + R_i -> Cl_i + C_o + R_o` | `(-1,0,+1,0,-1)` for Na, or `(0,-1,+1,0,-1)` for K | 0 |
| SR4 sequential carbonate, `1:1:1:2` | `Cl_o + Na_o + R_i + 2K_i -> Cl_i + Na_i + R_o + 2K_o` | `(+1,-2,+1,0,-1)` | 0 |
| SR5 | Inherits the source of its transported core; a nontransported regulatory ligand has zero net source | core-dependent | 0 if the core passes |
| SR6 | Phosphorylation changes transporter state, not ion stoichiometry | core-dependent | 0 if the core passes |

For the exact 2025 directions, the `1:1:1:1` charge check is

\[
 q_{1111}=+1-1-1-(-1)=0. \tag{9}
\]

For the sequential carbonate proposal it is

\[
 q_{1112}=+1-2-1-2(-1)=0. \tag{10}
\]

By contrast, interpreting `1:1:1:1` as Cl inward while Na, HCO3, and K all
move outward would give `q_s=-2`, not zero.  That alternative direction is
electrogenic and must be rejected under the voltage-independence evidence.

## 4. SR1: one-cation alternating-access `1:1:2`

### 4.1 Smallest explicit graph and exact reduction

For `C` equal to Na or K, a strict loaded-carrier alternating-access graph is

\[
 E_o
 \rightleftarrows E_oCl
 \rightleftarrows E_iCl
 \rightleftarrows E_i
 \rightleftarrows E_iCB_2
 \rightleftarrows E_oCB_2
 \rightleftarrows E_o. \tag{11}
\]

The binding edges carry the factors `a_Cl,o`, `a_Cl,i`,
`a_C,i a_B,i^2`, and `a_C,o a_B,o^2`.  Equation (4) gives

\[
 \boxed{
 \mathcal A_C^{112}
 =\log\frac{a_{Cl,o}a_{C,i}a_{B,i}^2}
 {a_{Cl,i}a_{C,o}a_{B,o}^2}
 }. \tag{12}
\]

The reversal surface is

\[
 a_{Cl,i}^{*}=a_{Cl,o}
 \frac{a_{C,i}}{a_{C,o}}
 \left(\frac{a_{B,i}}{a_{B,o}}\right)^2. \tag{13}
\]

Fast binding/release permits an exact two-conformation, two-channel reduction.
Let `O` and `I` be the outward- and inward-facing macrostates.  Let

- `a` be the effective `O -> I` Cl-loading rate;
- `b` the reverse Cl rate;
- `c` the `I -> O` `CB2`-export rate; and
- `d` its reverse.

Then

\[
 \dot p_I=(a+d)p_O-(b+c)p_I,
 \qquad p_O+p_I=1, \tag{14}
\]

and at stationarity

\[
 p_O=\frac{b+c}{a+b+c+d},\qquad
 p_I=\frac{a+d}{a+b+c+d}, \tag{15}
\]

\[
 \boxed{J=\frac{ac-bd}{a+b+c+d}},
 \qquad
 \boxed{\mathcal A_C^{112}=\log\frac{ac}{bd}}. \tag{16}
\]

One two-state edge pair is insufficient: every stationary two-state chain with
only one reversible edge has zero net current.  Two physically distinct edge
pairs are the minimal Markov representation of an exchanger cycle.

### 4.2 Audit disposition

- The Na and K variants are both electroneutral and reversible.
- At the salivary tuple, the Na-only cycle has the wrong direction, whereas
  the K cycle strongly favors chloride loading (Section 11).
- A pure Na SR1 is therefore not a physiological salivary chloride loader.
- A pure K SR1 is thermodynamically viable, but cannot alone explain the direct
  Na-transport experiment.
- A four- or six-state SR1 and a scalar reversible law with the same `J(u)` are
  exactly indistinguishable to steady whole-cell balances.  State resolution
  earns its complexity only through transient occupancy, saturation, mutation,
  or mixed-cation predictions.

## 5. SR2: Na and K cycles sharing conformational states

### 5.1 Minimal shared-state graph

The explicit graph shares the empty and Cl-loaded conformations in (11), but
has distinct `NaB2` and `KB2` loaded states.  After fast binding reduction, it
is a two-node multigraph with one Cl edge pair and two return edge pairs:

\[
 O\underset{b}{\overset{a}{\rightleftarrows}} I,
 \qquad
 I\underset{d_{Na}}{\overset{c_{Na}}{\rightleftarrows}} O,
 \qquad
 I\underset{d_K}{\overset{c_K}{\rightleftarrows}} O. \tag{17}
\]

Here the arrows in the two return channels are written in their physiological
`I -> O` direction.  Put

\[
 C=c_{Na}+c_K,\qquad D=d_{Na}+d_K,
 \qquad \Sigma=a+b+C+D.
\]

Then

\[
 p_O=\frac{b+C}{\Sigma},\qquad
 p_I=\frac{a+D}{\Sigma}. \tag{18}
\]

The branch and Cl currents are

\[
 J_s=c_sp_I-d_sp_O
 =\frac{c_s(a+D)-d_s(b+C)}{\Sigma}, \tag{19}
\]

\[
 J_{Cl}=ap_O-bp_I
 =\frac{aC-bD}{\Sigma}
 =J_{Na}+J_K. \tag{20}
\]

Each fundamental cycle obeys

\[
 \boxed{\mathcal A_s=\log\frac{ac_s}{bd_s}}
 =\log\frac{a_{Cl,o}a_{s,i}a_{B,i}^2}
 {a_{Cl,i}a_{s,o}a_{B,o}^2}. \tag{21}
\]

The third, cation-exchange cycle is not an extra free affinity:

\[
 \mathcal A_{Na}-\mathcal A_K
 =\log\frac{c_{Na}d_K}{d_{Na}c_K}
 =\log\frac{a_{Na,i}a_{K,o}}
 {a_{Na,o}a_{K,i}}. \tag{22}
\]

Thus independent arbitrary equilibrium constants for the Na and K loops would
violate a Wegscheider identity.  Binding affinities and barriers may be
cation-specific, but the completed passive transport loops still have their
equilibria fixed by the solution electrochemical potentials.

The exact entropy production is

\[
 \boxed{
 \frac{\dot S_{prod}}{R}
 =J_{Na}\mathcal A_{Na}+J_K\mathcal A_K\ge0
 }. \tag{23}
\]

Because the loops share states, one branch can run against its own isolated
affinity while another branch supplies enough dissipation to keep (23)
nonnegative.  Rejecting a shared-state SR2 merely because one `J_s A_s` is
negative would be a thermodynamic error.  This coupled operation is a genuine
difference from Task-12 independent passive scalar branches, which imposed the
sign on each branch separately.

### 5.2 Exact source and reduction boundary

Carrier stationarity gives (20), so the whole-cell source is

\[
 s_{SR2}=(-J_{Na},-J_K,J_{Na}+J_K,-2(J_{Na}+J_K)). \tag{24}
\]

If `J_Na=-J_K`, the anion source is zero and the transporter implements a pure
Na/K exchange slip cycle.  The same source subspace can be generated by two
independent scalar branches; shared conformations constrain *how the currents
vary jointly* with mixed Na/K activities, not the source vector at one state.
Pure-cation dose curves therefore do not identify sharing.  Mixed Na/K
competition, simultaneous Na/K/anion flux near reversal, or transient
occupancy is required.

SR2 must retain distinct Na and K channels.  Replacing both by `[Na]+[K]`
removes the cation-specific mutant prediction, ignores (22), and cannot
represent the T756A--T448I Na/K separation.

## 6. SR3: asymmetric sequential Na-in/K-out cycle

### 6.1 Exact graph, affinity, and reversal

The literal 2025 proposal has neutral loaded pairs and can be represented by

\[
 E_o
 \rightleftarrows E_o(Cl\,Na)
 \rightleftarrows E_i(Cl\,Na)
 \rightleftarrows E_i
 \rightleftarrows E_i(B\,K)
 \rightleftarrows E_o(B\,K)
 \rightleftarrows E_o. \tag{25}
\]

Because both `Cl-Na` and `B-K` complexes have net charge zero, this graph has a
natural microscopic explanation for weak voltage dependence.  Its completed
forward event is

\[
 Cl_o+Na_o+B_i+K_i
 \rightleftharpoons
 Cl_i+Na_i+B_o+K_o. \tag{26}
\]

Local detailed balance gives

\[
 \boxed{
 \mathcal A_{1111}
 =\log\frac{a_{Cl,o}a_{Na,o}a_{B,i}a_{K,i}}
 {a_{Cl,i}a_{Na,i}a_{B,o}a_{K,o}}
 }. \tag{27}
\]

The reversal condition is equality of the numerator and denominator.  For
example,

\[
 a_{Cl,i}^{*}=a_{Cl,o}
 \frac{a_{Na,o}}{a_{Na,i}}
 \frac{a_{B,i}}{a_{B,o}}
 \frac{a_{K,i}}{a_{K,o}}. \tag{28}
\]

The minimal fast-binding reduction again has two conformations and two edge
pairs, now a `Cl-Na` pair channel and a `B-K` pair channel.  Equations
(14)--(16) apply with `a proportional to a_Cl,o a_Na,o` and
`c proportional to a_B,i a_K,i`.

### 6.2 Exact decomposition and observational equivalence

The SR3 source has the exact decomposition

\[
 \boxed{
 (+1,-1,+1,-1)
 =(0,0,+1,-1)+(+1,-1,0,0)
 }. \tag{29}
\]

The first vector is an AE2-like `Cl/B` exchanger; the second is a `Na_o/K_i`
exchange.  The affinity decomposes in the same way:

\[
 \boxed{
 \mathcal A_{1111}
 =\log\frac{a_{Cl,o}a_{B,i}}{a_{Cl,i}a_{B,o}}
 +\log\frac{a_{Na,o}a_{K,i}}{a_{Na,i}a_{K,o}}
 }. \tag{30}
\]

Therefore a single coupled SR3 cycle and an AE2-like exchanger plus a Na/K
exchanger, constrained to equal event currents, are exactly indistinguishable
in steady species balances.  Anion-only observations cannot identify SR3 at
all: its Cl/HCO3 projection is collinear with AE2.  Demonstrating molecular
coupling requires simultaneous cation and anion fluxes while independently
changing the Na/K exchange affinity.

### 6.3 Thermodynamic feasibility versus experimental rejection

At the salivary tuple, (27) gives `A_1111=+6.643`, so the forward chloride
loading direction is strongly favorable.  Thermodynamics does **not** reject
the cycle.

The direct transport direction does.  In the 2016 low-external-Cl experiment,
AE4-expressing cells lost Cl, gained HCO3, and showed an AE4-associated rise in
intracellular Na.  Thus Na moved with HCO3 and opposite Cl.  The reverse of
(26) predicts Cl and Na both leaving the cell.  The recovery after restoring
external Cl could include Na/K-ATPase, as the authors noted, but the initial
AE4-associated Na rise under low Cl is the contradictory direction.

Accordingly:

> **SR3 disposition.**  Retain `1:1:1:1` as a thermodynamically valid 2025
> structural hypothesis, but reject it as the sole fixed salivary AE4 cycle.
> It can enter a source-complete family only together with a second transport
> mode and an independently testable condition/species switch.  A species
> difference between the directly tested mouse Na flux and human 2025
> construct is possible in principle but presently unsupported.

## 7. SR4: explicit carbonate variants

### 7.1 Required chemistry

The carbonate state is not interchangeable with bicarbonate.  The acid-base
reaction is

\[
 B\rightleftharpoons H^+ +R,
 \qquad
 K_{a2}=\frac{a_Ha_R}{a_B},
 \qquad
 a_R=K_{a2}\frac{a_B}{a_H}. \tag{31}
\]

At fast acid-base equilibrium,

\[
 \frac{a_{R,i}}{a_{R,o}}
 =\frac{a_{B,i}}{a_{B,o}}
 \frac{a_{H,o}}{a_{H,i}}. \tag{32}
\]

An implementation must either add `R_i` as a state with explicit reaction
flux, or derive it from both `B_i` and `H_i` and project carbonate transport
through the acid-base Jacobian.  Simply replacing `R` by `B`, or subtracting a
carbonate event directly from only the bicarbonate balance, violates carbon,
charge, and alkalinity accounting.

Useful conserved chemical coordinates are total inorganic carbon and carbonate
alkalinity contributions,

\[
 T_C=[CO_2]+[B]+[R],
 \qquad
 Alk_C=[B]+2[R]-[H^+]
 \quad\text{(before other buffers)}. \tag{33}
\]

Export of one carbonate removes one carbon unit but two negative-charge
equivalents.  This is not the source of one exported bicarbonate.

### 7.2 One-cation `1:1:1` carbonate mode

The alternative corresponding to a one-cation AE4 cycle is

\[
 Cl_o+C_i+R_i\rightleftharpoons Cl_i+C_o+R_o, \tag{34}
\]

with

\[
 \boxed{
 \mathcal A_C^{111,R}
 =\log\frac{a_{Cl,o}a_{C,i}a_{R,i}}
 {a_{Cl,i}a_{C,o}a_{R,o}}
 }. \tag{35}
\]

Using (32),

\[
 \mathcal A_C^{111,R}
 =\log\frac{a_{Cl,o}a_{C,i}a_{B,i}a_{H,o}}
 {a_{Cl,i}a_{C,o}a_{B,o}a_{H,i}}. \tag{36}
\]

Thus a carbonate model can look like a bicarbonate model with a fixed bias at
one fixed pH, but its reversal shifts by exactly the transmembrane proton
activity ratio.  A pH perturbation separates them.

### 7.3 Sequential `1:1:1:2` carbonate mode

The 2025 sequential carbonate proposal is

\[
 Cl_o+Na_o+R_i+2K_i
 \rightleftharpoons
 Cl_i+Na_i+R_o+2K_o. \tag{37}
\]

Its affinity is

\[
 \boxed{
 \mathcal A_{1112,R}
 =\log\frac{a_{Cl,o}a_{Na,o}a_{R,i}a_{K,i}^2}
 {a_{Cl,i}a_{Na,i}a_{R,o}a_{K,o}^2}
 }. \tag{38}
\]

or, under fast carbonate equilibrium,

\[
 \mathcal A_{1112,R}
 =\log\frac{a_{Cl,o}a_{Na,o}a_{B,i}a_{H,o}a_{K,i}^2}
 {a_{Cl,i}a_{Na,i}a_{B,o}a_{H,i}a_{K,o}^2}. \tag{39}
\]

The minimal loaded complexes, `Cl-Na` and `CO3-K2`, are each neutral.  The
mode is strongly favorable under the illustrative salivary tuple, but it
inherits the SR3 conflict with the directly measured Na direction.  In
addition, the 2025 work did not directly measure carbonate transport, D709
protonation, or a carbonate binding constant.  SR4 is therefore a
source-proposed, thermodynamically coherent hypothesis, not a calibratable
selected mechanism.

## 8. SR5: cation coordination plus allosteric regulation

The 2016 Hill fits constrain effective response curves, not molecule counts:
`EC50_Na=49 mM`, `nH_Na=2.0`, `EC50_K=62 mM`, and `nH_K=1.8`.  A noninteger
macroscopic Hill slope is not a transport stoichiometry.

A microscopic regulatory site may be represented by states
`E_{x,l,r}`, where `x` is conformation, `l` the transported load, and
`r` the regulatory occupancy.  For nontransported cation binding on side `x`,
local detailed balance requires

\[
 \frac{k_{on,x}^{(s)}a_{s,x}}{k_{off,x}^{(s)}}
 =\exp[-\beta(g_{bound,x}^{(s)}-g_{free,x})]. \tag{40}
\]

If occupancy changes a transition barrier but the regulatory ion is released
to the same side, a convenient admissible change is

\[
 k_{rs}^{+}\mapsto g_s(a)k_{rs}^{+},
 \qquad
 k_{sr}^{-}\mapsto g_s(a)k_{sr}^{-},
 \qquad g_s(a)>0. \tag{41}
\]

This changes kinetics but not the edge ratio, cycle affinity, or reversal.  If
the cation remains bound across the membrane and is released on the opposite
side, it is transported and must appear in the source vector and in (4).

The 2025 observations require cation-specific parameters: the T756A--T448I
double mutant loses Na-dependent activity while retaining substantial
K-dependent activity, and the simulations show different Na/K coordination
networks involving T448, D709, T713, and T756.  A single pooled-cation gate is
therefore rejected.  The data do not identify which microscopic edge changes,
however.  A Na-specific capacity reduction and a larger retained K capacity in
a simpler two-branch law can reproduce the steady mutant hierarchy.  Full SR5
states become identifiable only with mixed-cation competition, pre-steady-state
kinetics, or state/conformation-resolved measurements.

When regulatory binding is fast, SR5 reduces to a positive gate multiplying an
SR1/SR2 thermodynamic core.  It is then steady-state equivalent to the Task-12
allosteric scalar family and must not be credited as a uniquely state-resolved
explanation.

## 9. SR6: PKA-regulated transitions

A thermodynamically explicit minimal model duplicates the transport graph in
unphosphorylated (`U`) and phosphorylated (`P`) forms and connects matching
conformations:

\[
 E_{r,U}\rightleftharpoons E_{r,P}. \tag{42}
\]

Within each phosphoform, every ion-transport cycle must separately obey (4)
with the same transported-ion affinity.  In transition-state notation,

\[
 k_{rs}^{\sigma,+}=\kappa
 e^{-\beta(G_{rs}^{\ddagger,\sigma}-G_r^\sigma)},
 \qquad
 k_{sr}^{\sigma,-}=\kappa
 e^{-\beta(G_{rs}^{\ddagger,\sigma}-G_s^\sigma)}, \tag{43}
\]

so changing a barrier changes both directions consistently.  Phosphorylation
may alter state energies and redistribute individual edge ratios, but those
energy changes telescope around a completed transport loop; the total loop
affinity remains the ion chemical affinity.

If phosphorylation is conformation-independent, the exact minimal regulatory
reduction is

\[
 \dot p_P=k_{phos}u_{PKA}(1-p_P)-k_{deph}p_P, \tag{44}
\]

\[
 J_s=(1-p_P)J_s^U+p_PJ_s^P. \tag{45}
\]

At fixed or fast `p_P`, (45) is just an effective scalar activity law.  The
2021 experiments establish increased macroscopic activity, H89 sensitivity,
and S173 dependence; they do not identify the regulated edge or an activation
time constant.  A transition-specific SR6 therefore contains a new modelling
choice unless independently constrained.

Kinase/phosphatase cycling is chemically driven.  A fully closed model would
include ATP, ADP, and phosphate affinities in cycles that cross between `U`
and `P`.  There is no evidence that one ATP hydrolysis is tightly coupled to
each AE4 ion cycle.  Accordingly, PKA may change barriers, occupancy, active
fraction, or capacity, but it may not be used to shift AE4 reversal by changing
only one clockwise transport rate.  Such a model is rejected unless the
energetic coupling reaction and its stoichiometry are declared and supported.

## 10. Whole-cell exact invariants and source geometry

### 10.1 General current identity

Let

- `N` be NKCC1 cycle influx;
- `P` be Na/K-ATPase cycle flux;
- `K_out` be outward K-channel particle flux;
- `C` be apical Cl particle efflux;
- `J2` be forward AE2 exchange;
- `H` be NHE1 Na influx; and
- `R_B` be net intracellular bicarbonate production by the acid-base module.

For any bicarbonate AE4 source `S=(S_N,S_K,S_Cl,S_B)`, the steady balances are

\[
 0=N-3P+H+S_N, \tag{46}
\]

\[
 0=N+2P-K_{out}+S_K, \tag{47}
\]

\[
 0=2N+J_2-C+S_{Cl}, \tag{48}
\]

\[
 0=R_B-J_2+S_B. \tag{49}
\]

With the historical acid/base topology, the proton balance gives `H=R_B` at
steady state.  Using electroneutrality

\[
 S_N+S_K-S_{Cl}-S_B=0, \tag{50}
\]

equations (46)--(49) yield the exact, AE4-mechanism-independent identity

\[
 \boxed{C=P+K_{out}}. \tag{51}
\]

State resolution does not create a new steady electrical current.  It can
alter chloride delivery only by moving the state and thereby changing pump and
K-efflux support.  With explicit carbonate, (51) remains the charge/current
identity only if carbonate charge, proton chemistry, and membrane currents are
all retained consistently.

### 10.2 SR1/SR2 anion inverse

Define

\[
 A=C-2N.
\]

For SR1/SR2, put `J_Sigma=J_Na+J_K`.  Equations (48)--(49) give

\[
 A=J_2+J_\Sigma,
 \qquad
 R_B=J_2+2J_\Sigma. \tag{52}
\]

Hence

\[
 \boxed{J_\Sigma=R_B-A},
 \qquad
 \boxed{J_2=2A-R_B}. \tag{53}
\]

This exact inverse identifies only total SR2 current.  Na and K branch currents
remain equivalent under anion-only observation and require cation balance
information.

### 10.3 SR3 rank collapse

For SR3, `S_Cl=J_3` and `S_B=-J_3`, so

\[
 \boxed{A=R_B=J_2+J_3}. \tag{54}
\]

The two-row AE2/SR3 anion stoichiometric matrix has determinant zero.  No
anion-balance inverse exists.  The cation balances are

\[
 0=N-3P+H+J_3,
 \qquad
 0=N+2P-K_{out}-J_3, \tag{55}
\]

and therefore

\[
 K_{out}=N+2P-J_3,
 \qquad
 \boxed{C=N+3P-J_3}. \tag{56}
\]

The SR3 forward event dissipates the Na-in/K-out gradient while directly
displacing K-channel support by `J3`.  Whether its state feedback increases
`N+3P` enough to compensate is a whole-cell question; thermodynamics alone
does not decide it.

### 10.4 Carbonate coordinates

For SR4, the exact source must be applied to `R`, total inorganic carbon, and
alkalinity rather than to the bicarbonate balance alone.  Acid-base reaction
flux cancels from `T_C` but redistributes `B`, `R`, and `H`.  Consequently the
Task-10/12 two-anion inverse does not extend to SR4 without an observed or
modelled carbonate/pH coordinate.  A carbonate cycle and a bicarbonate cycle
may have the same apparent rate at one fixed pH, yet separate under pH
perturbation because of (32), (36), and (39).

## 11. Direction audit at salivary concentrations

The following calculation uses the 2016 thermodynamic tuple, supplemented by
the WT intracellular/extracellular pH values used in the source audit:

\[
\begin{aligned}
 &[Cl]_o=124.6,\quad [Cl]_i=50.1,\\
 &[B]_o=24.7,\quad [B]_i=19.0,\\
 &[Na]_o=151.6,\quad [Na]_i=15.5,\\
 &[K]_o=3.4,\quad [K]_i=139.5\quad {\rm mM},\\
 &pH_i=6.91,\qquad pH_o=7.40.
\end{aligned} \tag{57}
\]

For carbonate, (32) gives

\[
 \frac{[R]_i}{[R]_o}
 =\frac{19.0}{24.7}10^{6.91-7.40}
 =0.24892. \tag{58}
\]

At 310 K, `RT=2.57748 kJ/mol`.  Positive affinity favors chloride loading.

| Cycle | `A` | `Delta G=-RT A`, kJ/mol | Direction verdict |
| --- | ---: | ---: | --- |
| SR1 `1:1:2`, Na branch | `-1.89405` | `+4.8819` | Cl loading unfavorable; reject as the sole physiological loading branch |
| SR1 `1:1:2`, K branch | `+4.10065` | `-10.5694` | Cl loading favorable |
| SR3 `1:1:1:1`, Na in/K out | `+6.64342` | `-17.1233` | Strongly favorable, but fails the direct Na-direction constraint |
| SR4 `1:1:1` carbonate, Na-out branch | `-2.75995` | `+7.1137` | Cl loading unfavorable |
| SR4 `1:1:1` carbonate, K-out branch | `+3.23475` | `-8.3375` | Cl loading favorable |
| SR4 `1:1:1:2` carbonate, Na in/2K out | `+9.22944` | `-23.7887` | Strongly favorable, but unmeasured chemistry and the same Na-direction conflict |

Useful reversal checks are

\[
 [Cl]_i^{*,112,Na}=7.538\ {\rm mM},
 \qquad
 [Cl]_i^{*,112,K}=3025.0\ {\rm mM}, \tag{59}
\]

\[
 [Cl]_i^{*,1111}=3.846\times10^4\ {\rm mM},
 \qquad
 [K]_i^{*,1111}=0.1817\ {\rm mM}. \tag{60}
\]

The absurdly large SR3 chloride reversal value does not invalidate the cycle;
it shows that the imposed Na and K gradients make it a strongly dissipative
Na-in/K-out pathway.  Feasibility says nothing about whether the protein
actually uses this stoichiometry.

## 12. Minimal reductions and steady-state equivalence classes

| Family | Smallest reduction that preserves the claimed feature | What is lost by further reduction | Exact steady-state equivalence |
| --- | --- | --- | --- |
| SR1 | Two conformations with a Cl edge pair and one `CB2` edge pair | A single edge pair cannot carry stationary exchange current | Equivalent to a scalar reversible flux with the same `J(u)` |
| SR2 | Two conformations with shared Cl and distinct NaB2/KB2 edge pairs | Pooling removes Na/K competition, slip, mutant hierarchy, and multicycle coupling | At one state, equivalent to two scalar branch currents with the same `(J_Na,J_K)`; sharing is identified only across mixed-cation conditions/dynamics |
| SR3 | Two conformations with neutral `Cl-Na` and `B-K` edge pairs | Pooling destroys the opposite cation source signs | Exactly source-equivalent to AE2-like Cl/B exchange plus Na/K exchange at equal event currents |
| SR4 | SR1 or SR3 transport graph plus explicit carbonate/pH chemistry | A bicarbonate substitution destroys charge, carbon, and pH dependence | At one fixed pH it can mimic a biased bicarbonate law; pH/alkalinity perturbation separates it |
| SR5 | Base transport graph plus a binary or cation-specific regulatory occupancy only if occupancy memory is measured | Fast occupancy collapses to a positive scalar gate | Fast-gate SR5 is equivalent to a scalar allosteric core; mutant hierarchy alone identifies cation specificity, not the internal edge |
| SR6 | `U/P` copies of the minimum transport graph plus matched phosphorylation transitions | Fast/conformation-independent phosphorylation collapses to `p_P(t)` and an effective activity | At fixed `p_P`, equivalent to a weighted scalar law; transition-specific regulation requires kinetic data |

Two models with hidden Markov states are observationally equivalent at a given
steady condition whenever their resolved whole-cell sources agree,

\[
 S J(u)=\widetilde S\widetilde J(u), \tag{61}
\]

even if their internal occupancies differ.  Agreement at one state is weaker
than agreement as functions of `u`; mixed-ion and reversal perturbations are
needed to test functional equivalence.  Model selection must compare each SR
family with its reduction, not only with other large state graphs.

## 13. Explicit rejection and retention conditions

### Reject before whole-cell fitting

1. **Nonzero source charge.**  Reject any proposed cycle with `q_s != 0`
   unless a measured counterion/current is added and the voltage term in (2)
   is retained.
2. **Wegscheider violation.**  Reject a graph if any independent cycle fails
   (4), or if it sustains flux at equal electrochemical activities without a
   declared energy source.
3. **Wrong unicyclic direction.**  Reject a passive unicycle if its simulated
   current has the opposite sign from its affinity.
4. **Na-only physiological loader.**  Reject pure-Na SR1 under the stated
   salivary concentrations; its chloride-loading affinity is negative.
5. **SR3 as the sole fixed AE4 mode.**  Reject it because its reverse assay
   direction contradicts the directly measured AE4-associated Na influx.
   Thermodynamic favorability does not cure that contradiction.
6. **Carbonate without carbonate chemistry.**  Reject any SR4 implementation
   that aliases `CO3` to `HCO3`, omits pH/protonation, or applies the carbonate
   event only to the bicarbonate balance.
7. **Pooled Na/K state.**  Reject it for SR2/SR5 because it cannot encode the
   T756A--T448I Na/K mutant separation or the exact cation-exchange affinity
   (22).
8. **Hill slope as transported count.**  Reject identification of `nH=1.8` or
   `2.0` with a fractional or exact ion stoichiometry without independent flux
   ratios.
9. **One-way PKA bias.**  Reject changing only one clockwise transport rate if
   it shifts reversal, unless ATP coupling and its stoichiometry are explicit.
10. **Dynamic binding without bound-ion balance.**  Reject a state ODE that
    omits transient substrate storage while claiming exact conservation.

### Retain, but with qualification

- Retain SR1-K and thermodynamically consistent SR2 as viable chloride-loading
  mechanisms; direct Na transport requires an Na-capable mode somewhere in the
  source-complete family.
- Retain SR2 shared-state coupling because it can drive one branch uphill while
  satisfying total entropy production.  Do not call it identified by the
  existing pure-cation dose curves.
- Retain SR3 and sequential SR4 as explicit 2025 hypotheses, not accepted sole
  mechanisms.  Their decisive missing experiment is simultaneous AE4-specific
  Na, K, Cl, and HCO3/carbonate flux under one reversible protocol.
- Retain SR5 cation-specific coordination, because the mutant hierarchy rules
  out pooling; keep microscopic rate placement nonidentified.
- Retain SR6 only as a thermodynamically neutral rate/occupancy regulation
  unless a transport-coupled phosphorylation energy cycle is demonstrated.

## 14. Single most discriminating transporter experiment

The cleanest discriminator is a reversible, AE4-isolated experiment that
simultaneously measures Na, K, Cl, and pH/total inorganic carbon while imposing
independently variable Na and K gradients on both sides.  In particular, under
low external Cl:

- SR1/SR2 predicts the active cation moves with bicarbonate and opposite Cl;
- the literal SR3/SR4 sequential proposal predicts Na moves with Cl while K
  moves with bicarbonate;
- mixed-cation titration distinguishes a shared-state SR2 from independent
  branches through cross-saturation and the possible Na/K slip current; and
- repeating at two controlled pH values separates HCO3 from carbonate through
  the exact proton-ratio shift in (32).

The experiment must report simultaneous flux ratios, not only an
alkalinization slope.  That one protocol would decide the principal
thermodynamic equivalence classes above without reference to a whole-gland
secretion phenotype.

## Primary and project sources

- Peña-Münzenmayer et al. (2016), *J Gen Physiol*,
  [DOI 10.1085/jgp.201611571](https://doi.org/10.1085/jgp.201611571),
  [open full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC4845690/).
- Peña-Münzenmayer et al. (2021), *Am J Physiol Gastrointest Liver
  Physiol*,
  [DOI 10.1152/ajpgi.00145.2021](https://doi.org/10.1152/ajpgi.00145.2021),
  [open full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC8887885/).
- Catalán et al. (2025), *Am J Physiol Cell Physiol*,
  [DOI 10.1152/ajpcell.00346.2024](https://doi.org/10.1152/ajpcell.00346.2024),
  [open full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC12463136/).
- Task-10 exact balance reduction:
  `analysis/10_identifiability_discrimination/reduction.md` and
  `analysis/10_identifiability_discrimination/theory.md`.
- Task-11 model/Markov concordance:
  `analysis/11_forensic_reconstruction/paper_code_concordance.md`.
- Task-12 source vectors, thermodynamic calculations, and fixed-chassis
  identities:
  `analysis/12_ae4_mechanism_reconstruction/ae4_candidate_family.md`,
  `chassis.md`, and `structural_explanation.md`.

The 2018 printed four-state-to-polynomial AE4 reduction is not used as a
microscopic derivation here: Task 10 showed that its rational steady-state
denominator and printed state equation do not algebraically reduce to the
reported constant-activity polynomial.  Equations (11)--(45) are independent
thermodynamically constrained constructions.
