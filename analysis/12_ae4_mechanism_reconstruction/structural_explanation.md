# Structural explanation and mechanism localization

## Scope and evidence boundary

This note explains the candidate results using exact balances before using
numerical trajectories. It distinguishes four evidence levels:

- **Exact:** conservation, electroneutrality, steady balance algebra, and the
  fixed-chassis water/lumen reduction.
- **Primary experimental evidence:** Peña-Münzenmayer et al. (2015), DOI
  `10.1074/jbc.M114.612895`; Peña-Münzenmayer et al. (2016), DOI
  `10.1085/jgp.201611571`; and Peña-Münzenmayer et al. (2021), DOI
  `10.1152/ajpgi.00145.2021`. Together these constrain the knockout protocol,
  cAMP/PKA regulation, voltage independence, AE4-associated Na and K movement,
  cation dose responses, and the three stoichiometries the authors explicitly
  evaluated.
- **Historical implementation evidence:** the internally closed seven-state
  baseline and its fitted Na-only AE4 source vector established in Task 11.
- **New analysis:** the general balance inverse, the fixed-state forcing
  obstruction, the pooled-partition identity, and the chloride-to-flow
  monotonicity result below.

The three ratios considered in the primary paper were thermodynamic
possibilities, not directly measured stoichiometries. Hill slopes likewise do
not measure transported ion count. No knockout datum is used in any derivation
or calibration in this note; the observed knockout ratio enters only in the
post-freeze validation calculation.

The primary study also left a carbonate-containing mode and a low-chloride
bicarbonate self-exchange mode unresolved. The seven-state chassis has no
carbonate state, and the self-exchange coefficients are unknown. They are
recorded as source-viable but untestable chemistry; silently replacing
carbonate by bicarbonate or inventing the unknown coefficients would not be an
AE4-only reconstruction of this chassis.

## 1. General electroneutral AE4 cycle

Let one positive AE4 event import $c$ chloride ions and export $m$
monovalent cations and $b$ bicarbonate ions. Its intracellular source vector
in the order $(Cl,Na,K,HCO_3)$ is

\[
s_4=(c,-m_{Na},-m_K,-b),\qquad m_{Na}+m_K=m.
\]

The net charge moved into the cell is

\[
-c-m+b.
\]

Therefore the exact electroneutrality condition is

\[
\boxed{b=c+m}. \tag{1}
\]

This is also why the electrical-potential terms cancel from the cycle
affinity. The primary paper's three explicitly considered possibilities all
satisfy (1):

| $Cl:Cation:HCO_3$ | $HCO_3/Cl$ | cation/Cl | Intracellular source per event |
| --- | ---: | ---: | --- |
| 1:1:2 | 2 | 1 | $(+1,-m_{Na},-m_K,-2)$ |
| 1:2:3 | 3 | 2 | $(+1,-m_{Na},-m_K,-3)$ |
| 2:1:3 | 1.5 | 0.5 | $(+2,-m_{Na},-m_K,-3)$ |

A purported `1 Cl : 2 HCO3` cycle with no transported cation has charge
residual $+1$ per event. It cannot be inserted into this electroneutral,
voltage-independent chassis without an additional charge-carrying species or
current. A cation-gated but non-cation-transporting cycle can instead be
electroneutral only with equal anion counts, for example `1 Cl : 1 HCO3`.
That pure-gating construction conflicts with the AE4-associated Na and K
movement in the primary experiments. A transported-cation cycle with an
*additional* allosteric cation gate remains viable.

## 2. Exact steady exchanger decomposition

Use positive amount fluxes

- $C$: apical chloride efflux,
- $N$: NKCC1 cycle influx,
- $J_2$: forward AE2 cycle flux,
- $J_4$: forward AE4 event flux, and
- $B$: buffer production of intracellular bicarbonate equivalents.

The steady chloride and bicarbonate balances give

\[
A:=C-2N=J_2+cJ_4,
\qquad
B=J_2+bJ_4. \tag{2}
\]

Since $b-c=m>0$, this system has the exact inverse

\[
\boxed{J_4=\frac{B-A}{m}},
\qquad
\boxed{J_2=\frac{bA-cB}{m}}. \tag{3}
\]

For nonnegative forward exchanger rates, the exact feasible wedge is

\[
\boxed{A\le B\le \frac bc A}. \tag{4}
\]

The candidate-specific formulas are:

| Stoichiometry | $J_4$ | $J_2$ | Forward wedge |
| --- | --- | --- | --- |
| 1:1:2 | $B-A$ | $2A-B$ | $A\le B\le2A$ |
| 1:2:3 | $(B-A)/2$ | $(3A-B)/2$ | $A\le B\le3A$ |
| 2:1:3 | $B-A$ | $3A-2B$ | $A\le B\le1.5A$ |

Thus Na/K weighting, Hill cooperativity, saturation, and positive allosteric
gating cannot change a candidate's anion-balance ray. Only transported
stoichiometry can change $b/c$. Conversely, an anion-only measurement cannot
distinguish Na from K branches that share the same $c:m:b$ counts.

### Cation/current invariant

Let $P$ be Na/K-pump cycle flux and $K_{out}$ be outward Ca-activated K
particle flux. Adding the steady Na and K balances, using $B=J_2+bJ_4$,
$A=J_2+cJ_4$, and Eq. (1), gives

\[
\boxed{C=P+K_{out}}. \tag{5}
\]

This identity is independent of Na/K partition and of which of the three
electroneutral stoichiometries is used. AE4 changes chloride secretion only by
moving the state to a point at which the pump/K-efflux support and the
NKCC1/AE2 loading balances settle at different values. It does not create an
independent steady current source.

For parallel one-cation 1:1:2 branches, let $J_{Na}$ and $J_K$ be their
forward cycle rates and let $H$ be NHE1 influx. The separate steady cation
balances sharpen Eq. (5):

\[
0=N-3P+H-J_{Na},
\qquad
0=N+2P-K_{out}-J_K, \tag{5a}
\]

so

\[
\boxed{C=N+3P-J_K=2N+H-J_{Na}-J_K}. \tag{5b}
\]

Equation (5b) is counterintuitive but important. A K-coupled AE4 cycle imports
chloride while exporting K, yet it does not supply an independent steady
secretory current: at fixed NKCC1 and pump support, its outward K load reduces
the K-channel current available in Eq. (5) by the same amount. A recruited K
branch raises $C$ only if the induced increase in $N+3P$ exceeds $J_K$.

This penalty is not peculiar to 1:1:2. For any of the candidate
stoichiometries, let $X_K$ denote the outward AE4 **K-particle** flux, including
all $m$ cations per event. The K balance and (5) give

\[
K_{out}=N+2P-X_K,
\qquad
\boxed{C=N+3P-X_K}. \tag{5c}
\]

Changing from 1:1:2 to 2:1:3 changes the anion balance ray and chloride per
event, but it retains one exported cation per event and therefore does not
remove this direct K-current displacement. For 1:2:3, two exported K particles
per pure-K event enter $X_K$.

## 3. Fixed-state forcing obstruction

The historical WT state is closed by the non-AE4 chassis plus the Na-only
AE4 contribution

\[
J_0(Cl,Na,K,HCO_3)=J_0(+1,-1,0,-2),
\qquad J_0=0.03888333868
\]

in the historical flux unit. Holding the state and every non-AE4 term fixed,
the non-AE4 residual is exactly the negative of this vector. A replacement
AE4 law supplies only a scalar event rate $j$ times its source vector $s$.
The same state can close if and only if

\[
js=J_0(1,-1,0,-2). \tag{6}
\]

Equation (6) requires the two signatures to be collinear. Changing a positive
capacity, Hill coefficient, saturation denominator, or allosteric gate changes
only $j$, not $s$, and therefore cannot repair a species-signature
mismatch.

For the 2018 intracellular mole-fraction partition,

\[
\alpha=25/145,\qquad \beta=120/145,
\]

the chloride and bicarbonate rows force $j=J_0$. The resulting Na and K
residuals are then

\[
(r_{Na},r_K)=(+\beta J_0,-\beta J_0)
=(+0.03217931472,-0.03217931472). \tag{7}
\]

Alternatively, the K row would require $j=0$, while the chloride row requires
$j=J_0\ne0$. Exact fixed-state closure is impossible. This does **not**
reject Na/K-coupled biology; it proves that a Na-only calibrated chassis
cannot accept a different cation signature at precisely the same fitted state
without state relaxation or a declared non-AE4 recalibration. The primary
comparison correctly reports the residual rather than hiding it.

### Knockout trajectory invariance

Every candidate enters the fixed chassis only through its AE4 source term.
After setting AE4 activity to zero, that term is identically zero regardless
of its stoichiometry, affinity, Hill exponent, saturation, or partition. Hence,
for a common initial state and calcium protocol, all candidates have the same
knockout initial-value problem

\[
\dot u=F_0(u,t),\qquad u(0)=u_0.
\]

Whenever this problem has a unique solution, the entire AE4-knockout
trajectory is exactly candidate-independent. If the knockout steady state is
unique, it is candidate-independent even when reached from different initial
states. Candidate-to-candidate variation in a KO/WT ratio must therefore come
from the WT denominator, not from a different knockout mechanism.

This has two consequences. First, a mechanism cannot be credited with
"lowering the KO" on this chassis; it can only support a higher WT response.
Second, a KO/WT flow ratio cannot discriminate two kinetic laws that share the
same calibrated WT trajectory or WT equilibrium. A ratio obtained only by an
unconstrained inflation of stimulated WT flow is not a reconstruction unless
that WT response also satisfies independent WT evidence.

### Source-required stimulus gate and its exact threshold

The fixed-row theorem does **not** justify restricting AE4 activity to one
time-independent scale. The held-out 2015 experiment applied both `0.3 µM`
carbachol and `5 µM` isoproterenol for 10 min. It found total saliva
`35 ± 4.7%` lower in AE4-null glands, with similar secretion for the first
2--3 min followed by a sustained deficit. The same study tied the phenotype
to β-adrenergic/cAMP-dependent chloride reuptake. The 2021 primary study then
showed directly that isoproterenol activates native AE4 in AE2-null acini,
that H89 prevents this activation, and that forskolin or catalytic PKA
activates heterologous AE4 through a mechanism requiring AE4 residue S173.
The native assay response is of order a few-fold, rather than an arbitrary
free density change.

Represent a source-constrained stimulus input by $u(t)\in[0,1]$ and a gate

\[
g_F(t)=1+u(t)(F-1),\qquad F\ge1. \tag{7a}
\]

Normalizing $u=0$ at rest preserves the resting AE4 calibration. A common
positive gate multiplying a reversible branch preserves its equilibrium and
thermodynamic sign. It also preserves its species signature, so it cannot
repair the mixed-cation residual in Eq. (7) at the frozen rest row. A distinct
K branch that is silent at rest and PKA-recruited during stimulation can evade
that particular obstruction, but branch-specific recruitment is a new
modelling decision: the primary PKA study established activation of AE4, not
selective recruitment of K stoichiometry.

Let $R_0=H_{KO}/H_{WT}(1)$ be a frozen model's current knockout/WT ratio for
any positive endpoint or integrated observable $H$, and let $R_*=0.65$ be the
nominal held-out target. Since the KO trajectory is independent of the WT gate
for a common initial condition, any successful gate must satisfy exactly

\[
\frac{H_{WT}(F)}{H_{WT}(1)}\ge\frac{R_0}{R_*}. \tag{7b}
\]

For the static C1 lineage row, this requires a `53.624%` increase in endpoint
WT flow and a `54.714%` increase in the post-step WT integral. If $H_{WT}$ is
positive and monotone in $F$, the required mean logarithmic gate elasticity is

\[
\bar E_H
=\frac{\log[H_{WT}(F)/H_{WT}(1)]}{\log F}
\ge\frac{\log(R_0/R_*)}{\log F}. \tag{7c}
\]

At $F=3$, the endpoint and integral lower bounds are `0.391` and `0.397`;
at $F=1.25$, they are `1.924` and `1.955`. These are threshold diagnostics,
not proofs of response monotonicity. The nonlinear pump--leak network can
feed back on the AE4 driving force, so only a knockout-blind simulation of the
source-bounded gate can determine whether the threshold is crossed.

Treating the paper's literal `35 ± 4.7%` report as a validation band gives
$R_*\in[0.603,0.697]$ (not a reconstructed confidence interval). Across that
band, the required endpoint WT gain is `43.3--65.6%` and the integral gain is
`44.3--66.8%`. At $F=3$, the corresponding endpoint elasticity interval is
`0.327--0.459`; a gated prediction should be compared with this band rather
than judged against exactly `0.650`.

### State relaxation: what fixed-row incompatibility does not prove

Let $\lambda$ interpolate a candidate's cation structure and let the unchanged
non-AE4 chassis plus that AE4 law define $F(u,\lambda)=0$ at rest. If
$F_u(u_0,0)$ is nonsingular, the implicit-function theorem gives a local
candidate-specific equilibrium branch

\[
\frac{du^*}{d\lambda}
=-F_u(u^*,\lambda)^{-1}F_\lambda(u^*,\lambda). \tag{7d}
\]

Thus failure to close the *same* historical row does not exclude a nearby
steady state with shifted Na, K, Cl, bicarbonate, and volume. This escape is
real in the current implementation. Solving all seven resting equations for
the projected-scale C3 explicit Na/K law gives the positive root

\[
u^*=(118.76, 5.618, 27.55, 30.36, 114.71, 46.31, 10.21)
\tag{7e}
\]

in chassis state order, with pH `6.908` and maximum raw residual
$8.26\times10^{-15}$. Relative to the historical row, intracellular Na shifts
about `+21%`, K `-4.4%`, Cl `-7.4%`, and height `-4.0%`. This root proves that
the exact incompatibility is calibration-row-specific. It is not, however, an
admissible replacement WT row under the KO-blind two-SEM
resting-concentration screen, fixed before held-out evaluation:
$Cl_i=46.313$ mM is `-2.525` reported SEM from the WT target
$50.10\pm1.50$ mM, just outside the two-SEM gate. It also does not rescue the
phenotype: stimulation from the moved C3 root gives AE4-KO/WT endpoint ratio
about `1.0295`, the wrong direction, while AE2-KO/WT is about `0.9846`.

One KO-blind state-relaxation slice was applied to all three reversible C6
stoichiometries. For each discrete stoichiometry, an outer continuation
selected the AE4 capacity that recovers the independently published `1.3 pL`
resting volume, then all seven state equations were solved. These particular
volume-normalized roots close numerically and preserve branchwise
thermodynamic sign, but all fail the independent resting-$Cl_i$ gate:

| C6 stoichiometry | Fixed capacity | $Cl_i$ (mM; SEM) | pH | Na, K branch flux | Net flux | Max raw residual |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1:1:2 | $3.22984\times10^{-8}$ | 43.476; -4.416 | 7.0073 | -0.02725, +0.05093 | +0.02369 | $1.32\times10^{-13}$ |
| 1:2:3 | $8.64183\times10^{-12}$ | 40.675; -6.283 | 7.0415 | -0.02646, +0.03857 | +0.01211 | $7.49\times10^{-14}$ |
| 2:1:3 | $2.01702\times10^{-11}$ | 45.739; -2.907 | 6.9761 | -0.01686, +0.03571 | +0.01885 | $1.13\times10^{-13}$ |

Here negative Na-branch flux means that the Na branch runs in reverse while
the K branch runs forward; the branch flux and log-forward/reverse affinity
signs agree. These rows test the additional *exact-volume* convention, not the
entire candidate family. A subsequent KO-blind capacity continuation found
connected rest-screen-passing 2:1:3 bands away from that exact-volume point.
The table gives the frozen band and the state range over its five KO-blind
prediction quantiles, fixed before held-out evaluation:

| Kinetic core | Connected capacity band | Volume (pL) | Resting $Cl_i$ SEM | Resting pH SEM |
| --- | ---: | ---: | ---: | ---: |
| C6 2:1:3 mass action | $[3.5105005\times10^{-12},5.9165257\times10^{-12}]$ | 1.37486--1.41303 | -1.99977 to -1.81374 | +1.63977 to +1.99962 |
| C7 2:1:3 affinity law | [0.0033498713, 0.0072388250] | 1.33246--1.41460 | -1.99963 to -1.74714 | +1.12010 to +1.99905 |
| C8 saturated 2:1:3 | [0.0344779857, 0.0378955063] | 1.39480--1.41074 | -1.99990 to -1.90849 | +1.85804 to +1.99986 |

All capacities were selected without knockout information. The 15 frozen
roots meet the two-SEM chloride/pH gate and an explicit post-hoc $\pm10\%$
volume sensitivity band while closing the resting equations to at most
$4.82\times10^{-13}$ in raw chassis units. That volume band is a modelling
sensitivity, not a measured uncertainty interval or a predeclared primary-data
gate. Failure of the exact-volume roots therefore cannot support a family-wide
exclusion. All three bands are boundary-sensitive: among the prediction
points the nearest margins to a two-SEM boundary are only `0.00038` SEM for
C6, `0.00037` SEM for C7, and `0.00010` SEM for C8.

Within this declared continuation, every surviving thermodynamic interval
uses 2:1:3. That makes 2:1:3 numerically necessary for *rest compatibility in
this tested source-default capacity slice*, not an experimentally identified
AE4 stoichiometry and not a sufficient property for the knockout phenotype.
The continuation fixes the assay-default Na:K weights `1.5:1.6`,
$K_{eq}=1$, saturation parameter $\sigma=1$, and bath bicarbonate at `21 mM`;
those axes are not generally scalar-equivalent to capacity. It therefore does
not exclude untested combinations along those axes or constitute a global
dynamic-family theorem.

A rigorous candidate-specific-rest test should continue roots from the
historical row, solve all seven unscaled equations in positivity-preserving
coordinates, use multiple starts to detect other branches, require physical
$H_i$, $CO_{2,i}$, and flow, check Jacobian stability, and compare every moved
state against independent WT concentration uncertainty and explicitly labelled
volume sensitivity before stimulating it. Its capacity and stimulus fold must
remain fixed by WT/assay information; the knockout phenotype cannot select the
root, branch, or gate.

## 4. Thermodynamic direction and the necessary role of K

For the positive cycle defined above, ideal-solution thermodynamics gives

\[
\frac{\Delta G}{RT}
=c\ln\frac{Cl_i}{Cl_o}
+m\ln\frac{C_o}{C_i}
+b\ln\frac{HCO_{3,o}}{HCO_{3,i}}. \tag{8}
\]

Negative $\Delta G$ favors chloride influx. At equilibrium,

\[
Cl_i^*=Cl_o
\left[
\left(\frac{C_i}{C_o}\right)^m
\left(\frac{HCO_{3,i}}{HCO_{3,o}}\right)^b
\right]^{1/c}. \tag{9}
\]

Using the exact physiological tuple used by Peña-Münzenmayer et al. (2016),
$Cl_o=124.6$, $Cl_i=50.1$, $HCO_{3,o}=24.7$,
$HCO_{3,i}=19$, $Na_o=151.6$, and $Na_i=15.5$ mM, gives:

| Stoichiometry | $\Delta G/RT$, Na only | $Cl_i^*$, mM | $Na_i^*$, mM | $HCO_{3,i}^*$, mM |
| --- | ---: | ---: | ---: | ---: |
| 1:1:2 | +1.8940 | 7.538 | 103.016 | 48.982 |
| 1:2:3 | +4.4368 | 0.593 | 142.487 | 83.377 |
| 2:1:3 | +1.2453 | 26.879 | 53.848 | 28.776 |

All three Na-only cycles run in the wrong thermodynamic direction at the
physiological point. The small differences from the paper's rounded thresholds
(`7/0.5/26`, `104/143/54`, and `49/84/29` mM) are rounding only.

For a nonselective 1:1:2 calculation with
$(Na+K)_i=(Na+K)_o=155$ mM, the cation term cancels:

\[
\Delta G/RT=-0.38636,
\qquad Cl_i^*=73.728\ {\rm mM},
\qquad HCO_{3,i}^*=15.662\ {\rm mM}. \tag{10}
\]

The observed state therefore favors chloride influx. This yields a sharp
structural statement:

> **Necessary thermodynamic statement.** Under passive reversible transport,
> the primary-paper concentrations, one of the three source-left-viable
> electroneutral stoichiometries, and no external energy source, Na-only AE4
> cannot mediate physiological chloride loading. A sufficiently strong K
> branch (or another explicitly modelled energy source) is necessary.

For parallel 1:1:2 Na/K mass-action branches with equal standard constants,
the inferred $K_i=139.5$, $K_o=3.4$ mM values make the Na branch negative
and the K branch positive. If their unit-capacity rates are $j_{Na}<0$ and
$j_K>0$, and nonnegative branch capacities are normalized to sum to one, then

\[
J=(1-f_K)j_{Na}+f_Kj_K>0
\quad\Longleftrightarrow\quad
f_K>\frac{-j_{Na}}{j_K-j_{Na}}=0.3895. \tag{10a}
\]

This is a new threshold calculation, not a measured selectivity. Equal branch
capacities, which reproduce the simple pooled propensity for a one-cation
cycle, satisfy it.

### Why the published pooled partition is not a reversible branch law

For a one-cation cycle, write identical branch propensities as

\[
J_{Na}=FNa_i-RNa_o,\qquad
J_K=FK_i-RK_o.
\]

Let $f_i=Na_i/(Na_i+K_i)$. Allocating the pooled *net* flux by $f_i$, as
in the 2018 balance equations, differs from the explicit Na branch by

\[
J_{Na}-f_i(J_{Na}+J_K)
=R\,\frac{Na_iK_o-Na_oK_i}{Na_i+K_i}. \tag{11}
\]

This error vanishes only when the reverse term is zero or the intracellular
and extracellular Na fractions are equal. Neither is generally true. The
printed C2 law can be reproduced as a historical/published empirical
candidate, but it is not a species-resolved thermodynamic formulation.

## 5. What kinetics can and cannot change

A reversible passive law can be written generally as

\[
J_4=L(u)\,[1-\exp(\Delta G/RT)],\qquad L(u)>0, \tag{12}
\]

or as a forward-minus-reverse mass-action numerator divided by a positive
saturation denominator. Positive occupancy, saturation, and allosteric factors
change magnitude but do not change the zero or sign imposed by $\Delta G$.
Consequently:

- saturation cannot rescue the wrong direction of Na-only transport;
- a Hill factor cannot rescue it unless the candidate abandons thermodynamic
  sign control;
- Hill slopes near 2 can represent multiple transported cations (the 1:2:3
  cycle), cooperative occupancy with one transported cation, or transport plus
  an allosteric gate; and
- knockout flow alone cannot distinguish those placements because they can
  share the same zero-flux knockout and, for 1:1:2 variants, the same balance
  signature.

The primary AE4-associated cation-flux experiments make a *pure* gating-only
candidate source-inadmissible, but do not establish that a cation accompanies
every cycle or separate one transported cation plus a gate from a genuinely
two-cation cycle. That distinction requires simultaneous AE4-specific Na/K
flux and anion-flux measurements near reversal, not another flow-only
knockout.

### Why a capacity fold need not produce a steady-flux fold

Equation (3) makes the network-demand constraint explicit. Define the
candidate's steady event-flux demand from the non-AE4 balance aggregates as

\[
D(u):=\frac{B(u)-A(u)}{m}. \tag{12a}
\]

If a common activity factor $G>0$ multiplies a kinetic kernel $\phi(u)$, every
steady state on a smooth candidate-specific branch must satisfy both

\[
J_4=G\phi(u^*)
\quad\text{and}\quad
\boxed{G\phi(u^*)=D(u^*)}. \tag{12b}
\]

Thus a threefold capacity factor is not a threefold steady event-flux factor.
The state moves until the kinetic propensity supplies the event flux demanded
by the other balances. Differentiating (12b) along a regular steady branch
with respect to $\log G$ gives

\[
\frac{dJ_4}{d\log G}
=\nabla D\cdot\frac{du^*}{d\log G}
=J_4+G\nabla\phi\cdot\frac{du^*}{d\log G}. \tag{12c}
\]

The direct $+J_4$ capacity effect is therefore cancelled by the state-induced
change in propensity except to the extent that the non-AE4 network demand
$D(u^*)$ itself changes. If $D$ is locally constant, the cancellation is exact:
$dJ_4/d\log G=0$ and $d\log|\phi|/d\log G=-1$. More generally, a nearly
inelastic demand produces a nearly capacity-invariant steady flux even though
the microscopic turnover capacity changed substantially.

For a passive reversible kernel
$\phi=L_0(u)[1-\exp(\Delta G/RT)]$, any bounded family of steady states with
bounded $D$ satisfies $\phi=D/G\to0$ as $G\to\infty$. If $L_0$ is also bounded
away from zero along that branch, the transporter moves toward its reversal
condition while carrying the demand-set flux. Conversely, if $G\to0$ and
$\phi$ remains bounded, a nonzero demand cannot close: the state must move
enough to change $D$ or the steady branch terminates. These are capacity
limits, not claims that the full transient response is monotone.

The common-gate and K-recruitment diagnostics make the signature consequence
especially transparent. For a forward basal C1 flux $j>0$ and full-stimulus
fold $F$, their intracellular sources in $(Na,K,Cl,HCO_3)$ order are

\[
s_{common}=(-Fj,0,+Fj,-2Fj),
\qquad
s_{Krecruit}=(-j,-(F-1)j,+Fj,-2Fj). \tag{12d}
\]

They have the same instantaneous anion source at the same state but differ by

\[
s_{Krecruit}-s_{common}
=((F-1)j,-(F-1)j,0,0), \tag{12e}
\]

the same Na-to-K substitution direction that causes the fixed-row obstruction.
At a steady state, Eq. (5b) gives $C=N+3P-J_K$. Therefore recruited forward K
transport directly displaces K-channel current and lowers supported chloride
delivery unless its induced state change increases $N+3P$ by more than
$J_K$. This is why adding a K branch can fix passive AE4's thermodynamic
direction yet worsen secretion on this particular pump/channel chassis: the
thermodynamic and whole-network criteria are different constraints.

## 6. Exact chloride-to-flow map on the fixed chassis

The seven-state chassis has, at steady cell volume,

\[
q_a=b_1(\Omega_l-\Omega_i),\qquad
q_b=b_2(\Omega_i-I_e),\qquad q_a=q_b,
\]

and

\[
q_t=b_3(\Omega_l-I_e),\qquad
Q=q_a+q_t.
\]

With luminal electroneutrality $Cl_l=Na_l+K_l$, define

\[
\kappa=\frac{b_1b_2}{b_1+b_2}+b_3>0,
\qquad d=x_l-I_e.
\]

Then

\[
\boxed{Q=\kappa(2Cl_l+d)},
\qquad
\boxed{C=QCl_l}. \tag{13}
\]

Eliminating $Cl_l$ gives the positive branch

\[
\boxed{
Q(C)=\frac{\kappa d+\sqrt{(\kappa d)^2+8\kappa C}}{2}
}, \tag{14}
\]

with

\[
\boxed{
\frac{dQ}{dC}=\frac{2\kappa}{2Q-\kappa d}>0
}. \tag{15}
\]

Thus steady flow is strictly monotone in apical chloride delivery on this
fixed water chassis. Voltage and cation partition matter only insofar as they
change the attainable $C$; Eq. (5) shows that this occurs through pump and K
efflux support.

For the frozen parameters,

\[
\kappa=5.71406152\times10^{-5},\quad
d=-239.2998643,
\]

and the historical WT row has

\[
Q=5.31415477\times10^{-4},\quad
C=0.06605494376
\]

in code-native units. A nominal held-out flow ratio $r_Q=0.65$ requires,
from (14),

\[
\frac{C_{KO}}{C_{WT}}
=r_Q\frac{r_QQ_{WT}/\kappa-d}{Q_{WT}/\kappa-d}
=0.641489, \tag{16}
\]

or a `35.851%` loss of apical chloride flux. By contrast, the historical
seven-state endpoint ratio `0.998556` corresponds to only about `0.150%`
chloride-flux loss under the same steady closure.

Combining (2), (5), and (13) gives an exact necessary-and-sufficient target
condition under the fixed steady chassis:

\[
\boxed{
C_{WT}-C_{KO}
=cJ_{4,WT}+2(N_{WT}-N_{KO})+(J_{2,WT}-J_{2,KO})
}
\tag{17}
\]

must equal the chloride loss implied by (16), equivalently

\[
(P+K_{out})_{WT}-(P+K_{out})_{KO}
\]

must equal that same amount. Equation (17), not a rate-law label, is the
precise threshold a successful mechanism must cross.

Equations (13)--(17) are steady-state identities. The reported numerical
phenotype is also evaluated at the fixed simulation endpoint and over a
post-stimulus time integral. Those transient observables inherit this exact
map only if the relevant volume and lumen derivatives have become negligible;
the map is used here to localize the limiting steady pathway, not to relabel
every endpoint value as a steady state.

## 7. Static mechanism-family result and localization

The numerical values in this section are frozen machine rows from
`results/12_ae4_mechanism_reconstruction/model_comparison.csv`,
`resting_equilibria.csv`, `genotype_resting_equilibria.csv`,
`rest_capacity_envelope.csv`, and `pka_activation_envelope.csv`. Solver
crosschecks are recorded in `solver_crosscheck.json`; `summary.json` joins the
held-out value only after the KO-blind candidate and grid freeze.

### Frozen Round-1 result

Calibration used the WT balance vector and independent assay constants, then
froze every AE4-specific quantity before the knockout simulations. The common
AE4-knockout endpoint flow was `0.001358119885421` in code-native units, as
predicted by the knockout-invariance result above. The discriminating results
were therefore all changes in the WT denominator:

| Candidate class | WT baseline result | Held-out AE4-KO result | Structural interpretation |
| --- | --- | --- | --- |
| Historical C1, Na-only 1:1:2 | Closes to a raw residual of $1.78\times10^{-12}$ | endpoint KO/WT `0.9985561402`, a `0.144386%` loss; stimulated-integral KO/WT `1.0056381561`, a `0.563816%` *gain* | Lineage control only: primary-assay inadmissible and not thermodynamic |
| Published pooled C2 | AE4-vector residual `0.0328154`; raw WT residual `0.001143` | endpoint `1.029056`, integral `1.036809` | Ratios use a drifting, nonclosed WT and are not predictions |
| Explicit parallel C3 | AE4-vector residual `0.0331963`; raw WT residual `0.001157` | endpoint `1.029524`, integral `1.037283` | Correct species branching exposes rather than repairs the cation mismatch |
| Cooperative C4 | Nonnegative WT projection selects zero AE4 capacity; raw residual `0.00270964` | endpoint and integral `1` | Cooperativity cannot rotate the source signature; this occupancy law is also not affinity-certified |
| C5 and reversible C6--C8 defaults | Nonnegative WT projection selects zero AE4 capacity; raw residual `0.00270964` | endpoint and integral `1` | Thermodynamic sign is valid, but at the fixed row the fitted direction cannot supply the historical residual |

The apparent C2/C3 increases do not establish a biological inverse phenotype:
both rows fail the protocol-defined $10^{-8}$ WT closure gate by many orders of
magnitude. Among rows that actually preserve the baseline, C1 is the
closure-only lineage control, not a member of the jointly WT-closing and
primary-assay-admissible **static fixed-row** set. That set is empty. It is
therefore incorrect to call `0.144386%` the maximum of an admissible set. It is the
largest endpoint decrement in the closure-only lineage comparison and is
about 240 times smaller than the nominal 35% target; its integrated metric
even changes in the opposite direction.

### Stimulus-gated diagnostic status

The frozen production integrations give the four diagnostics below. None is
used to fit a candidate or to derive an analytical statement.

| Candidate, $F=3$ | Endpoint KO/WT | Post-step integral KO/WT | Interpretation |
| --- | ---: | ---: | --- |
| C1 common gate | `0.9985561417` | `1.0056381227` | Essentially identical to static C1; Na-only remains assay-inadmissible |
| C1 K recruitment | `1.0340323883` | `1.0411256187` | Novel and thermodynamically uncertified; K recruitment worsens both metrics |
| C3 common gate from moved root | `1.0384804790` | `1.0461011968` | Direct Na/K plus PKA evidence, but wrong phenotype direction and resting $Cl_i$ gate failure |
| C6 1:1:2 common gate from exact-volume root | `1.1221686891` | `1.1301267548` | Thermodynamic Na/K diagnostic; wrong direction, but this particular root fails the resting $Cl_i$ gate |

For C1 common gating, the endpoint ratio is invariant across the assay folds
`1.25`, `1.6`, and `3` to the shown precision. A separate endpoint flux check
finds $J_4=0.03780358$ at both $F=1$ and $F=3$, while the ungated kernel falls
almost exactly threefold. This is the capacity compensation in (12b)--(12c),
not an implementation that ignored the PKA input. In contrast, K recruitment
rotates the cation signature as in (12e); the induced $N+3P$ response does not
outgrow the direct $J_K$ penalty in (5b). The C3 endpoint is already above one
at $F=1.25$ and becomes more adverse as the gate increases.

The subsequent gated screen closes that loophole numerically within the
finite source-default slice rather than by rejecting the resting states.
Every frozen capacity-quantile point in the rest-admissible 2:1:3 intervals
and every tested gate fold has KO/WT greater than one.

The stimulated evidence brackets must remain distinct: $F=1.25$ is the
recombinant forskolin result, whereas approximately $F=3$ is the native
AE4-specific isoproterenol response. The $F=1.6$ total-WT exchange result is a
sensitivity rather than an AE4-specific fold. $F=1$ is only a no-activation
reference and is excluded from the source-complete stimulated domain.

The minimum frozen $F=3$ rows in each kinetic family are:

| Rest-passing 2:1:3 family | Endpoint KO/WT | Post-step integral KO/WT |
| --- | ---: | ---: |
| C6 mass action | `1.0201802845` | `1.0272280979` |
| C7 affinity | `1.0317095637` | `1.0384748907` |
| C8 saturation | `1.0709456182` | `1.0774352787` |

Across the complete positive AE4-specific source domain
$F\in\{1.25,3\}$, the minimum is the C6 near-lower-boundary $F=1.25$ row:
endpoint `1.0137777584` and integral `1.0209390021`. Thus even the most
favorable activated row predicts increases of `1.378%` and `2.094%`, not a
reduction. The lower no-activation reference (`1.0125760086`) is excluded from
that source-complete stimulated domain. The independent BDF values at the
source-domain minima are `1.0137777580` and `1.0209390004`, preserving both
the sign and the printed precision.

An independent C8 audit brackets its connected rest-pass capacity band as
$G\in[0.0344779857,0.0378955063]$ and finds the $F=3$ endpoint ratio remains within
`1.070946--1.079120` across the frozen interval points. Thus the inverse sign
is not an isolated capacity selection at the narrow resting-data boundary.

The production capacity-by-fold envelope contains 75 Radau rows. All 30
positive AE4-specific rows (five capacity quantiles by three cores by folds
$F=1.25,3$) were rerun independently with BDF: every ratio remains above one,
with maximum BDF/Radau differences $1.36\times10^{-9}$ for the endpoint and
$3.37\times10^{-9}$ for the integral. A separate crosscheck covers 15
representative rows (one in-band capacity per core across all five gate folds).
The two crosscheck sets contain 39 unique rows; no claim of dual-solver
coverage over all 75 is made.

Stronger source-supported gates move the result monotonically farther from the
observed reduction. Equation (5c)
explains the shared sign: 2:1:3 changes the anion ray, but its forward K branch
still displaces K-channel support. Across all three tested source-default rate
laws, state feedback does not increase $N+3P$ enough to overcome $X_K$.

The full knockout residual independently challenges the chassis, not just an
AE4 rate law. Solving the AE4-null resting equations gives a common screened
root

\[
u^*_{KO}=(118.5653,5.55935,65.6607,21.1323,123.7075,48.7822,58.9044).
\]

with maximum raw residual $5.2\times10^{-14}$, pH `7.3785`, and cell volume
`2.974 pL`. It is algebraically closed but biologically poor: pH is `+24.4`
reported SEM above the AE4-null target $6.89\pm0.02$, and chloride is `+7.68`
SEM above the AE4-null target $36.50\pm1.60$ mM. Equivalently, chloride is only
`2.63%` below the measured WT reference rather than the primary knockout's
`27 ± 3.1%` reduction. Relative to the moved C3 and exact-volume C6 WT
roots, knockout chloride is actually higher, reversing the observed
phenotype. An initialization sensitivity that started each genotype at its
own exact screened root preserved the endpoint ratios and their qualitative
signs at the quoted precision: approximately `0.99856`, `1.03403`, `1.03848`,
and `1.12217` for common C1, K-recruited C1, C3, and exact-volume C6 1:1:2.
This does **not** prove identical prestimulus trajectories or exact convergence
at the stimulus time: integrated ratios changed by about
$3.9\times10^{-4}$. The defensible claim is endpoint/sign robustness for
those four diagnostics. The candidate-independent KO acid-base/volume failure
remains a separate chassis falsification signal. It is not used as a substitute
for the rest-passing 2:1:3 tests; those independently produce the adverse
ratios above.

The corresponding **C1-specific** AE2-null root does not rescue the
resting-homeostasis audit. AE4 remains active in that genotype, so unlike the
AE4-null root this state is not candidate-independent. Its chloride is
`49.5249` mM versus the direct AE2-null target
$54.50\pm1.80$ mM (`-2.764` SEM, outside the two-SEM gate), whereas its pH
`6.9141` is within the direct $6.95\pm0.05$ target (`-0.718` SEM). Even the
AE2 study's own-control chloride reference, `53.40` mM, leaves the root
`-2.153` SEM low. Thus the evidence supports *near-normal AE2-null secretion*,
not an unqualified claim that this chassis reproduces AE2-null resting
chloride. The failure of both genotype-specific chloride audits reinforces the
missing-homeostasis diagnosis.

The same frozen audit solved candidate-specific AE2-null resting roots for all
15 rest-admissible 2:1:3 capacity points. Every root fails the direct
$54.50\pm1.80$ mM chloride target: $Cl_i=45.9749$--`46.5775 mM`, or
`-4.736` to `-4.401` SEM. Their pH residuals span `+1.216` to `+2.926` SEM,
so some also fail pH. Yet their stimulated AE2-KO/WT endpoint-flow ratios span
only `0.98894--0.99893`. This is the same distinction: a small acute secretion
effect does not validate chronic resting ion homeostasis.

### Closure theorem and Round-2 boundary

The preceding result is not a failed parameter search. Under the assumptions

1. the non-AE4 vector field and historical WT state are fixed;
2. an AE4 candidate enters only through its declared species source vectors;
3. its kinetic, saturation, occupancy, and gating factors are scalar and
   nonnegative; and
4. no second pathway is silently recalibrated,

the frozen historical WT state is a steady state **if and only if** the
candidate's net AE4 contribution at that state equals
$J_0(+1,-1,0,-2)$ in $(Cl,Na,K,HCO_3)$ order. For one fixed-partition cycle,
this occurs if and only if its source is a positive multiple of the historical
1:1:2, Na-only source and its scalar rate supplies the corresponding $J_0$.
For explicit branches, a nonzero K capacity is compatible only if its *net K
source at that row* is zero or is cancelled by another declared AE4 branch.
Neither escape occurs in the tested passive source-supported branches. This
is a theorem about the frozen calibration row, not a universal claim about a
model allowed to relax to another WT state.

At the numerical amount-balance tolerance $\epsilon=10^{-8}$, the allowed
baseline net transported-K fraction is bounded by

\[
f_{K,net}:=\frac{|J_{K,AE4}(u_0)|}{J_0}
\le \frac{\epsilon}{J_0}
=2.5718\times10^{-7}. \tag{18}
\]

The minimal Round-2 interpolation found ten WT-closing rows, all at the
exactly pure-Na boundary; none was primary-assay admissible. Thus the numerical
grid reaches the analytical closure boundary only by collapsing to the
excluded Na-only limit. For scale, assigning branch capacity in proportion to
the reported assay maxima would give
$f_K=1.6/(1.5+1.6)=0.516$. In a **fixed-partition event law** such as C2,
that is also the event-wise K source fraction and is about $2.0\times10^6$
times the tolerance bound. This comparison does not apply directly to
explicit C3/C6--C8 branch capacities because their Na and K driving
propensities differ. Assigning assay maximum to physiological capacity is
itself a new modelling assumption, and direct K capacity does not by itself
require nonzero net K flux at $u_0$. The tested explicit branch laws did,
however, produce a nonzero K source there and failed closure; the historical
baseline has effectively no room for baseline net K transport.

Sweeps of an equilibrium constant away from the passive convention $K_{eq}=1$
are energetic-bias diagnostics, not source-supported fits. They cannot be
promoted as AE4 reconstructions merely because a pure-Na boundary row closes.
Even the most generous such closure-only endpoint diagnostic, a saturated
1:1:2 pure-Na law with $K_{eq}=100$, gives only about a `0.1476%` decrement
and gives no integrated secretion decrement. This is an envelope of an
unsupported relaxation, not an attainable maximum of the empty admissible
**static fixed-row** set.

This gives the requested necessary/sufficient distinction:

- K participation is necessary, and over the threshold in Section 4
  sufficient, to repair the **thermodynamic direction** of passive AE4 at the
  primary-paper physiological concentrations.
- A source collinear with the Na-only 1:1:2 historical residual is necessary
  and sufficient for **frozen-row WT closure** under the assumptions above.
- Those properties are mutually incompatible in the tested static AE4-only
  family.
  Cooperativity, saturation, gating, and the two alternative stoichiometries
  cannot reconcile them because none creates a new cation-balance pathway.

### First missing pathway

Write source corrections in $(Na,K,Cl,HCO_3)$ order. Replacing the historical
Na-only event by the published intracellular partition requires the additional
fixed-row correction

\[
d=(-\beta J_0,+\beta J_0,0,0),
\qquad \beta=120/145. \tag{19}
\]

Among individual fixed-chassis pathway signatures, the Na/K-ATPase cycle
$(-3,+2,0,0)$ has the greatest alignment with $d$:

\[
|\cos\theta|=\frac{5}{\sqrt{26}}=0.98058,
\qquad
\frac{\|d-d_{proj}\|}{\|d\|}=\frac{1}{\sqrt{26}}=0.19612. \tag{20}
\]

A pump change alone therefore cannot exactly repair the two cation balances.
The smallest exact positive two-direction decomposition using existing chassis
signatures is

\[
d=\frac{\beta J_0}{2}(-3,+2,0,0)
 +\frac{\beta J_0}{2}(+1,0,0,0), \tag{21}
\]

namely equal code-native cycle-flux increments of Na/K pump and NHE1 (or an
equivalent Na-loading direction). This is a **flux-space localization**, not a
validated parameter retuning and not evidence that either protein changes in
the knockout. It says that the first missing model property needed to admit a
mixed-cation AE4 at the frozen WT row is Na/K-pump-led cation-homeostasis
closure, with an accompanying Na-loading route if that row is to be preserved
exactly.

The primary knockout study reports unchanged NHE-dependent stimulated
alkalinization and unchanged resting pH. Therefore the NHE1 vector in (21)
must not be reinterpreted as evidence for a genotype-dependent NHE1 retuning;
it is the algebraic Na-loading direction needed to preserve this particular
fitted row. The same study's slow alkalinization is mechanistically compatible
with bicarbonate supply to AE4, but its normalized fluorescence rate does not
provide the buffer-capacity conversion required for a quantitative flux
retuning on this time-uncertified chassis.

This resolves the minimal Round-4 one-component test analytically. Scaling
NHE1 alone cannot cancel any K residual because its K component is zero.
Scaling the pump alone cannot match both cations because its fixed Na:K ratio
is $-3:2$, not $-1:1$; its optimal projection necessarily leaves the
`19.612%` residual in (20). No nonlinear parameterization whose fixed-row
effect remains a scalar multiple of either signature can evade those facts.
The two-source equality (21) is the first exact repair, but the missing flux
conversion and unchanged NHE readout make a quantitative joint retuning
unsupported. Thus Round 4 is a one-component rejection plus flux-space
localization, not merely an omitted simulation.

Equation (5) separately explains how a cation-homeostasis error can propagate
to secretion: steady apical chloride current equals pump plus K-efflux
support. The 35% phenotype cannot occur on the fixed water branch unless that
support changes enough to cross the chloride-flux threshold in (16)--(17).
Neither (20) nor (21) proves that a pump retuning is sufficient to create the
phenotype; a joint correction still requires a knockout-blind test. A
water-gain change could alter Eq. (14), but it cannot repair the explicit Na/K
source mismatch in (19).

The exact frozen-row result remains: cation balance blocks every tested
source-supported mixed-cation law at the historical state, the C1 lineage
decrement is `0.144386%`, and the unsupported static closure-only envelope is
about `0.1476%` with no integrated decrement. If a mixed-cation candidate were
required to preserve that row, the first additional balance direction would
be Na/K-pump-led cation homeostasis, subject to the limitations above.

KO-blind continuation did produce rest-screen-passing C6, C7, and saturated
C8 2:1:3 roots, so state relaxation closed the resting-state loophole in the
tested slice. Their held-out dynamics then closed it in the opposite way:
every rest-admissible capacity selection and stimulated evidence-band gate
($F\in\{1.25,3\}$) produced an inverse phenotype, not a secretion decrement.
Consequently the structural Task-12 classification for the fixed chassis and
finite source-default domain is

> **AE4 ALONE INSUFFICIENT; MISSING PATHWAY LOCALIZED.**

Over the jointly rest-admissible, assay-compatible, thermodynamic gated rows
in that finite domain, the maximum predicted *positive* secretion reduction
is `0%`; even the closest rows have KO/WT greater than one. This is not a
universal exclusion of AE4 biology or of untested $K_{eq}$, saturation,
Na:K-weight, and bath combinations. The `0.144386%` C1 decrement and `0.1476%`
static diagnostic envelope are not counterexamples because they exclude the
direct K evidence and, for the latter, use unsupported energetic bias. The
required additional model property is stimulus-dependent cation-current
support sufficient to raise $(P+K_{out})_{WT}-(P+K_{out})_{KO}$ to the exact
threshold in (16)--(17). At the frozen historical row the closest missing
direction is Na/K-pump-led cation homeostasis; the full dynamics also implicate
K-channel and acid-base/volume support through (5c) and the KO residual. This
is a localization, not proof that retuning any one of those components is
sufficient.

## Reproducibility

The exact identities and primary-paper threshold calculations are implemented
independently in `src/ae4_mechanism_reconstruction/structure.py` and checked in
`tests/test_ae4_structure.py`. They do not read or execute files under
`archive/`.
