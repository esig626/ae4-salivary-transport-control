# AE4 candidate mechanism family

## Scope and evidence boundary

This report derives the AE4 candidate family used in Task 12. It starts from
the primary transport experiments and applies conservation, electroneutrality,
dimensions, and thermodynamics before using the 2018 model or the historical
implementation. The held-out knockout secretion ratio is not used anywhere in
the derivation or in candidate parameterization.

The evidence labels are:

- **Primary experimental evidence:** Peña-Münzenmayer et al. (2016),
  [DOI 10.1085/jgp.201611571](https://doi.org/10.1085/jgp.201611571).
- **Published:** Vera-Sigüenza et al. (2018),
  [DOI 10.1007/s11538-017-0370-6](https://doi.org/10.1007/s11538-017-0370-6).
- **Historical implementation evidence:** the immutable Na-only branch audited
  in Task 11.
- **New modelling decision:** a mathematical closure not uniquely selected by
  those sources.

The implementation is
`src/ae4_mechanism_reconstruction/mechanisms.py`. The registry in that module
is machine-readable and the corresponding exact checks are in
`tests/test_ae4_mechanisms.py`.

The source-constrained stimulation extension also uses Peña-Münzenmayer et
al. (2021), [DOI 10.1152/ajpgi.00145.2021](https://doi.org/10.1152/ajpgi.00145.2021).

## Constraints from the primary study

The following are measurements, not conclusions imported from the 2018
model.

| Constraint | Direct observation | What it permits or excludes |
| --- | --- | --- |
| Coupled substrates | Removing external Cl produced HCO3-linked alkalinization and Cl loss; removal of external Na or HCO3 suppressed it (Figs. 1-2) | AE4 is cation-dependent and couples Cl to bicarbonate chemistry |
| Na transport | AE4-transfected cells accumulated Na when the exchange direction was imposed, and the change reversed when external Cl was restored (Fig. 3) | Na is transported; a pure Na-gated anion exchanger is insufficient |
| Electroneutrality | Exchange changed membrane voltage by less than 2 mV; rates at -100 and 0 mV were essentially identical; macroscopic current did not differ from controls (Fig. 4) | A physiologically used cycle must have zero net transported charge |
| K transport | Under Na-free imposed gradients, AE4-transfected cells showed K-dependent exchange and sustained K efflux (Fig. 8) | K is transported, not merely a substitute allosteric ligand |
| Cation selectivity | Na, K, Li, Rb, and Cs all supported robust exchange (Fig. 8) | A strictly Na-selective physiological mechanism is excluded |
| Na/K dose response | Na: EC50 49 mM, Hill 2.0, Rmin 0.3, Rmax 1.5; K: EC50 62 mM, Hill 1.8, Rmin 0.4, Rmax 1.6 (Fig. 9; response rates in the paper's assay convention) | Na and K have similar apparent affinity and capacity; cooperative response is measured |
| Meaning of cooperativity | The authors explicitly leave two interpretations open: more than one transported cation per cycle, or transported cation plus allosteric regulation | A Hill exponent must not be silently identified with stoichiometry |
| Candidate stoichiometries | The thermodynamic analysis explicitly examined Cl:cation:HCO3 ratios 1:1:2, 1:2:3, and 2:1:3 (Fig. 7) | All three integer, electroneutral ratios must be retained until tested |
| Reversal | Restoring external Cl reversed the Na response; the study also demonstrated exchange-mode substrate substitution under low-Cl conditions | An irreversible one-way uptake law is not source adequate |
| Physiological cation | Na-only thermodynamic calculations required nonphysiological concentrations; K-supported 1:1:2 transport gave a physiological inward-Cl region (Figs. 7 and 10) | Historical Na-only C1 is a negative control, not a primary-supported physiological mechanism |

At the primary study's pooled-cation 1:1:2 calculation, inward Cl was favored
for intracellular Cl below about 73 mM or intracellular HCO3 above about
16 mM when the other reported concentrations were held fixed. These are
thermodynamic calculations based on an assumed stoichiometry, not direct
measurements of a reversal point.

Two source-level possibilities cannot be represented honestly on the current
seven-state chassis:

1. The authors state that Na-dependent Cl/CO3 exchange cannot be ruled out.
   They do not identify its stoichiometry, and the chassis has no carbonate
   state. `C6_CARBONATE_UNTESTABLE` records this rather than substituting HCO3
   silently.
2. Under nonphysiological Cl-free conditions, AE4 supports a
   bicarbonate/Na-bicarbonate exchange mode with unknown bicarbonate
   coefficients `x:y` (Figs. 5-6). This establishes substrate substitution and
   reversibility but not a defined physiological Cl-cycle stoichiometry.
   `C6_HCO3_SELF_EXCHANGE_UNTESTABLE` records the unresolved mode.

## Beta/cAMP/PKA activation evidence and boundary

The 2021 primary study establishes a stimulus dependence that is absent from
all static C1-C8 cores:

- In isolated submandibular acinar cells lacking Ae2, and hence isolating
  native Ae4 exchanger activity, 5 micromolar isoproterenol increased the
  initial chloride-uptake rate by approximately threefold in Fig. 1C.
  Pretreatment with 10 micromolar H89 prevented that increase.
- In wild-type acinar cells, total Ae2-plus-Ae4 exchanger activity increased
  by approximately 1.6-fold in Fig. 1B; this value is not AE4-specific.
- In Ae4-expressing CHO-K1 cells, 10 micromolar forskolin increased the
  alkalinization rate by approximately 25%, an effect prevented by H89.
  Constitutively active PKAc also increased both HCO3-linked alkalinization
  and directly measured Cl flux.
- The S173A Ae4 mutant was not activated by PKAc, whereas S273A retained PKA
  dependence. This supports regulation associated with the cytosolic S173
  site, although the study did not directly demonstrate its phosphorylation.
- The authors explicitly state that their macroscopic assay cannot
  distinguish higher turnover, a larger number of active plasma-membrane
  transporters, or a PKA-induced stoichiometric change.
- The same discussion reports that muscarinic stimulation does not activate
  Ae4. Therefore the beta/cAMP/PKA input is independent of the model's calcium
  input; calcium is not used as a proxy for PKA activation.

The full-activation fold values fixed before held-out simulation are consequently
`(1.25, 1.6, 3.0)`: recombinant forskolin, approximate wild-type total
exchange, and approximate native Ae4-only isoproterenol, respectively. The
last two are figure-level readings rather than published fitted constants.
They define an assay-observed sensitivity envelope, not parameters estimated
from secretion or knockout output. Mapping one of these initial-rate assay
folds to a sustained whole-cell capacity multiplier is a new modelling
assumption.

## Common sign, balances, and thermodynamics

Let a positive event rate `J` denote the physiological direction: `a` Cl ions
move from bath to cell while `b` monovalent cations and `c` bicarbonate ions
move from cell to bath:

\[
a\,Cl_e+b\,C_i+c\,B_i
\rightleftharpoons
a\,Cl_i+b\,C_e+c\,B_e,
\qquad B=HCO_3^-.
\]

The exact intracellular source vector in `(Na,K,Cl,HCO3)` order is

\[
sJ=(-b_{Na}J,-b_KJ,+aJ,-cJ),
\qquad b_{Na}+b_K=b.
\]

The transported charge entering the cell per event is

\[
q/F=-a-b+c.
\]

The direct electroneutrality evidence therefore requires

\[
\boxed{c=a+b}.
\]

This identity is satisfied by each explicitly considered ratio:
`1:1:2`, `1:2:3`, and `2:1:3`. It also makes the voltage contribution to the
electrochemical affinity cancel exactly.

For a branch carrying cation species `s`, define the dimensionless forward
affinity

\[
\mathcal A_s=
\log\!\left[
K_{eq}\,
\frac{Cl_e^a C_{s,i}^b B_i^c}
     {Cl_i^a C_{s,e}^b B_e^c}
\right].
\]

Concentrations are divided by the same 1 mM reference before products are
formed in code. Thus `A_s` is dimensionless. Positive affinity favors Cl
entry, the forward free-energy change is `-RT A_s`, and equilibrium is

\[
\boxed{\mathcal A_s=0}.
\]

For a passive transporter between chemically identical aqueous compartments,
`K_eq=1` is the primary thermodynamic convention. A different value must be
independently justified; it is not a knob to force the knockout phenotype.

Two independently implemented sign-controlled laws are used:

\[
J_s=V_s\left(K_{eq}^{1/2}P_s^+-K_{eq}^{-1/2}P_s^-\right)
\tag{mass action}
\]

and

\[
J_s=2V_s\sinh(\mathcal A_s/2).
\tag{affinity law}
\]

Here `P_s+` and `P_s-` are the normalized forward and reverse substrate
products. A symmetric saturating alternative is

\[
J_s=V_s\frac{K_{eq}^{1/2}P_s^+-K_{eq}^{-1/2}P_s^-}
{1+\sigma(K_{eq}^{1/2}P_s^++K_{eq}^{-1/2}P_s^-)}.
\]

Its denominator is strictly positive, so saturation can reduce magnitude but
cannot change the thermodynamic direction or equilibrium.

## Candidate table

`AE4-specific free parameters` counts quantities not fixed by the stated
independent evidence in the primary comparison. Every active candidate has one
capacity normalized from WT closure. The held-out KO result is never a target.

| Family / implementation IDs | Transported species and balance insertion | Selectivity / cooperativity | Reversal and saturation | AE4-specific free parameters in primary comparison | Exact frozen-source match with one capacity? | Evidence and new assumptions |
| --- | --- | --- | --- | ---: | --- | --- |
| **C1** | Historical `1 Cl in : 1 Na + 2 HCO3 out`; no K term | Strict Na-only | Historical forward-minus-reverse polynomial; not detailed-balance certified | 1 capacity; historical default is also retained as a lineage check | **Yes.** It is the source of the required historical vector. | **Historical implementation evidence:** exact Na-only law, fixed 21 mM reverse bicarbonate, folded reverse coefficient. Contradicted as a physiological class by direct K transport and primary thermodynamics |
| **C2** | Printed pooled `(Na+K)` 1:1:2 rate; the net cation term is divided by intracellular Na/K mole fractions exactly as printed | Nonselective pooled cation; no measured Hill dependence | Empirical polynomial with printed `k+`, `k-`; not thermodynamically certified | 1 capacity; printed rates fixed | **No.** Its nonzero K partition conflicts with the required zero K source. | **Published:** Eq. 31 and cell balances. **New assumption:** pooled rate is meaningful despite the Markov reduction problem. Exact printed partition is not a reversible species-resolved cycle |
| **C3** | Separate Na-coupled and K-coupled 1:1:2 branches; each branch enters its own cation balance | Relative Na/K capacity explicit | Branchwise empirical mass action using the C2 constants | **1 in Round 1:** capacity is fitted and Na:K weights are fixed at the independent Rmax ratio `1.5:1.6`. A fitted branch ratio would raise the count to 2. | **No under the source-supported mixed-cation convention.** It becomes shape-compatible only by disabling the K branch. | **Primary evidence:** both ions transported, similar EC50/Rmax. **New decisions:** independent parallel cycles with common kinetic constants, and Rmax ratio treated as a capacity ratio |
| **C4** | One transported cation in parallel 1:1:2 branches | Measured Hill occupancy replaces linear donor-cation dependence: Na `(49,2.0)`, K `(62,1.8)` | Empirical bidirectional occupancy difference; its zero need not equal chemical-affinity zero | 1 capacity; Hill values fixed | **No.** Active K transport creates a K-balance residual. | **Primary evidence:** Hill fits. **New decision:** cooperativity belongs to occupancy while only one cation is transported. Fails the strict C7 thermodynamic criterion unless a microscopic scheme is supplied |
| **C5** | Parallel transported-cation 1:1:2 core plus an additional positive gate | Direct cation transport retained; Na/K Hill gates and Rmin/Rmax floors fixed from assay | Thermodynamic mass-action sign multiplied by donor-side gate; gate cannot move equilibrium | 1 capacity | **No.** Active K transport creates a K-balance residual. | **Primary evidence:** direct cation flux plus added allostery explicitly left viable. **New decision:** phenomenological donor-side gate; no microscopic allosteric states |
| **C6_112**, **C6_123**, **C6_213** | Branch-resolved ratios `1:1:2`, `1:2:3`, `2:1:3`; exact `(Na,K,Cl,HCO3)` signatures are applied | Equal Na/K branch capacity is the primary minimal convention; primary Rmax values are similar | Reversible dimensionless mass action | 1 capacity per discrete stoichiometry | **112:** no with active K (only if K is disabled). **123/213:** no for any cation partition because their anion ratios differ. | **Primary evidence:** all three were explicitly considered, not measured. **New decision:** pure-Na and pure-K branches in a mixed-cation environment |
| **C7_112**, **C7_123**, **C7_213** | Same three C6 balance signatures | Same branch convention | `2 sinh(A/2)`; exact sign and equilibrium controlled by each branch affinity | 1 capacity | Same structural verdict as C6. | **Conservation/thermodynamics:** source-justified stoichiometries with a new symmetric rate law |
| **C8_*_MA**, **C8_*_SAT** | Each viable stoichiometry is held fixed within a mass-action/saturating pair | Same branch convention | Unbounded mass action versus symmetric substrate saturation; both have identical equilibrium | 1 for MA; 2 for SAT (capacity plus saturation strength) | Same structural verdict as C6; saturation cannot repair a source-vector mismatch. | **Alias warning:** each `C8_*_MA` is exactly the corresponding `C6_*`, not independent evidence. **Primary evidence:** saturating cation response. **New decision for SAT:** common full-substrate occupancy denominator; `sigma=1` is an unconstrained Round-1 default and the cation assay does not identify it |

The three `C8_*_MA` IDs exist only to make each requested saturation-versus-
mass-action pair explicit. Algebraically and in the implementation,
`C8_112_MA == C6_112`, `C8_123_MA == C6_123`, and
`C8_213_MA == C6_213` for every shared context and parameter set. They are
aliases of the C6 hypotheses and must be collapsed to their C6 canonical IDs
for ranking, complexity accounting, or claims about independent mechanistic
support. Only the `C8_*_SAT` members add a distinct rate law, together with the
extra unconstrained saturation parameter.

## Stimulus-dependent extension fixed before held-out validation

Let `u_PKA` be a dimensionless beta/cAMP/PKA input on `[0,1]`, separate from
calcium, and let `F_PKA` be one of the assay folds fixed without access to
knockout outputs. The positive
macroscopic gate is

\[
g_{PKA}(u)=1+u_{PKA}(F_{PKA}-1).
\]

All stimulus candidates default to `u_PKA=0`; therefore their resting source
is bit-for-bit the corresponding static core and the original WT calibration
is preserved. A full beta/cAMP stimulus uses `u_PKA=1`. An instantaneous step
aligned with the model's combined stimulation epoch is a protocol assumption:
the source study measured initial transport rates under treated conditions and
does not identify an activation time constant. A future measured cAMP/PKA time
course can replace this input without changing any AE4 balance law.

| Candidate | Resting source | Full-stimulus source | Round-1 parameter count and exact baseline status | Evidence versus new assumptions |
| --- | --- | --- | --- | --- |
| **C1_PKA_COMMON** | Exact C1 source | `g_PKA` times every C1 source term | 1 fitted resting capacity; `F_PKA` is independently selected from the assay envelope. Exact resting closure. | **Round-2 diagnostic; assay-inadmissible as a complete mechanism. Source:** PKA increases macroscopic AE4 activity. **New:** interpret that increase as common turnover gating and map assay activation to the model step. It remains Na-only and conflicts with direct K transport. |
| **C1_PKA_K_RECRUIT** | Exact C1 source; no K term | With basal C1 flux `J0`, recruited `J_K=u_PKA(F_PKA-1)J0`; source is `(-J0,-J_K,J0+J_K,-2(J0+J_K))` in `(Na,K,Cl,HCO3)` order | 1 fitted resting capacity; no additional continuous fit because all recruited cycles are discretely assigned to K. Exact resting closure. | **Explicitly novel Round-2 diagnostic; never a Round-1 source-justified candidate.** PKA activation and K transport are separately supported, but selective K recruitment is not measured. Resting Na-only transport and K direction-locking to C1 are also new assumptions. This is electroneutral but not a thermodynamically validated K law. |
| **C3_PKA_COMMON** | Exact C3 parallel Na/K source | `g_PKA` times both branch sources | 1 fitted capacity with Na:K weights fixed `1.5:1.6`; its resting fixed-vector residual is exactly C3's and is not repaired by the gate. | **Direct-source combination:** experiments establish Na/K transport and PKA-dependent macroscopic activation. **New:** a common factor acts on both branch sources; the operational law deliberately does not claim whether turnover or active transporter number changed. |
| **C6_112_PKA_COMMON**, **C6_123_PKA_COMMON**, **C6_213_PKA_COMMON** | Exact corresponding reversible C6 source for `1:1:2`, `1:2:3`, or `2:1:3` | `g_PKA` times both thermodynamic branch sources; each branch affinity and `K_eq=1` are unchanged | 1 fitted resting capacity per discrete stoichiometry, Na:K weights fixed `1.5:1.6`, and `F_PKA` independently selected. Resting fixed-vector status is exactly the core's: 112 cannot match with active K; 123/213 cannot match the anion ratios. | **Simplest assay-plus-thermodynamics source combination (112) and its two source-viable stoichiometric alternatives:** direct Na/K transport, reversibility, the three considered ratios, and PKA activation are experimental constraints. **New:** pure-cation branches, a common positive macroscopic gate, `K_eq=1`, and mapping an assay fold to the model stimulus. The primary study considered but did not measure any of the three stoichiometries. |

The Round-3 moving-root computation also applies this same common positive
wrapper to the source-defined `C7_213` affinity core and `C8_213_SAT`
saturating core. The result rows are labelled
`C7_213_PKA_COMMON` and `C8_213_SAT_PKA_COMMON` for traceability, but these are
**sensitivity adapters, not additional registry mechanisms or additional
experimental evidence**. Each inherits its core's one fitted resting capacity,
balance signature, affinity/reversal condition, and kinetic assumptions; the
only added discrete input is a separately reported assay fold. The wrapper
adds no phenotype-fitted parameter and must scale total and branch diagnostic
fluxes consistently.

For a common-gate candidate, `g_PKA>0` preserves electroneutrality, source
direction, every branch-affinity sign, and any equilibrium already supplied by
its core. It cannot repair a resting source-vector mismatch. The three C6 PKA
members share the same wrapper but are not aliases: their distinct
stoichiometries produce distinct balance rays and equilibria. Registering all
three avoids treating PKA activation as evidence for the otherwise unmeasured
1:1:2 ratio; `C6_112_PKA_COMMON` is the minimal representative, not a uniquely
supported stoichiometry. For the recruited-K diagnostic,
electroneutrality follows exactly from each 1:1:2 branch, but no claim of
detailed balance is made because `J_K` is not controlled by a separate K
chemical affinity.

The genotype multiplier remains outside the stimulus law. Thus WT activity is
`m=1`, AE4 knockout is `m=0`, and every source term is exactly zero in knockout
for every value of `u_PKA` and `F_PKA`. Neither the `0.65` held-out secretion
ratio nor any other knockout output enters resting calibration, the fold
envelope, or the stimulus protocol.

### C1 exact historical law

The negative control is reproduced as

\[
J_1=G_1\left[
k_+Cl_eB_i^2Na_i
-(k_-\,118.7219)Cl_i(21)^2Na_e
\right].
\]

At the historical WT row, the bracket is
`0.00205186299263`; multiplying by the historical capacity
`18.9502607221547` gives `0.03888333867648` in chassis flux units. These
values are a regression of historical evidence, not independent biological
calibration.

### C2 exact printed law and its balance limitation

The printed 2018 constitutive law is

\[
J_2=G_2\left[
k_+Cl_eB_i^2(Na_i+K_i)
-k_-Cl_iB_e^2(Na_e+K_e)
\right].

The published intracellular cation terms are

\[
-\frac{Na_i}{Na_i+K_i}J_2,
\qquad
-\frac{K_i}{Na_i+K_i}J_2.

The implementation deliberately preserves those expressions so C2 is an
exact test of the printed mechanism. They are not equivalent to reversible
Na and K branches when `J_2<0`: a reverse cycle draws cation from the bath, so
its Na/K identity should depend on external availability or explicit branch
kinetics, not the intracellular fractions. C2 can therefore be judged in its
forward operating regime, but it cannot pass a general reversal audit merely
because its scalar polynomial changes sign.

### C3 branch resolution

For cation `s` in `{Na,K}`,

\[
J_{3,s}=G_3w_s
\left[k_+Cl_eB_i^2C_{s,i}-k_-Cl_iB_e^2C_{s,e}\right],
\qquad w_{Na}+w_K=1.

The Na branch contributes only to the Na balance and the K branch only to the
K balance. This removes C2's ambiguous pooled balance insertion. A freely
fitted `w_K/w_Na` costs one additional parameter; fixing it from the similar
primary `Rmax` values avoids that cost.

### C4 occupancy versus C5 allostery

The measured Hill function is

\[
h_s(C)=\frac{C^{n_s}}{EC50_s^{n_s}+C^{n_s}}.

C4 places `h_s` in the directional transported-cation occupancy. That is a
literal empirical test of the observed cooperative response, but it removes
the cation chemical activity from the rate product. Its flux zero can
therefore differ from `A_s=0`. C4 is retained as the requested empirical
occupancy hypothesis and explicitly marked thermodynamically uncertified.

C5 instead retains the linear transported-cation activity in the reversible
1:1:2 core and multiplies it by a positive gate

\[
g_s(C)=f_s+(1-f_s)h_s(C),
\]

where the source-fit ratios give `f_Na=0.3/1.5=0.20` and
`f_K=0.4/1.6=0.25`. The donor-side concentration is used according to the
direction of a branch. Because `g_s>0`,

\[
\operatorname{sign}(g_sJ_s)=\operatorname{sign}(J_s).

\]

Thus extra allostery can alter rate magnitude but cannot rescue a cycle whose
thermodynamic direction is wrong. The gate is a new phenomenological decision,
not a recovered microscopic scheme.

## Fixed-row compatibility and moved-root scope

The fixed historical chassis was calibrated around a Na-only 1:1:2 AE4
source. At its WT row, the AE4 contribution required to leave every other
equation unchanged is

\[
r_* = J_*(-1,0,+1,-2),
\qquad J_*=0.03888333867648,

\]

in `(Na,K,Cl,HCO3)` order.

For a candidate with one nonnegative scalar capacity, stoichiometry `(a,b,c)`,
and forward cation fractions `(p_Na,p_K)`, exact preservation of that same
row requires

\[
(-bp_{Na},-bp_K,a,-c)\lambda=(-1,0,1,-2)J_*.

\]

For nonzero flux this implies all of the following:

\[
p_K=0,\qquad a=b,\qquad c=2a.

\]

Therefore only a Na-only multiple of `1:1:2` can reproduce the exact
historical source vector with one capacity. Any candidate carrying direct K,
or either alternative primary stoichiometry, must leave a nonzero WT closure
residual at that frozen row unless the state moves or another balance pathway
changes. This is an exact **chassis-compatibility result**, not a biochemical
rejection of the primary evidence. The common calibration code reports that
residual explicitly and does not retune a non-AE4 parameter to erase it.

A separately identified state-relaxation branch confirms the distinction. C3
has a moved resting root in chassis order `(Na_l,K_l,height,Na_i,K_i,Cl_i,HCO3_i)`
of approximately `(118.76,5.618,27.55,30.36,114.71,46.31,10.21)`, with
`pH_i=6.908` and maximum raw residual `8.5e-14`. This is a valid C3 state root,
not a pass of the frozen historical-row calibration. `C3_PKA_COMMON` should be
tested from both branches with the distinction retained: its positive PKA gate
does not alter either resting root at `u_PKA=0`.

## Analytical exclusions and non-exclusions

### Rejected before simulation

1. **Electrogenic ratios with `c != a+b`.** They contradict the voltage-clamp
   and current-clamp evidence and would require a membrane-potential term in
   the affinity. None is admitted as an electroneutral AE4 candidate.
2. **Pure cation gating with no cation transport.** A `1 Cl:1 HCO3` anion
   exchange multiplied by a Na/K gate is electroneutral, but direct Na and K
   flux measurements exclude it as the full AE4 mechanism. It is registered
   as `C5_GATE_ONLY_REJECTED`; additional gating on a transported-cation core
   remains viable as C5.
3. **A fractional transported stoichiometry equal to Hill 1.8 or 2.0.** A Hill
   coefficient is a response slope, not a molecule count. The primary authors
   explicitly state the alternative allosteric interpretation. C4, C5, and
   integer C6 candidates test the possibilities separately.
4. **Irreversible uptake-only laws.** Direct reversal and exchange-mode data
   require a reverse path. A one-way flux cannot satisfy the observed behavior
   or a finite equilibrium condition.
5. **Using a positive gate or saturation denominator to reverse affinity.** A
   positive multiplier cannot change a sign. Such a result would be a coding
   or sign error, not a mechanism.
6. **Silently replacing carbonate by bicarbonate.** The carbonate alternative
   is source-viable but not identified. The current state system cannot test it
   without a carbonate/acid-base extension.

### Retained but qualified

- C1 is retained as the historical negative control even though primary K
  transport and the Na-only thermodynamic calculation oppose it as the
  physiological mechanism.
- C2 is retained exactly because it is the published claim being tested. Its
  pooled balance insertion and printed Markov-to-polynomial derivation remain
  unresolved.
- C4 is retained as a data-matched occupancy interpolation but must be marked
  thermodynamically invalid if its flux sign disagrees with chemical affinity.
- The three C6 stoichiometries remain discrete hypotheses. The primary paper
  evaluated them thermodynamically; it did not measure which ratio operates.
- C5, C7, and C8 are new mathematical closures constrained to preserve direct
  cation transport, electroneutrality, and the correct reversal condition.

## Parameter discipline and anti-leakage rule

The primary comparison uses the same protocol for every active candidate:

1. Fix every non-AE4 chassis equation and parameter.
2. Fix EC50, Hill slopes, and minimal Na/K capacity assumptions from the
   primary assay where applicable.
3. Fit only the candidate's nonnegative AE4 capacity to WT closure. For C3,
   either fix relative Na/K capacity from the assay or pay one extra parameter.
   For C8-SAT, report the extra saturation parameter and do not tune it to KO.
   For stimulus candidates, keep `u_PKA=0` during this resting fit and select
   `F_PKA` only from the 2021 assay envelope fixed before held-out simulation.
4. Freeze the candidate.
5. Only then evaluate WT, AE4 KO, AE2 KO, thermodynamic direction, and all
   residual phenotypes.

The knockout flow observation does not appear in the mechanism API, the
candidate registry, or the capacity-fit function. It can only be joined to a
frozen prediction table after calibration.

## Implementation contract and checks

`evaluate_ae4(candidate_id, context, parameters, activity_scale=m,
stimulus_activation=u)` accepts a duck-typed context containing the eight
intracellular/bath concentrations. The stored `parameters.activity_scale` is
the calibrated capacity; the supplied `m` is a multiplicative genotype factor
(`1` WT, `0` KO). The optional `u` is the independent beta/cAMP/PKA input; if
it is omitted, the evaluator reads an optional context field
`pka_activation`, otherwise defaults to resting zero. The output contains
explicit intracellular Na, K, Cl, and HCO3 balance terms, not merely a pooled
scalar flux.

The test suite checks:

- exact C1 turnover and historical capacity;
- exact C2 pooled equation and printed cation partition;
- explicit C3 branch insertion;
- all three primary stoichiometric balance signatures;
- electroneutrality for every active candidate;
- zero flux at equilibrium and sign reversal across it for C6-C8;
- identical signs for mass-action and saturating pairs;
- rejection of pure gating without cation transport; and
- preservation of calibrated WT capacity and exact zeroing at AE4 KO;
- bit-for-bit restoration of each static core when `u_PKA=0`;
- common source scaling at full PKA activation; and
- preservation of equilibrium and affinity-controlled reversal by all three
  common-gated C6 stoichiometries; and
- the exact electroneutral source vector of the novel recruited-K diagnostic.

These tests establish algebra and software behavior. Whether a candidate
reproduces the held-out phenotype is a separate result reported by the common
comparison and adversarial audit.
