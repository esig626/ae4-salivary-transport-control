# AE4 transport module: evidence, thermodynamics, and whole-cell interface

## Decision

The smallest defensible whole-cell AE4 module is a reversible, shared-carrier
`SR2-112` core with separate Na and K branches and a common chloride path. It
is electroneutral branch by branch, obeys local detailed balance, exposes the
unmeasured mixed-bath Na/K routing, and returns separate Na, K, Cl, total
inorganic carbon, and alkalinity amount sources.

An `SR5` fast cooperative gate is retained only as a nested extension. The
2016 pure-cation dose curves show an apparent cooperative response, but they do
not distinguish a second transported cation, a nontransported catalytic
cation, heterogeneous transporter populations, or ordinary state saturation.
The 2025 work requires cation-specific mutation logic but does not identify a
complete atomic transport cycle.

The default regulation coupling is a positive common active-capacity gain from
the external beta/cAMP/PKA subsystem. State-specific rate modulation remains a
nested hypothesis. No AE4 carrier or PKA time constant is selected from the
held-out AE4-null secretion curve.

## Scope and evidence labels

This module uses no AE4-null secretion magnitude or timing target. It audits:

- [Peña-Münzenmayer et al. 2016, JGP 147:423-436](https://pmc.ncbi.nlm.nih.gov/articles/PMC4845690/), DOI `10.1085/jgp.201611571`;
- [Catalán et al. 2025, AJP Cell Physiology 328:C2070](https://pmc.ncbi.nlm.nih.gov/articles/PMC12463136/), DOI `10.1152/ajpcell.00346.2024`;
- the frozen Task 13 code, tests, and machine-readable artifacts.

Claims below are labelled as:

| Label | Meaning |
| --- | --- |
| `MEASURED16` | Direct 2016 ion-sensitive dye, pH, current, or voltage result in its stated assay |
| `INTERPRET16` | The authors' mechanistic or thermodynamic interpretation, not a measured stoichiometry |
| `MEASURED25` | Direct 2025 functional mutagenesis, expression control, pH, or voltage-sensitivity result |
| `MD25` | Homology-model or molecular-dynamics result |
| `HYPOTHESIS25` | Proposed cycle, ion species, binding site, or sequence not established by transport measurement |
| `EXACT` | Algebraic consequence of a declared reaction, conservation law, or passive Markov graph |
| `MODEL` | New modeling choice that must be stress-tested |

## Task 13 SR2/SR5 reproduction

The complete frozen Task 13 regression command was rerun on 2026-08-27 UTC:

```bash
PYTHONPATH=src python -m unittest discover -s tests -p 'test_state_resolved_ae4_*.py' -v
```

All **72 tests passed** in 10.6 s. This reproduced carrier closure, branch
charge, local detailed balance, affinity/current signs, reversal, the
T448I-T756A Na-collapse/K-retention limit, frozen evidence hashes, CTMC
snapshots, the failed inherited-chassis WT calibration, and independent
null-field checks.

The frozen independent reproduction contains 17 SR2/SR5/SR6 conditions, 288
state rows, and 365 edge rows. Independent SVD and matrix-tree stationary
solutions reproduce the saved occupancies and currents. Across those snapshots
the largest reported SVD occupancy difference from the saved solution is
`5.71e-12`, the largest matrix-tree occupancy difference is `2.90e-13`, and
the largest independent edge-current difference is `1.32e-12`.

The decisive transporter-level behavior is:

| Family/gauge | States/edges | Na apparent EC50/Hill | K apparent EC50/Hill | Double mutant Na/K | Task 13 disposition |
| --- | ---: | ---: | ---: | ---: | --- |
| SR2-112, slow/unit/fast | `8/9` | `48.64-48.80 / 1.14-1.16` | `53.69-59.17 / 1.17` | Na `0`, K `0.67-0.76` of WT | Shared-pool and mutant logic pass; cooperative dose shape fails |
| SR5-112, unit | `12/15` | `49.000 / 2.040` | `60.816 / 2.052` | Na `0`, K `0.853` | Conditional assay-expressivity pass |
| SR5-112, fast | `12/15` | `49.000 / 2.061` | `61.278 / 2.065` | Na `0`, K `0.837` | Conditional assay-expressivity pass |

These are conditional observation-map results, not identified microscopic
parameters. In particular:

- the 2016 intracellular concentrations during the external cation sweep were
  not measured and Task 13 used a declared sensitivity tuple;
- the reported `Rmin` values were retained as assay backgrounds, not converted
  into cation-independent AE4 current;
- pure-cation data do not measure the shared-pool competition coefficient;
- all transporter-compatible SR5 barrier-grid points on the fixed Task 13
  reference state gave negative net intracellular Cl source. That is a
  mixed-bath routing/gauge result, not evidence against AE4 biology;
- no conditional transporter survivor produced an accepted WT whole-cell
  calibration on the inherited seven-state chassis.

Task 13 therefore licenses the transporter structure and negative-control
regressions, but not its kinetic gauges or chassis as a modern parameter set.

## Audit of the 2016 primary evidence

### Direct findings and assay contexts

| Source location | Direct result | Assay context | Constraint on this module |
| --- | --- | --- | --- |
| Fig. 1 | AE4-dependent alkalinization in native `Ae2-/-` mouse SMG acinar cells required external Na; activity was nearly absent in `Ae2-/-;Ae4-/-` cells | Isolated SMG acini; external Cl `128.3 -> 4 mM`; BCECF; 37 C | AE4 has a Na-dependent anion-exchange mode in native salivary acini (`MEASURED16`) |
| Figs. 2-3 | In AE4-expressing CHO cells, low external Cl produced HCO3 entry, Cl loss, and an AE4-associated intracellular Na increase | Recombinant mouse/human AE4; SBFI/SPQ/BCECF; imaging solutions at 37 C | Na is transported with HCO3 and opposite Cl in this protocol (`MEASURED16`) |
| Fig. 3 caveat | Restoring external Cl returned the Na signal toward baseline, but the authors explicitly note Na/K-ATPase could contribute | Same CHO protocol | Do not treat the recovery trace alone as a particle-stoichiometry measurement |
| Fig. 4 | AE4-associated transport changed membrane voltage by less than 2 mV; alkalinization was essentially the same at `-100` and `0 mV`; AE4 and control currents did not differ significantly | CHO whole-cell recording at about 22 C, simultaneous imaging | Net AE4 cycle is macroscopically electroneutral/weakly voltage-dependent over the tested range (`MEASURED16`) |
| Figs. 5-6, Table 4 | Under Cl-free/low-Cl conditions, HCO3 influx required intracellular Cl or HCO3 and was voltage-insensitive | Artificial whole-cell dialysis; 10 uM EIPA in the imaging protocol | AE4 can support HCO3/HCO3 exchange under nonphysiological low-Cl conditions; simple Na-HCO3 cotransport is not the only explanation (`MEASURED16`) |
| Fig. 8 | K as the main cation supported exchange; in Na-free cells an imposed inward Cl/outward K gradient caused AE4-associated K loss | CHO cells; PBFI has limited K/Na selectivity, so Na was removed | K is transported, not merely a Na substitute in an empirical gate (`MEASURED16`) |
| Fig. 8 | Li, Rb, and Cs also supported alkalinization | Pure-cation recombinant assay | The site is broadly monovalent-cation permissive in that assay; this does not imply equal native branch rates |
| Fig. 9 | Na: `EC50=49 mM`, `nH=2.0`, `Rmin=0.3`, `Rmax=1.5`; K: `62 mM`, `1.8`, `0.4`, `1.6` | External cation values `5,25,50,100,125,150 mM`; external Cl `4 mM`, HCO3 `25 mM`; at least four experiments/condition | Preserve apparent cation cooperativity, but keep it in the assay observation context (`MEASURED16`) |

The paper reports no fit covariance or uncertainty for the Hill summaries, and
the intracellular concentrations during the dose series were not measured.
`nH=2.0` or `1.8` is therefore not a transported-cation count.

### Reversal and stoichiometry discipline

The low-Cl and restored-Cl traces show reversible exchange behavior, and the
Na and K experiments establish their directions in those protocols. They do
not measure a complete simultaneous Na:K:Cl:HCO3 particle ratio.

The paper discusses electroneutral candidate counts `Cl:cation:HCO3` of
`1:1:2`, `1:2:3`, and `2:1:3` and calculates their free energies. Those are
`INTERPRET16`, not measured stoichiometries. It also explicitly leaves
Cl/CO3 exchange and OH substitution unresolved. The modern implementation
must therefore treat `1:1:2` as the smallest discrete working model, not the
identified molecular count.

## Audit of the 2025 primary evidence

### Direct functional and control results

The 2025 experiments used transiently transfected human AE4 in HEK-293 cells,
BCECF alkalinization after external Cl was changed from `128.3` to `4 mM`, and
room-temperature recordings. The output is a normalized pH-fluorescence rate,
not a simultaneous ion-flux ratio.

| Source location | Direct result (`MEASURED25`) | What it does and does not constrain |
| --- | --- | --- |
| Fig. 3 | WT human AE4 produced reversible, Na-dependent alkalinization; nontransfected cells did not | Confirms cation-dependent exchange in this expression system; does not identify stoichiometry |
| Fig. 4 | T448I and T448G reduced the Na-supported initial rate by about 50% | The T448 hydroxyl-bearing position contributes to function; this does not identify a binding edge |
| Fig. 4 | D709A, S446A, and T756A reduced activity by about 30%; D709N, S447A, T754A, I758G, and I758R did not significantly reduce it | Requires residue-specific effects, but not a unique coordination sequence |
| Fig. 5 | S446A-T448I and D709A-T448I resembled T448I; T756A-T448I was comparable with nontransfected cells in Na bath | The combined T448/T756 perturbation collapses detectable Na-supported exchange |
| Fig. 6 | WT, S446A, T448I, and T448G retained similar Na- and K-supported behavior; D709A/T756A were cation-sensitive; T756A-T448I regained substantial K-supported transport while Na-supported transport remained absent | A single pooled `[Na]+[K]` branch with one common mutation multiplier is rejected |
| Figs. S2-S3 | Mutants retained membrane colocalization and total protein expression comparable with WT | The main functional reductions are not explained by gross expression/mistrafficking in this assay |
| Fig. S4 | Resting calculated HCO3 did not increase in mutants, and valinomycin depolarization did not reveal voltage sensitivity | Supports an electroneutral macroscopic cycle; does not select bicarbonate count or carbonate |

The robust mutation rule for modeling is therefore cation-resolved
conductance: the T448I-T756A Na path can approach assay background while the K
path remains nonzero. Changing both directions of a branch by the same factor
implements this hierarchy without inventing a mutation-dependent equilibrium.
The approximate `0.5` and `0.7` Task 13 branch factors are qualitative assay
encodings, not microscopic rate estimates.

### MD and proposed cycles remain hypotheses

The structural work used an NDCBE-based homology model, manually placed one
Na/HCO3 or K/HCO3 pair at a putative site, and selected replicas from 120-ns
MD simulations. It was not an AE4 structure and did not simulate Cl in the
proposed transport step.

- D709 and T713 coordination of Na/K, HCO3 interactions with G449 and K879,
  and the loss/preservation of T713 interactions in the double mutant are
  `MD25`.
- G449, P450, T713, and K879 were not functionally mutated in this paper.
- MD resolved one modeled HCO3 at the site. A second HCO3 site, one transported
  carbonate, D709 protonation control, and S1/S2 ordering are
  `HYPOTHESIS25`.
- The authors explicitly state that further experiments are required.

The sequential proposals are thermodynamically coherent hypotheses, not the
accepted AE4 cycle:

| Proposed positive Cl-loading event | Intracellular source `(Na,K,Cl,HCO3,CO3)` | Charge | Status |
| --- | --- | ---: | --- |
| `Cl_o + Na_o + HCO3_i + K_i -> Cl_i + Na_i + HCO3_o + K_o` | `(+1,-1,+1,-1,0)` | `0` | `HYPOTHESIS25`; Na and Cl always move together, conflicting with the direct 2016 low-Cl Na direction if used as the sole mode |
| `Cl_o + C_i + CO3_i -> Cl_i + C_o + CO3_o` | Na: `(-1,0,+1,0,-1)`; K: `(0,-1,+1,0,-1)` | `0` | Smallest carbonate alternative; carbonate transport was not measured |
| `Cl_o + Na_o + CO3_i + 2K_i -> Cl_i + Na_i + CO3_o + 2K_o` | `(+1,-2,+1,0,-1)` | `0` | `HYPOTHESIS25`; unmeasured carbonate and the same Na/Cl direction conflict |

At the 2016 salivary concentration tuple, the sequential HCO3 proposal has
affinity `+6.643` and is thermodynamically favorable. That does not rescue it
as the sole mode: its fixed source direction contradicts direct Na transport
in the 2016 low-Cl experiment. It may be tested only as a parallel or
condition-dependent cycle with an independently measurable switch.

## Minimal SR2-112 derivation

Let `o` denote the interstitial/bath side, `i` the cytosol, and
`B=HCO3-`. Positive `J_s` is chloride loading through cation branch
`s in {Na,K}`:

\[
Cl_o+s_i+2B_i\rightleftharpoons Cl_i+s_o+2B_o.
\]

### Source matrix, carbon, and charge

In intracellular species order `(Na,K,Cl,B)`, the source matrix is

\[
S_{AE4}=
\begin{pmatrix}
-1&0\\
0&-1\\
+1&+1\\
-2&-2
\end{pmatrix},
\qquad
\dot n_i^{AE4}=S_{AE4}
\begin{pmatrix}J_{Na}\\J_K\end{pmatrix}.
\]

With charge row `z=(+1,+1,-1,-1)`, `z^T S_AE4=(0,0)` exactly. In the
modern conserved acid-base coordinates, each exported bicarbonate removes one
TIC unit and one alkalinity equivalent, so

\[
\begin{aligned}
\dot n_{Na,i}&=-J_{Na}, & \dot n_{K,i}&=-J_K,\\
\dot n_{Cl,i}&=J_\Sigma, & \dot n_{TIC,i}&=-2J_\Sigma,\\
\dot n_{Alk,i}&=-2J_\Sigma, & J_\Sigma&=J_{Na}+J_K.
\end{aligned}
\]

This gives the conserved-coordinate charge identity
`d(Na+K-Cl-Alk)/dt=0` for AE4 alone.

### Affinity and reversal

Using dimensionless activities, the branch affinity is

\[
\mathcal A_s=
\log\frac{a_{Cl,o}a_{s,i}a_{B,i}^2}
{a_{Cl,i}a_{s,o}a_{B,o}^2}.
\]

The exact reversal surface is

\[
[Cl]_i^*=[Cl]_o\frac{[s]_i}{[s]_o}
\left(\frac{[B]_i}{[B]_o}\right)^2
\]

under the ideal-activity approximation. With the 2016 illustrative tuple

`Cl_o/Cl_i=124.6/50.1`, `B_i/B_o=19/24.7`,
`Na_i/Na_o=15.5/151.6`, and `K_i/K_o=139.5/3.4`,

| Branch | Affinity | `Cl_i*` | Interpretation |
| --- | ---: | ---: | --- |
| Na | `-1.89405` | `7.538 mM` | Isolated Na-112 branch runs opposite physiological Cl loading at this tuple |
| K | `+4.10065` | `3025.0 mM` | K-112 branch favors Cl loading |

These are `EXACT` consequences of the declared event and an ideal-activity
calculation, not measured branch fractions.

### Shared carrier and Na/K competition

The smallest exact QSS reduction has outward and inward conformations `O,I`:

\[
O\underset{b}{\overset{a}{\rightleftarrows}}I,
\qquad
I\underset{d_s}{\overset{c_s}{\rightleftarrows}}O,
\quad s\in\{Na,K\}.
\]

The common edge is the Cl path. The two return edge pairs are cation-resolved.
A zero internal-energy gauge and symmetric barriers give

\[
\frac{a}{b}=\frac{a_{Cl,o}}{a_{Cl,i}},
\qquad
\frac{c_s}{d_s}=\frac{a_{s,i}a_{B,i}^2}{a_{s,o}a_{B,o}^2},
\qquad
\log\frac{ac_s}{bd_s}=\mathcal A_s.
\]

Define `C=c_Na+c_K`, `D=d_Na+d_K`, and
`Sigma=a+b+C+D`. Then

\[
p_O=\frac{b+C}{\Sigma},\qquad
p_I=\frac{a+D}{\Sigma},
\]

\[
J_s=T_{AE4}G_R(t)(c_sp_I-d_sp_O).
\]

`T_AE4` is in fmol, rates are in `s^-1`, and `G_R(t)>=0` is the external
regulatory capacity gain; hence `J_s` is in fmol/s. The two branches compete
through the same `p_O,p_I`. Their relative attempt-rate fraction

\[
f_{Na}^{mix}=\frac{k_{Na}}{k_{Na}+k_K}
\]

is explicitly exposed because no mixed Na/K competition experiment identifies
it. Equal Na/K attempts are a named model, not a measured equality.

The shared graph also permits cation slip. If `J_Na=-J_K`, net Cl, TIC, and
alkalinity sources vanish while the transporter still exchanges Na and K. The
slip affinity is

\[
\mathcal A_{Na}-\mathcal A_K
=\log\frac{a_{Na,i}a_{K,o}}{a_{Na,o}a_{K,i}}.
\]

The correct multicycle thermodynamic check is total entropy production,

\[
\frac{\dot S_{prod}}{R}
=10^{-15}\left(J_{Na}\mathcal A_{Na}+J_K\mathcal A_K\right)\ge0,
\]

where the factor converts fmol/s to mol/s. A shared branch may not be judged
from net Cl current alone.

## SR5 cooperative extension

The exact Task 13 SR5 graph adds a same-reservoir catalytic-cation detour to
each loaded branch. The extra cation appears on both sides of that detour and
has zero net source. It may change conductance but not the SR2 affinity or
reversal.

For whole-cell integration, the optional fast-gate limit is

\[
g_s(x)=f_s+(1-f_s)\frac{x^{n_s}}{K_s^{n_s}+x^{n_s}},
\qquad
c_s\mapsto g_sc_s,\quad d_s\mapsto g_sd_s.
\]

This symmetric scaling preserves local detailed balance. Using the 2016
apparent `EC50/Hill` values in `g_s` is only a shape sensitivity. The reported
`Rmin` values remain assay backgrounds and are not converted into AE4 current.
A literal single catalytic binding site would have its own integer occupancy
law; the noninteger apparent Hill summaries do not identify that law.

Thus the retained hierarchy is:

1. `SR2_SHARED_112_QSS`: smallest source/thermodynamic skeleton;
2. `SR5_FAST_GATE_112_QSS`: nested cooperative sensitivity;
3. explicit SR5 carrier states only if mixed-cation or pre-steady-state data
   require occupancy memory.

## Mutation logic

A mutation changes branch conductance by scaling forward and reverse rates
together:

\[
(c_s,d_s)\mapsto \rho_s^{mut}(c_s,d_s),\qquad \rho_s^{mut}\ge0.
\]

This preserves `c_s/d_s`, affinity, and reversal. The only robust current
hierarchy required by the 2025 experiment is

\[
\rho_{Na}^{T448I-T756A}\approx0,
\qquad
\rho_K^{T448I-T756A}>0.
\]

The data do not decide whether T448, T756, D709, or the MD-proposed T713
interaction changes binding, release, or conformational flipping. Assigning a
mutation to one microscopic edge is therefore rejected unless the reverse
edge/state energy is changed consistently and new kinetic data select that
placement.

## QSS versus explicit carrier-state dynamics

Task 13 used carrier occupancies algebraically. Its slow/unit/fast common
barriers were gauge choices, not physical seconds. Neither the 2016 initial
rate traces nor the 2025 pH traces report an independently isolated AE4
relaxation spectrum, transporter abundance, or bound-substrate inventory.
There is therefore no evidence-based physical carrier time constant to insert.

For the reduced two-state graph,

\[
\dot p_I=(a+d_{Na}+d_K)(1-p_I)
-(b+c_{Na}+c_K)p_I,
\]

and its exact relaxation time is

\[
\tau_{carrier}=\frac{1}{a+b+c_{Na}+c_K+d_{Na}+d_K}.
\]

The implementation reports this value but requires all attempt rates from the
caller; it supplies no source-claimed default. QSS is acceptable only
conditionally if an independently calibrated rate set makes carrier
relaxation well separated from ionic, water, and beta/cAMP/PKA dynamics, bound
substrate is negligible compared with compartment amounts, and explicit-state
sensitivity checks do not change the accepted whole-cell result. An
order-of-magnitude separation is a useful declared numerical gate, not a
measured fact.

If QSS fails, the full explicit graph must track edge-level reservoir exchange
and bound-ion inventory. Integrating carrier occupancy while applying only
completed-cycle sources would lose transient substrate storage and is not
conservation-correct.

## Coupling to beta/cAMP/PKA regulation

The 2021 evidence requires a PKA-dependent AE4 response and S173 dependence,
but it does not identify direct S173 phosphorylation, an AE4 phosphorylation
rate, or the regulated carrier edge. The transporter module therefore accepts
the output of the separate dynamic regulatory subsystem.

| Coupling | Thermodynamic treatment | Status |
| --- | --- | --- |
| Common active capacity | `J_s(t)=G_R(t)J_s^0`, `G_R>=0` | Minimal retained coupling; preserves occupancy, affinity, reversal, and Na/K routing |
| Regulated active fraction | `J_s=(1-p_R)J_s^U+p_RJ_s^R` | Equivalent to common capacity if the two layers differ only by uniform turnover |
| Matched state-specific barrier | Multiply both directions of one edge by the same positive factor | Admissible nested hypothesis; changes occupancy/relaxation but not cycle affinity |
| One-way clockwise rate change | Changes a forward/reverse product without an explicit ATP-coupled cycle | Rejected; it creates an unsupported reversal shift |
| Stoichiometry switch under PKA | Changes source vector | Rejected unless simultaneous ion-flux data establish it |

S173A should retain basal transport but remove the regulated increment in a
qualitative mutant simulation; S273A should retain regulation. Neither rule
proves that S173 itself is phosphorylated.

## Implementation and tests

The implementation is in `src/modern_full_model/transporters.py`.

```python
flux = evaluate_ae4_qss(
    environment,
    parameters,
    regulation_gain=regulatory_output,
)
cell_sources = flux.intracellular_conserved_sources_fmol_s
```

The exact conserved-source keys are `na_i`, `k_i`, `cl_i`, `tic_i`, and
`alkalinity_i`, all in fmol/s. The result also reports explicit `hco3_i`,
separate Na/K branch currents, branch affinities, entropy production,
transported charge/current, QSS residual, relaxation time, and the mixed-bath
Na attempt fraction.

Focused tests in `tests/test_modern_ae4.py` cover:

- exact source charge, carbon, and alkalinity accounting;
- the Task 13 salivary affinities and reversal values;
- zero current at global equilibrium;
- carrier QSS closure and nonnegative total entropy production;
- common-capacity regulation without a reversal shift;
- double-mutant Na collapse with retained K transport;
- SR5 gating without an affinity change;
- explicit reporting of the unmeasured mixed-bath fraction and carrier time;
- output units and source-map keys.

## Unresolved mechanism choices

The modern full model must carry these as explicit uncertainty rather than
silently selecting them:

1. net count/species: `1:1:2 HCO3` versus `1:2:3`, `2:1:3`, or explicit
   carbonate;
2. mixed-bath Na/K attempt-rate ratio and possible slip;
3. origin of the apparent cooperative cation response;
4. the microscopic edge affected by T448/T756/D709 changes;
5. common-capacity versus state-specific PKA modulation;
6. physical AE4 carrier relaxation time and transporter amount;
7. whether any 2025 sequential mode exists in parallel with the directly
   supported 2016 branch directions.

The most discriminating transporter experiment is a reversible AE4-isolated
protocol in voltage-clamped AE4-expressing cells with independently controlled
intracellular and extracellular Na, K, Cl, pH, and inorganic carbon. It should
measure AE4-specific Na, K, Cl, and total-carbon flux simultaneously at
subsecond-to-seconds resolution across mixed-cation and reversal conditions,
using WT and T448I-T756A with and without independently measured PKA
activation. Particle ratios separate one transported cation from a catalytic
cation, mixed-bath cross-saturation identifies shared-pool routing, pH changes
separate HCO3 from CO3, and relaxation spectra test common-capacity versus
state-specific regulation.

Until that experiment exists, the correct whole-cell choice is the minimal
SR2-112 QSS core with explicit uncertainty, not a hard-coded 2025 cycle.
