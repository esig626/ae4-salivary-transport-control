# Round 5: minimal chassis reconstruction

## Controlling conclusion

The smallest independently supported whole-cell topology extension was built
and tested, not merely projected: total Na/K-ATPase capacity and total
Ca-activated K conductance were redistributed between apical and basolateral
membranes, local extracellular K and local membrane voltage were used, lumen
Na/K terms were added, and both membrane-current equations were re-solved.
The extension is exactly nested in the fixed chassis at zero apical fractions.

It does **not** reconstruct or identify the AE4-null phenotype on the available
evidence.

- With the smallest conditional transporter-level survivor, every WT-only
  capacity attempt for M0, pump-only, K-only, and coupled split controls fails
  exact balance closure at the derived `height=200 um` domain boundary. There
  is no valid WT capacity and therefore no admissible WT denominator.
- In the candidate-independent AE4-null field, re-solved split-topology roots
  move the released chloride and pH errors only modestly. The best predeclared
  20--40% sensitivity tested is still `6.969` Cl SEM and `23.871` pH SEM from
  the released means.
- At both named conditional Stage-B target slices, the cation split leaves the
  nonzero bicarbonate equation unchanged. No choice of pump/K localization
  fractions can make either full target slice stationary without an acid-base
  addition.
- A conserved carbon/buffer plus genotype-invariant NHE1 block has three exact
  independent reaction directions but rank only one when projected directly
  onto the released Cl/pH pair. At least a two-dimensional reaction-flux
  nullspace remains before unknown carbon totals, buffer pools, and boundary
  fluxes are counted.

Thus the source-supported apical cation topology is a real omission in the
historical comparison chassis, but it is not a calibrated repair. The current
evidence cannot select the split fractions, parameterize the carbon/NHE block,
or produce a WT-validated model for a held-out secretion prediction. This is a
chassis/evidence bottleneck, not a conclusion that AE4 biology is unimportant.

Machine-readable results are in
`results/13_state_resolved_ae4/chassis_reconstruction_scan.csv`; implementation
and tests are in `src/state_resolved_ae4/chassis.py`,
`src/state_resolved_ae4/run_chassis_reconstruction.py`, and
`tests/test_state_resolved_ae4_chassis.py`.

## Stage gate and information firewall

Round 5 opened only after `stage_a_summary.json` recorded:

- five conditional transporter-assay survivors;
- zero valid WT capacity calibrations;
- zero fixed-chassis WT Gate-2 passes; and
- `strict_secretion_used_in_any_fit=false`.

Only `Cl_i=36.50 +/- 1.60 mM` and `pH_i=6.89 +/- 0.02` were released for
Stage-B localization. The reconstruction script does not read
`heldout_targets.csv`, and every output row has
`secretion_target_read=False`. No saliva magnitude, time-course value, null
volume, null Na/K concentration, or null lumen composition is used for a
fraction, capacity, root, or path choice.

The M1 fraction values were declared before calculation:

- `0` is the exact historical nested control;
- `0.2` and `0.4` are published-model/parotid-informed sensitivity values, not
  measurements of submandibular functional fractions; and
- `0.398` is the historical area-ratio initialization, not morphometric or
  flux evidence.

The strongest direct apical pump/K evidence comes from mouse parotid acini
([Almassy et al. 2012](https://doi.org/10.1085/jgp.201110718),
[Almassy et al. 2018](https://doi.org/10.1007/s00424-018-2109-0)). The
AE4-null phenotype is from mouse submandibular gland. In gland-matched
submandibular acini, local apical Ca activated Cl but did not produce a local K
current until Ca spread more broadly (Harmer et al. 2005, PMID 15308468).
Therefore none of the tested `f_P` or `f_K` values is called calibrated.

## M1 equations: split cation-current architecture

Let `P_a,P_b` be pump-cycle rates and `J_Ka,J_Kb` positive outward K particle
rates. The implementation conserves total capacity/conductance:

\[
P_a=f_P P_{max}\,v_{pump}(Na_i,K_l),\qquad
P_b=(1-f_P)P_{max}\,v_{pump}(Na_i,K_e),
\tag{1}
\]

\[
g_{K,a}=f_K g_{K,total},\qquad
g_{K,b}=(1-f_K)g_{K,total}.
\tag{2}
\]

The historical pump law is evaluated with the local external K:

\[
v_{pump}(Na_i,K_o)=k_{pump}
\frac{(10^{-3}K_o)^2(10^{-3}Na_i)^3}
{(10^{-3}K_o)^2+\alpha(10^{-3}Na_i)^3}.
\tag{3}
\]

With `V_a=phi_i-phi_l`, `V_b=phi_i-phi_e`, and `V_t=V_a-V_b`,

\[
J_{K,a}=\frac{g_{K,a}}F(V_a-E_{K,a}),\quad
E_{K,a}=V_T\log(K_l/K_i),
\tag{4}
\]

\[
J_{K,b}=\frac{g_{K,b}}F(V_b-E_{K,b}),\quad
E_{K,b}=V_T\log(K_e/K_i).
\tag{5}
\]

The two QSS current balances are solved simultaneously:

\[
J_{Cl,a}+J_{K,a}+P_a+J_{t,Na}+J_{t,K}=0,
\tag{6}
\]

\[
P_b+J_{K,b}-J_{t,Na}-J_{t,K}=0.
\tag{7}
\]

The cell cation contributions are

\[
\dot n_{Na,i}|_{pump}=-3(P_a+P_b),
\qquad
\dot n_{K,i}|_{pump,K}=2(P_a+P_b)-J_{K,a}-J_{K,b},
\tag{8}
\]

and the newly complete lumen terms are

\[
\dot n_{Na,l}|_{M1}=J_{t,Na}+3P_a-q_{out}Na_l,
\tag{9}
\]

\[
\dot n_{K,l}|_{M1}=J_{t,K}+J_{K,a}-2P_a-q_{out}K_l.
\tag{10}
\]

At `f_P=f_K=0`, Eqs. (1)--(10) reproduce the fixed-chassis RHS exactly to
floating-point precision. For all tested nonzero splits, both current
residuals, the differentiated proton identity, and the redundant lumen-Cl
balance close below `2e-12` in unit tests. Thus the result is not caused by
omitting the cation terms that make membrane location meaningful.

## WT-only nested calibration attempts

The WT test uses SR2 shared-state `1:2:3` with the frozen slow-common gauge,
the simplest of the five conditional transporter-assay survivors. The loaded
energy lumps, Na/K attempt ratio, and common barrier were read from the frozen
pre-Stage-B ledger. Exactly one whole-cell AE4 capacity was attempted against
WT resting chloride; pH was retained as an independent gate, not included in
the objective. The seven balance equations and WT Cl give eight residuals for
seven reduced steady coordinates plus log capacity. Five log-capacity starts
were retained per topology.

| Nested topology | `f_P` | `f_K` | Attempted capacity | Derived height (um) | Max raw RHS | WT pH diagnostic | Valid capacity / WT gate |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| M0 fixed basal | 0 | 0 | `6.6883e-6` | `199.999995` | `3.0505e-5` | `7.52055` | no / no |
| pump only | 0.2 | 0 | `1.69169e-5` | `199.999989` | `3.8655e-5` | `7.52127` | no / no |
| K only | 0 | 0.2 | `6.14451e-6` | `199.999987` | `3.1472e-5` | `7.52022` | no / no |
| coupled split | 0.2 | 0.2 | `4.46740e-3` | `199.999995` | `4.3730e-5` | `7.52293` | no / no |
| coupled split | 0.4 | 0.4 | `2.25735e-2` | `199.999984` | `5.4024e-5` | `7.53393` | no / no |

Acceptance required maximum raw RHS `<=1e-8` on the declared derived-height
domain `2 < height < 200 um`. Every attempt instead approaches the upper
domain boundary and misses closure by more than three orders of magnitude.
The capacities are failed optimizer coordinates, not parameter estimates.
The M0 value differs modestly from the earlier eight-start failed attempt,
which is further evidence that an invalid boundary minimum must not be
promoted. No valid WT row exists, so the AE2 dynamic secretion comparison is
`N/A at upstream WT gate`, not silently passed or deleted.

## AE4-null re-solved topology sensitivity

After AE4 deletion every state-resolved candidate short-circuits to the same
zero source. The following roots therefore test only the non-AE4 topology.
They were followed by `0.005` fraction continuation and required maximum raw
RHS `<=1e-8`; reported accepted roots close at `2.5e-12` or better.

| Topology sensitivity | `f_P` | `f_K` | Cl (mM) | pH | Height (um) | Normalized Cl/pH distance |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| M0 | 0 | 0 | `48.78222` | `7.378474` | `65.66070` | `25.6016` |
| pump only | 0.2 | 0 | `48.26484` | `7.372056` | `62.99661` | `25.1994` |
| K only | 0 | 0.2 | `48.73701` | `7.378026` | `65.42914` | `25.5718` |
| coupled | 0.2 | 0.2 | `48.30658` | `7.374313` | `63.35390` | `25.3150` |
| coupled | 0.4 | 0.4 | `48.33636` | `7.375324` | `63.55970` | `25.3688` |
| pump-high/K-low | 0.4 | 0.2 | `47.65116` | `7.367421` | `60.34311` | **`24.8677`** |
| pump-low/K-high | 0.2 | 0.4 | `48.95392` | `7.380977` | `66.63866` | `25.7533` |
| historical area sensitivity | 0.398 | 0.398 | `48.33420` | `7.375296` | `63.54864` | `25.3670` |

The distance is only a transparent display metric,

\[
D=\sqrt{[(Cl-36.5)/1.6]^2+[(pH-6.89)/0.02]^2},
\tag{11}
\]

not a fit objective. Fractions were not optimized. The best listed sensitivity
still predicts `Cl=47.6512 mM` and `pH=7.36742`, or `+6.969` and `+23.871`
reported SEM from the released values. It reduces `D` by only `2.87%` from M0.

Pump-only continuation toward 40% ceased closing at `f_P=0.295`, and K-only
continuation ceased at `f_K=0.290`. These failures are saved. They are **not**
claims that no remote branch exists; they prevent silently substituting a
failed numerical point for the requested controls. Coupled paths remain
solvable because the apical pump and apical K routes support each other's lumen
K/current closure, but that structural coherence does not repair the ionic
phenotype.

## The two conditional target slices remain nonstationary

Round 4 proved that released Cl/pH define a five-dimensional target-state
manifold. Only its two declared diagnostic slices were tested here:

1. `model_anchored`: unmeasured lumen Na/K, height, and intracellular Na/K held
   at the common fixed-chassis null root; and
2. `wt_reference_cations`: the same root lumen values with historical
   `height=28.7`, `Na_i=25`, and `K_i=120`.

Neither slice is an observed knockout state. Every M0/M1 row has a nonzero
seven-state RHS at both slices. More decisively, at a fixed state M1 changes
only cation/current/lumen terms. The bicarbonate RHS is invariant across all
fractions:

| Slice | Invariant `dHCO3_i/dt` raw row | Exact stationary under any tested M1 split? |
| --- | ---: | --- |
| model anchored | `-0.00654501062243` | no |
| WT-reference cations | `-0.00180648185387` | no |

Thus no pump/K fraction can make either full target slice stationary without
an acid-base source or a different unmeasured target-state member. The result
does not identify which of those two explanations is correct, because the
slice coordinates are assumptions and the required source rotates with rank
five across the target manifold.

## M2 acid-base family and why it is not identifiable

The smallest conservation-correct reaction block tested structurally is

\[
CO_2+H_2O\rightleftharpoons H^++HCO_3^-,
\qquad
BH\rightleftharpoons B^-+H^+,
\tag{12}
\]

with genotype-invariant basolateral NHE1

\[
Na_e+H_i\rightleftharpoons Na_i+H_e.
\tag{13}
\]

For rates `r_C`, `r_B`, and `J_NHE`, the exact source signatures in
`(Na_i,HCO3_i,H_i,CO2_i,BH,B-)` are

\[
\sigma_{NHE}=(1,0,-1,0,0,0),
\]

\[
\sigma_C=(0,1,1,-1,0,0),\qquad
\sigma_B=(0,0,1,0,-1,1).
\tag{14}
\]

The implementation verifies exact conservation of inorganic carbon, buffer
sites, and charge for mass-action forward/reverse rates. The three signatures
have exact full rank `3`. In the direct released `(Cl_i,H_i)` observation
projection, however, they have rank `1` and flux nullity `2`. Unknown total
inorganic carbon, total buffer, CO2 boundary flux, and any apical bicarbonate
loss add further degrees of freedom.

This is why an acid-base-only nonlinear fit was not manufactured: the 2015
experiment did not measure the matched carbon totals, buffer pool, carbonic-
anhydrase flux, or absolute NHE flux. It found no detectable AE4-genotype
change in its NHE-dependent alkalinization assay, so an AE4-null-specific NHE1
multiplier is prohibited. Likewise, unchanged NKCC1 assay activity prohibits
an arbitrary knockout capacity increase. The carbonate AE4 variants also
remain an unresolved chemistry interface: the fixed seven-state adapter
correctly rejects them because it lacks distinct carbonate/proton/carbon
states; this is not a biological rejection of those families.

## M3 smallest coupled reconstruction

The smallest source-supported coupled candidate is M1 plus conserved
carbon/buffer chemistry and shared NHE1. Its topology can span lumen/current,
Na/H, and carbon directions that M1 alone cannot. It cannot presently be
calibrated without introducing at least:

- submandibular `f_P` and `f_K` or membrane-partitioned functional fluxes;
- total intracellular inorganic carbon and buffer pool/capacity;
- a matched CO2/carbonic-anhydrase boundary flux or kinetic scale; and
- sufficient AE4-null state coordinates to select one member of the
  five-dimensional Cl/pH-compatible manifold.

Adding arbitrary values for these quantities would create enough degrees of
freedom to fit two released outputs by dimension and would not constitute a
minimal source-constrained reconstruction. M3 is therefore identified as the
next mechanistic chassis, but not accepted as a predictive model.

## What was repaired, what was not, and what would decide it

M1 repairs a precise structural defect: apical pump/K events now have their
correct lumen Na/K and membrane-specific current consequences, while total
capacity and conductance remain fixed. It does not repair WT balance closure,
the common-null Cl/pH mismatch, or the acid-base balance at either conditional
target slice. M2 repairs conservation structure conceptually and in unit
tests, but its rates/pools are not identifiable from the released evidence.

The decisive missing measurement is a matched native mouse-submandibular
balance experiment under the same WT/AE4-null stimulation protocol. At minimum
it must measure `Na_i`, `K_i`, `Cl_i`, pH, and cell volume simultaneously;
adding luminal/effluent Na and K plus membrane-sided ouabain/K-channel current
partition would directly identify the two M1 extended directions

\[
(3P_a,\ J_{K,a}-2P_a)
\]

instead of inferring them from intracellular Cl/pH. Parallel total inorganic
carbon or CO2-flux/NHE measurements would then determine whether the remaining
stationarity error is carbon/NHE support rather than another unmeasured member
of the target manifold. Without these measurements, there is no WT-validated
chassis row to send to the held-out secretion round.
