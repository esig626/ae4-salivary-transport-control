# Task 13B WT calibration and resting-root geometry

## Result and scope

The corrected-water G2b model produced **nine historical WT resting
branches** under its then-declared profile panel.  They span three unresolved
Na/K routing choices and three apical-topology choices.  Every branch closes on exact
cell and lumen electroneutral manifolds, has a full-column-rank normalized
Jacobian, avoids every core-state and fitted-parameter bound, and passes the
declared numerical and WT-resting gates.  They are no longer production
candidates: the old `G=10` common conductance multiplier fails the native
whole-cell source map, and their downstream G5 gland scale fails the absolute
geometry gate.  `g5_supersession.json` controls this rejection while the old
gate and manifest remain immutable historical evidence.

This is a historical resting-model result, not a validated full-model
classification.
R0 is retained only as a static negative control and root scaffold.  The
mandatory dynamic beta/cAMP/PKA subsystem and stimulated WT flow-shape gate
must also pass before any genotype secretion prediction can be evaluated.
The first dynamic staging of the declared reference is numerically robust but
does not pass the WT minute-flow-shape envelope, so the resting result is not
promoted to a final model.

The calibration code does not import or open the strict AE4-null secretion
ledger.  No knockout secretion magnitude, ratio, timing, or digitized curve
was used in a target, score, branch rule, parameter bound, or model choice.

## Permitted targets and parameter tiers

The fit uses the smallest successful WT-only closure set in this generation.
The target roles are deliberately narrower than the search bounds.

| Quantity | Role | Value or gate | Tier and provenance |
|---|---|---:|---|
| WT intracellular Cl | fitted target | 50.10 +/- 1.50 mM SEM | Tier 2; direct matched WT measurement, P15 Table 1 |
| WT intracellular pH | fitted target | 6.91 +/- 0.07 SEM | Tier 2; direct matched WT measurement, P15 Table 1 |
| WT cell volume | fitted closure convention | 1.30 pL | Tier 3; published-model lineage, not a primary uncertainty |
| WT intracellular Na | out-of-objective gate | 10--35 mM | Tier 3 historical-lineage physiology screen |
| WT intracellular K | out-of-objective gate | 100--160 mM | Tier 3 historical-lineage physiology screen |
| apical voltage | out-of-objective gate | -80 to -30 mV | Tier 3 historical-lineage physiology screen |
| basolateral voltage | out-of-objective gate | -90 to -35 mV | Tier 3 historical-lineage physiology screen |

Tier meanings follow the Task 13B prompt: Tier 0 is exact/source-known, Tier 1
is directly estimated from a primary assay, Tier 2 is calibrated to WT resting
physiology, Tier 3 is a weakly constrained bounded nuisance or profile axis,
and Tier 4 is unsupported and should be fixed or removed.  No Tier 4 parameter
was newly fitted by this calibration layer.

The selected reference has the following fitted or profiled values.  These
values define a model construction; they are not claims of molecular
identification.

| Parameter | Value | Unit | Tier | Interpretation |
|---|---:|---|---:|---|
| AE4 carrier amount at the declared attempt-rate gauge | 0.1501264265070657 | fmol carrier | 2 | effective WT-rest capacity; all declared attempt rates retain the 1 s^-1 gauge, so carrier amount and attempt rate are not separately identified |
| NHE1 capacity | 0.007440980590414482 | fmol s^-1 | 2 | WT acid-base closure nuisance; no knockout pH target used |
| common membrane conductance scale | 10 | dimensionless | 3 | historical profile, now rejected: it does not represent the native 31.4 nS CaCC and 14 nS total-K literals |
| OTHER cell-impermeant osmoles | 124.71681310732373 | fmol osmoles | 3 | excludes the separate finite 30-fmol buffer-site pool |
| native mixed-bath Na attempt fraction | 0.10 | fraction | 3 | unresolved routing axis |
| apical pump fraction | 0.10 | fraction | 3 | historical panel value, replaced by the area-density-derived transfer `0.075075` with `0.056878--0.094586` sensitivity |
| apical K fraction | 0.10 | fraction | 3 | historical panel value; native reconstruction treats `0.20/0.30/0.40` as explicit modeling assumptions, not measurements |

## Exact charge-manifold root system

The twelve-state core contains two exact left-null charge directions.  A
least-squares solve of all twelve raw balances is therefore rank deficient.
The production calibration does not use `ModernFullModel.solve_resting_root`.
It instead solves in ten physical coordinates

\[
x=(Na_i,K_i,TIC_i,pH_i,V_i,Na_l,K_l,TIC_l,pH_l,V_l).
\]

Total alkalinity is reconstructed from finite carbon/buffer speciation.  Cell
and lumen chloride are eliminated exactly:

\[
Cl_i=Na_i+K_i-TA_i-X_i/V_i,
\qquad
Cl_l=Na_l+K_l-TA_l.
\]

Here `X_i` is the registered fixed cell-anion charge, which is not varied by
the calibration.  The independent core rows are
`(0,1,2,3,5,6,7,8,9,11)`.  Alkalinity rows 4 and 10 are omitted during the
solve and verified afterwards.  Their exact identities are

\[
R_{Na,i}+R_{K,i}-R_{Cl,i}-R_{TA,i}=0,
\qquad
R_{Na,l}+R_{K,l}-R_{Cl,l}-R_{TA,l}=0.
\]

Amount, volume, and regulatory rows use separate declared numerical scales.
No maximum over heterogeneous fmol/s, pL/s, fmol-equivalent/s, and s^-1 rows
is used as a physical norm.  Numerical rank is computed from a dimensionless
Jacobian: every residual row is scaled, and every physical column is
multiplied by its declared bound span before the SVD.

## Deterministic multistart and profile design

The predeclared panel crossed:

- apical pump/K fractions `(0.03,0.03)`, `(0.10,0.10)`, and `(0.10,0.30)`;
- native mixed-bath Na attempt fractions `0.05`, `0.10`, `0.20`, `0.50`, and
  `0.90`;
- common membrane conductance scales `10`, `30`, and `100`.

All 45 slices received a three-start deterministic screen: one preferred
anchor and two non-collinear quartile anchors.  Every slice with a passing
root was then rerun with seven starts using the same named anchors plus a
seeded SciPy Latin-hypercube design.  All distinct normalized-coordinate
clusters were retained.

| Conductance | Routing 0.05 | Routing 0.10 | Routing 0.20 | Routing 0.50 | Routing 0.90 |
|---:|---|---|---|---|---|
| 10 | pass in all 3 topologies | pass in all 3 topologies | pass in all 3 topologies | no exact interior WT root | no exact interior WT root |
| 30 | no passing slice | no passing slice | no passing slice | no passing slice | no passing slice |
| 100 | no passing slice | no passing slice | no passing slice | no passing slice | no passing slice |

Thus nine of 45 slices pass.  Each passing slice contains one detected branch,
has normalized rank/nullity 13/0, has no disallowed boundary hit, and was
recovered from five to seven of seven confirmation starts.  The complete
failure attempts and raw, unit-specific residual diagnostics remain in the
machine-readable artifacts; failed slices were not silently dropped.

## Declared reference and fixed-parameter geometry

All nine resting branches are retained as historical reconstruction evidence,
not current source-compatible candidates.  For the old downstream WT dynamic
staging, the predeclared balanced-apical, K-biased, G=10 slice was the
reference; ties were then broken by the largest converged-start cluster and
smallest dimensionless residual.  This rule was fixed without consulting any
knockout secretion result.

Reference root `G2_BALANCED_APICAL_K_BIASED_G10:M00` has:

| Observable | Value |
|---|---:|
| intracellular Na | 14.7567734958 mM |
| intracellular K | 111.1037463712 mM |
| intracellular Cl | 50.1000000000 mM |
| intracellular pH | 6.9100000000 |
| cell volume | 1.3000000000 pL |
| apical voltage | -62.2174607678 mV |
| basolateral voltage | -72.6637953276 mV |
| lumen Na/K/Cl | 146.1356841095 / 4.0873611205 / 138.1242234065 mM |
| lumen pH / volume | 7.2764536489 / 0.10164106694 pL |
| resting local outflow | 0.000328213388 pL s^-1 |

Its largest dimensionless independent residual is
`3.64e-13`; normalized fitted-system rank/nullity is 13/0; charge, current,
omitted-row, positivity, and boundary gates pass.

With these parameters fixed, a separate 33-attempt search (the exact
calibrated root plus 32 named/Latin-hypercube starts) found one distinct branch.
Thirty starts converged to it, with normalized rank/nullity 10/0 and no bound
hit.  Two initial Latin points were charge-inadmissible and explicitly
reported; one remote start terminated at the intracellular-pH lower bound with
a dimensionless residual near 1.96.  None supplied a second valid branch.

The same root was independently transcribed and refined without calling
`ModernFullModel.evaluate`: the independent solve recovered rank/nullity 10/0,
the same Cl, pH, volume, and voltages, and native-unit RHS no larger than
approximately `2.52e-16`.

## Water-scale correction and generation decision

The pre-audit G2a root used placeholder hydraulic coefficients and is rejected.
The Palk/2018 published-model values convert as
`P_model = P_source * 10^9` because numeric mM contributes `10^-3 mol/L` and
one litre is `10^12 pL`.  G2b therefore uses

- apical `0.00432 pL s^-1 (mOsm/L)^-1`;
- basolateral `0.0515 pL s^-1 (mOsm/L)^-1`;
- paracellular `0.000260 pL s^-1 (mOsm/L)^-1`.

These are absolute, area-folded historical-lineage coefficients in the source
equation `q=P Delta osm`; no extra membrane-area factor is applied.  The old
root is recorded only as a rejected generation and is not a candidate.

## Minimality, dynamic regulation, and unresolved degrees

The successful square system fits exactly three positive WT-only quantities:
AE4 effective capacity, NHE1 capacity, and OTHER impermeant osmoles.  The
historical conductance scale was profiled and required to pass both individual
voltage gates; that profiling is now rejected by the native source map.  Na,
K, and both voltages remain out of the old objective. Generation-local numerical
ablations fix one fitted quantity at a time to its declared pre-fit preferred
value and reoptimize the other two.  Their exact results are recorded in
`wt_parameter_ablations.csv/json`; this tests minimality within G2b but is not
a structural-identifiability claim.

| Quantity removed from fit | Fixed pre-fit value | Best max dimensionless residual | Boundary/rank result | Exact root |
|---|---:|---:|---|---|
| AE4 effective carrier gauge | 0.05 fmol carrier | 0.6232 | NHE1 driven to its lower bound; normalized rank 11/12 | no |
| NHE1 capacity | 0.005 fmol s^-1 | 0.01266 | no bound; normalized rank 12/12 | no |
| OTHER impermeant osmoles | 110 fmol | 0.5755 | no bound; normalized rank 12/12 | no |

The exact-root tolerance is `1e-7`; therefore none of the three two-parameter
ablations closes the selected G2b system.  This supports the three-quantity
set as generation-minimal relative to the declared pre-fit values.  It does
not prove that a different added mechanism could not compensate.

The dynamic G3 confirmation attaches the mandatory R1 effective
beta/cAMP/PKA/AE4 activation state.  Its beta-free equilibrium is exactly zero.
That state boundary is reported and accepted only because it is the analytic
source-supported basal equilibrium; no core coordinate or fitted parameter
receives this exemption.  The 30-s R1 time constant, its matching numerical
RHS scale, and the 0.25 activated-capacity sensitivity are unmeasured modeling
choices and are not fitted to secretion.  R1 resting closure alone cannot
repair or certify stimulated WT dynamics.

The frozen-reference baseline G5 WT-only R0--R3 gain/timing ensemble is
recorded in `dynamic_frozen_manifest.json`, `regulatory_fit_summary.csv`, and
`wt_dynamic_gate.json`.  Both stiff methods pass the numerical crosscheck, but
the predeclared scale-free minute-flow-shape gate fails.  The later old joint
repair passed that shape gate but required `294,706,403` or `260,704,116`
modeled cells at the selected root.  Both exceed the generous Kondo-derived
one-SMG bound of `51,469,231` cells and are hard failures, not unassessed
warnings.  The native dynamic contract instead freezes R1 time constants
`10/30/90 s`, R2 pairs `(3,10)/(10,30)/(30,90) s`, and their symmetric R3
isomorphs; this contract remains an ensemble because the 2021 evidence has no
acute kinetic trace.

A separate exact-root sustainment audit then ran all 18 dynamic R1--R3
gain/timing members.  Zero of 18 repaired the shape gate.  The best max/min
ratio was `1.732366` against the allowed `1.111111`, corresponding to a
`42.275%` minute-1-to-minute-10 decline; larger regulatory gain worsened the
shape.  These WT-only results are in
`regulatory_sustainment_sensitivity.csv` and
`regulatory_sustainment_summary.json`.  They localize the defect away
from merely choosing another predeclared acute beta/cAMP/PKA time constant.

### Current native reroot, hydraulic diagnostics, and absolute-source screen

The native-source reroot now has a complete three-start WT-rest checkpoint
through the absolute N/AN source classes.  Canonical production mode `H1`
comprises 45 panels and 135 deterministic attempts; 35 numerical roots close,
but zero pass
the independent WT gates because all fail intracellular Cl and K.  Scale-one
HW and HWQ constructions are exact alias controls of H1, not extra panel rows.
Every nonunit HW/HWQ row is explicitly diagnostic-only: HW scales the three
hydraulic coefficients, while HWQ also scales local outflow.

The completed nonunit screens contain 225 HW panels and 225 HWQ panels, with
675 attempts and 175 numerical roots in each mode.  Neither mode has a
WT-passing root.  The combined H1/HW/HWQ checkpoint is therefore 495 panels,
1,485 attempts, 385 numerical roots, and zero WT-passing roots or panels.  For
every tested `H=2.5/3.5/5/6.75/8`, the 35-root intracellular ranges agree
between HW and HWQ to displayed precision: Na `8.107--29.441`, K
`70.840--92.515`, Cl `23.589--29.243 mM`, and pH `6.7522--6.9440` over the H
panel.  This is a numerical diagnostic, not an analytic invariance claim.  It
shows that scaling outflow together with water permeabilities does not repair
the decisive steady-rest Cl/K defect.  Neither diagnostic mode is eligible to
rescue production, regardless of its numerical-root count.

The subsequent production `N_ABS_NKCC`/`AN_ABS_AE4_NHE` tier adds 630 panels
and 1,890 attempts.  Optimizer success is reported for 1,847 attempts, all
1,890 attempts remain admissible, and 981 meet the numerical-root criterion.
After within-panel clustering, 354 panels contain one numerical root and 276
are rootless.  The source-class landscape is:

| Source scale | 0.5 | 2 | 4 | 5 | 6 | 7 | 8 | total |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `N_ABS_NKCC` numerical roots | 40 | 35 | 33 | 31 | 30 | 30 | 30 | 229 |
| `AN_ABS_AE4_NHE` numerical roots | 45 | 25 | 20 | 10 | 10 | 10 | 5 | 125 |
| `N_ABS_NKCC` preliminary WT passes | 0 | 0 | 10 | 0 | 0 | 0 | 0 | 10 |
| `AN_ABS_AE4_NHE` preliminary WT passes | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

All ten preliminary passes are `N_ABS_NKCC` scale 4.  Every one of the five
declared pump/K topologies passes in two contexts only: pump-capacity scale 1
with `AE4NA05`, and pump-capacity scale 2 with `AE4NA20`.  Their WT-rest ranges
are Cl `51.0172--53.0228 mM`, pH `6.83860--7.04827`, Na
`14.3306--18.9814 mM`, K `105.935--118.116 mM`, fixed cell volume `1.3 pL`,
apical voltage `-36.1641` to `-34.5019 mV`, basolateral voltage `-73.5186` to
`-66.5896 mV`, and fitted OTHER `110.024--126.008 fmol`.  The smallest Cl and
pH gate margins are `0.07723 mM` and `0.001728`, respectively.  Maximum scaled
residuals span `1.01e-12--2.41e-11`; all ten normalized Jacobians have rank
11/nullity 0, with no boundary hits.  These are three-start screen results,
not confirmed production roots.

The combined H1/HW/HWQ/N/AN ledger therefore contains 1,125 panels, 3,375
attempts, and 739 numerical roots.  A hierarchy-code audit paused promotion
before confirmation and hardened panel-atomic persistence, both-quartile start
coverage, duplicate-ID/eligibility checks, exact scale-one aliases, complete
upstream-grid validation, and fail-closed 33-start intake.  This engineering
gate changes no model equation, source scalar, WT acceptance interval,
parameter bound, or scientific calibration.  All ten rows remain queued for
seven-start confirmation and, only if retained, explicit 33-start geometry.
They are not yet eligible for WT dynamics, a native pre-reveal freeze, or
held-out reveal.

The exact unresolved degrees after resting calibration are:

1. three surviving native Na/K routing choices (`0.05`, `0.10`, `0.20`);
2. the pump transfer `0.075075` and its `0.056878--0.094586` area-density
   sensitivity, plus K fractions `0.20/0.30/0.40` that remain assumptions;
3. literal native whole-cell maxima of `31.4 nS` CaCC and `14 nS` total K,
   with the Palk gate; no free common conductance multiplier is licensed;
4. the carrier-amount/attempt-rate gauge split;
5. independent physiological values for NHE1 capacity and OTHER impermeant
   osmoles;
6. local-lumen dead volume/outflow geometry and compatibility with the
   Kondo-derived upper cell-count bound when mapping pL/s to gland uL/min;
7. the regulatory family, gain, and activation/deactivation time constants;
8. the WT stimulus/sustainment mechanism needed to pass the minute-flow-shape
   gate.

Accordingly, no old resting branch is claimed to identify topology, routing,
AE4 molecule count, or regulatory kinetics. `STRICT_NATIVE`, both hydraulic
diagnostics, and the complete AN absolute tier fail the WT-rest gate.  The N
absolute tier supplies ten preliminary scale-4 rows, but the next active
iterations are the predeclared seven-start confirmation and 33-start root
geometry, using WT gates only.  Model selection by AE4-null secretion remains
prohibited until a replacement WT dynamic/absolute-scale gate and pre-reveal
manifest are frozen.

## Machine-readable artifacts

- `wt_root_attempts.csv`: every screen/confirmation attempt, boundary hit, and
  unit-specific residual maximum;
- `wt_root_table.csv`: all 45 profile decisions and every retained branch;
- `wt_root_states.csv`: full conserved state vectors;
- `wt_parameter_candidates.csv`: fitted/profiled values with Tier 0--4 roles;
- `wt_solver_diagnostics.json`: raw residual vectors, normalized SVDs, charges,
  currents, conservation checks, and branch membership;
- `wt_selected_candidate.json` and `wt_selected_root_geometry.json`: declared
  reference construction and deep fixed-parameter geometry;
- `wt_parameter_ablations.csv/json`: generation-local minimality checks;
- `wt_dynamic_rest_calibration.json` and
  `wt_dynamic_reference_candidate.json`: mandatory R1 resting confirmation;
- `regulatory_sustainment_sensitivity.csv` and
  `regulatory_sustainment_summary.json`: exact-root R1--R3 WT-only sustainment
  audit;
- `wt_model_generations.json`: rejected G2a, retained G2b scaffold, and G3
  dynamic-regulation decision.
- `native_source_panel.csv`, `native_source_root_attempts.csv`,
  `native_source_wt_roots.csv`, and `native_source_wt_summary.json`: canonical
  H1, diagnostic-only HW/HWQ, and the completed three-start production N/AN
  ledger.  Ten N-scale-4 roots remain pre-confirmation; these files are not a
  final pre-reveal freeze.
