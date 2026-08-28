# Task 13B model-generation log

This log is append-only at the generation level. A retained generation must state its parent, defect, equation/topology change, new parameters and provenance, calibration/validation data, gate movement, numerical status, and retain/revert decision. The strict AE4-null secretion magnitude and time course are excluded from every calibration and model-selection objective.

## G0 — Task 13 inherited-chassis reproduction

- **Parent:** Task 13 frozen state-resolved package.
- **Defect addressed:** establish that Task 13 decisive checks reproduce before extension.
- **Equations/topology:** unchanged seven-state historical-lineage comparison chassis; SR2/SR5/SR6 transporter registry; nested M1 cation split.
- **New states/parameters:** none.
- **Independent evidence:** existing frozen Task 13 evidence ledger and artifacts.
- **Calibration data:** none newly fitted.
- **Validation data:** frozen hashes, invariant tests, independent roots/ranks, and dual-solver diagnostic only.
- **Gate movement:** Gate A/C/F controls reproduce; WT Gate B/E remains failed on this chassis.
- **Numerical status:** 72 Task 13 state-resolved tests passed on 2026-08-27 UTC.
- **Decision:** retain as negative control only; do not patch or promote the seven-state chassis.

## G1 — conserved acid-base/carbon reconstruction

- **Parent:** G0 negative control.
- **Defect addressed:** the inherited single-HCO3-row closure did not conserve
  finite carbon, alkalinity, buffer charge, and water-facing osmoles.
- **Equations/topology:** explicit cell/lumen total inorganic carbon and total
  alkalinity; finite buffer; neutral CO2 exchange; NHE1, AE2, and AE4 sources
  enter TIC/TA with charge-correct signatures.
- **New states/parameters:** cell/lumen TIC and TA amounts, finite buffer pools,
  CO2 exchange coefficients, and other-impermeant osmoles.
- **Evidence/provenance:** conservation and charge are derived constraints;
  numerical buffer/exchange values remain declared modeling decisions.
- **Calibration/validation:** WT Cl/pH only; no knockout secretion input.
- **Gate movement:** carbon, buffer, and charge identities pass analytically and
  in tests; physiological rate values remain uncertified.
- **Decision:** retain the conservation architecture.

## G2a — finite lumen and modern cation topology with placeholder hydraulics

- **Parent:** G1.
- **Defect addressed:** missing finite luminal ions, membrane-specific current
  closure, and nonzero apical pump/K topology.
- **Equations/topology:** finite luminal Na/K/Cl/TIC/TA/volume; two algebraic
  membrane potentials; apical/basolateral pump and K partitions; charge-complete
  paracellular paths and water/outflow.
- **New parameters:** pump/K routing fractions and placeholder water
  coefficients.
- **Evidence/provenance:** Almássy topology supports nonzero apical pathways but
  not the provisional fractions.
- **Gate movement:** algebraic charge/current and positivity tests improved;
  the provisional root failed after the published hydraulic-unit conversion.
- **Decision:** reject G2a; preserve only as a water-unit regression.

## G2b — corrected-water WT resting-root ensemble

- **Parent:** G2a.
- **Defect addressed:** replace placeholder hydraulics with the exact Palk/2018
  published-model conversion.
- **Equations/topology:** unchanged G2a sources; `P_a=0.00432`, `P_b=0.0515`,
  `P_t=0.000260 pL s^-1 (mOsm L^-1)^-1`.
- **New fitted quantities:** AE4 effective capacity gauge, NHE1 capacity, and
  OTHER impermeant osmoles; the old common membrane conductance multiplier was
  a profiled nuisance axis.
- **Calibration/validation:** WT resting Cl/pH/volume construction and
  out-of-objective Na/K/voltage gates only.
- **Gate movement:** 9/45 profiled slices produced interior charge-manifold
  roots; all nine occurred at the old common conductance scale `10`.
- **Decision:** retain the 12-state/root machinery, but supersede every G2b
  root as a production root after the native conductance source-map audit.

## G3 — mandatory dynamic beta/cAMP/PKA regulation

- **Parent:** G2b scaffold.
- **Defect addressed:** static AE4 activation did not represent the mandatory
  beta -> cAMP/PKA -> AE4 regulatory chain.
- **Equations/topology:** R0 static control; R1 one effective activation state;
  R2 cAMP plus activation; R3 PKA plus regulated fraction. Capacity modulation
  composes outside the thermodynamic AE4 cycle.
- **New parameters:** source-unmeasured gains `1.25/1.70` and predeclared kinetic
  sensitivities; S173A removes the increment and S273A retains it.
- **Evidence/provenance:** 2021 direction/hierarchy is measured; gains are
  context sensitivities and acute time constants are modeling assumptions.
- **Gate movement:** positivity, beta-free recovery, construct hierarchy, and
  R2/R3 exact isomorphism pass. R0 remains a negative control; R4
  transition-specific regulation is rejected as unsupported.
- **Decision:** retain R1 as the smallest production-eligible dynamic family
  and R2/R3 as frozen structural sensitivities; do not claim kinetic
  identification or direct S173 phosphorylation.

## G4 — state-resolved AE4 integrated with conserved whole-cell balances

- **Parent:** G2b + G3.
- **Defect addressed:** embed direct Na/K, electroneutral, reversible AE4
  transport rather than a pooled scalar event law.
- **Equations/topology:** SR2 shared-carrier Na/K branches in QSS, with local
  detailed balance and common-capacity regulatory coupling.
- **New parameters:** carrier-amount/attempt-rate gauge and assumed native
  Na-routing fractions `0.05/0.10/0.20`.
- **Evidence/provenance:** 2016 transport and 2025 mutation asymmetry constrain
  the module; routing fractions and QSS rate gauge remain assumptions.
- **Gate movement:** transporter thermodynamics, mutation logic, conservation,
  and reversal tests pass; absolute carrier kinetics remain unmeasured.
- **Decision:** retain the source-consistent module; do not identify a unique
  microscopic cycle.

## G5a — first WT 600-second stimulated ensemble

- **Parent:** G4.
- **Defect addressed:** establish a physical-seconds WT trajectory with distinct
  CCh/calcium and IPR/beta inputs before any genotype test.
- **Equations/topology:** old fixed `0.058 -> 0.55 uM` Ca input and frozen
  R0--R3 ensemble.
- **Calibration/validation:** WT-only graphical `9--10 uL/min` envelope;
  knockout saliva remained sealed.
- **Gate movement:** solvers/conservation passed, but `q_600/q_60` shape failed
  (`max/min=1.7334`, allowed at most `1.1111`).
- **Decision:** reject as G5; retain failure.

## G5b/G5c — targeted regulatory and NKCC sustainment attempts

- **Parent:** G5a.
- **Defect addressed:** test whether supported AE4 regulation or a CCh/Ca-linked
  NKCC source can repair WT sustainment.
- **Equations/topology:** complete R1--R3 gain/timing sweep; then N0 static,
  N1 algebraic, and N2 one-state NKCC-capacity families preserving NKCC
  stoichiometry and reversal.
- **New parameters:** declared WT-only NKCC gains/times; no genotype-specific
  compensation.
- **Gate movement:** 0/18 regulatory-only members repaired the shape. Nominal
  NKCC profiles failed; only deliberately out-of-scale slow stress probes
  passed and were production-ineligible.
- **Decision:** reject the regulatory-only and stress-test repairs; continue.

## G5d — superseded joint effective-Ca/NKCC ensemble

- **Parent:** G5c.
- **Defect addressed:** profile WT effective Ca with algebraic N1 capacity
  without using AE4-null data.
- **Equations/topology:** old roots plus `CA010/M1.75`, `CA012/M2.0`, and failed
  `CA015/M1.75` profiles across 20 regulatory members and two solvers.
- **Numerical result:** 18 root/profile groups passed the old within-profile
  flow-shape and shared-scale tests; AE2 deletion changed modeled flow by only
  about `+0.060%` on the reference profile.
- **Subsequent source audit:** the ensemble required approximately
  `2.485e8--3.401e8` cells (`323--442 uL` of 1.3-pL modeled cells), exceeding
  the generous Kondo one-SMG acinar upper bound of `51,469,231` cells. Its
  common conductance multiplier also did not implement the literal Palk-lineage
  `31.4 nS` CaCC and `14 nS` total-K maxima with the Palk gate.
- **Decision:** `REJECTED_SUPERSEDED` with reason codes
  `ABSOLUTE_GLAND_GEOMETRY_SCALE_FAIL` and
  `NATIVE_CONDUCTANCE_SOURCE_MAP_FAIL`. The immutable old manifest/gate retain
  their original bytes and hashes; `g5_supersession.json` is controlling.
  Held-out reveal is prohibited.

## G5N0 — strict native-source reroot

- **Parent:** conserved G1/G4 equations, not the rejected G5d roots.
- **Defect addressed:** replace the old conductance/common-scale and gland-scale
  assumptions before any new stimulated calibration.
- **Equations/topology:** literal `g_CaCC,max=31.4 nS`, total
  `g_K,max=14 nS`, Palk gate `K_Ca=0.26 uM`, Hill `1.46`, zero K-background
  regularizers, total-K partition only, pump fraction
  `0.056878/0.075075/0.094586`, K routing `0.20/0.30/0.40`, and AE4 routing
  `0.05/0.10/0.20`.
- **Calibration/validation:** deterministic charge-manifold native reroot using
  WT rest only; primary dynamic profile is `C_eff=0.10 uM` with the frozen
  reference R3/gain sensitivity. `C_eff` is in uM and is an assumed spatially
  averaged WT input, not measured free Ca and not a knockout fit.
- **Gate rule:** strict hydraulics `H=1` is production. `H>1` may be rerooted
  only as a named diagnostic and cannot rescue strict failure.
- **Mode/identity controls:** `H1` is the canonical production spelling.
  Scale-one `HW` and `HWQ` are exact alias controls of H1 and are not counted
  as new diagnostic panels; every nonunit HW/HWQ member is explicitly marked
  diagnostic-only. HW scales `Pa/Pb/Pt`; HWQ scales the same three hydraulic
  coefficients and the local outflow rate together.
- **Completed H1 result:** the 45 H1 panels (15 native-pump plus
  30 conditional pump-capacity `0.5/2` sensitivities) produced 35 numerical
  roots and **zero WT-passing roots**.  Every numerical root failed the
  independent WT Cl and K gates; some also failed Na and/or pH.
- **Completed hydraulic diagnostics:** HW and HWQ each contain 225 panels
  (`675` three-start fits), each produce 175 numerical roots, and each produce
  **zero WT-passing roots**.  Across H1+HW+HWQ the completed ledger therefore
  contains 495 panels, 1,485 attempts, 385 numerical roots, and zero passing
  panels or roots.  At each nonunit H (`2.5/3.5/5/6.75/8`), HWQ reproduces the
  same 35-root intracellular ranges as HW to displayed precision: Na
  `8.107--29.441`, K `70.840--92.515`, Cl `23.589--29.243 mM`, and pH
  `6.7522--6.9440` over the H panel.  This is a numerical diagnostic, not a
  proved analytic identity; it shows that coupled outflow scaling does not
  repair the decisive steady-rest Cl/K defect.
- **Decision:** strict native is not promoted, and neither hydraulic mode can
  rescue it because both are diagnostic-only and fail the WT gates anyway.
  The production N/AN absolute-source hierarchy was therefore run as the next
  predeclared tier and is reported below; no reveal authorization can be
  inferred from this log.

## G5N/G5AN — absolute source-class hierarchy after strict failure

- **Parent:** G5N0, only if strict native geometry/sustainment fails.
- **Defect addressed:** test the smallest physically interpretable WT loading
  alternatives without global area scaling.
- **Equations/topology:** `N_abs` scales whole-cell NKCC capacity while
  preserving stoichiometry/reversal; `AN_abs` scales the native AE4+NHE cycle
  together while preserving their ratio and thermodynamics. AE4 regulation
  composes outside `AN_abs`.
- **New parameters:** production scalars `{0.5,1,2,4,5,6,7,8}`; larger
  `{12,16,24,32}` values are diagnostic stress tests only. Scalar one is an
  exact nest. Every passing row is retained; no first-passer selection.
- **Calibration/validation:** WT source-class calibration only; never AE4-null
  magnitude or timing.
- **Completed three-start scope:** the production N/AN grid contains 630 panels
  and 1,890 deterministic attempts across pump contexts and nonunit scalars
  `0.5/2/4/5/6/7/8`.  Of those attempts, 1,847 reported optimizer success, all
  1,890 remained admissible, and 981 met the numerical-root criterion.  The
  latter collapse to one numerical cluster in each of 354 panels; 276 panels
  are rootless.  This is the preliminary screen, not the required seven-start
  confirmation or 33-start root-geometry result.
- **Numerical landscape:** `N_ABS_NKCC` produces
  `40/35/33/31/30/30/30` numerical roots at scales
  `0.5/2/4/5/6/7/8`, respectively (229 total).  `AN_ABS_AE4_NHE` produces
  `45/25/20/10/10/10/5` (125 total).  AN has zero WT-passing roots.  N has ten
  preliminary WT-passing roots, all and only at scale 4.
- **Preliminary N-scale-4 passes:** all five declared topologies
  (`PHIGH_KHIGH`, `PHIGH_KLOW`, `PLOW_KHIGH`, `PLOW_KLOW`, `PNOM_KMID`) pass in
  each of exactly two contexts: pump-capacity scale 1 with `AE4NA05`, and
  pump-capacity scale 2 with `AE4NA20`.  No first-passer pruning was applied.
  Across these ten roots, Cl is `51.0172--53.0228 mM` (minimum gate margin
  `0.07723 mM`), pH is `6.83860--7.04827` (minimum margin `0.001728`), Na is
  `14.3306--18.9814 mM`, K is `105.935--118.116 mM`, and cell volume is
  `1.3 pL`.  Apical voltage ranges from `-36.1641` to `-34.5019 mV`, and
  basolateral voltage from `-73.5186` to `-66.5896 mV`; fitted OTHER spans
  `110.024--126.008 fmol`.  Maximum scaled residuals are
  `1.01e-12--2.41e-11`, every normalized Jacobian has rank 11/nullity 0, and
  no boundary hit is recorded.
- **Combined checkpoint:** H1, HW, HWQ, N, and AN now total 1,125 panels,
  3,375 attempts, and 739 numerical roots, with the same ten preliminary WT
  passes.  HW/HWQ remain diagnostic-only and contribute none of them.
- **Engineering gate:** the hierarchy was paused before confirmation while a
  code audit hardened panel-atomic rerun replacement, preservation of both
  quartile starts in the three-/seven-start designs, fail-closed duplicate-ID
  and eligibility handling, exact scale-one alias coverage, complete upstream
  grid checks, and 33-start intake validation.  These are persistence and
  execution-integrity fixes, not changes to equations, source scalars, WT gate
  thresholds, parameter bounds, or biological calibration.
- **Decision:** retain all ten preliminary N-scale-4 rows for the mandatory
  seven-start confirmation, followed by explicit 33-start geometry if they
  survive.  None is dynamics-eligible or a frozen production model yet.
  Nonunit values remain explicit context hypotheses, not identified
  primary-source capacities; AN remains a completed negative tier.  Held-out
  reveal remains prohibited.

## G5V — beta/volume-sensitive anion-exit diagnostic

- **Parent:** native-source architecture.
- **Defect addressed:** Catalán 2015 requires a beta-associated,
  volume-sensitive anion-exit topology distinct from CaCC and AE4 loading.
- **Equations/topology:** V0 exact control and diagnostic
  `V1 = g_V beta max(V_i/V_rest-1,0)`; no new kinetic state and no epsilon gate.
- **Numerical diagnostic:** on the existing exact rest root, IPR-only R1 with
  `Ca=0.058 uM` shrank `V_i` from `1.3000` to `1.273730 pL`; maximum swelling
  was zero, so V1 deadlocked identically for every `g_V`.
- **Evidence/provenance:** the source establishes the pathway/topology but does
  not supply a valid whole-cell `g_V` map. Even an exogenously imposed `+12.5%`
  swelling gate would leave the capacity unlicensed.
- **Decision:** reject V1 as a rescue and retain V0 only as an exact control.
  The missing beta/volume-sensitive absolute flux remains a substantive
  full-model limitation.

## G5-native final reconciliation (2026-08-28)

- The ten preliminary `N_ABS_NKCC` scale-4 roots all survived the declared
  seven-start confirmation and 33-start geometry. Each retained one branch;
  the final production lineage is complete and independently hash-matched.
- The native dynamic contract propagated all ten roots through 20 frozen
  R0--R3 members, two stimulus arms, and two stiff solvers (800 production
  cases), plus 320 nearby-state cases. Numerical, conservation, solver,
  sustainment, regulatory-direction, and co-stimulation gates passed.
- All ten failed the one-SMG absolute WT scale. Required observation scales
  were `6362.31--8210.72`, versus the frozen ceiling `3088.15386`.
- G6 targeted repair used only the predeclared high-capacity N and AN stress
  directions (`12/16/24/32`). It produced zero additional WT-rest passers and
  was rejected. No holdout datum was accessed.
- Final decision: retain the equation, root, and negative diagnostic records;
  reject every model as a validated full model; close as Task 13B Outcome 4.
