# Agent I — independent numerical reproduction

## Scope and firewall state

This is a pre-reveal independent implementation audit. I read `AGENTS.md` and the full Task 13B prompt before inspecting the modern model. I did not open, parse, print, or hash `heldout_targets.csv`, and I did not use an AE4-null secretion magnitude or time course. No held-out numerical value appears in the implementation or tests.

The independent transcription is in `src/modern_full_model/independent.py`. It imports the immutable parameter container but does not import the production whole-cell assembly and does not call its `evaluate` or `rhs` methods. The following equations were transcribed separately:

- TIC/TA carbonate speciation with a finite monoprotic buffer;
- reversible NKCC1, NHE1, and AE2 affinities and bounded rates;
- the SR2-112 AE4 shared-carrier QSS cycle, including common-capacity regulation;
- the coupled apical/basolateral current equations, solved with an explicit two-by-two determinant rather than the production linear solver;
- pump and channel source stoichiometry and current-to-fmol/s conversion;
- osmotic water fluxes, local luminal outflow, CO2 exchange, and all 12 conserved-coordinate source balances;
- exact R0--R3 regulatory families, including R1 effective activation, R2 cAMP-plus-activation, and R3 PKA-plus-regulated-fraction dynamics with a common AE4 capacity map.

## Charge-manifold root method

The independent root routine does not optimize all 12 core states. It removes intracellular and luminal chloride amounts as dependent coordinates:

\[
n_{Cl,i}=n_{Na,i}+n_{K,i}-n_{TA,i}-X_i,
\qquad
n_{Cl,l}=n_{Na,l}+n_{K,l}-n_{TA,l}.
\]

It solves the ten independent source rows `(Na_i, K_i, TIC_i, TA_i, V_i, Na_l, K_l, TIC_l, TA_l, V_l)`. This prevents a conserved bulk-charge leaf from masquerading as an additional root degree of freedom. Residual rows are scaled by a declared 100 s source-to-state scale. Jacobian singular values are computed only after both rows and free-coordinate columns have been scaled, so fmol and pL coordinates are not compared as raw numbers.

## Equation-level pre-reveal results

At the executable, explicitly uncalibrated starting parameter set:

- all 12 independently transcribed RHS rows match the production assembly at the reference state to a maximum absolute difference of `9.03e-17` in their native source units;
- independently reconstructed cell pH agrees exactly at printed precision;
- the apical voltage differs by `4.17e-17 V`;
- independent current-closure residuals are below `2e-26 A`;
- charge/carbon source-accounting residuals are below `1.5e-16 fmol/s`, and water accounting is below `2.1e-27 pL/s`;
- the independent charge-manifold search does **not** find a WT root: maximum scaled residual `0.5839`, maximum amount/equivalent-source debug residual `0.1339` in its registered fmol/s row unit, and maximum volume-source residual `2.22e-5 pL/s`;
- the scaled ten-by-ten local Jacobian is full rank (`rank 10`, `nullity 0`) at that failed least-squares point.

This agreement rules out a simple production-assembly transcription error as the cause of the starting parameter set's root failure. It does not validate the starting parameter set and is not a claim about AE4 biology.

Agent F supplied a provisional WT-only candidate `G2_BALANCED_APICAL_K_BIASED:M00`. It was internally reproducible only while the executable defaults still contained placeholder hydraulic coefficients. The dimensional audit then corrected those values to the Palk-lineage physical conversion (`P_a=4.32e-3`, `P_b=5.15e-2`, `P_t=2.60e-4 pL s^-1 (mOsm/L)^-1`). Under the corrected coefficients, the stale state immediately has `dV_i/dt = 0.3373 pL/s` and `dV_l/dt = 0.03431 pL/s`. An independent 1% charge-manifold refinement fails with maximum scaled residual `377.5`. I reject that provisional root and preserve it only as a regression showing why water units must be frozen before calibration.

## Frozen corrected-water resting roots

The corrected-water calibration was frozen in `wt_selected_candidate.json` with SHA256 `41b91c09a69b6f50f1ea96ec49df9c56681bcc21ada8bc68b515f1b34a0eea44`. The selected branch is `G2_BALANCED_APICAL_K_BIASED_G10:M00`; its R0 state has cell chloride `50.1 mM`, pH `6.91`, cell volume `1.3 pL`, and apical/basolateral voltages `-62.2175/-72.6638 mV`.

I reconstructed parameters and states from the frozen CSV records and evaluated/refined every row marked `WT_REST_PASS_R0_CONTROL_ONLY` using only the independent charge-manifold equations. All nine retained roots reproduce:

- all nine have a successful independent refinement, scaled Jacobian rank `10`, and nullity `0`;
- maximum state discrepancy from the supplied frozen state is `5.07e-12` across the panel; every other root is within `8.53e-14`, as recorded at full precision in `independent_retained_roots.csv`;
- the largest independently evaluated native amount/equivalent source magnitude is `2.92e-14 fmol/s`, and the largest volume-source magnitude is `1.87e-15 pL/s` before refinement;
- no production whole-cell evaluator is called for these decisive root checks.

These results independently endorse the numerical existence and local full-rank character of the frozen corrected-water R0 roots. They do not by themselves establish a valid stimulated WT model.

## Frozen R1 dynamics and solver checks

The finalized R0--R3 production equations were transcribed without importing the production regulator. At deterministic off-equilibrium states, every independent regulatory RHS and capacity output agrees with the production API to numerical roundoff. The default R2 and symmetric R3 equations are independently confirmed to be dynamically isomorphic. An independently written R1 state remains positive and increases AE4 capacity after a beta-input step. Short independent Radau and BDF trajectories reproduce one another within the declared test tolerances while preserving charge, carbon, and water accounting. A scalar absolute tolerance is interpreted per state in its native unit; the final frozen crosscheck must also use a state-scale-derived vector tolerance.

The frozen dynamic manifest SHA256 is `6283c8aa4ee3449592a1228d82333ba8b06da280460a8fdd562c479e19f9804a`; its companion hash table SHA256 is `cba21c627e8ce91f6bb8d87ae4351d6835ebaec5f43225e4384feb9c18226729`. I reproduced preferred member `R1_G125_P21_PROSE_TREFERENCE` from the selected root with the independently transcribed whole-cell and R1 equations. The physical protocol is 600 s; CCh changes only calcium from `0.05` to `0.55 uM`, IPR changes only normalized beta input from `0` to `1`, the exact `t=0` coordinate remains basal, and the stimulated right limit is `t=1e-6 s`. R1 uses `tau=30 s` and changes the common AE4 capacity multiplier from `1` toward `1.25`.

With the frozen vector absolute tolerances, `rtol=1e-7`, and maximum step `2 s`:

- independent Radau gives cumulative flow `0.499430642730 pL`, endpoint chloride `20.0380117092 mM`, endpoint pH `6.81025670564`, endpoint flow `0.000636064675009 pL/s`, and endpoint AE4 capacity multiplier `1.24999999948`;
- independent BDF gives cumulative flow `0.499430611422 pL`, endpoint chloride `20.0380117252 mM`, endpoint pH `6.81025670612`, and endpoint flow `0.000636064675565 pL/s`;
- solver disagreement is `3.13e-8 pL` in cumulative flow, `1.91e-8` in the largest endpoint state coordinate, and `3.81e-7` on the declared normalized state comparison;
- independent amount/equivalent conservation residuals remain below `1.37e-15 fmol/s`, and water accounting below `2.17e-19 pL/s`.

Thus the numerical dynamic gate passes independently. The fixed WT minute-wise flow-shape gate does not: minute samples decline from `0.00110255` to `0.000636065 pL/s`, giving max/min ratio `1.73340`, above the frozen envelope ratio `1.11111`. The corresponding scalar-multiplier interval is empty (`lower 14149.5 > upper 9069.85 uL min^-1 per pL s^-1`). No gland scaling can repair that shape discrepancy. The inferred effective cell count (`1.902e8`) remains an unassessed warning because no certified primary geometry range is available; it is not promoted to an additional pass/fail gate.

## WT-only NKCC1 sustainment profiles

The next reconstruction generation adds one effective calcium-to-NKCC1 state without changing the reversible NKCC1 affinity or its `1 Na : 1 K : 2 Cl` stoichiometry:

\[
u_{Ca}=\operatorname{clip}\!\left(\frac{Ca-0.058}{0.55-0.058},0,1\right),
\quad
\dot n=\frac{u_{Ca}-n}{\tau_N},
\quad
M_N=1+(M_{N,\max}-1)n.
\]

I transcribed this N2 adjustment directly into the independent homeostasis equation and reproduced two frozen WT-only cases with the same R1 AE4 regulation:

- the source-centered `N2_M2_T120` case is numerically stable and solver-reproducible, but still fails the fixed WT shape gate: Radau cumulative flow `0.575729623588 pL`, minute max/min ratio `1.275461`, and endpoint Na/K/Cl `25.7255/94.5849/31.9556 mM`;
- the `N2_M6_T60` stress test crosses the numerical threshold with cumulative flow `0.671013481820 pL` and ratio `1.099181`, but requires an almost sixfold endpoint NKCC1 capacity and is explicitly not accepted as a physiologically credible model;
- Radau/BDF normalized state disagreements are below `4.68e-7` in both cases, charge/carbon residuals below `1.26e-15 fmol/s`, and water residuals below `2.17e-19 pL/s`.

This preserves the crucial distinction between a mathematical shape pass and a source-disciplined physiological reconstruction. Full-precision values and hashes are in `independent_n2_sustainment.json`.

## Current native-source hierarchy intake

The old corrected-water roots and joint G5 trajectories above are preserved
numerical reproductions, not current production candidates.  In particular,
the earlier `1.902e8`-cell scale warning is superseded by the binding
Kondo-derived one-SMG upper bound and is now a hard geometry failure.

The replacement native-source root ledger has completed H1, HW, HWQ, and the
three-start production N/AN absolute-source grid.  H1 is the canonical
production hydraulic mode.  Scale-one HW and HWQ are exact alias
controls of H1; they do not create duplicate candidate panels.  At nonunit H,
HW scales `Pa/Pb/Pt`, whereas HWQ also scales the local outflow rate.  Both are
explicitly diagnostic-only and cannot rescue a failed production model.

The completed counts are:

| Mode | Panels | Attempts | Numerical roots | WT-passing roots |
|---|---:|---:|---:|---:|
| canonical H1 | 45 | 135 | 35 | 0 |
| diagnostic HW | 225 | 675 | 175 | 0 |
| diagnostic HWQ | 225 | 675 | 175 | 0 |
| production `N_ABS_NKCC` | 315 | 945 | 229 | 10 preliminary |
| production `AN_ABS_AE4_NHE` | 315 | 945 | 125 | 0 |
| combined | 1,125 | 3,375 | 739 | 10 preliminary |

For each `H=2.5/3.5/5/6.75/8`, HWQ reproduces the corresponding HW 35-root
intracellular ranges to displayed precision: Na `8.107--29.441`, K
`70.840--92.515`, Cl `23.589--29.243 mM`, and pH `6.7522--6.9440` over the H
panel.  This is a numerical diagnostic, not an independently proved analytic
identity.  It confirms that adding proportional outflow scaling does not
alter the decisive steady-rest Cl/K failure in the tested panel.

Within the absolute grid, N numerical-root counts at scales
`0.5/2/4/5/6/7/8` are `40/35/33/31/30/30/30`; AN counts are
`45/25/20/10/10/10/5`.  All ten WT passes are N scale 4 and occur in two
contexts across all five topologies: pump-capacity scale 1 with `AE4NA05`, and
pump-capacity scale 2 with `AE4NA20`.  Their Cl and pH ranges are
`51.0172--53.0228 mM` and `6.83860--7.04827`, with minimum gate margins
`0.07723 mM` and `0.001728`; Na and K span `14.3306--18.9814 mM` and
`105.935--118.116 mM`.  Cell volume is `1.3 pL`, all normalized Jacobians have
rank 11/nullity 0, scaled residuals are `1.01e-12--2.41e-11`, and no boundary
hit is recorded.

These counts are an intake audit of the production reroot ledger; they are not
presented as 3,375 independently re-solved fits.  The ten passing rows are
three-start results only.  A hierarchy-code audit paused promotion before
confirmation and added panel-atomic persistence, both-quartile start coverage,
fail-closed ID/eligibility/upstream checks, exact scale-one aliases, and strict
33-start intake checks.  Those are engineering invariants, not scientific
retuning.  Independent reproduction remains pending until the rows survive
seven-start confirmation and explicit 33-start geometry.

## Open independent-reproduction items

1. Complete the seven-start confirmation and explicit 33-start geometry for
   every one of the ten preliminary N-scale-4 rows, with no first-passer
   pruning.
2. Independently reproduce every retained confirmed production root, then its
   WT dynamics, dual-solver/conservation checks, and one-SMG geometry gate.
3. Keep HW/HWQ excluded from production confirmation even if a future
   diagnostic rerun closes additional numerical roots.
4. Only after explicit lead authorization, run any licensed held-out
   comparison.

## Preserved disagreement

The present code-level agreement, nine old corrected-water resting roots, and
two-solver R1/N2 reproductions are strong. However, the independently computed
baseline and source-centered N2 trajectories reproduce decisive fixed WT
flow-shape failures; the only N2 numerical shape pass is a rejected high-gain
stress test.  The replacement native hierarchy has zero WT passes through H1,
HW, HWQ, and AN, while N supplies ten preliminary scale-4 passes that have not
yet undergone the required seven-/33-start checks or an independent reroot.
Therefore Agent I endorses the implementation's existing numerical
reproducibility but does **not** endorse these generations as a physiologically
validated G5 full model or as eligible for WT dynamics or held-out reveal.
This is not a final classification of AE4 biology or the full transporter
class.

## Final root-table reproduction

After the ten retained panels completed 33-start geometry, the independent
equation transcription parsed the exact final panel, attempt, and root byte
snapshots. It reproduced all 739 numerical roots, all ten production WT-gate
classifications, full rank/nullity-zero Jacobians for the production passers,
and reported zero gate or hierarchy disagreement. The native summary then
validated the commit ledger and marked the root freeze ready for dynamics.

The subsequent production WT dynamic contract had zero scientific group
passers because of the shared absolute one-SMG scale failure. Independent
genotype reproduction and holdout comparison were therefore correctly not
run; they cannot turn a failed WT denominator into a validated model.
