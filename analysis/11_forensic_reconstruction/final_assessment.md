# Final forensic assessment

## Decision

| Question | Decision |
| --- | --- |
| Identity of the archived MATLAB snapshot as the implementation that generated the 2018 article | **`UNLIKELY`** |
| Relationship to the published project | Directly related calibration/reduction branch; shared ancestry is strongly supported |
| Conditional clean-reconstruction gate | **Closed**: identity is below `VERY LIKELY` and the central phenotype does not reproduce |
| Trustworthy physiological parameter-to-observation map | **No** |
| Identifiability/discrimination gate | **Not reopened** |
| Project classification | **`STOP`** remains in force |

`UNLIKELY` is a grade of exact publication-generating identity, not a claim that
the files are unrelated. Geometry and several fitted parameter transformations
show unusually strong lineage. The decisive contrary evidence is equation
non-identity, a different calcium/time protocol, missing and stale dependencies,
and failure of both available dynamical variants to reproduce the published
AE4-knockout phenotype. The evidence therefore does not meet the prompt's
`VERY LIKELY` threshold.

No clean 2018 implementation was created, no physiological identifiability
calculation was resumed, and no manuscript prose was drafted.

## Identity evidence matrix

| Dimension | Evidence of common lineage | Evidence against exact 2018 identity | Assessment |
| --- | --- | --- | --- |
| Equation identity | Same compartments and transporter families; conservation-consistent volume sign; H, luminal Cl, and voltage can be eliminated algebraically under the code's charge/QSS assumptions | Every archived variant makes AE4 Na-specific, whereas the paper uses Na+K and partitions the cycle between both cation balances; CO2 is forced to QSS instead of remaining dynamic; NHE1 proton powers and channel gates differ | Strong ancestry; strong non-identity |
| Parameter identity | The recovered geometry gives exactly 1.3 pL; height/area transformations recover published tight-junction, AE2, NHE1, and AE4 activities to about 0.3% or better; corrected NKCC coefficients agree | Resting HCO3, bath HCO3/CO2, channel gates, AE4 reverse coefficient, water coefficients, buffer scaling, and several family-specific powers of ten differ | Strong calibration lineage; not one transparent parameter set |
| Baseline | The Original calibration closes its 12 displayed residuals to `1.39e-17`; the seven-state WT raw RHS is below `5.24e-12`; many displayed resting values are reproduced | Closure is constructed from the same target row and fitted conductances; HCO3 is 10 rather than 12.1 mM; absolute flow units are undocumented; the paper's literal tables still fail their own closures | Historical baseline reproduced, published baseline not independently reproduced |
| Phenotype/output | Translated equations reproduce the qualitative negligible AE2 effect | Seven-state AE4 KO gives endpoint flow ratio `0.998556` and first post-step ratio `1.101877`; five-state AE4 KO gives endpoint ratio `1.059019` and is not settled. None approaches the published approximately 24% reduction | Decisive evidence against exact generating identity |
| Figure correspondence | Historical continuation figures reproduce the later manuscript's claimed roughly 25% reduction as activity tends to zero | The archived drivers do not save those figures; their three-state/MATCONT source is absent. The 20 analysis/sensitivity figures do not duplicate the ten published figures, whose calcium and time-course shapes differ | Same scientific question, different figure-generating workflows |
| Chronology | MATLAB EPS exports and the Golubitsky manuscript are from 2017, during the online-first publication interval | `Par.mat` was serialized on 2018-10-10; standalone manuscript panels are from 2018; chronology cannot make this snapshot a preserved prepublication release | Mixed snapshot, not a frozen release |
| Dependency structure | `Par.m` exactly serializes all 45 `Par.mat` fields; the Original calibration reproduces the outer parameters | Original driver calls a missing four-argument `Salivary`; the archived function requires six arguments; `Conductances.mat` is absent; the five-state call also has arity drift; the three-state continuation code is absent | Research working directory with missing/stale components |

The detailed evidence is in
[`code_lineage.md`](code_lineage.md),
[`paper_code_concordance.md`](paper_code_concordance.md),
[`reproduction.md`](reproduction.md), and
[`figure_and_manuscript_provenance.md`](figure_and_manuscript_provenance.md).

## Direct answers

### 1. What are these historical MATLAB files?

They are three connected layers of a research branch:

1. `Original/Parameters.m` and the first part of
   `Original/Saliva_Ae4.m` form a sequential baseline-calibration worksheet.
   Its 12 entries are a residual audit, not a 12-state simulator.
2. `Par.m` and `Par.mat` hand the calibrated values to later functions.
   `Par.mat` contains one `1 x 1` struct with 45 scalar-double fields, exactly
   matching `Par.m`; it contains no trajectories, solver state, or hidden
   knockout results.
3. `Salivary.m`, `Salivary2.m`, and `Salivary_ex.m` are a seven-state
   algebraically reduced engine, a five-state fixed-lumen fork, and an
   exploratory driver. They are closely related to the unpublished
   Golubitsky mathematical-analysis project but do not implement its final
   three-state continuation system.

This conclusion follows from equations, parameter transformations, state
dimensions, and function interfaces—not from the directory name `Original`.

### 2. Which files are most likely to have generated the 2018 results?

No available file can be attributed as the generator of the published results.
The `Original/` pair is the strongest candidate for an **ancestral calibration
workflow**, but it has no complete compatible solver in the archive and already
uses Na-only AE4. The outer seven- and five-state functions are later or parallel
reductions with materially different dynamics. The archived continuation
figures came from an absent three-state/MATCONT workflow, not from the included
driver.

Thus the exact generator is missing from the repository evidence. Calling any
included variant the published implementation would require inferring across
both missing dependencies and substantive model changes.

### 3. What exact evidence supports that conclusion?

- The Original calibration reproduces the fitted outer parameters with maximum
  relative difference `2.51e-15`, and every one of the 45 binary fields exactly
  matches `Par.m`. This proves direct lineage.
- `Salivary.m` integrates seven states and eliminates H, CO2, luminal Cl, and
  both membrane potentials; the paper reports ten differential states plus two
  quasi-steady potentials.
- Archived AE4 uses intracellular/extracellular Na only. It places the complete
  exchanger flux in the Na balance and none in K. Restoring the paper's Na/K
  partition at the historical resting row creates equal-and-opposite residuals
  of `0.0321793147` in historical flux units.
- The archived calcium input is one step from 0.05 to 0.55 at `t>100`, with
  hidden RHS multipliers of `10000` and `1000`. The published input is a
  multistage minute-scale curve.
- The seven-state AE4-KO endpoint flow change is only `-0.144%`; the first
  post-step sample is `+10.19%` relative to WT. The five-state endpoint is
  `+5.90%` and still moving. The published result is approximately `-24%`.
- The archived drivers have stale arities, a missing `Conductances.mat`, and no
  code for the three-state continuation figures.
- The 20 EPS analysis plots are 2017 MATLAB exports used by the Golubitsky
  manuscript; they are determinant, continuation, and parameter-sensitivity
  plots, not duplicates of the published figures. `Par.mat` itself is dated
  2018-10-10.

### 4. Which published inconsistencies are resolved by the implementation?

The historical evidence resolves or strongly corroborates limited conventions:

- `d omega_i/dt = q_b-q_a` is the conservation-consistent volume sign.
- The Palk-corrigendum NKCC coefficients must be converted for mM numerical
  states; the archive contains exactly that conversion and gives turnover
  `0.5980945 s^-1` at the printed ion state.
- Signed channel/tight-junction terms are consistently `I/F` particle-flux
  equivalents, with NaK cycle flux multiplied by `F` in the voltage equations.
- The recovered cell height/area geometry explains several fitted
  activity/conductance transformations and gives exactly 1.3 pL.
- `ode15s` tolerances are `RelTol=AbsTol=1e-6` for the historical runs.
- Consistent use of `K_B=1e-4` in both code layers supports a missing minus sign
  in the paper's printed AE2 value.

The code also removes its own 59.6 water mismatch and 67.8 pA electrical
residual, but only by fitting different water coefficients, gates, and
conductances to the selected row. Those are implementation-specific
calibrations, not demonstrated corrections to the published parameter tables.

### 5. Which inconsistencies remain real?

- The published AE4 Na/K mechanism and absolute density/activity units are not
  recovered by the Na-only historical law.
- The historical code's own chemical-equilibrium CO2 value is `2.907908` mM,
  not its algebraically forced `6.600124` mM; the aggregate CO2 sign remains
  conservation-inconsistent with the two individually printed influxes.
- Absolute flow, current, transporter-density, CO2-permeability, and physical
  time conversions do not share a documented unit ledger.
- The published calcium input, complete runnable dependency set, and
  publication-generating numeric trajectories are absent.
- Historical CO2 QSS, first-power NHE1 law, altered channel gates, bath
  compositions, and family-specific scale factors are substantive differences
  rather than corrections. Eliminating H by differentiating electroneutrality
  is algebraically defensible under the code's other balances and QSS charge
  closure, but it does not rescue the other transient and parameter differences.

Every Task 10 blocker and its required disposition is tabulated in
[`blocker_resolution.md`](blocker_resolution.md).

### 6. Does the implementation reproduce the published baseline and AE2/AE4 phenotype?

It reproduces **its own fitted historical baseline**, including many displayed
published values, but not the published specification as an independent
parameter set. It qualitatively reproduces the negligible AE2 effect. It does
**not** reproduce the central AE4-knockout flow reduction in either available
dynamical variant. Literal MATLAB/Octave execution was unavailable, so the
trajectory evidence is a direct equation transcription integrated with SciPy
BDF; the exact baseline/calibration conclusions do not depend on that solver.

### 7. Is there now a trustworthy physiological parameter-to-observation map?

No. A trustworthy map would require a common equation set, physical unit
ledger, reproduced baseline, correct stimulation protocol, and reproduced WT
and knockout outputs. The archive provides none of those as one validated
executable system. Its fitted resting point cannot certify derivatives or
observation geometry for the published model.

### 8. Does `STOP` change?

No. The classification remains **`STOP`**. The clean-reconstruction and
identifiability conditions are explicitly false, so `QUICK PAPER` and
`PUSH HARDER` are not available outcomes on this evidence.

### 9. What is the strongest mathematical result now available?

The strongest independently checked result remains the Task 10 exact
steady-balance coordinate inversion

\[
J_4=B-C+2N,\qquad J_2=2C-4N-B,
\]

equivalently `J2=2A-B`, `J4=B-A` for `A=C-2N`. It gives an exact conditional
decomposition of AE2/AE4 cycle fluxes and regular full-state forcing geometry.
The forensic work strengthens the negative gate: the historical code cannot
turn that structural identity into a certified physiological activity map.
The unpublished three-state continuation/sensitivity results are provenance
evidence and are not adopted as new results.

### 10. What are the next three steps?

1. Obtain the actual publication-generating source, parameter file, calcium
   time course, and unit/scaling ledger from the original authors, or obtain an
   explicit correction that uniquely specifies them.
2. Execute that source in a compatible MATLAB environment and require separate
   regression against the published resting row, WT trace, AE2 knockout, and
   approximately 24% AE4-knockout flow phenotype without retuning.
3. Only after those gates pass, build an independent clean implementation with
   convention-by-convention provenance and then recompute physiological
   implicit sensitivities, panels, perturbation geometry, and global overlap.

## Archive integrity and conditional-work record

The file-level audit accounts for all 67 archive items. Sixty-five locally
materialized files retain their authoritative SHA-256 values; the two connector
transport exceptions retain their authoritative remote Git blob identities and
sizes. No path under `archive/` was edited. See
[`archive_audit.md`](archive_audit.md) and
`results/11_forensic_reconstruction/archive_hashes.csv`.

Because exact identity is graded `UNLIKELY` and the central phenotype failed,
the prompt's conditional work was intentionally not performed:

- no `src/reconstructed_2018/` namespace;
- no clean-model regression suite;
- no physiological `F_u`, `F_theta`, panel ranking, global-equivalence,
  singularity, or mechanism-class computation; and
- no manuscript draft.
