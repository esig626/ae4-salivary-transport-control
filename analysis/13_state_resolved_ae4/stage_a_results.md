# Fixed-chassis Stage A results

## Controlling result

**Every conditionally transporter-compatible state-resolved bicarbonate family
fails before a predictive AE4-null secretion test on the fixed seven-state
chassis.** The failure is at the WT whole-cell gate, not evidence that AE4
biology is dispensable or that state-resolved transport is biologically
irrelevant.

Five model/gauge rows can conditionally reproduce the frozen transporter-level
summary after an explicit assay observation map:

- SR2 shared-state `1:2:3`, slow-common-transition gauge;
- SR5 explicit catalytic-cation `1:1:2`, unit and fast gauges;
- SR6 explicit U/P layers on SR5 with the loaded-transition PKA placement,
  unit and fast gauges.

These are **conditional assay-expressivity survivors**, not identified
microscopic mechanisms. Each uses two fitted loaded-state energy combinations,
one Na/K attempt-rate ratio, a common assay-output scale, the reported
`Rmin_Na=0.3` and `Rmin_K=0.4` as fixed observation-map nuisances, a discrete
common-barrier gauge, and—under SR6—one fitted regulated-layer barrier factor.
The dose-response context has source bath concentrations but unmeasured
intracellular assay concentrations fixed to a declared sensitivity tuple.
Raw curve errors and parameter covariances are unavailable, so the comparison
is not AIC and establishes no formal complexity winner.

A KO-blind 17-point log-spaced SR5/SR6-core barrier grid over `10^-4` through
`10^4` adds a sharper source-direction result. At every grid point that passes the
conditional transporter summary, the physiological baseline unit-capacity
source has `Cl_i < 0`, with the broad accepted-gauge limit approaching
`Na/Cl=-1.0205`, `K/Cl=+0.0205`, and `HCO3/Cl=-2`. The fixed historical WT row
instead demands `Cl_i=+0.0388833`, `Na_i=-0.0388833`, `K_i approximately 0`,
and `HCO3_i=-0.0777667`. The few tiny-gauge points that can flip the chloride
direction fail the transporter summary or are dominated by an unmeasured,
nearly cancelling Na/K slip. Thus the WT insertion failure is not an artifact
of choosing only the three headline gauges. See
`results/13_state_resolved_ae4/barrier_gauge_source_scan.csv`.

No carbonate row is counted as a biological transporter failure here. SR4
requires an alkalinity/pH observation map and an explicit carbonate/proton
whole-cell state; the fixed bicarbonate-only chassis cannot accept its source
vector without changing non-AE4 chemistry. SR3 and SR4b also retain the direct
Na-direction and K-only stand-alone conflicts documented in
`state_cycle_models.md`.

## Anti-leak freeze

Before calibration, code verified the curator-frozen files:

- `evidence_freeze.md`: `48022d5f69c4459edf6c9cceb914dd5ce092ffe73ffaffdbf9ceedc578005631`
- `question_tree.md`: `34343e076d46b00ce6eadc7a419688747b13917afce576611a2149c79aba019e`

Transporter and WT candidate tables were then hashed in
`results/13_state_resolved_ae4/stage_a_reveal_marker.json` before any AE4-null
target value was joined. No AE4-null ion, secretion magnitude, or secretion
time-course value entered a fit, a gauge selection, a root selection, or a
bound.

After the original reveal marker, an adversarial review found that the first
transporter summary projected BCECF onto HCO3 source alone for SR4. Amendment
`A01` corrected the transporter-only observation map to inorganic-carbon
alkalinity `HCO3 + 2 CO3` and added the finite barrier grid. The original
bytes and hash remain at `transporter_fit_summary_pre_reveal.csv`; the current
required summary and chronology are recorded in `amendment_01.json`. The
amendment read no knockout ionic or secretion value, retained exactly five
conditional survivors, and did not alter the WT/null Stage-A decision.

## WT capacity/root result

For each of the five conditional transporter rows, one chassis-native AE4
capacity was *attempted* against WT resting chloride only. None is accepted as
calibrated:

| Family/gauge | Attempted capacity | Derived height at joint attempt (um) | Raw closure error | Valid calibration? |
| --- | ---: | ---: | ---: | --- |
| SR2-123 / slow | `6.0e-6` | `199.999994` | `3.05e-5` | no |
| SR5 / unit | `3.94e-5` | `199.999991` | `3.05e-5` | no |
| SR5 / fast | `1.9e-5` | `199.999983` | `3.08e-5` | no |
| SR6 loaded-edge / unit | `5.1e-5` | `199.999986` | `3.06e-5` | no |
| SR6 loaded-edge / fast | `3.4e-5` | `199.999984` | `3.07e-5` | no |

The attempts approach the declared derived-height domain boundary and do not
close to the required `1e-8`; they are preserved as failed rows, not promoted
as fitted capacities. Freezing each attempted capacity and independently
re-solving the seven resting balances from 10 valid reduced-coordinate starts
finds one distinct bounded exact root per row. All five roots are nearly the
same:

- `Cl_i = 48.7820--48.7822 mM`;
- `pH_i = 7.37844--7.37846`;
- volume `2.9737--2.9740 pL`;
- maximum raw residual `5.94e-13` or smaller.

The chloride is within the deliberately broad WT chloride gate, but pH is
about `+6.69` reported WT SEM from `6.91 +/- 0.07`. Therefore every row fails
Gate 2. The root search covers the documented reduced domain corresponding to
positive states and `2 < height < 200 um`; it is not a global uniqueness
proof. An adversarial solve found an additional exact positive null root at
height about `1.25e14 um`, which is recorded as a remote, physiologically
inadmissible branch outside this domain.

## Stage-A AE4-null reveal

Deleting AE4 sets every candidate capacity to exactly zero, so the AE4-null
vector field is candidate-, state-graph-, PKA-, and barrier-gauge-independent.
A 20-start bounded solve gives one distinct local root:

| Quantity | Prediction | Stage-A target | Residual |
| --- | ---: | ---: | ---: |
| `Cl_i` | `48.7822 mM` | `36.50 +/- 1.60 mM` | `+7.676` reported SEM |
| pH | `7.37847` | `6.89 +/- 0.02` | `+24.424` reported SEM |
| volume | `2.974 pL` | no direct AE4-null target exists | diagnostic only |

The exact root residual is `5.77e-14`. This result is structurally immune to
the unidentifiable Na/K microscopic barrier gauge because no AE4 edge remains
after deletion.

## Secretion and timing disposition

No row passes the WT whole-cell gate, so none supplies an admissible frozen WT
denominator for a predictive AE4-null secretion ratio. The secretion magnitude
and digitized early/late curve remain strict holdouts; reporting a ratio from
the failed-capacity rows would falsely promote a chassis-invalid diagnostic.

The current adapter also solves transporter occupancy quasi-steadily on every
whole-cell RHS call. It has explicit state competition and edge-specific PKA
regulation, but no carrier-state memory. It therefore cannot attribute a
physical 2--3 minute delay to AE4 state relaxation. The immediate PKA step is
only a code-time sensitivity because no PKA activation time constant and no
historical code-time-to-minute map are source constrained.

## Stage-A decision

Stage A fails for all implemented, conditionally transporter-compatible
bicarbonate candidates on the unchanged fixed chassis. The reason is precise:
the fixed non-AE4 network cannot jointly close an admissible WT pH state and,
after AE4 deletion, its common bounded ionic root has the wrong chloride and
acid-base state. Per the frozen question tree this triggers Stage B
missing-balance localization; it does not justify stopping and does not support
the statement “AE4 itself is insufficient.”

Machine-readable evidence is in:

- `results/13_state_resolved_ae4/transporter_fit_summary.csv`
- `results/13_state_resolved_ae4/wt_capacity_attempts.csv`
- `results/13_state_resolved_ae4/stage_a_roots.csv`
- `results/13_state_resolved_ae4/stage_a_null_roots.csv`
- `results/13_state_resolved_ae4/model_scores.csv`
- `results/13_state_resolved_ae4/rejected_models.csv`
