# Transporter-level constraint and calibration

## Result

The smallest state graph compatible with direct Na transport, direct K
transport, electroneutral reversal, and Na/K-specific mutant behavior is the
shared-state SR2 skeleton. The available experiments do not identify its
microscopic barriers, shared-pool competition, or stoichiometric member.

Under one declared low-Cl assay context, two distinct cooperative mechanisms
can express the reported Na/K dose summaries:

1. SR2 `1:2:3`, in which two cations are transported per cycle; and
2. SR5 on SR2 `1:1:2`, in which one cation is transported and a second
   same-reservoir cation occupies an explicit catalytic state but cancels from
   net stoichiometry.

The present dose-response data cannot discriminate those alternatives. SR6
adds an explicit U/P regulatory layer. A loaded-transition placement can match
the approximate recombinant forskolin fold in the unit and fast gauges; the
common-Cl-transition placement cannot do so in those same gauges. This is an
edge-placement discrimination within a declared kinetic gauge, not evidence
of direct S173 phosphorylation.

## Freeze and permitted targets

Calibration began only after exact verification of the frozen evidence and
question-tree hashes. The calibration API accepts only `CAL-T` and `CAL-WT`
IDs and raises on Stage-A/STRICT targets. Transporter rows use E16-03 through
E16-07, E21-03, and E25-04. No 2015 AE4-null value is available to these
functions.

## Observation map

The reported Fig. 9 quantity is BCECF alkalinization, not a selected
microscopic branch current. The fitted observable is therefore positive total
inorganic-carbon alkalinity source

\[
S_{Alk}=S_{HCO_3}+2S_{CO_3},
\]

under the reverse low-external-Cl assay. A wrong-sign net source is retained
as zero; it is never made positive with `abs(J_branch)`. For bicarbonate
families this reduces to total HCO3 source summed over every Na/K branch. For
SR4 it is only a fast-equilibrium alkalinity projection; it does not turn
carbonate into a bicarbonate chassis state.

E16-07 did not include a zero-cation point and did not report fit covariance.
The reported `Rmin_Na=0.3` and `Rmin_K=0.4` are therefore kept as fixed assay
observation-map nuisances, separate from AE4 current. One common assay scale
maps model amplitudes to `Rmax-Rmin=1.2` for both cations. This exactly records
the full reported `Rmin/Rmax` values without inventing a cation-independent
AE4 transport cycle.

The external dose grid and HCO3/Cl baths are source values. Intracellular
concentrations were not measured in this experiment; the diagnostic fixes
`Na_i/K_i/Cl_i/HCO3_i = 25/120/50/10 mM`. Consequently the fitted loaded-state
energy combinations are context-conditional lumps, not microscopic Kd.

## Conditional fit results

The hard numerical gate requires both transported cations, exact branch
electroneutrality, local detailed balance, the encoded double-mutant
Na-collapse/K-retention limit, EC50 within 30% of `49/62 mM`, Hill slopes
within 0.4 of `2.0/1.8`, common-scale `Rmax` errors no larger than `0.10`, and
the PKA fold when a regulatory layer is claimed.

| Model/gauge | Na EC50 / Hill | K EC50 / Hill | Mapped Rmax Na/K | PKA fold | Conditional result |
| --- | --- | --- | --- | ---: | --- |
| SR2-123 / slow | `49.11 / 2.11` | `59.47 / 2.13` | `1.50 / 1.60` | n/a | pass |
| SR5-112 / unit | `49.00 / 2.04` | `60.82 / 2.05` | `1.50 / 1.60` | core only | pass |
| SR5-112 / fast | `49.00 / 2.06` | `61.28 / 2.06` | `1.50 / 1.60` | core only | pass |
| SR6 loaded-edge / unit | `49.00 / 2.04` | `60.82 / 2.05` | `1.50 / 1.60` | `1.2500` | pass |
| SR6 loaded-edge / fast | `49.00 / 2.06` | `61.28 / 2.06` | `1.50 / 1.60` | `1.2500` | pass |

“Pass” means conditional expressivity under the observation map and gauge. It
does not identify a unique transporter. The double-mutant hierarchy is a
constructed branch-limit test (`Na` near assay background, K retained), not an
out-of-sample numerical prediction of a mutation rate.

## Rejected and unresolved families

- SR1 pure-Na and pure-K are assay submodels, not a one-protein mixed-cation
  mechanism.
- SR2-112 yields apparent slopes near one in this implementation and misses
  the cooperative summary.
- SR2-213 also yields near-unit cation slopes in the declared context.
- SR3 and SR4b are thermodynamically/electrically admissible hypotheses but
  conflict as stand-alone models with the direct Na direction and K-only
  operation.
- SR4a, after the explicit alkalinity correction, can place apparent EC50 near
  `49/62 mM` but yields Hill slopes about `1.11/1.13`; it remains a discrete
  carbon-species alternative and cannot enter the bicarbonate-only fixed
  chassis.
- Several rejected gauge fits end at extreme loaded-energy lumps or at the
  PKA `gamma<=10^6` diagnostic bound. They remain in the CSV and are not
  silently discarded.

## Complexity and nonidentifiability

Rows enumerate scalar dimensions rather than counting prose labels: two
loaded-state energy combinations, one Na/K attempt ratio, one assay scale, and
one additional PKA barrier factor for SR6. Two Rmin values are fixed nuisance
observations; the common-barrier gauge is a discrete unmeasured choice. Raw
curve data, errors, fit covariance, intracellular conditions, mixed-cation
competition, and relaxation time courses are unavailable. Therefore no AICc
or formal likelihood winner is reported.

A 17-point log-spaced common-barrier grid over eight decades shows that the
mixed-bath physiological source direction is gauge sensitive outside the
transporter-pass region, but every passing SR5/SR6-core grid point has negative
chloride source at the fixed baseline. The broad accepted-gauge limit is close to
`Na/Cl=-1.0205`, `K/Cl=+0.0205`, `HCO3/Cl=-2`. This opposite sign relative to
the fixed WT chloride-loading demand is why no arbitrary unit-rate gauge is
allowed to drive the whole-cell conclusion.

The chronological carbonate correction is recorded in
`results/13_state_resolved_ae4/amendment_01.json`; the original pre-reveal
summary is byte-preserved, and the amendment used no knockout evidence.
