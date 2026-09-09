# Regulatory sustainment sensitivity on the corrected-water WT staging root

## Decision

`COMMON_CAPACITY_REGULATION_IS_NOT_THE_G5_SUSTAINMENT_REPAIR_FOR_THIS_PROVISIONAL_ROOT`

All 18 predeclared dynamic `R1`--`R3` gain/timing members integrate
successfully and retain positive core states, but none passes the scale-free
10-minute WT flow-shape gate. Dynamic beta/cAMP/PKA regulation remains a
mandatory biological module. The narrower numerical result is that tuning its
unmeasured common-capacity kinetics cannot repair sustained WT secretion on
root `G2_BALANCED_APICAL_K_BIASED_G10:M00`.

This is a decisive nested-model result for the frozen staging root and the
predeclared ensemble, not a claim about AE4 biology and not final full-model
validation. The root artifact is explicitly labelled
`R0_WT_RESTING_REFERENCE_ONLY_NOT_FINAL_FULL_MODEL`.

This preserved old-generation calculation used the original `TSLOW=(30,120)
s` sensitivity shown below.  It is not silently recomputed.  The replacement
native dynamic contract corrects TSLOW to `(30,90) s` and R1 slow to `90 s`;
neither contract identifies physiological kinetics.

## Frozen scope

The calculation loads the exact state, variant, and WT-only calibration tuple
from `results/13B_modern_full_model/wt_selected_candidate.json`. It holds fixed

- the 12-state conserved whole-cell root;
- corrected water coefficients and every other whole-cell parameter;
- the SR2-112 AE4 transport module and carrier abundance;
- mixed-bath routing and apical pump/K fractions;
- the separate CCh/calcium and IPR/beta inputs;
- the physical 600 s protocol and 5 s output grid;
- the production Radau tolerances; and
- the graphical WT minute-flow envelope of `9--10 uL/min`.

Only the already-declared regulator varies:

| Axis | Values | Evidential status |
|---|---|---|
| common-capacity gain | `1.25`, `1.70` fully activated multiplier | alternative approximate readings of the internally discrepant 2021 report; sensitivity only |
| family | `R1`, `R2`, `R3` | nested dynamic alternatives |
| timing | fast `(3,10 s)`, reference `(10,30 s)`, slow `(30,120 s)` upstream/downstream | predeclared unmeasured sensitivities, never fitted |

This gives 18 dynamic members. Two R0 static-gain members are carried as
controls. No phenotype target file is loaded by the runner.

## Quantitative result

For a single scale factor to place all ten minute landmarks within
`9--10 uL/min`, the model flow max/min ratio must not exceed

\[
10/9=1.111111,
\]

or equivalently minute-10/minute-1 flow must be at least `0.9`.

| Ensemble subset | Members | Best max/min ratio | Worst max/min ratio | Minute-1-to-10 decline | Passing members |
|---|---:|---:|---:|---:|---:|
| `R1`--`R3`, 1.25 gain | 9 | `1.732366` | `1.733590` | `42.275--42.316%` | `0` |
| `R1`--`R3`, 1.70 gain | 9 | `1.769549` | `1.773417` | `43.488--43.612%` | `0` |
| complete dynamic ensemble | 18 | `1.732366` | `1.773417` | `42.275--43.612%` | `0` |

The best equivalence class is low-gain slow `R2/R3`; their matched equations
are analytically isomorphic and their computed shape ratios differ by less
than `1.7e-11`. Its minute flow falls from `0.001103586 pL/s` at 60 s to
`0.000637040 pL/s` at 600 s:

- late/early ratio: `0.577245`, versus at least `0.9` required;
- decline: `42.275%`;
- max/min ratio: `1.732366`, or `55.91%` above the allowed ratio;
- feasible scale interval: lower `14127.84`, upper `9061.37
  uL min^-1 per (pL s^-1)`—an empty intersection;
- after mean-flow scaling: `7.250--12.560 uL/min` across minute landmarks;
- minute landmarks within `9--10 uL/min`: `1/10`.

The slow low-gain cascade improves the low-gain R0 static control's shape
ratio by only `0.0192%`. Increasing the gain from the prose sensitivity to the
visual-figure sensitivity makes the best ratio `2.15%` worse. Thus neither a
delay nor a larger common AE4 capacity increment supplies the missing
sustainment.

All 20 trajectories are numerically successful and positive. Separate BDF
method checks of the best and worst dynamic members agree with Radau to
maximum coordinate-relative differences below `4.7e-7`; cumulative-flow
differences are below `7.1e-8` relative. This is a distinct stiff algorithm,
not the independent equation transcription. The failure is therefore not
caused by the stiff solver choice. The separate sensitivity runner also
matches the subsequently frozen production `wt_trajectory_summary.csv`
exactly for all 20 Radau endpoint flow, cumulative flow, chloride, pH, and
volume fields.

## Nested-model interpretation

1. `R2` and symmetric `R3` remain exactly non-discriminable, as designed;
   matched cumulative flows differ by at most `3.3e-12 pL/cell`.
2. `R1` versus `R2/R3` timing changes the scale-free shape ratio by at most a
   few parts in `10^-3`, far short of the required change.
3. The larger 2021-compatible gain moves the model away from the WT shape
   envelope rather than toward it.
4. Selecting an unmeasured time constant outside the declared set merely to
   flatten this trace would violate the calibration firewall and is not a
   source-supported reconstruction step.

The regulatory chain should therefore be retained with its uncertainty, but
the next G5 iteration must address an independently defensible nonregulatory
stimulated-flux, water, outflow, input, or calibration axis. This conclusion
is conditional on the current staging root; a materially different valid WT
root must receive the same frozen ensemble check.

## Reproduction artifacts

- `analysis/13B_modern_full_model/run_regulatory_sustainment_sensitivity.py`
- `results/13B_modern_full_model/regulatory_sustainment_sensitivity.csv`
- `results/13B_modern_full_model/regulatory_sustainment_summary.json`

Reproduce with:

```bash
PYTHONPATH=src python analysis/13B_modern_full_model/run_regulatory_sustainment_sensitivity.py
```
