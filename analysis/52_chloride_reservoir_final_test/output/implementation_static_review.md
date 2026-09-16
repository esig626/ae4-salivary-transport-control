# Task52 independent static implementation review

Scope: source inspection only, after independently verified 52A. The reviewer
read the complete binding Task52 stack, including R51H/R51I, the results ledger,
analysis index and relevant historical reports before reviewing the exact
Task37 model, acid-base, water, membrane, NBC and regulation code. No scientific
module was imported, no model or projector evaluated, and no test or trajectory
executed by this reviewer.

Reviewed Task52 files: `task52_model.py`, `project_onset.py`, and
`PROJECTION_DERIVATION.md`. No blocking correctness defect was identified in
the inspected implementation. This is a static software review, not numerical
verification or a mechanistic result.

## Shared auxiliary current

The active template replaces the stored maximal chloride conductance by
`g_parent + G_aux/h`, where `h` is the inherited calcium Hill gate. Multiplying
by that same gate in the unchanged membrane code gives `g_parent*h + G_aux`.
This implements the prescribed additive current for the frozen protocols
because calcium is constant throughout each active interval, `h>0`, beta is
restricted to exactly zero or one, and beta zero dispatches to the untouched
parent. It is not a general implementation for arbitrary time-varying calcium,
zero calcium gates, or fractional beta inputs; those are outside Task52 scope.

Both inherited electrical pathways consume the active template: the baseline
closure used for IPR-only and the NBC-inclusive closure used under active CCh.
Thus the zero-NBC shortcut does not suppress the auxiliary current. The native
K, pump, NBC and paracellular laws are retained; their currents can change
through the newly solved local voltage. No genotype expression enters the
auxiliary conductance calculation.

`cl_apical` and `cl_apical_total` are total currents. Native TMEM16A and
auxiliary components are explicitly named `cl_tmem16a` and `cl_auxiliary`.
Reports and quadrature must respect these labels to avoid counting the total
as TMEM16A or adding the auxiliary current twice. The supplied component-sum
residual allows numerical verification of the floating-point implementation.

## Paired supply

The override is applied after the final `MinimalNbcModel` evaluation, so its
active-CCh homeostasis recomputation cannot replace the imposed supply later.
The WT is evaluated first; the KO receives its instantaneous signed NKCC and
AE2 cycle values, with only the declared common sensitivity multiplier.

For differences `dN` and `dE` from the local KO values, the correction is
`(dN,dN,2*dN+dE,-dE,-dE)` in Na/K/Cl/TIC/TA amount coordinates. This is the
complete inherited electroneutral source vector. The implementation refreshes
homeostasis, cell charge-rate and carbon-accounting diagnostics. KO NHE1,
NBC, pump, K, water, CO2, membrane voltage and lumen remain local. The exact
KO genotype selects the inherited early zero-AE4 return.

The shared source idealisation retains signed AE2 flux; multiplying a negative
AE2 value increases its outward magnitude. It must not be described as a
separately positive AE2 chloride influx or silently clipped.

## Onset projection

The projector derives the exact fixed-pH relation
`TA=a*TIC+b/V+w`, including finite buffer charge, and eliminates K using the
unchanged bulk charge and bath-osmolarity equations. Central K/TIC and alternate
K/volume solutions match the prescribed equations. Alternate Na and TIC are
explicitly retained as concentrations. Finite buffer, fixed charge and
impermeant osmole amounts remain fixed and dilute with the altered volume.
Finite buffer osmoles are included once; fixed charge is not counted as a
second osmotic pool. TA is finally constructed with the inherited function.

Lumen and regulatory coordinates remain unchanged. Bath-isotonicity does not
imply a stationary water or transport state because the lumen is retained.
The implementation rejects nonpositive or singular solutions. Independent
recovered pH/speciation, charge/osmotic tolerances and broad physical gates
remain the responsibility of the root's prescribed 52B tests.

## Numerical evidence still required

Static inspection does not replace the 52B tests of total/component current
agreement and signs, cell/lumen chloride conservation, exact imposed supply,
WT off-state parent nesting, projector residuals and recovered pH/physiology.
The mandatory production reservoir identity must subsequently use signed
total apical chloride export and independently integrated N/A/E terms.
