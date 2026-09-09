# Paper viability audit

## Decision

**Classification: `STOP`.**

This is a decision about the proposed paper under the currently published
specification, not a claim that AE2/AE4 discrimination is scientifically
unimportant. The exact structural results are valid, but they do not support
the model-specific minimal physiological measurement result requested by the
prompt. No manuscript should be drafted from this branch.

## Audit questions

### What is the strongest new theorem or proposition?

The strongest model-specific proposition is a conditional steady-flux
decomposition. Let

\[
A=C-2N=J_2+J_4,\qquad B=J_2+2J_4,
\]

where \(C\) is outward CaCC chloride particle flux, \(N\) is NKCC1 cycle flux,
and \(B\) is the signed buffer/NHE1 aggregate under the published balances.
Then

\[
J_2=2A-B,\qquad J_4=B-A.
\]

For nonnegative forward cycle rates, the exact feasible wedge is
\(A\le B\le2A\). Nonzero published turnover factors then give the activities
by division. This decomposition is useful and model-specific in its mapping of
published fluxes, but it is an instance of an elementary two-signature linear
inverse problem. The broader nuisance-flux rank criterion is standard linear
algebra rather than a new general theory.

### Is it generic or specific to this model?

The constant-rank fiber theorem, the implicit-function sensitivity formula,
the two-output determinant condition, and the two-condition gradient-stacking
condition are generic. The coefficients \((1,1)\) and \((1,2)\), the forward
wedge, and the translation to CaCC/NKCC1/buffer fluxes are specific to the
printed AE2/AE4 balance structure.

### What is the strongest new physiological conclusion?

Only a conditional conclusion is supported: flow-related chloride transport
alone cannot assign the shared chloride-cycle aggregate to AE2 versus AE4;
an independent bicarbonate-equivalent aggregate separates their cycle fluxes
if nuisance fluxes are known or independently constrained. The current record
does not certify that any proposed pair such as \((Q,[Na]_i)\),
\((Q,[K]_i)\), or \((Q,[Cl]_i)\) has full rank at the physiological baseline.

### What is merely a numerical illustration?

The sampled \((J_2,J_4)\mapsto(A,B)\) CSV and its figure illustrate an exact
linear map on normalized knockout-to-WT cycle-rate coordinates. They are not a
simulation of the published parameter-to-physiology map. The water, current,
acid/base, and kinetic calculations are direct inconsistency diagnostics, not
calibrated predictions.

### Did escalation uncover anything stronger than local identifiability?

It produced two exact conditional statements: full-state activity uniqueness
when the two forcing columns have rank two, and the aggregate flux inverse
above. It also produced a sharp regular two-condition criterion. None becomes
a new model-specific global partial-observation or perturbation theorem without
an executable equilibrium map. The general forms have close prior precedent.

### Does the novelty audit support publication?

No, not in the current state. The search found no formal AE2/AE4 salivary
identifiability study, but absence of that exact application is insufficient.
Minimal-output selection, steady-state rank methods, multi-condition
identifiability, perturbation design, stoichiometric criteria, and pump-leak
inverse inference all have close precedents. The potentially publishable niche
would be a certified model-specific panel classification or perturbation
result, which this specification cannot presently support.

### Is the result enough for a short mathematical-biology paper?

No. A paper built only from a scalar-output dimension argument and a 2-by-2
stoichiometric inverse would be routine. Presenting the structural figure as
physiological observation geometry would overstate the evidence.

## Acceptance and blocker status

| Requested element | Status | Reason |
| --- | --- | --- |
| Independent symbolic reduction | Completed conditionally | Every elimination and extra convention is documented |
| Numerical verification against full published baseline | Blocked | Printed water, current, acid/base, and unit conventions do not close |
| Correct local identifiability mathematics | Completed | Regular and singular cases are separated; proofs are recorded |
| Minimal physiological observation panel | Not certified | Requires the unavailable physiological equilibrium Jacobian |
| Reproducible geometry | Completed only in exact flux coordinates | It is not the requested physiological observation image |
| Global/perturbation/design/testing scan | Completed | No stronger model-specific result survived the evidence gate |
| Targeted novelty audit | Completed | Generic ingredients have close precedents |

## Strictly necessary next work

Before reconsidering a paper, obtain either the authors' executable model or a
corrected public specification that resolves all of the following together:

1. concentration and whole-cell versus density units, including the AE4
   fourth-order rate scaling;
2. water-permeability/calibration values that reproduce both volume balance
   and reported flow;
3. channel-current, particle-flux, and transporter-flux conversions that close
   the membrane and mass balances;
4. the CO2 aggregate sign and acid/base closure; and
5. the exact baseline calcium input and fitted baseline parameters.

After reproducing one steady state, the minimum resumption sequence is:

1. verify \(F_u\), \(F_\theta\), and implicit derivatives against re-solved
   finite differences;
2. certify scaled ranks for all biologically feasible two-output panels and
   compute their observation images on an explicitly justified activity
   domain; and
3. test one published-model-supported perturbation for rank rotation and
   robustness under measurement noise.

Until then, further parameter sweeps or manuscript prose would not remove the
scientific blocker.
