# Task 29: exact AE4-null algebraic feasibility audit

Repository: `esig626/ae4-salivary-transport-control`
Branch: `codex/task-29-exact-null-feasibility-audit`

Read `AGENTS.md` first. Use maximum/deep reasoning if available, but remain a
single agent. Do not compensate for difficulty by spawning subagents or by
running a larger numerical search.

## Central question

Before changing another transporter or fitting any phenotype, determine whether
the CURRENT transport architecture, together with the single Task 28 NHE1
pH-regulation candidate, can in principle support a physiologically reasonable
exact AE4-null stationary state near the experimental resting observations.

The task is to distinguish four possibilities:

1. **direction/topology conflict**: required steady fluxes have signs that the
   existing constitutive laws cannot simultaneously produce;
2. **finite-capacity conflict**: the required directions are possible, but the
   frozen capacities cannot achieve the required equalities;
3. **electrical/cation/water closure conflict**: acid-base closure is possible,
   but pump/K/CaCC/current or osmotic volume closure makes the target state
   impossible under the existing architecture;
4. **no algebraic contradiction found**: the architecture can in principle
   support the state, so Task 28 failure is numerical/calibration-related rather
   than a proved structural impossibility.

A fifth valid outcome is **unresolved with the permitted analysis**. Do not
force a structural conclusion.

This task does NOT seek the 30-35% secretion deficit. It does NOT allocate 20%
AE4 loading. It does NOT run stimulation. It does NOT fit anything.

## Inputs and provenance

Start from the clean Task 28 seed already present on this branch. Production
model modules are frozen for this audit.

Read:

- `reference/task28_partial_handoff.md`
- `reference/native_wt_contract.json`
- `reference/R09_wt_rest.json`
- `reference/R10_wt_rest.json`
- `reference/sources.md`
- the relevant source files in `src/modern_full_model/`, especially
  `model.py`, `membranes.py`, `water.py`, `acid_base.py`, `transporters.py`,
  `parameters.py`, and `ae4_routing_only.py`.

The Task 28 handoff was supplied after an unpushed run. Treat it as provenance
and a description of the candidate NHE1 law, not as independently verified
repository evidence. Re-derive every equation needed for Task 29 directly from
code.

At task start, test ordinary Git authentication exactly once with a harmless
read such as `git ls-remote origin HEAD`. Record success/failure. Git failure
must not stop the scientific audit and must not trigger authentication work.

## Experimental state to interrogate

Use the 2015 AE4-null resting observations as context, not optimisation targets:

- nominal intracellular chloride: 36.50 mM;
- reported spread: 36.50 +/- 1.60 mM;
- nominal intracellular pH: 6.89;
- reported spread: 6.89 +/- 0.02.

WT context is approximately Cl_i = 50.10 +/- 1.50 mM and pH = 6.91 +/- 0.07.

Do NOT require knockout acidification. Near-unchanged pH is the relevant
observation.

Na_i, K_i, TIC_i, HCO3_i and volume are not experimentally fixed here. Do not
invent narrow physiological boxes merely to prove impossibility. Instead derive
what ranges/values would be REQUIRED by the equations, then compare those
requirements with the saved R09/R10 WT states and ordinary physical
admissibility. State clearly where biological evidence is absent.

The inherited volume admissibility ceiling of 3 pL may be used as an existing
model gate. Do not introduce a tighter one for this audit.

## Task 28 NHE1 candidate to audit

Do not edit the production transporter law. In a standalone audit helper, if
needed, evaluate the Task 28 candidate exactly as supplied:

    J_H = C_H * tanh(log((Na_e * H_i)/(Na_i * H_e)) / w)
          * g(pH_i) / g(pH_WT_ref)

    g(pH) = 1 / (1 + 10**(1.40 * (pH - 6.57)))

Use each background's fixed WT reference pH for the normalization exactly as
Task 28 did. Apply the same law to WT and exact AE4 null. This is a candidate
state-dependence law, not established acinar biology.

Verify its source stoichiometry and sign conventions, but do not rerun the
Task 28 constitutive test suite unless needed for a specific audit check.

## 1. Derive the exact stationary identities

Work in conserved AMOUNTS and verify signs from code. Do not reason from
concentration ODEs while silently dropping dilution terms.

Using the model's positive directions, define:

- `J_N`: NKCC1 cycle influx;
- `J_H`: NHE1 Na-in/H-out exchange rate;
- `J_2`: AE2 chloride-in/bicarbonate-out cycle rate;
- `J_4`: AE4 chloride-loading rate;
- `J_a`: apical CaCC chloride cell-to-lumen flux;
- `J_P`: total Na/K-pump cycle rate;
- `J_K`: total K channel cell-to-outside flux across both membranes;
- `J_CO2`: net CO2 entry to the cell, basolateral plus apical, under the
  code's sign convention.

Verify or correct the following intracellular amount balances from code:

    d nCl_i/dt = 2 J_N + J_2 + J_4 - J_a
    d A_i/dt   = J_H - J_2 - 2 J_4
    d T_i/dt   = J_CO2 - J_2 - 2 J_4
    d nNa_i/dt = J_N + J_H + AE4_Na_source - 3 J_P
    d nK_i/dt  = J_N + AE4_K_source + 2 J_P - J_K

At exact AE4 null, derive the complete necessary stationary identities. In
particular verify whether they reduce to:

    J_2 = J_H
    J_CO2 = J_H
    J_a = 2 J_N + J_H
    J_P = (J_N + J_H) / 3
    J_K = J_N + 2 J_P = (5 J_N + 2 J_H) / 3

If any identity is wrong because of an omitted implemented source, correct it
and document the source. Do not proceed with an incorrect simplified ledger.

Also derive the charge/current redundancy or extra constraint implied by these
identities and the algebraic membrane-current closure.

Write the derivation before numerical evaluation.

## 2. Derive transporter direction conditions analytically

For exact AE4 null near pH 6.89 and Cl_i 36.5 mM, derive from the actual
constitutive laws the sign conditions for:

### NHE1

From its reversible chemical affinity, derive the exact inequality on Na_i and
pH_i required for positive proton extrusion. Then include the Task 28 pH
modifier and prove what it can and cannot change. In particular, distinguish
changing RATE MAGNITUDE from changing the thermodynamic REVERSAL SIGN.

### AE2

From its actual affinity, derive the exact inequality involving Cl_i, external
Cl, intracellular HCO3 and external HCO3 required for positive chloride influx
/ bicarbonate extrusion. Convert this into a required intracellular HCO3
threshold at the nominal knockout Cl_i, then use the existing acid-base
speciation equations to express the corresponding TIC requirement at pH 6.89.

Do this at the nominal point and at the reported Cl/pH interval endpoints by
direct algebra, not a grid.

### NKCC1

Derive its sign condition and the dependence of its chemical drive on Na_i,
K_i and Cl_i. Show how lowering Cl_i changes the drive if Na_i/K_i do not
change, and derive what cation-product change would be needed to hold its
chemical affinity fixed.

### CO2

Write the exact implemented net CO2 source and determine what cell/lumen CO2
relationship is required for `J_CO2 = J_H` at the null stationary state.
Separate a sign conflict from a magnitude/capacity conflict.

Produce a compact `direction_conditions.md` and a machine-readable table of
these inequalities/thresholds.

## 3. Fixed-capacity acid-base feasibility

Now ask the narrow question:

> At experimental-like pH and Cl, can the corrected NHE1 candidate and current
> AE2/CO2 laws satisfy J_H = J_2 = J_CO2 with the existing R09 and R10 frozen
> capacities for ANY physically admissible values of the remaining state
> variables?

Do not run a multidimensional solver.

Use monotonicity, capacity ceilings and algebraic bounds first. Derive useful
necessary bounds such as:

- the maximum possible positive AE2 rate under its frozen capacity;
- the Na_i value/range needed to reduce positive NHE1 extrusion to a rate that
  AE2 could match at pH 6.89;
- the HCO3/TIC needed to put AE2 in the correct direction;
- the CO2 gradient needed to match that same rate.

If one-dimensional inversion is genuinely needed, use a bracketed scalar root
only. No Newton multivariate root finder, least squares, differential evolution,
multistart or continuation.

Evaluate R09 and R10 separately because their WT-normalized NHE1 modifier
references differ. Report whether the obstruction, if any, is:

- wrong sign;
- insufficient finite capacity;
- requirement for implausible Na/TIC/CO2;
- or no obstruction at this layer.

Do not call an extreme required value biologically impossible unless supported
by the inherited model gates or cited evidence. "Far from retained WT states"
is a valid weaker statement.

## 4. Cation and electrical closure

Only if acid-base sign/capacity feasibility is not already contradictory,
propagate the necessary null flux identities into the cation/electrical system.

Using the current NKCC1 law and the derived/required `J_H`, compute the REQUIRED:

    J_P = (J_N + J_H)/3
    J_K = (5 J_N + 2 J_H)/3
    J_a = 2 J_N + J_H

Compare these requirements with:

- the frozen Na/K pump capacity and substrate gating;
- the total K-channel conductance and its Ca-independent REST state;
- the apical CaCC conductance at REST calcium;
- the algebraic apical/basolateral/paracellular current closure.

Do not fit conductances. Determine whether a necessary rate exceeds a hard
implemented maximum, requires an impossible current direction, or merely differs
from the saved WT operating point.

If voltage is needed, evaluate the existing algebraic closure at a SMALL number
of explicitly declared fixed states. Do not solve for a full stationary state.

## 5. Water and cell-volume closure

Re-derive from `water.py`:

    dV_i/dt = q_b - q_a
    dV_l/dt = q_a + q_t - q_out

and the stationary relationship between q_out and lumen-minus-bath osmolarity.

At the nominal pH/Cl target, determine what the ion/carbon requirements derived
above imply for intracellular osmolarity and fixed-solute dilution as volume
changes. Ask whether `dV_i = 0` can be compatible with a volume <= 3 pL without
using a numerical volume search.

Do not interpret swelling as increased secretion. Do not calculate a knockout
secretion phenotype or integrate 600 s.

## 6. Classification

The final report must classify the current architecture using the strongest
statement actually established:

### A. DIRECTION/TOPOLOGY CONTRADICTION
Use only if required stationary flux directions cannot coexist for any state
within the explicitly analysed experimental/admissible domain, independent of
finite capacity tuning.

### B. FIXED-CAPACITY CONTRADICTION AT THE EXPERIMENTAL-LIKE STATE
Use if signs can coexist but the frozen R09/R10 capacities cannot satisfy the
required equalities. This is NOT a proof that the architecture can never work
under another WT calibration.

### C. ELECTRICAL/CATION OR WATER CLOSURE CONTRADICTION
Use only if a necessary derived flux exceeds an implemented hard maximum or
violates the algebraic closure under the analysed state domain.

### D. NO ALGEBRAIC CONTRADICTION FOUND
Use if all necessary conditions can overlap. Then state clearly that Task 28's
failed nonlinear solve does not establish structural impossibility.

### E. UNRESOLVED
Use if the permitted reasoning/evaluations do not decide the question.

If a contradiction is found, identify the MINIMAL missing degree of freedom or
mis-specified relationship suggested by the equations, for example NHE1
magnitude regulation, AE2 direction/capacity, carbon handling, cation closure,
or chloride allocation. Do NOT implement the repair in this task.

Do not propose cAMP or other signalling merely because the task remains
unresolved.

## 7. Numerical limits

This is deliberately tiny.

Absolute limits for the entire task:

- zero multidimensional stationary solves;
- zero parameter optimisations;
- zero ODE/stimulated integrations;
- at most 50 fixed-state transporter/full-RHS evaluations;
- at most 12 one-dimensional bracketed scalar roots;
- one scientific worker/process;
- numerical library threads = 1;
- at most 10 minutes total numerical execution.

Most of the task should be pencil-and-paper-style derivation checked against
code. If a conclusion can be obtained analytically, do not spend the numerical
budget confirming it repeatedly.

No parameter grids. Evaluating the two experimental interval endpoints and the
nominal point is not a grid and is allowed.

## 8. Outputs

Create only compact UTF-8 text outputs:

`analysis/29_exact_null_feasibility_audit/`

- `final_answer.md`
- `stationary_balance_derivation.md`
- `direction_conditions.md`
- `closure_reasoning.md`

`results/29_exact_null_feasibility_audit/`

- `contract.json`
- `fixed_state_evaluations.csv`
- `direction_thresholds.csv`
- `required_flux_ledger.csv`
- `capacity_bounds.csv`
- `verification.json`
- `budget.json`

A standalone audit script and focused tests may be added if useful, but do not
edit production model files. No binary arrays, archives or large traces.

The final answer must directly answer:

1. Can NHE1 and AE2 run in the simultaneously required directions near the
   experimental null pH/Cl?
2. With frozen capacities, can their rates be equal there?
3. Can net CO2 flux equal that same rate?
4. What Na_i, TIC/HCO3 and CO2 conditions are mathematically required?
5. Given those conditions, can NKCC1, pump, K and CaCC/current closure coexist?
6. Can cell volume close below the inherited 3 pL ceiling?
7. Is there a proven structural contradiction, only a fixed-parameter
   contradiction, or no contradiction?
8. What is the single most informative model relationship to repair next?

## 9. Publication

Do not contact GitHub during the scientific work. At the end, if ordinary Git
authentication passed the single startup check, commit the compact text outputs
and any standalone audit helper/tests to this branch and push once. Do not
merge and do not open a PR unless explicitly instructed.

If Git authentication is unavailable, preserve the local outputs and report
that publication was unavailable. Do not retry authentication, create manual
GitHub objects, encode Base64, create an archive, or rerun science because push
failed.
