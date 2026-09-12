# Task 28: repair transport relationships in a clean phase

Repository: `esig626/ae4-salivary-transport-control`.
Branch: `codex/task-28-transport-relationships`.
Read `AGENTS.md` first. This task supersedes the abandoned inverse search; do
not run Task 27 or import its search framework. Do not merge into main.

## Purpose

Test two specific hypotheses, rather than searching for a combination of
parameters that forces a knockout phenotype:

1. The NHE1 law needs intrinsic intracellular pH regulation in addition to
   the chemical driving force, so removal of AE4 does not leave an excessive
   proton extrusion demand that can be balanced only through alkalinisation,
   compensatory AE2 chloride uptake or swelling.
2. The WT model must give AE4 a meaningful role in productive chloride
   delivery during stimulation. An approximately 1% positive loading share
   cannot, without a demonstrated indirect mechanism, explain a large flow
   effect. WT secretion alone does not determine transporter allocation.

These are hypotheses to test, not established repairs. A constitutive change
must apply identically in WT and AE4 null. Do not add a genotype dependent
rescue, new transporter, acid production term, bicarbonate sink, signalling
cascade or arbitrary suppression of the pump or NKCC1.

## 0. Clean input and execution discipline

The old work is preserved at the references in `ARCHIVE_INDEX.md`. Do not
recreate the archive, remove old branches, rewrite history, or claim that
unpushed Task 27 results have been backed up. Leave other workspaces alone.

The essential scientific modules in `src/modern_full_model/` are unchanged
seed copies. `reference/native_wt_contract.json` provides the ten original
parameter payloads. `reference/R09_wt_rest.json` and `reference/R10_wt_rest.json`
provide routed WT amount states. The reference CSV preserves the larger
historical panel without requiring its trajectories.

Run `python scripts/check_phase28_inputs.py --check-original-code` and the
three setup tests before editing. Then implement a small input loader that
constructs dataclasses directly from `roots[*].whole_cell_parameters` and
`roots[*].ae4_parameters`. Do not import omitted calibration/task runners.
Retain the input protocol and R1 plus N1 settings from the reference manifest:
use `protocols`, `stimulated_nkcc1`, and
`regulatory_members.R1_G125_P21_PROSE_TREFERENCE`. The physiological model must
be wrapped with the imported routed AE4 evaluator. Do not accidentally use
the independent cation branch sources or the abandoned pooled mechanism.

First replay the two saved WT amount states through the unmodified RHS under
REST and reproduce their fluxes and residuals. No ODE run is needed for this
check. Record Python, NumPy and SciPy versions. The old reference metadata
lists files and tasks that are intentionally absent here; those are not
execution prerequisites. Never substitute an older model to make imports work.

## 1. Absolute compute limits

New calculations are limited to the two nominal backgrounds R09 and R10,
chosen before outcomes and representing the two inherited routing families.
The remaining eight backgrounds are preserved for later validation only.

Before starting numerical work create a persistent `results/budget.json` and
write counters after each case, including failed attempts and test runs:

- At most 100 stationary solver calls in total, including calibration,
  continuation, retries and any stationary stimulated calculations.
- At most 20,000 stationary residual evaluations in total, including numerical
  differentiation evaluations; each individual solver call is also bounded.
- At most four local WT calibration solver calls in total: one per background
  for NHE1 only and one per background for the combined repair. Prefer direct
  algebra; no multistart or retries disguised as another optimiser.
- At most 24 integrations of a stimulated trajectory, including WT controls,
  trial validations, final contrasts, BDF checks and tests.
- At most two worker processes, with numerical library threads set to one.
- At most 60 minutes of numerical execution, accumulated across restarts.

Stop at whichever limit is reached first. No automatic budget extension. The
scientific question may remain unresolved. Preserve partial outputs and report
what is missing. An individual long solve/integration must have a time limit
so the overall stop can actually be enforced.

Prohibited: Cartesian grids, exhaustive pairs, global or stochastic search,
Bayesian/evolutionary optimisation, automatic hyperparameter tuning, high
dimensional phenotype optimisation, Shapley enumeration, and any escalation
from a failed hypothesis to a larger search. Small fixed-state evaluations of
one pH response curve are allowed; they are not a parameter sweep.

## 2. Reason through the balances before changing the model

Write a short `analysis/balance_reasoning.md` with exact signs and units.
Use amounts, not concentration derivatives without dilution terms.

For the imported stoichiometry, denote positive chloride uptake by J4 (AE4),
J2 (AE2), JN (NKCC1 cycles), proton extrusion by JH (NHE1), and apical chloride
export by Ja. Verify from code:

    d nCl_i/dt = 2 JN + J2 + J4 - Ja
    d A_i/dt   = JH - J2 - 2 J4
    d T_i/dt   = JCO2_b + JCO2_a - J2 - 2 J4

At stationary state this gives Ja = 2 JN + JH - J4. This is a balance identity,
not a partial derivative proving that AE4 inhibits secretion. It also does not
by itself distinguish this model from the published 2018 balances. Explain
which rate relationships determine compensation.

Check NHE1 versus AE2 possible rates and signs in the intended pH/ion domain.
Distinguish insufficient capacity from an impossible direction at a specified
state. Do not spend numerical solves on an analytically contradicted state.
A fixed-state obstruction or failed continuation is not a global absence proof.

Verify water fluxes, cell/lumen volume balances and charge/current closure.
For stationary volumes derive the dependence of q_out on lumen minus bath
osmolarity using the existing hydraulic coefficients. Separate apical delivery,
paracellular chloride return, ductal output and temporary lumen storage.
Do not conclude that a larger cell necessarily secretes more water.

Read `reference/sources.md` and the 2018 paper/code comparison. Check units and
conditions against the papers rather than treating the audit as primary data.
Do not force acidification: the relevant resting pH observation is close to
unchanged. NHE assay activity, protein/capacity and net flux at altered gradients
must be distinguished. Do not equate an NKCC1 knockout flow reduction to its
WT chloride fraction.

## 3. One NHE1 constitutive correction, shared by both genotypes

Replace or extend the coarse NHE1 law with ONE parsimonious, source motivated
intracellular pH response. Prefer direct proton sensing, not a signalling ODE.
The existing chemical driving force already depends on H; explain precisely
what additional regulatory dependence the correction represents.

A permissible minimal form is a positive, bounded proton activation factor
multiplying the reversible exchange rate, for example a Hill modifier. This
is an option, not a measured acinar law. Establish the functional choice and
one set of response parameters from appropriate primary evidence before
knockout evaluation. Prefer salivary/acinar evidence; explicitly label any
borrowed cell type values and any quantitatively unmeasured assumption.
Do not infer transported proton stoichiometry from a Hill coefficient.

The law must retain 1 Na inward per 1 H outward, zero charge source, no carbon
source, proper reversal under reversed chemical driving, finite capacity,
and a declining proton extrusion response as the cell alkalinises at fixed
other inputs. No hard pH clamp or genotype dependence.

A fixed WT normalisation of capacity is allowed to isolate a change in state
dependence. If used, its reference is fixed once from WT and then applied to
both genotypes. It must not renormalise independently in the knockout and
must not undo finite saturation. Do not assume this correction can solve the
null acid/base balance merely because it has the desired qualitative slope.

Add inexpensive tests of stoichiometry, charge/carbon sources, reversal,
bounds, pH response and genotype independence before full model calculations.
Compare the isolated corrected response with the old law and the published
2018 NHE1 law. Do not silently restore unresolved old sign or unit conventions.

## 4. One WT chloride allocation repair, not a search

Compare exactly three model definitions:

A. Imported routed reference, unchanged.
B. Corrected NHE1 state dependence, without intentional chloride repartition.
C. The same corrected NHE1 law plus ONE WT chloride allocation hypothesis.

Keep A as the historical control. Construct B and C from WT information only;
never choose or revise a WT solution because of its null secretion result.
Finish and hash both WT constructions before running their genotype contrasts.

For C, use ONE declared exploratory target: AE4 provides 20% of the positive
basolateral chloride loading at the stationary WT stimulation condition
Ca = 0.25 uM. This is an explicit modelling assumption, not a measured fraction
and not an inference from the 30 to 35% knockout flow deficit. Do not additionally
try 10%, 30%, or any other allocation if it fails. Report nonidentification of
this share by WT data rather than presenting the imposed allocation as a result.

The positive pool contains max(J4,0), max(J2,0), 2 max(JN,0), and any other
actually implemented positive basal chloride loader. Report every signed
counterflux separately. Also report net basal loading, apical delivery and
paracellular return; none is interchangeable with the positive pool.

Set the share using stationary flux algebra and, only where necessary, one
bounded local WT closure solve per background, not an optimiser that repeatedly
runs 600 s trajectories. Rebalance the complete existing transporter sources,
not chloride in isolation. Reduce the contribution of other positive loaders
consistently while preserving a stated WT loading/physiology constraint.

The existing uncertain AE4/NKCC1/NHE1/AE2 capacities, pump capacity, K and CaCC
conductances and existing membrane allocation fractions may participate in
WT closure if necessary. Predeclare the exact adjustable fields and objective
before solving. Use algebra for linear capacity relations and a modest
reference-departure criterion for remaining nonuniqueness. Old 14 nS or other
fitted values are not immutable laws; do not impose an inherited arbitrary
capacity box as a biological impossibility proof. Report substantial changes.

Re-equilibrate the WT state when the constitutive law changes. Do not require
the old fitted state to remain exact. Use measured WT pH/Cl and reasonable
Na, K and volumes to state the physiology check; distinguish data uncertainty
from diagnostic acceptance margins. Exact trace or 1:1 time-signature matching
is not required. All current, charge, carbon and water balances still close.

Keep CO2 permeabilities, paracellular selectivity/conductance, hydraulic
permeabilities, outflow law and the stimulus/regulatory parameters unchanged
in this task. Do not lower water permeability or apical chloride conductance
only in the knockout to manufacture a secretion deficit. A WT CaCC rebalance
needed for electrical closure must be declared and retained in both genotypes.

Evaluate the ACTUAL stimulated AE4 share subsequently, including the integrated
positive flux over 0 to 600 s and the late interval. A stationary target or one
endpoint is not the integrated contribution. Do not add an ODE fitting loop if
the measured trajectory share differs from the imposed stationary share.

The inherited routed scalar AE4 law remains a provisional coarse description.
No opposing cation sources or charge error is allowed, but those checks alone
do not validate its entire thermodynamics. Do not claim that NHE1 repair also
establishes that unrelated AE4 modelling assumption.

## 5. Small matched evaluation

Use WT expression 1.0 and exact AE4 null 0.0. No expression dose response and
no new 5% run. Changing genotype sets only the complete AE4 source to zero.

First use Ca = 0.25 uM for the two backgrounds and the B/C definitions. Every
arm starts from its own correctly equilibrated REST state. Use the saved WT
reference states for A; do not rerun already available old knockout failures.
Report any missing reference as missing, not as a fabricated comparison.

Assess admissibility per case. Do not demand that all ten archived backgrounds
pass before computing a valid selected case. A failed null rest is not zero
flow and must not enter the phenotype statistics. Do not run production
secretion from a physiologically rejected rest. Retain failures and equations.
Check local REST stability in the independent conserved coordinates or perform
one small, budgeted perturbation check; an algebraic root need not be stable.

Use the inherited seconds axis, REST calcium 0.058 uM, CCh + IPR protocol,
R1/N1 settings, production Radau tolerances and conservation checks. Integrate
0 to 600 s. Audit the whole accepted trajectory, not merely the endpoints.
No arbitrary duration resting integration is required: use bounded direct
solutions with actual RHS/residual convergence tests.

Only if the central cases are numerically valid, validate C at Ca = 0.10 and
0.50 uM with no further tuning. Reserve at most four BDF confirmation runs
within the total 24 trajectory budget. This is a two-background diagnostic,
not a claim of ten-background robustness or experimental validation.

The known approximately 30 to 35% decrease is considered after the matched
results are frozen. Report direction and magnitude, but never continue tuning
to reach 0.65. A result outside that target is a valid outcome of this task.

## 6. Required explanation and outputs

Record for each model/background/genotype: Na, K, Cl, TIC, HCO3, pH, cell and
lumen volumes, both potentials, NHE1, AE2, AE4, NKCC1 cycles and chloride flux,
pump and K fluxes, CaCC efflux, both CO2 fluxes, each paracellular ionic flux,
cell/lumen/bath osmolarity, q_b, q_a, q_t, q_out and cumulative secretion.
Include stationary ledgers, selected time points and integrated transport.

Show whether reduced AE4 export is accommodated by less NHE1 extrusion,
AE2 reversal/increased chloride uptake, NKCC1 loading, altered cation balances,
carbon storage, fixed-solute dilution or luminal return transport. Retain
complete charge and water accounting, including the lumen storage correction
between integrated apical inflow and integrated ductal outflow.

Describe changes as model effects, not established biological causation.
B versus C can include coupled WT recalibration, so that difference alone is
not a clean single-parameter causal effect. No Shapley or combinatorial ablation.

Write only compact UTF-8 outputs under `results/28_transport_relationships/`
and `analysis/28_transport_relationships/`: contract and budget JSON; source
and parameter ledgers; WT freeze; resting states; flux ledger; chloride
allocation table; matched secretion comparisons; numerical verification;
focused test log; and a concise report. Store local arrays only if useful,
untracked. No binary arrays or archives in the publication deliverable.

Answer: Did pH regulation improve? Was an admissible exact null obtained? Did
actual WT chloride delivery acquire a meaningful AE4 contribution? What
replaced the lost flux? How did luminal salt, osmolarity and water change?
Which conclusion remains unproven? Was the budget obeyed? If unsuccessful,
identify the next algebraic question without launching another search.

## 7. Checkpoint and publication

Save the small text checkpoint after each case. Test only the focused setup,
conservation and new-law checks; do not rediscover the historical broad suite.

Do not contact GitHub repeatedly during numerical execution. If an ordinary
Git checkout and credentials are available, commit explicitly listed code,
tests and compact text outputs, then push this branch once and open a draft
PR for review. Do not merge or update main. Warn that the clean-tree phase
must not be merged as an accidental deletion of historical working-tree files.

If Git is unavailable or push fails, preserve local text outputs and report
them directly. Do not reconstruct Git history, attempt credentials, create
manual Git objects, encode Base64, package archives, or upload trajectory
files. Publication failure must not cause simulation reruns.
