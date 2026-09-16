# Task 50 decision log

Every scientific decision must be recorded here before dependent compute. Each decision must cite the relevant master-ledger result/section and state why it does not reopen a closed axis.

## D50-00 — task scope

**Decision:** Task 50 will test/repackage exactly one reduced effective law:

`g_Cl,eff = g_parent * [1 - lambda * a_Ca * beta * (1 - e_AE4)]`

with fixed `lambda = 0.89488127156712` inherited from Task 41.

**Ledger basis:** master ledger R41, R42, R48, R49; sections 5-7.

**Classification:** TARGET-CALIBRATED CONSTRUCTION.

**Why not a reopened search:** no new mechanism family, no parameter optimisation, no stoichiometric/routing/NKCC/NHE/pump/K/VRAC change. The task only conditions Task 41's already selected constructive direction on the existing beta protocol input.

## D50-01 — timing and signalling

**Decision:** exact experimental time course is not a success/failure target. No cAMP, PKA, phosphorylation, beta-receptor or calcium signalling ODE may be added. `beta` is an algebraic protocol input only.

**Ledger basis:** `docs/PHENOTYPE_TARGET_CONVENTION.md`; master ledger sections 0, 6 and 7.

**Why not a reopened search:** this reduces complexity rather than adding it.

## D50-02 — reuse before recomputation

**Decision:** reuse frozen Task 40/41 results whenever exact vector-field, initial-state and protocol equivalence proves identity. Do not rerun a production trajectory merely to reproduce an already frozen result. If one required equivalence cannot be proven, run only the minimum specific case needed to resolve that equivalence.

**Ledger basis:** master ledger R40-R41 and remaining-shot budget.

**Reason:** only three funded scientific Codex runs remain; duplicate trajectories provide no new scientific information.

## D50-03 — interpretation boundary

**Decision:** Task 50 must not call the new scalar a molecular AE4-to-TMEM16A recruitment fraction. The permitted interpretation is an unresolved beta-conditioned effective network coupling controlling stimulated apical chloride secretory capacity.

**Ledger basis:** master ledger R41, R42 and section 9.

**Classification:** TARGET-CALIBRATED CONSTRUCTION, not SOURCE FACT.

## D50-04 — audit the stored coefficient before any implementation

**Decision:** complete a read only provenance and algebra audit of the Task 41
selected coefficient. The public driver, candidate 08 design and frozen input
record store `0.10511872843288446`, whereas the Task 50 instructions use the
shorter decimal `0.10511872843288`. Preserve both literals. Do not change the
mandated `lambda`, round the frozen coefficient, import a production runner,
evaluate the scientific RHS or integrate a trajectory during 50A.

**Ledger basis:** R41, sections 0, 7, 8 and 9; phenotype convention; D50-02.
These sections and the phenotype convention were consulted again before this
decision. The historical reports and both ledgers were read. The novelty audit
will search all fetched remote branch heads, as well as the current source,
prompts, analyses and results, for the coupling and its synonyms.

**Reason:** exact inheritance and exact identity must refer to the actual frozen
coefficient rather than a rounded report value. Machine precision comparison
and exact algebraic identity are different assertions. A finite numerical
trajectory cannot establish equality of unequal declared coefficients.

**Why not a reopened search:** this is the prescribed 50A audit of one existing
construction. No parameter selection, mechanism test or production calculation
is authorised by this entry. If the required exact identity fails, apply the
task's stop rule and publish the contradiction before scientific source edits.

**Classification:** provenance audit / FORMAL DEDUCTION concerning a
TARGET-CALIBRATED CONSTRUCTION; not a biological rejection.

## D50-05 — apply the exact identity stop rule and publish 50A

**Decision:** stop before 50B. Publish the completed 50A audit with the failed
exact inheritance gate, the algebraic proof, source provenance, novelty record,
status and master ledger qualification. No dependent scientific compute is
authorised. Do not describe 50B, 50C or 50D as completed.

**Ledger basis:** R41 and sections 7, 8 and 9, consulted again with the phenotype
convention immediately before this decision; Task 50 prompt's exact equivalence
and stop instructions; D50-04. The actual frozen parameter is the longer decimal
in the driver, design and manifest, all preserved since Task 41 commit
`86f135e075aa91a70a752b650d7bb894353fee04`.

**Evidence:** `EQUIVALENCE_PROOF.md` and `output/coefficient_audit.json` establish
the nonzero exact factor difference `-4.46e-15*a*(1-e)` at beta one. All seven
frozen execution hashes match. The four parent nesting identities remain valid
formally. Initial state and standard protocol records are consistent across
the three saved cases. This is a precision/provenance conflict, not a claim of
physiologically meaningful error or a rejection of the effective coupling.

**Why no resolving trajectory:** a numerical agreement tolerance could address
approximation, but cannot prove the required exact algebraic equality between
the specified and frozen coefficients. Silently changing lambda, rounding the
historical coefficient, or weakening exact identity would alter a binding
instruction. No such change is made.

**Why not a reopened axis:** only the existing coefficient and reuse premise
are audited. No parameter is chosen, no alternative model is introduced, and
Task 51 is not started. The run remains a TARGET-CALIBRATED CONSTRUCTION audit.

## D50-06 — authorised continuation with the specified decimal

**Authorisation:** after the published 50A audit the user answered "Yes please"
to retaining exactly `lambda=0.89488127156712` and treating the stored coefficient
rounding difference through numerical equivalence. This explicitly supersedes
D50-05's stop for this precision issue only. It does not change the parameter,
permit fitting or open another scientific task. This is continuation of Task 50.

**Ledger basis:** R41, R50A, sections 5 to 10 and the phenotype convention,
read again before this decision. The original audit, full frozen coefficient
and exact difference remain valid historical facts. All other nesting and scope
requirements remain binding. Classification: TARGET-CALIBRATED CONSTRUCTION.

**Decision:** publish this amended 50A boundary before implementation. Keep the
exact algebra for the displayed coefficient distinct from numerical comparison
with the full frozen Task 41 coefficient. Reuse frozen mutant results only with
that explicit qualification after the prescribed software checks pass.

## D50-07 — minimal implementation and predeclared software comparisons

**Decision before compute:** add one Task 50 wrapper that inherits the unchanged
Task 41 conductance reconstruction, RHS dispatch and public integration helper.
Override only the evaluation of the multiplier and its diagnostic labels.
Use the existing normalised calcium activation, protocol beta input and genotype
expression, with one fixed literal lambda and no configurable new parameter.
Do not construct or store a second recruitment coefficient in the new wrapper.

**Ledger basis:** R40, R41, R48, R50A, sections 2, 5, 7, 8 and D50-06. No closed
mechanism axis is reopened. The inherited 12 core states and existing AE4
regulatory coordinate remain unchanged; no signalling ODE is added.

**Verification fixtures:** use the frozen Task 40 WT rest and saved final states
of the three Task 41 candidate 08 cases. These are fixed software fixtures, not
a parameter panel. Check the required factor boundaries and representative
activation/expression values; exact parent equality in the off cases; complete
RHS, currents, voltages, chloride sources and charge bookkeeping at beta one;
unchanged fixed state water, AE4 and homeostasis evaluations; parameter identity
apart from the single apical chloride conductance; and mocked integration
dispatch. No integration or equilibrium solve is authorised by these tests.

**Numerical tolerances fixed before seeing comparison output:** off case RHS and
physical diagnostics must be bit identical. Factor comparison to the full frozen
coefficient has absolute allowance `8e-15`, covering the known `4.46e-15`
decimal discrepancy plus floating point evaluation. For active full model
comparisons use relative tolerance `2e-12` and absolute tolerances `2e-13` for
RHS and amount flux/source quantities, `2e-14 V` for voltages, `2e-25 A` for
currents, and `1e-24 S` for conductances. Near zero RHS components require an
absolute allowance because charge cancellation and the inherited voltage solve
are not exact arithmetic. These are software comparisons, not biological
equivalence bands. Apply the unchanged parent conservation gates as well.

Report actual differences. Passing fixed state comparisons does not establish a
global trajectory error bound or bit identity of the rounded mutant solutions.
The existing production integrator and voltage closure tolerances remain
unchanged. No threshold adjustment or parameter tuning is authorised.

## D50-08 — retain current instructions while loading frozen scientific inputs

**Decision before the next verification attempt:** the original Task 40 fixture
loader rejects today's `AGENTS.md` against its historical governance hash before
any test or model evaluation. The read only audit checks all 150 manifest files:
149 match and only `AGENTS.md` differs. Its actual Git blob
`94465907d246d530e78f18cf1237fad5088cfa51` is the binding Task 50 start version.

Use a task local loader that verifies the 149 unchanged inputs and explicitly
pins this current governance blob; also verify Task 40's own execution hashes,
rest vector and parameter hashes. Then construct the same parent with the
unchanged factory. Do not edit AGENTS, historical manifests, parent source or
production runners. This corrects a fixture's historical governance assumption,
not a scientific equation or an acceptance tolerance. Preserve attempt 01 and
`output/legacy_manifest_audit.json`.

**Ledger basis:** R40, R41, R50A continuation and section 8, checked again with
the phenotype convention. Current instructions take precedence over archived
instruction hashes. No closed axis is reopened. Attempt 01 performed zero
model evaluations and zero tests, so no failed scientific comparison is hidden.

## D50-09 — accept implementation checks and preserve the precision scope

**Decision:** accept the ten passing tests in verification attempt 02 and publish
50B before extracting phenotype summaries. The 21 parent comparisons are bit
identical, and the six complete active evaluations agree within the unchanged
D50-07 allowances. Largest factor difference is `4.5102810375396984e-15`;
largest RHS component difference is `8.090750291955828e-15` in native units.
The successful run used 80 core evaluations including constructor checks, zero
production integrations, zero stationary solves and zero fits. Attempt 01 was
the separately recorded governance hash setup issue, with zero evaluations.

**Ledger basis:** R40, R41, R50A continuation, R50B setup qualification and the
phenotype convention were checked again before accepting this result. The
parent scientific files remain unchanged. This is software verification of a
TARGET-CALIBRATED CONSTRUCTION, not independent biological validation.

**Dependent work after remote verification:** verify/hash the frozen Task 40/41
source, states, protocols and outputs and read the existing summary fields.
Do not integrate any trajectory. Use the already reported cumulative and
endpoint deficits and all saved 60 s cumulative observations to describe the
finite window persistence without selecting a favourable time. Numerical reuse
of mutant results carries the authorised coefficient precision qualification;
it does not establish an exact solution identity or a global error bound.
No new CCh only or AE2 trajectory is needed for the exact parent IVP statement.
Do not transfer Task 48's different bath/genotype rest outputs to this parent.
