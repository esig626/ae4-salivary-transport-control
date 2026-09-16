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
