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
