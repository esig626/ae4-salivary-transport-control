# Task 52 decision log

## D52-01 — final-shot scope

Task52 is the last funded/scientific Codex shot. It is numerical only. No mechanism-family search is permitted.

## D52-02 — controlling mechanism

The only candidate is the chloride-reservoir/supply mechanism established in R51H and supported by the R51I reduced diagnostic:

- chronic AE4 loss lowers the measured intracellular chloride reservoir;
- WT retains ordinary beta/PKA-activated AE4 chloride supply;
- non-AE4 chloride supply is experimentally controlled/matched across genotypes in the central model;
- shared apical chloride demand converts the reservoir/supply difference into differential export;
- no AE4-dependent channel multiplier is allowed.

## D52-03 — Palk disposition

Palk/Benjamin NKCC is excluded from the Task52 candidate. It may be cited only as a historical negative diagnostic explaining excessive compensation in Tasks39–51.

## D52-04 — initial-state projection

Central measurement-constrained initialisation:

1. load accepted Task37 WT rest and active Task37 parameters;
2. retain Task37 intracellular Na and cell volume for both genotypes;
3. impose source central values WT Cl=50.10 mM, pH=6.91 and AE4-KO Cl=36.50 mM, pH=6.89;
4. solve only intracellular K and TIC such that exact cell bulk electroneutrality and bath-isotonicity hold using unchanged fixed charge, finite buffer, impermeant osmoles and acid-base equations;
5. derive TA from the existing speciation function;
6. keep lumen state and geometry unchanged;
7. reject before trajectories if either projected state is nonpositive or outside inherited broad physical limits.

Predeclared projection sensitivity: retain Task37 intracellular Na and TIC and solve K plus cell volume from the same two constraints. No other projection may be invented after outputs are seen.

These projected states are stimulation onset states, not claimed stationary solutions of the incomplete chronic KO transport model.

## D52-05 — matched non-AE4 supply

Use a paired WT/KO ODE.

At each RHS evaluation:

- evaluate pre-Palk NKCC and AE2 cycles from the WT state;
- use those WT cycles in WT as usual;
- impose the same cycle values and exact stoichiometric sources on KO;
- KO NHE1, NBC, pump, K channels, CO2, water, lumen and electrical closure remain its own local functions of its own state;
- WT AE4 is ordinary inherited AE4; KO AE4 is exact zero.

Fixed KO supply sensitivities multiply the shared WT NKCC and AE2 cycles together by 1.05 or 1.10. These are diagnostics, not fits.

## D52-06 — auxiliary beta demand

Shared effective beta/IPR apical chloride current:

`I_aux = beta * G_aux * (V_a - E_Cl)`.

Exact same `G_aux` and law in WT and KO. Insert before final electrical closure, include in cell/lumen chloride sources and current diagnostics, and do not multiply secretion directly.

Fixed values only:

- 0 S;
- 2.32e-9 S (DCPIB-sensitive source-scale equivalent);
- 4.49e-9 S (total IPR-current source-scale upper sensitivity).

No swelling gate is used in Task52. This is a protocol-level effective representation of the observed IPR-associated anion conductance, not a claim of direct molecular cAMP gating or complete apical attribution.

## D52-07 — case matrix

Combined CCh+IPR, 600 s:

1. matched supply, central projection, G=0;
2. matched supply, central projection, G=2.32e-9 S;
3. matched supply, central projection, G=4.49e-9 S;
4. central projection, G=2.32e-9 S, KO supply x1.05;
5. central projection, G=2.32e-9 S, KO supply x1.10;
6. alternate projection, matched supply, G=2.32e-9 S.

Controls:

7. CCh-only, central projection, matched supply, beta=0;
8. IPR-only, central projection, matched supply, G=2.32e-9 S.

No 5%-AE4 trajectory because no measured chronic 5%-AE4 chloride starting state exists.

## D52-08 — no-reveal tuning rule

All equations, projections, conductance values, supply sensitivities, solver settings, case matrix and acceptance diagnostics are frozen before any Task52 phenotype comparison.

No output may trigger parameter changes or additional cases.

## D52-09 — reservoir accounting

For every paired case compute and numerically verify:

`delta(t)=nCl_WT(t)-nCl_KO(t)`

and

`∫(J_WT-J_KO)dt = delta(0) + 2∫(N_WT-N_KO)dt + ∫(A_WT-A_KO)dt + ∫(E_WT-E_KO)dt - delta(T)`.

Central matched-supply NKCC and AE2 differences should be zero to numerical tolerance by construction.

A fluid deficit without this mass accounting is not sufficient to claim the reservoir mechanism.

## D52-10 — interpretation

Success means independently fixed/shared physiology plus the measured reservoir produces a substantial, persistent fluid deficit with all physical/conservation gates passing.

Failure is final. Do not add another mechanism. Task50 remains the target-calibrated quantitative proof of missing beta-conditioned coupling if Task52 does not succeed.

## D52-11 — exact Task37 source isolation and numerical audit

Before model evaluation, the complete binding stack and overlapping Task37–51
reports were read. R51H/R51I and the Task52-start ledger control over the older
R51G recommendation. Materialise exact Task37 source dependencies in a task-local
snapshot only after verified 52A: current optional NKCC selectors can import the
forbidden Palk module. Use Task37's original donor-weighted AE4 allocation;
do not substitute Task40's 50:50 routing. The required chloride identity is
routing-independent; R51G's equal-routing stationary identity is not a Task37
constraint. No closed scientific axis is reopened.

Inherited solver/physical gates and the prospective independent quadrature,
onset handling, failure handling, exact eight-case matrix and reporting windows
are fixed in ARCHITECTURE_AND_NUMERICS_FREEZE.md and output/case_matrix_52A.json.
The projection is a deterministic equation solve only; no resting-model solve,
parameter fit or production trajectory occurs before its authorised boundary.

## D52-12 — conservative implementation and fixed-state verification

After independently verified 52A (`acf3e132ce698d64a220b52a966c8472e86405a3`),
the exact historical snapshot is materialised without modifying its bytes.
R51H/R51I, Task52-start sections2–6 and D52-04–09 remain controlling.
For each fixed active protocol, adding `G_aux/gate(Ca_active)` to a temporary
copy of the maximal chloride conductance adds exactly `beta*G_aux` to the
effective conductance. Dispatch to the original parent at beta zero preserves
off-state identity. Both genotypes share the same two templates. This inserts
the auxiliary current into the unchanged full NBC and IPR-only current
closures; no biological conductance or water parameter is fitted. Diagnostics
retain explicitly named native TMEM16A, auxiliary and total currents.

The KO NKCC/AE2 override is applied after the final local NBC/homeostasis
evaluation with source increment `(dN,dN,2*dN+dE,-dE,-dE)`. Charge/carbon
diagnostics are rebuilt; NHE1/NBC/pump/K/current/water remain local. At fixed
state those two overridden pathways are electroneutral and do not enter the
electrical closure. Alternate projection retains Na and TIC concentrations,
with unchanged finite-buffer, fixed-charge and fixed-osmole amounts.

Verification is confined to prescribed projected onsets, the original accepted
state and a synthetic reversed chloride-driving fixture. All physical and
conservation gates retain52A tolerances; source-vector/equality tests use exact
identity where algebra requires it and round-off tolerances otherwise.
No scientific axis is reopened, and no trajectory is authorised at52B.

## D52-13 — immutable52C numerical freeze

After independently verified52B, freeze all four projected onset states,50 source/input identities and exactly eight cases. The production runner directly advances the26-coordinate paired RHS, preserves accepted endpoints plus integer-second diagnostics, and stops on the first sampled physical/conservation failure. Independent8-point and4-point Gauss quadratures include explicit2DN mass-term checks and separately bounded onset epsilon accounting. Failure records preserve partial coverage and unvalidated accepted endpoints. Static review corrections concern bookkeeping only; the tested scientific implementation is byte-unchanged. No scientific evaluation, integration, parameter change or tuning occurred at52C.

## D52-14 — fresh 52D execution from immutable 52C on the recovery branch

The current user instruction authorises only `analysis/task-52D-recovery-from-52C`
starting exactly at `36f7c3d0bb2353bf55ac0ed9e22582720097c829`, and only the
three central combined-stimulus 52D cases. It supersedes the old operational
branch name and automatic continuation through 52E/52F. No earlier unpublished
52D work is recovered, compared, reconstructed or reused.

Before numerical execution, fresh clone/remote fetch and independent GitHub
ref/commit reads agree on HEAD and tree. All 52 recorded 52C artifact SHA256
checks, all 50 immutable source/input SHA256 checks, the freeze SHA and ledger
prefix pass. `output/publication_52C.json` records this validation of the
existing freeze; no 52A, 52B or 52C scientific work is repeated. Runtime Python,
NumPy and SciPy versions equal the freeze. An independent read-only audit
confirms these checks and the runner's required output coverage.

The frozen runner contains the old branch literal. The auditable operational
launcher `output/launch_52D_recovery.py` changes only that literal in the
input-verification function in memory to the explicitly authorised branch.
Every other integrity gate and every scientific function remains unchanged.
For receipt generation only, `checkpoint.BRANCH` will likewise use the actual
recovery branch. Existing scientific source/input files remain byte-identical.

The scientific scope is exactly D52-07 cases 1–3, D52-09 reservoir accounting,
and D52-13 frozen numerics, grounded in master ledger R52A/R52B and the
R51H/R51I/Task52-start continuations. No closed axis is reopened. Use one
scientific process and one BLAS thread, the stored central onsets, unchanged
600 s target, Radau settings, quadrature and failure gates. Preserve every
attempt. Publish the three-case 52D summaries/checkpoint and new ledger result,
commit/push immediately, independently verify remote hashes, then stop.
No 52E, 52F, tuning, additional case, mechanism search or final classification.

## D52-15 — stop on frozen output-schema failure; no scientific replay

The first fresh52D process exited1 at run_frozen_cases.py's diagnostic-column
assertion before saving case01 arrays/summary. Static inspection identifies
active NBC current keys absent at exact-zero onset. This is not an explicit
Radau solver-status failure and therefore does not license BDF or a replay under
D52-13/NUMERICAL_FREEZE. No valid terminal time, numerical result or mass budget
is established. Cases02/03 remain unlaunched. Preserve the failure evidence,
append the new operational qualification to master-ledger R52D, and publish an
incomplete status. Do not create a completed52D checkpoint or run52E/52F.
Master-ledger R52A/R52B and R51H/R51I scientific scope is unchanged; no closed
axis is reopened. All50 immutable source/input hashes still pass. A bounded
serialisation-only remedy is described, not applied, in SOFTWARE_BLOCKER_52D.md.

## D52-16 — authorised serialisation-only repair and one fresh 52D attempt

Controlling authorisation: AUTHORIZATION_52D_REPORTER_RECOVERY.md at verified
branch HEAD f2b56f903d33a0ac5f5c292db92cead8b6157f1f. It expressly supersedes
D52-15's operational no-replay stop for this bounded reporting correction and
one fresh execution of cases01–03. The earlier failure at9c0fde42 remains
preserved byte-for-byte. Master-ledger R52A/R52B and R51H/R51I scientific
architecture, D52-07 case definitions and D52-13 numerics remain controlling;
no scientific decision, changed hypothesis or reopened axis is introduced.

The new external run_reporter_recovery_52D.py leaves every frozen source/input
file untouched. Four exact-match reporting substitutions adapt run_attempt
in memory: name-based union encoding, NPZ boolean presence mask with NaN for
absent fields, JSON endpoint nulls for absent fields, and finiteness metadata
that checks original present values/states while accepting labelled missing
entries. A genuine reported zero remains distinct from missing. Original
named diagnostic definitions and values are preserved. All integration,
observe/diagnose, physical/conservation, quadrature and summary formulae remain
unchanged. Frozen input verification retains its50 hashes,52C receipt/freeze
SHA and ledger-prefix gates; operational branch/authorisation-parent checks
are explicit. New stage/case records go under output/reporter_recovery_52D/
so the original attempted stage and case directory cannot be overwritten.

File-only and synthetic-dictionary checks pass with zero scientific evaluations.
Independent review confirms exactly four reporting edits and unchanged solver,
gates, reservoir and headline calculations. The integration try/except,
observe and window ASTs match the frozen function. Validation is recorded in
output/reporter_recovery_validation.json. The authorisation HEAD, all52
historical52C receipt hashes, all50 present scientific/input hashes and exact
freeze SHA9cebb8fcf85557643f6c7838ff57f9fb9603f5a7479cbe6fb63ba8c40b417ba2
were verified before compute. Runtime remains one scientific process/one BLAS
thread. Execute only52D, preserve any inherited gate failure honestly, publish
all three outcomes and stop without52E or52F.

### D52-16 setup-path qualification before scientific execution

The first reporter-recovery launch stopped before creating a stage record,
case directory or model: the frozen main's receipt hash lookup followed the
new output namespace rather than the unchanged input directory. Its traceback,
launcher version and pre-run receipts are preserved in
output/reporter_recovery_setup_failure_01/. Under the same authorised stage-
replay bookkeeping scope, exactly that hash lookup is redirected to SOURCE_OUT.
No scientific evaluation or production attempt occurred in this setup failure;
the one authorised fresh three-case numerical execution remains pending.

## D52-17 — publish the three frozen 52D numerical outcomes and stop

All three authorised central cases completed 600 s on their first Radau
scientific attempt after the reporting-only repair. Every original sampled
physical/conservation gate and independent reservoir audit passes; no fallback
or additional case occurred. All 50 immutable inputs remain byte-identical.
Independent file-only review reconstructs headline integrals, checks complete
integer grids, delta(t), the Gauss8/onset budget and Gauss4 cross-check, and
verifies the presence mask and original failure evidence. It calls no model.

Following master-ledger R52A/R52B, R51H/R51I and the current authorisation,
report NUMERICAL PREDICTIONS under the frozen MODEL IDEALISATIONS only. The
cumulative KO fluid deficits are 20.0394286%, 20.0757527% and 20.1015822%.
This records the requested 52D summaries without experimental/Task50 comparison,
mechanism search, sensitivity claims or 52F classification. Publish the normal
52D checkpoint, summaries and complete compact outputs immediately, verify
remote SHA/tree/output hashes independently, record publication, then stop.
No scientific decision reopens a closed axis; 52E and 52F remain unexecuted.
