# Task 52 — chloride-reservoir final numerical test

## Objective

Use the single remaining scientific/Codex shot to test one frozen mechanism in the full conservation-explicit model:

**measured chronic chloride reservoir + ordinary WT AE4 supply + matched non-AE4 chloride supply + shared beta/IPR apical chloride demand.**

This is the final numerical test. Do not search for another mechanism.

## Repository and branch

Repository: `esig626/ae4-salivary-transport-control`

Work only on:

`analysis/task-52-chloride-reservoir-final-test`

Read `AGENTS.md` first, then `analysis/52_chloride_reservoir_final_test/CODEX_START_HERE.md` and every binding file listed there.

The exact prepared scientific parent is canonical main merge:

`d2b14c2821dd7e38703620389097cc8f1052977a`

The Task50 preservation branch must remain untouched:

`archive/task-50-working-effective-coupling-benchmark`

at `b16c30094b95f61a73d8f1977cd58e79c7bb50f6`.

## Absolute prohibitions

Do not:

- use or import the Palk/Benjamin NKCC model in the Task52 candidate;
- use the Task41/Task50 AE4-dependent CaCC multiplier;
- add any AE4-expression term to an apical channel;
- tune any value against the approximately 35% secretion phenotype or Task50 outputs;
- use a grid, sweep, stochastic search, global optimiser, evolutionary method, mechanism enumeration or subset search;
- add a new transporter or signalling state;
- change AE4 stoichiometry/routing;
- alter NHE1, NBC, pump, K, CO2, water, geometry or paracellular laws for phenotype rescue;
- solve a new chronic resting-state model and call it experimentally identified;
- invent a 5%-AE4 chronic chloride value;
- retune after seeing Task52 outputs;
- start any Task53.

## Authoritative scientific parent

Use the accepted **pre-Palk Task37** conservation-explicit model as the transport chassis. The Task37 model contains the accepted mechanistic NHE1, stimulus-recruited NBC, AE4 regulation, electrical closure, water/lumen model and the pre-Palk NKCC implementation.

Historical later models may be read for diagnostics and tested utility code only.

## Stage 52A — source, architecture and case freeze

Before any new model evaluation:

1. Reconstruct and verify the exact Task37 scientific parent and all required source/result hashes.
2. Verify that the active Task52 dependency plan contains no `nkcc1_palk2010.py`, Task41 recruitment or Task50 effective-coupling module.
3. Verify the Task50 archive ref independently and record it.
4. Freeze the exact case matrix from D52-07.
5. Freeze solver settings before any production trajectory. Use the inherited accepted stiff-solver settings unless a purely software compatibility issue requires an explicitly documented equivalent setting. Do not loosen physical/conservation gates.
6. Freeze numerical quadrature and observation windows before output.
7. Create machine-readable `output/checkpoints/52A.json`.
8. Update `CURRENT_STATUS.md` and ledger only for genuinely new source/provenance facts.
9. Commit, push and independently verify 52A remotely before 52B.

No production trajectory at 52A.

## Stage 52B — measurement-constrained initial states and paired model implementation

### Central onset projection

Load the accepted Task37 WT resting state and active parameter set.

For both WT and AE4 KO:

- retain Task37 intracellular Na concentration;
- retain Task37 cell volume;
- impose source central resting chloride and pH:
  - WT: `[Cl]_i=50.10 mM`, pH=6.91;
  - AE4 KO: `[Cl]_i=36.50 mM`, pH=6.89;
- solve only intracellular K and TIC from:
  1. exact cell bulk electroneutrality;
  2. cell osmolarity exactly equal to the active bath osmolarity;
- derive TA through the unchanged existing acid-base/speciation implementation;
- retain unchanged fixed-cell charge, buffer amount, impermeant osmoles, lumen state and geometry;
- regulatory state begins at the inherited rest value.

This must be a deterministic two-equation projection, not a secretion fit or a stationary-root optimisation.

Require positive K/TIC, exact charge/osmotic closure and inherited broad physical ranges. If central projection fails these predeclared conditions, publish the failure and stop without inventing another central state.

### Alternate projection sensitivity

Also implement but do not yet run production:

- retain Task37 intracellular Na and TIC;
- impose source Cl/pH;
- solve K plus cell volume from the same exact charge and bath-isotonicity equations.

This is the only allowed initial-state projection sensitivity.

### Paired matched-supply model

Implement a paired WT/KO RHS so the two cells are advanced on one common time grid.

At each evaluation:

1. Evaluate the inherited **pre-Palk** WT NKCC and AE2 cycles from the WT state.
2. WT uses those cycle values.
3. KO receives the same instantaneous NKCC and AE2 cycle values, with the exact inherited cycle stoichiometric sources applied to the KO conserved coordinates.
4. KO NHE1, NBC, pump, K currents, CO2, lumen, water and electrical closure are evaluated from the KO state normally.
5. WT AE4 remains the inherited WT transporter with its existing beta/PKA regulation.
6. KO AE4 expression is exactly zero.

Do not share NHE1/NBC/pump/water/voltage between genotypes.

Implement one explicit KO supply sensitivity argument, default `1.0`, permitted only at fixed `1.05` and `1.10`; it scales the imposed WT NKCC and AE2 cycles together for KO. No other value is permitted.

### Shared auxiliary beta/IPR chloride demand

Implement the identical current in both genotypes:

`I_aux = beta * G_aux * (V_a - E_Cl)`.

It must be inserted into the existing apical current closure and chloride cell/lumen source balances exactly, with conservation diagnostics.

Allowed `G_aux` values only:

- `0.0` S;
- `2.32e-9` S;
- `4.49e-9` S.

No swelling gate. No AE4 input. No secretion multiplier.

### 52B tests

At minimum test:

- no Palk import/dependency;
- exact KO AE4 zero;
- central/alternate projector charge and osmotic closure;
- auxiliary current zero at beta=0;
- auxiliary current identical law/parameter in WT and KO;
- current sign and Cl source sign;
- total cell+lumen Cl conservation for the auxiliary current;
- paired WT supply is exactly imposed on KO at supply multiplier 1;
- x1.05/x1.10 multipliers are the only accepted sensitivity values;
- WT parent behaviour is unchanged when the paired override and auxiliary current are disabled and the original Task37 state is supplied;
- all inherited conservation diagnostics remain valid.

Run software/fixed-state evaluations only. No production trajectory until 52B is published and independently verified.

Publish `52B` immediately.

## Stage 52C — immutable numerical freeze

After verified 52B and before production:

Freeze machine-readable:

- central projected WT and KO onset states, including exact source values and projection residuals;
- alternate projected states;
- case matrix;
- solver settings;
- `G_aux` values and source interpretation;
- KO supply multipliers;
- quadrature windows;
- required diagnostics and acceptance gates;
- exact model/source/input hashes.

Run a dependency audit proving Task52 candidate code has no Palk/Task41/Task50 scientific dependency.

Create `output/parameter_and_case_freeze.json` and checkpoint `52C`.

No phenotype-informed change after remotely verified 52C.

## Stage 52D — primary combined-stimulus production cases

Run exactly the first three central paired cases under combined CCh+IPR for 600 s:

1. matched supply, central projection, `G_aux=0`;
2. matched supply, central projection, `G_aux=2.32e-9 S`;
3. matched supply, central projection, `G_aux=4.49e-9 S`.

Use one scientific process / one BLAS thread. Preserve complete compact trajectories/diagnostics without printing them to terminal.

For each case record at integer-second resolution or the inherited denser production standard:

- all conserved states;
- pH/speciation;
- water flows and cumulative secretion;
- voltages and reversal potentials;
- TMEM16A and auxiliary chloride currents/fluxes;
- NKCC, AE2, AE4, NHE1, NBC, pump and K fluxes;
- charge/current/carbon/water/speciation residuals.

Compute cumulative/broad-window summaries at the binding phenotype observation windows.

### Exact reservoir audit

For each case compute:

`delta(t)=nCl_WT(t)-nCl_KO(t)`

and the numerical identity

`D_J = delta(0) + 2D_N + D_A + D_E - delta(600)`

where

- `D_J = ∫(J_WT-J_KO)dt`;
- `D_N = ∫(N_WT-N_KO)dt` in NKCC cycle units;
- `D_A = ∫(A_WT-A_KO)dt`;
- `D_E = ∫(E_WT-E_KO)dt`.

The independently quadrature-computed `D_J` must agree with the RHS within a predeclared numerical quadrature tolerance. Do not loosen conservation tolerances to obtain agreement.

Publish 52D before running sensitivities.

## Stage 52E — predeclared sensitivities and controls

Run exactly:

4. central projection, combined CCh+IPR, `G_aux=2.32e-9 S`, KO matched-supply multiplier `1.05`;
5. same with multiplier `1.10`;
6. alternate initial-state projection, matched supply, combined CCh+IPR, `G_aux=2.32e-9 S`;
7. CCh-only, central projection, matched supply, beta=0, auxiliary current therefore exactly zero;
8. IPR-only, central projection, matched supply, `G_aux=2.32e-9 S`.

No other cases.

Report whether the full-model effect is robust to modest non-AE4 compensation and initial-state projection uncertainty.

For IPR-only, do not equate model chloride concentration slope to SPQ fluorescence recovery without an identified observation map. Direction/positive response may be discussed separately from quantitative assay matching.

Publish 52E before final interpretation.

## Stage 52F — final project lock

After verified 52E:

1. Audit every Task52 input, checkpoint and output hash.
2. Re-read the preserved Task50 benchmark and verify its archive ref has not moved.
3. Compare Task52 outputs with the experimental phenotype and Task50 only now; do not retune.
4. Produce final report:

`analysis/52_chloride_reservoir_final_test/CHLORIDE_RESERVOIR_FINAL_REPORT.md`

The report must include:

- architecture and provenance;
- exact measured onset projection values;
- case table;
- full secretion outcomes;
- chloride export outcomes;
- exact reservoir mass-budget decomposition;
- electrochemical driving-force trajectories/summary;
- supply-compensation sensitivity;
- projection sensitivity;
- protocol controls;
- physical/conservation status;
- comparison with Task50;
- explicit limitations;
- final classification.

### Allowed final classifications

Use one of:

**A. SOURCE-CONSTRAINED MECHANISTIC SUCCESS**  
Full model gives a substantial robust AE4-null secretion deficit of the experimental order without genotype-specific channel coupling or phenotype tuning; reservoir accounting closes and physical gates pass.

**B. PARTIAL MECHANISTIC SUPPORT**  
Reservoir mechanism creates a substantial chloride-export/fluid deficit but falls materially short of the experimental order or is sensitive to the fixed projection/supply uncertainty. State exactly what is and is not supported.

**C. FULL-MODEL MECHANISTIC FAILURE**  
The reduced R51I effect does not survive full coupled dynamics at useful magnitude. State which feedback erases it. Task50 remains the quantitative proof that additional unidentified beta-conditioned coupling is required.

Do not manufacture a more favourable category.

Append the controlling result to the mandatory ledger stack, update `CURRENT_STATUS.md`, commit/push, independently verify remote `52F`, and stop.

No further scientific execution after verified 52F.
