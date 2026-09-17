# Mandatory research ledger continuation — Task52 start

**Binding status:** append-only continuation of the project research ledger. Task52 must read this file together with, in order:

1. `docs/MANDATORY_RESEARCH_LEDGER.md`
2. `docs/MANDATORY_RESEARCH_LEDGER_PRE_TASK52_APPENDIX.md`
3. `docs/MANDATORY_RESEARCH_LEDGER_R51H.md`
4. `docs/MANDATORY_RESEARCH_LEDGER_R51I.md`
5. `analysis/pre_task52_algebraic_gate/CHLORIDE_RESERVOIR_SOLUTION.md`
6. `analysis/pre_task52_algebraic_gate/REDUCED_RESERVOIR_SIMULATION.md`

The repository, not conversational memory, is authoritative.

## R52A — final numerical shot frozen architecture

**Budget:** exactly one funded/scientific Codex shot remains and Task52 consumes it. No later mechanism search is authorised.

**Scientific question:**

> In the full conservation-explicit pre-Palk model, can the experimentally measured chronic AE4-KO chloride depletion plus loss of ordinary beta-activated AE4 supply generate a substantial stimulated secretion deficit when non-AE4 chloride supply and apical secretory conductances are shared across genotypes?

The Task50 multiplier is a frozen benchmark only. It must not enter the Task52 equations.

### 1. Scientific parent

Task52 must use the **pre-Palk Task37 lineage** as its transport chassis. Palk/Benjamin is excluded from the candidate model.

Do not import `src/modern_full_model/nkcc1_palk2010.py` into the Task52 model dependency chain.

Historical Palk Tasks39–51 remain valid negative diagnostics showing how state-driven NKCC compensation erased AE4 loss; they are not the final candidate architecture.

### 2. Measured genotype onset states

Primary source values from Peña-Münzenmayer et al. JBC 2015:

- WT resting `[Cl]_i = 50.10 +/- 1.50 mM`, pH `6.91 +/- 0.07`;
- AE4 KO resting `[Cl]_i = 36.50 +/- 1.60 mM`, pH `6.89 +/- 0.02`.

The old shared-WT-rest acute perturbation is forbidden as the central Task52 initialisation because it deletes the measured chronic chloride reservoir difference.

The central initial-state projection is deterministic and phenotype-free:

- start from the accepted Task37 WT resting state and active Task37 parameter set;
- retain its intracellular Na concentration and cell volume for both genotypes;
- impose the measured central genotype `[Cl]_i` and pH values above;
- solve only intracellular K and TIC so that (i) exact cell bulk electroneutrality and (ii) cell osmolarity equal to the active bath osmolarity both hold;
- derive TA from the existing acid-base speciation at the imposed pH/TIC and unchanged finite buffer;
- keep fixed-cell charge, fixed osmoles, geometry and lumen state unchanged;
- require positive K/TIC, exact charge closure, osmotic closure, normal speciation, and inherited broad physiological ranges before any trajectory.

This is an initial-state projection, not a resting-state model or stationary solve. Task52 must not claim that the parent transport laws maintain these chronic states before stimulation.

One predeclared initialisation sensitivity is allowed, fixed before production: retain Task37 intracellular Na and TIC instead, and solve K plus cell volume from the same exact charge and osmotic constraints. It is a projection sensitivity only, not a second mechanism.

No optimisation against secretion is allowed in either projection.

### 3. Shared non-AE4 chloride supply

The primary 2015 isolated NKCC assay did not detect increased functional NKCC activity in AE4 KO and the authors did not attribute the phenotype to NKCC up-regulation. This does not prove exact physiological flux equality.

Task52 therefore uses a transparent matched-supply experiment:

- evaluate the inherited **pre-Palk** NKCC and AE2 fluxes on the WT state;
- feed those same instantaneous NKCC and AE2 cycle fluxes, with their exact stoichiometric source vectors, to the matched KO system at the same time;
- WT retains its own ordinary pre-Palk NKCC/AE2 evaluation;
- NHE1, NBC, pump, K channels, CO2, water, lumen, voltage and all other laws remain genotype-local and fully dynamic;
- WT retains ordinary AE4 transport/regulation;
- KO has exact zero AE4 expression.

The central model therefore asks specifically whether the measured chloride reservoir plus AE4 supply is sufficient **without large state-driven KO chloride compensation**.

This paired shared-supply construction is an explicit experimental-control idealisation. It is not a claim that native WT and KO NKCC fluxes are exactly equal at every instant.

Predeclared compensation sensitivities, using the same WT supply trajectory and no fitting:

- KO NKCC+AE2 supply multiplier `1.05`;
- KO NKCC+AE2 supply multiplier `1.10`.

Do not run larger rescue values or tune this multiplier to phenotype.

### 4. Shared beta/IPR auxiliary anion demand

Catalán et al. PNAS 2015 independently establish an IPR/beta-associated TMEM16A-independent, DCPIB/NPPB-sensitive chloride conductance. Task52 treats this as a **shared effective apical demand pathway**, not as an AE4-specific gate and not as a claim of direct molecular cAMP gating.

Use the same effective law in WT and KO:

`I_aux = beta * G_aux * (V_a - E_Cl)`.

It must participate in the exact apical electrical closure and cell/lumen chloride amount balances. It contains no AE4 expression term and no water-output multiplier.

No swelling gate is used in Task52. The experimentally observed IPR-induced conductance is represented directly at protocol level because Task51 showed that trying to infer conductance from swelling is sign-inappropriate: outward chloride current reduces swelling.

Predeclared fixed `G_aux` values only:

- `0 S` — no-auxiliary control;
- `2.32e-9 S` — approximate DCPIB-sensitive whole-cell equivalent from the recovered Catalán SI/source extraction;
- `4.49e-9 S` — approximate total IPR-induced whole-cell equivalent, an explicit upper sensitivity.

These are source-scale effective sensitivities, not identified apical conductances. The unknown apical fraction and rectification limitation must remain explicit. Do not search or fit `G_aux`.

### 5. Primary and control protocols

Primary numerical test: combined CCh+IPR, 600 s, using the inherited Task37 stimulus doses/inputs except for the added shared `G_aux` term.

Required central cases:

1. matched-supply WT/KO, central initial-state projection, `G_aux=0`;
2. matched-supply WT/KO, central projection, `G_aux=2.32e-9 S`;
3. matched-supply WT/KO, central projection, `G_aux=4.49e-9 S`.

Required fixed sensitivities:

4. central projection, `G_aux=2.32e-9 S`, KO shared supply x`1.05`;
5. central projection, `G_aux=2.32e-9 S`, KO shared supply x`1.10`;
6. alternate charge/osmotic initial-state projection, matched supply, `G_aux=2.32e-9 S`.

Required protocol controls at the central projection / matched supply:

7. CCh-only, beta=0, auxiliary term therefore exactly off;
8. IPR-only with `G_aux=2.32e-9 S`.

Do not invent a 5%-AE4 chronic initial chloride state. The primary Task52 mechanistic test is WT versus complete AE4 KO. The Task50 5% result remains a benchmark only.

### 6. Frozen equations and diagnostics

For each paired production case report at minimum:

- cumulative secretion 0–600 s and broad-window 60–600 s mean flow;
- endpoint flow;
- cumulative total apical chloride export, with TMEM16A and auxiliary components separately;
- WT/KO `[Cl]_i`, Cl amount, `E_Cl`, `V_a` and `E_Cl-V_a` over time;
- NKCC cycles, AE2 cycles and their exact WT/KO imposed difference;
- AE4 loading;
- NHE1, NBC and pump cycles;
- Na, K, TIC, TA, pH, volume and lumen states;
- charge/current/carbon/water/speciation conservation diagnostics.

The exact reservoir accounting must be reported:

`delta(t) = n_Cl,WT(t) - n_Cl,KO(t)`

and

`Integral(J_WT-J_KO) dt = delta(0) + 2 Integral(N_WT-N_KO) dt + Integral(A_WT-A_KO) dt + Integral(E_WT-E_KO) dt - delta(T)`.

For the central matched-supply cases, the imposed NKCC/AE2 differences should be numerically zero, so the identity reduces to the R51H/R51I reservoir form up to numerical integration error.

No result may be called a reservoir mechanism success unless this accounting closes and the deficit is not created by a hidden genotype-dependent channel multiplier.

### 7. Success / failure interpretation

A successful mechanistic result requires all of:

- substantial AE4-null cumulative secretion reduction of the experimentally observed order, judged under `docs/PHENOTYPE_TARGET_CONVENTION.md` rather than exact 35.000% fitting;
- persistence over a broad stimulation window;
- correct physical/conservation checks;
- no Task50 multiplier;
- identical apical conductance laws across genotypes;
- no phenotype-tuned KO supply;
- reservoir identity closes quantitatively.

The reduced R51I diagnostic found ~`21.16–25.64%` less KO apical chloride export across the unfitted voltage-closure interval, central `23.5147%`, and discharged ~81% of the available reservoir-plus-AE4 advantage. Task52 tests whether the full coupled water/ion model preserves enough of that effect to generate the fluid phenotype.

If the full result is weak, state why: remaining reservoir storage, electrical accommodation, NHE/NBC/pump feedback, water/lumen feedback, or modest allowed supply compensation. Do not rescue it with another mechanism.

If Task52 fails, the final project conclusion is that Task50 remains the quantitative proof-of-sufficiency for an additional unidentified beta-conditioned network coupling.

### 8. Publication rule

Task52 must publish incrementally and stop after its final verified report. Suggested checkpoints:

- `52A`: ledger/source/architecture and exact case matrix freeze; no production trajectories;
- `52B`: initial-state projector + paired shared-supply/auxiliary-current implementation tests;
- `52C`: source-scale parameter/case freeze and exact no-Palk/Task50 dependency audit;
- `52D`: central combined-stimulus production cases;
- `52E`: predeclared sensitivities and protocol controls;
- `52F`: final comparison, uncertainty/limitations, article-facing lock and ledger update.

After each checkpoint: update status/decision/ledger if warranted, commit, push, independently verify the remote SHA, then proceed.

No Task53. Stop after remotely verified 52F.
