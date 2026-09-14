# Task 30: replace the coarse NHE1 law with a transporter-specific model

Repository: `esig626/ae4-salivary-transport-control`
Branch: `codex/task-30-nhe1-model-repair`

Read `AGENTS.md` first. Use deep reasoning, but this is an implementation task, not another diagnosis project.

## Objective

The current generic reversible `tanh` NHE1 law drives too much proton extrusion near physiological pH. Task 29 showed that exact-null acid/base closure then requires intracellular Na far above the retained WT states, and R10 subsequently demands NKCC1 influx above its frozen capacity.

Repair NHE1 itself.

Implement the transporter-specific NHE1 law used in Vera-Sigüenza et al. (2018), Eq. 26, as the active Task 30 NHE1 model. Keep the existing generic law selectable as a historical comparator, but do not use it for the repaired model.

Do not modify AE4, AE2, NKCC1, pump, K, CaCC, CO2, paracellular or water laws in this task.

## Published NHE1 law

The 2018 model gives

    J_NHE1 = G_NHE1 * [
        (H_i / (H_i + K_H))**2 * (Na_e / (Na_e + K_Na))
        - (Na_i / (Na_i + K_Na)) * (H_e / (H_e + K_H))**2
    ]

with positive `J_NHE1` meaning Na entry to the cell and H extrusion from the cell.

Use the published values:

- `G_NHE1 = 0.0305 fmol/s`
- `K_H = 4.5e-4 mM`
- `K_Na = 15 mM`

Source: Vera-Sigüenza et al., *A Mathematical Model Supports a Key Role for Ae4 (Slc4a9) in Salivary Gland Secretion*, Bulletin of Mathematical Biology 80, 255-282 (2018), Eq. 26, PMCID PMC5792321. The paper states that the functional form is based on Falkenberg and Jakobsson (2010); `G_NHE1` is model-derived.

Use the equation as printed in the 2018 paper, including the squared proton saturation factors. Do not silently substitute the historical archived implementation, which used a different proton exponent.

## 1. Implementation

Add a clean NHE1 model selection in the homeostasis layer.

Requirements:

- Preserve the existing generic thermodynamic/tanh NHE1 evaluator for legacy reproduction.
- Add a named published-2018 NHE1 evaluator with the equation and parameters above.
- The Task 30 model must explicitly select the published-2018 evaluator.
- Do not add the Task 28 extra pH gate on top of the 2018 law.
- Preserve 1 Na inward : 1 H outward stoichiometry.
- Positive NHE1 contributes `+J` to intracellular Na amount and `+J` to total alkalinity, with zero direct TIC source and zero net transported charge.
- Do not make NHE1 genotype dependent.
- Keep units explicit. `H_i` and `H_e` in the equation are mM, consistent with `K_H`.
- Maintain backward compatibility with the clean seed parameter payloads. New published-law parameters may be added with explicit provenance, but loading the inherited roots must not break.

Write focused unit tests before whole-cell runs. Independently calculate at least three fixed-state values of Eq. 26 and require machine-precision agreement. Also test source stoichiometry, finite bounds, increased forward activity under acidification at fixed Na, and genotype independence.

## 2. No model-family search

This task implements ONE repair law: published 2018 Eq. 26.

Do not compare a menu of NHE1 models, do not fit Hill exponents, do not invent a second family if this one fails, and do not perform parameter sweeps.

The only optional adjustment allowed is a single WT-only scalar adjustment of `G_NHE1` for a background if and only if the published value `0.0305 fmol/s` cannot produce an admissible WT rest. If needed:

- vary only `G_NHE1`;
- use no knockout information;
- target WT resting pH 6.91 only, while requiring full stationary closure;
- use at most four scalar `G_NHE1` evaluations for that background;
- do not tune WT chloride or secretion;
- freeze the resulting `G_NHE1` before exact-null evaluation;
- report clearly that this is a WT calibration, not a published value.

Do not invoke this calibration merely because the re-equilibrated WT differs slightly from the saved state. First attempt the published value exactly.

## 3. Whole-cell repair test

Use only R09 and R10.

For each background:

1. Construct the routed whole-cell model with the published-2018 NHE1 law and all other inherited parameters unchanged.
2. Re-equilibrate WT REST directly from the saved WT state.
3. Verify the full amount RHS, current closure, charge, carbon, positivity and volume gates.
4. Freeze the WT law/parameters.
5. Set only AE4 expression from 1.0 to exact 0.0.
6. Solve exact-null REST starting from the repaired WT rest.
7. If needed, permit exactly one additional null warm start biased toward the experimental resting observations (`Cl_i = 36.5 mM`, `pH = 6.89`) while preserving positive amounts and the inherited model constraints. This is a solver start, not a target fit.
8. Do not alter any parameter after seeing the null result.

Report WT and null:

- Na_i, K_i, Cl_i, TIC_i, HCO3_i, pH_i, cell volume;
- lumen concentrations and volume;
- NHE1, AE2, AE4 and NKCC1 fluxes;
- pump, K and CaCC fluxes;
- CO2 fluxes;
- membrane voltages;
- osmolarities and water fluxes;
- all numerical/conservation residuals.

The experimental AE4-null pH `6.89 +/- 0.02` and Cl `36.50 +/- 1.60 mM` are validation context, not fit targets.

For project-level interpretation, explicitly flag whether the repaired null avoids the previous pathological state:

- intracellular Na should remain of ordinary cellular magnitude rather than the Task 29 90-230+ mM requirement; use `Na_i <= 30 mM` as a declared diagnostic screen, not as a claimed measured confidence interval;
- cell volume must satisfy the inherited `<= 3 pL` gate;
- pH should remain near the observed value; use `|pH_i - 6.89| <= 0.10` as a diagnostic success band, not as a fitted objective;
- no transporter may require flux above its implemented finite capacity.

Do not reject a numerically valid state solely because chloride does not yet reproduce 36.5 mM. Chloride allocation is a separate later problem.

## 4. Small dynamic confirmation only after REST repair

If an admissible exact-null REST is obtained for a background, run exactly one matched WT/null production stimulation at Ca = 0.25 uM using the inherited 600 s CCh + IPR protocol and production tolerances.

Maximum dynamic runs: four total, one WT and one null for each of R09/R10.

This dynamic check answers only whether the repaired NHE1 model remains numerically and physiologically stable under stimulation and what secretion direction emerges without further tuning.

Report cumulative secretion and the actual AE4/NKCC1/AE2 chloride fluxes, but DO NOT tune anything to the known ~35% secretion reduction. Do not run Ca 0.10 or 0.50 in this task.

If a null REST is inadmissible, do not run its trajectory and do not convert missing flow to zero.

## 5. Compute limits

Hard global limits:

- 12 stationary solver calls total;
- 10,000 stationary residual evaluations total;
- 8 scalar `G_NHE1` evaluations total across both backgrounds;
- 4 stimulated integrations total;
- one numerical worker;
- one numerical library thread;
- 30 minutes numerical execution.

Stop at the first exhausted limit and preserve partial results. No automatic expansion.

Forbidden: combinatorial searches, allocation searches, grids, random/global optimisation, phenotype fitting, multi-parameter calibration, subagents, Shapley analysis, or changing another transporter because NHE1 repair fails.

## 6. Deliverables

Write compact UTF-8 outputs under:

- `analysis/30_nhe1_model_repair/`
- `results/30_nhe1_model_repair/`

Required:

- `final_answer.md`
- `nhe1_law_comparison.md`
- `resting_states.csv`
- `flux_ledger.csv`
- `matched_secretion.csv` if any trajectories are valid
- `contract.json`
- `verification.json`
- `budget.json`
- focused test log

The final answer must state plainly:

1. Was the published NHE1 law implemented correctly?
2. Did WT REST remain/recover physiological?
3. Did exact AE4-null REST become admissible?
4. Did Na_i remain ordinary rather than rising to the Task 29 extreme values?
5. Did pH stay near the observed knockout pH?
6. What did AE2 and NKCC1 do after AE4 removal?
7. If dynamics ran, what secretion direction emerged without tuning?
8. What remains for the later chloride-allocation phase?

## 7. Git/publication

Test ordinary Git authentication once at task start. If it works, commit and push this branch once after completion. Do not merge to main.

If authentication fails, continue the scientific task, preserve the local text outputs, and report the failure. Do not retry authentication, construct manual Git objects, encode Base64, create archives, or use alternate publication machinery.