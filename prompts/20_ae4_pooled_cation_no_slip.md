# Task 20: AE4 pooled-cation no-slip mechanism test

## Mission

Tasks 14B, 18 and 19 establish that the current physiological AE4 implementation gives the wrong genotype direction even after chloride-load reallocation and Na/K-pump impairment diagnostics. Task 15 showed why the present transporter is suspicious: its independent Na and K branches run in opposite directions and cancel by roughly 97–99% in WT.

Task 20 tests a narrower mechanism hypothesis:

> Does replacing the unsupported independent opposing Na/K AE4 cycles with one reversible pooled-cation net cycle restore a biologically sensible AE4-loss phenotype while retaining the directly supported transport facts?

Work on branch:

`codex/task-20-ae4-pooled-cation-no-slip`

Tracking issue: #19.

Repository inheritance is merged Task 19 on `main`. Do not edit `archive/`. Do not merge to `main` until Task 20 is complete and reviewed.

## Evidence boundary

Use the 2016 Peña-Münzenmayer et al. JGP paper and the existing repository evidence ledger.

The experiments support:

- AE4-dependent Cl/HCO3 exchange;
- Na transport capability;
- K transport capability under imposed gradients;
- broad monovalent-cation support;
- macroscopic electroneutrality;
- reversibility under manipulated gradients.

They do **not** establish that, in a physiological mixed Na/K bath, AE4 contains two independent net cycles that can run in opposite directions simultaneously.

The paper's physiological Figure 10 is explicitly a working model. It assumes the cation pathway does not distinguish strongly between Na and K and reasons using the combined cation pool. The 1:1:2 Cl:cation:HCO3 stoichiometry is a working-model assumption, not a measured microscopic law. Hill values do not determine transported-ion count.

Task 20 must preserve that evidence discipline.

## Legacy model remains available

Do not delete or silently rewrite the current `SR2_SHARED_112_QSS` model.

Implement the pooled-cation model as a separate selectable AE4 mechanism so that:

- all legacy results remain reproducible;
- unit tests can compare legacy and pooled-cation fluxes at identical states;
- Task 20 can be interpreted as a nested mechanism comparison.

Do not alter unrelated production equations.

## New physiological mechanism: one pooled-cation net cycle

The primary Task 20 AE4 model must have exactly **one scalar net cycle flux** `J_AE4`.

There must not be independent physiological `J_Na` and `J_K` cycles whose net signs can oppose one another.

Use the existing 1:1:2 working stoichiometry for this diagnostic:

`Cl_o + C_i + 2 HCO3_i <-> Cl_i + C_o + 2 HCO3_o`

where `C` is a pooled monovalent cation carried by the same net cycle.

Positive `J_AE4` means intracellular chloride loading.

### Pooled thermodynamic affinity

For the primary WT model use equal Na/K selectivity weights, reflecting the 2016 nonselective-cation interpretation:

`w_Na = 1`, `w_K = 1`.

Define effective pooled cation activities

`C_i = w_Na * Na_i + w_K * K_i`

`C_o = w_Na * Na_o + w_K * K_o`.

Use one scalar cycle affinity

`A = log(Cl_o / Cl_i) + log(C_i / C_o) + 2 * log(HCO3_i / HCO3_o)`.

Construct reversible forward and reverse effective rates so their ratio obeys local detailed balance:

`log(k_forward / k_reverse) = A`.

A symmetric barrier gauge such as `k0*exp(+A/2)` and `k0*exp(-A/2)` is acceptable. The absolute attempt-rate/carrier product is an unmeasured kinetic lump and may be absorbed into the fitted AE4 capacity.

Require

`J_AE4 * A >= 0`

within numerical tolerance.

Do not import the legacy branch affinities into the pooled model.

### Na/K bookkeeping without slip

The total transported cation source is tied to the same scalar cycle:

`source_cation_i = -J_AE4`.

Partition that source between Na and K using a bounded **donor-side nonselective availability rule**.

For positive chloride-loading flux (`J_AE4 > 0`), the cation donor side is intracellular:

`p_Na = Na_i / (Na_i + K_i)`

`p_K = 1 - p_Na`.

For reverse flux (`J_AE4 < 0`), use the extracellular donor pool:

`p_Na = Na_o / (Na_o + K_o)`

`p_K = 1 - p_Na`.

At exactly zero net flux the cation sources are zero; the partition itself is irrelevant.

Then

`source_Na_i = -p_Na * J_AE4`

`source_K_i = -p_K * J_AE4`.

Thus Na and K can both be transported, but their **net** transported amounts cannot run in opposite directions through AE4 at the same state.

If a smooth donor-side interpolation is required solely for numerical differentiation near `J_AE4 = 0`, use the smallest smooth approximation that converges to the same two donor-side limits. Document it. Do not reintroduce an independent slip degree of freedom.

### Remaining source vector

Keep the existing working stoichiometric source coupling:

- `source_Cl_i = +J_AE4`;
- `source_HCO3_i = -2 J_AE4`;
- `source_TIC_i = -2 J_AE4`;
- `source_alkalinity_i = -2 J_AE4`.

Verify electroneutrality exactly.

### Regulation

The existing beta/cAMP/PKA `regulation_gain` must scale active AE4 capacity only, as before.

Do not change regulation kinetics or calcium inputs.

## Assay-level tests before whole-cell work

Before any WT rebalancing or genotype simulation, add focused tests independent of phenotype outcome.

Required tests:

1. global equilibrium gives zero pooled AE4 flux;
2. reversing the complete thermodynamic gradient reverses `J_AE4`;
3. positive `J_AE4` gives Cl loading and cation/HCO3 efflux under the declared working cycle;
4. Na and K intracellular sources always have the same sign or one is zero;
5. `abs(source_Na_i) + abs(source_K_i) = abs(J_AE4)`;
6. transported charge is zero;
7. `J_AE4 * A >= 0`;
8. pure/near-pure Na availability produces Na-supported transport;
9. pure/near-pure K availability produces K-supported transport;
10. removing one cation renormalizes the carried-cation partition rather than disabling the entire exchanger;
11. the pooled model cannot generate an independent Na/K exchange loop with zero net Cl flux;
12. legacy `SR2_SHARED_112_QSS` outputs remain unchanged when that model is selected.

The tests must not use the AE4-knockout secretion phenotype.

## Whole-cell comparison design

Use all ten inherited accepted WT conserved states as fixed physiological targets.

Do not discard the states simply because they originated from two legacy routing families. The new pooled model has one primary equal-selectivity cation rule; retain the original root labels only for provenance and robustness accounting.

Use exactly **two new-architecture WT loading conditions** to control computational cost:

### Condition P0: inherited AE4 chloride flux

At each inherited WT root, require the pooled model's net AE4 chloride source to equal the inherited legacy WT net AE4 chloride source at that state.

This isolates the effect of changing AE4 cation architecture while keeping AE4 productive chloride loading unchanged.

### Condition P10: AE4 share 0.10

At each inherited WT root, require pooled AE4 to carry exactly 10% of the inherited positive basolateral chloride-loading pool.

Reduce the other positive basolateral chloride loaders proportionally so the total positive loading pool remains equal to the inherited WT value.

Do **not** run a new 0.30 condition in Task 20.

Legacy comparator results for inherited allocation and legacy 0.10 share already exist and should be reused rather than recomputed when hashes match.

## Fixed-WT inverse rebalancing

For both P0 and P10, keep every inherited WT conserved state exactly fixed.

Use the Task 18 inverse-rebalancing philosophy:

- solve for existing uncertain conductances/capacities needed to restore the complete WT steady-state equations;
- allow total K conductance, CaCC, pumps, NHE1, CO2 exchange and paracellular capacities to move as reference parameters rather than immutable laws;
- preserve transporter topology, bath, chemistry, geometry and all hard conservation laws;
- require nonnegative physical capacities and valid membrane fractions;
- allow algebraic voltages to reclose self-consistently;
- do not use genotype data in feasibility or optimization.

Among multiple exact WT solutions, use the same Task 18 canonical minimal-deviation objective:

1. minimize the largest absolute log-fold change from inherited capacities;
2. subject to that, minimize the sum of squared log-fold changes.

Record all parameter changes.

## WT freeze before genotype evaluation

Before any new AE4-loss or AE2-loss result is accessed:

1. verify every fixed WT conserved state byte-for-byte;
2. verify full resting RHS, charge, current, carbon, alkalinity and water closure;
3. verify the pooled AE4 target flux/share;
4. verify absence of opposing Na/K AE4 net sources;
5. save all feasible and infeasible P0/P10 decisions;
6. write and hash the complete WT candidate manifest;
7. commit and push the WT checkpoint.

No genotype result may revise a WT parameterization after this freeze.

## WT dynamics

For every frozen feasible pooled-cation WT model, run the inherited CCh + IPR protocol at calcium

- 0.10 uM;
- 0.25 uM;
- 0.50 uM.

Use production Radau and inherited conservation/numerical gates.

Absolute gland-scale mapping remains nonblocking.

Do not use AE4-loss time-course shape to accept or reject WT models.

## Post-freeze genotype evaluation

Only after the WT checkpoint is pushed:

### AE4 near-loss

Reduce AE4 expression from 1.0 to 0.05 with all other Task 20 parameters frozen.

Do not refit after AE4 reduction.

Do not require exact zero.

For every valid connected resting state, run 600 s dynamics at all three calcium values and report

`R_AE4 = total_AE4_0.05 / total_WT`.

### AE2 loss

Evaluate matched AE2 deletion with all Task 20 parameters frozen and report

`R_AE2 = total_AE2_loss / total_WT`.

### Secondary native-state checks

Without fitting to them, report whether AE4 near-loss moves:

- resting intracellular Cl in the experimentally observed lower direction;
- resting pH substantially or not;
- intracellular Na and K;
- NKCC1, NHE1, AE2 and pump fluxes.

The experimental approximately 35% secretion reduction is context only after results are frozen. It is not a target or acceptance band.

## Primary scientific comparisons

Compare four conceptual conditions using new results plus reused legacy results:

1. legacy architecture at inherited AE4 loading;
2. pooled-cation architecture at matched inherited AE4 loading (P0);
3. legacy architecture at 10% AE4 loading;
4. pooled-cation architecture at 10% AE4 loading (P10).

The key question is whether removing the independent Na/K slip architecture changes the direction of the AE4-loss phenotype.

Do not claim success based on one root or one calcium value.

## Interpretation discipline

A successful mechanism signal requires a robust decrease in secretion after AE4 near-loss across the pooled-cation ensemble, not merely a numerical difference from the legacy model.

If P0 remains wrong-direction but P10 becomes right-direction, conclude that both architecture and productive loading allocation matter.

If both P0 and P10 remain wrong-direction, the no-slip pooled-cation replacement is insufficient.

If P0 alone flips direction, architecture is sufficient without increasing AE4's WT chloride share.

Do not infer that the pooled model is the unique true microscopic mechanism. It is a deliberately less committal coarse-grained alternative consistent with the experimental uncertainty.

## Required classifications

Choose exactly one primary classification:

1. `POOLED-CATION AE4 RESTORES ROBUST AE4-LOSS SECRETION DEFICIT`
2. `POOLED-CATION AE4 RESTORES DEFICIT ONLY WITH 10-PERCENT WT LOADING`
3. `POOLED-CATION AE4 REMAINS WRONG-DIRECTION`
4. `POOLED-CATION AE4 GIVES MIXED ROOT OR CALCIUM DEPENDENCE`
5. `POOLED-CATION FIXED-WT REBALANCING REQUIRES EXTREME SHIFTS`
6. `TASK 20 NUMERICALLY INCONCLUSIVE`

## Required artifacts

Create at minimum:

- `results/20_ae4_pooled_cation/contract.json`
- `results/20_ae4_pooled_cation/assay_tests.json`
- `results/20_ae4_pooled_cation/legacy_vs_pooled_reference.csv`
- `results/20_ae4_pooled_cation/wt_inverse_solutions.csv`
- `results/20_ae4_pooled_cation/wt_parameter_payloads/`
- `results/20_ae4_pooled_cation/frozen_wt_manifest.json`
- `results/20_ae4_pooled_cation/wt_dynamic_results.csv`
- `results/20_ae4_pooled_cation/ae4_5pct_results.csv`
- `results/20_ae4_pooled_cation/ae2_results.csv`
- `results/20_ae4_pooled_cation/state_flux_summary.csv`
- `results/20_ae4_pooled_cation/final_comparison.csv`
- `results/20_ae4_pooled_cation/test_results.txt`
- `analysis/20_ae4_pooled_cation/evidence_boundary.md`
- `analysis/20_ae4_pooled_cation/mechanism.md`
- `analysis/20_ae4_pooled_cation/wt_rebalancing.md`
- `analysis/20_ae4_pooled_cation/genotype_interpretation.md`
- `analysis/20_ae4_pooled_cation/final_answer.md`

The final answer must state plainly:

- whether the new AE4 law eliminates opposing Na/K net transport by construction and in every simulated state;
- whether all ten inherited WT states can be supported under P0 and P10;
- how much capacity rebalancing is required;
- whether AE4 5% lowers or raises secretion under P0 and P10;
- whether resting Cl moves in the observed direction;
- whether AE2 remains comparatively neutral;
- whether results are robust across all roots and calcium values;
- whether the pooled-cation mechanism is sufficient or further AE4 mechanism work is required.

## Final actions

Run focused and inherited tests needed to prove the new mechanism is isolated and the legacy model remains reproducible.

Commit and push completed Task 20 work to `codex/task-20-ae4-pooled-cation-no-slip`.

Do not edit `archive/`.

Do not merge to `main` automatically. Open a PR after completion and report the result for review before merge.
