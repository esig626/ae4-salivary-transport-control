# Task 18: fixed-WT inverse transporter rebalancing

## Mission

Task 17 kept each inherited WT state fixed and asked whether the chloride-loading burden could be shifted toward AE4 while preserving the rest of the model under inherited capacity restrictions. It found that the 0.10 and 0.30 AE4-share targets were blocked by K-source constraints, especially the inherited 14 nS total K conductance assumption.

Task 18 asks the intended scientific question without treating inherited calibration numbers as physical laws:

> Can the exact same accepted WT conserved state be supported when AE4 carries a larger fraction of positive basolateral chloride loading, provided the existing uncertain transporter capacities and conductances are allowed to rebalance?

This is a fixed-WT inverse problem. It is not a new WT state search and it is not a fit to the AE4-loss phenotype.

Work on branch:

`codex/task-18-fixed-wt-inverse-rebalance`

Tracking issue: #13.

Use completed Task 17 commit:

`fb2afcf3e23af1d19b297280ad67953370f7c9a9`

as the repository inheritance. The unchanged scientific model baseline remains completed Task 14B commit:

`4d4403f279fc36aec940c9e4759486d02f07e3f6`.

Do not edit `archive/`. Do not merge to `main`.

## Central fixed-state rule

For every one of the ten accepted WT roots, copy the complete inherited conserved WT state exactly.

Do not optimize, continue, round, reconstruct, or otherwise change any WT conserved state coordinate.

The fixed quantities include all intracellular and luminal conserved amounts and both volumes, plus the inherited resting regulatory state where applicable.

The intervention is a reparameterization of the flux decomposition that supports the same WT state.

Algebraic membrane voltages may reclose self-consistently if required by the modified conductances and pumps. They are not conserved state coordinates. Report their deviations from the inherited values, but do not reject a solution merely because an old voltage target or historical-lineage value moved. Reject only if the production electrical closure fails, the voltage is nonfinite, or the resulting parameterization is physically undefined.

## Exactly three chloride-allocation conditions

Use exactly:

1. inherited baseline;
2. AE4 positive basolateral chloride-loading share `0.10`;
3. AE4 positive basolateral chloride-loading share `0.30`.

No other share values, midpoint refinements, threshold searches, or phenotype-guided values are allowed.

For each root, calculate the inherited positive basolateral chloride-loading pool

`Lplus = sum of all positive basolateral intracellular Cl sources`.

At the inherited roots this is expected to contain AE4 and NKCC1, but determine the set from the actual signed production flux ledger rather than assuming it.

For target share `f`, require

`J_AE4_Cl = f * Lplus`

and require all other positive basolateral chloride loaders to share the remaining

`(1-f) * Lplus`

in proportion to their inherited positive contributions.

Thus the total positive basolateral chloride-loading pool is unchanged at the fixed WT state.

Scale AE4 through the complete coupled transporter. Never scale the chloride source alone. Preserve the AE4 Na/K/Cl/HCO3/TIC/alkalinity source stoichiometry, routing, affinities, reversibility and thermodynamic ratios.

Scale each other positive chloride loader through its complete coupled transporter law.

Signed negative chloride counterfluxes such as AE2 are not part of the positive-loading denominator and are not changed merely to manufacture the target share.

## Preserve the chloride balance, not old conductance numbers

The purpose of this task is to alter who carries the WT chloride load, not to change the WT chloride state.

At the exact fixed WT state, the complete chloride RHS must remain equal to the inherited value, zero within inherited tolerance.

The actual WT fluxes of chloride pathways outside the positive basolateral loading pool should remain at their inherited values wherever this can be achieved by adjusting their capacity under the reclosed voltage. This applies in particular to AE2, apical CaCC and paracellular Cl.

If maintaining an identical nonloading chloride flux is algebraically incompatible with electrical closure, retain the exact total chloride balance as the hard requirement and report the smallest unavoidable redistribution among those nonloading chloride pathways. Do not use genotype information to choose that redistribution.

## Hard physics versus reference parameters

Task 18 must explicitly separate structural constraints from inherited parameter values.

### Hard constraints

Treat as hard only:

- exact fixed WT conserved state;
- complete production steady-state RHS equal to zero within inherited numerical tolerances;
- exact target AE4 positive-loading share for the declared condition;
- preservation of the total positive basolateral chloride-loading pool;
- transporter stoichiometry and source signs;
- existing transporter topology;
- AE4 Na/K routing and thermodynamic ratios;
- bath composition;
- acid-base chemistry and buffer bookkeeping;
- geometry;
- charge and current closure;
- water and lumen accounting;
- nonnegative physical capacities and permeabilities;
- membrane fractions in `[0,1]` where a fraction parameterization is used;
- finite model quantities.

### Reference values, not laws

Do **not** use any of the following as an automatic feasibility wall:

- inherited total K conductance 14 nS;
- inherited CaCC maximum conductance 31.4 nS;
- old NHE1 calibration bounds;
- old AE4 carrier calibration boxes;
- old common-membrane-conductance boxes;
- previous pump capacity values;
- previous paracellular capacity values;
- previous CO2-exchange capacity values;
- any conductance or capacity imported from parotid or another preparation and then carried into the SMG model as a modeling assumption.

Those values are reference points and provenance diagnostics. They may inform how surprising a solution is, but they must not reject an otherwise exact fixed-WT solution.

If the repository contains a genuinely direct same-system measurement that should constrain a capacity, report it separately. Do not silently convert it into a hard bound unless the evidence ledger already establishes a measurement interval that is directly applicable to the modeled SMG preparation.

## Adjustable parameter set

Before solving, construct a parameter provenance table for every existing capacity-like parameter that can affect the resting balances at fixed concentrations.

At minimum consider the current production parameters controlling:

- total K conductance and its apical/basolateral allocation;
- CaCC conductance/capacity;
- apical and basolateral Na/K pump capacities and their allocation;
- NHE1 capacity;
- neutral CO2 exchange capacities;
- non-Cl paracellular Na, K and HCO3 capacities/permeabilities;
- paracellular Cl capacity if needed only to preserve the fixed WT chloride balance under voltage reclosure;
- any other already-existing active capacity multiplier needed to obtain full source rank.

Do not activate a transporter or background pathway that is absent from the existing topology. Do not invent a new mechanism.

Keep AE2 mechanism and genotype identity fixed during the WT inverse solve. AE2 is evaluated only after the WT checkpoint.

The target-share rule determines the AE4 carrier scaling and the proportional scaling of the other positive basolateral Cl loaders. Those are not free phenotype-fitting axes.

## Stage A: source and rank audit

At each inherited WT root:

1. reproduce the complete resting production flux ledger;
2. identify the positive basolateral chloride-loading pool;
3. construct the fixed-state source/sensitivity matrix for all eligible capacity parameters;
4. include the electrical closure equations and any required algebraic voltage variables;
5. report rank, nullity and structurally dependent rows;
6. demonstrate that the inherited baseline parameterization reproduces the exact WT state.

Do not declare infeasibility from one chosen basis. Preserve the full nullspace of the eligible parameter set.

## Stage B: solve the fixed-WT inverse problem

For each root and target share 0.10 and 0.30, solve simultaneously for a physically valid parameterization satisfying all hard constraints.

Because capacities span orders of magnitude, parameterize strictly positive capacities in log space where practical. Fractions should be represented in a way that respects `[0,1]` without clipping accepted solutions.

Do not impose inherited upper bounds on adjustable capacities.

Use robust multistart or an equivalent global-enough strategy appropriate to the problem dimension. If the fixed-state equations reduce exactly to an affine or convex feasibility problem after a suitable flux parameterization, exploit that structure rather than launching an expensive nonlinear search for sport.

A failed local optimization is not an infeasibility certificate.

If no solution is found, distinguish:

- structural rank incompatibility;
- positivity/sign incompatibility;
- true algebraic infeasibility under the existing topology;
- numerical nonconvergence.

Only call a target structurally infeasible if the conclusion no longer depends on inherited numerical capacity bounds such as 14 nS or 31.4 nS.

## Stage C: canonical minimal-deviation solution

When the fixed-WT constraints admit multiple solutions, choose one canonical solution for downstream simulation using a predeclared parameter-only objective that does not contain genotype information.

Use this lexicographic rule:

1. minimize the largest absolute log-fold change among adjustable positive capacities relative to the inherited parameterization;
2. subject to that optimum, minimize the sum of squared log-fold changes;
3. for partition variables, use the corresponding independent positive membrane capacities rather than an arbitrary penalty on a fraction itself whenever possible.

This gives the smallest overall rebalancing of the existing network without treating any inherited capacity as sacred.

For any baseline-zero parameter that cannot be represented by a log ratio, do not activate it unless it is already an existing active degree of freedom with a declared nonzero reference scale. Do not create hidden background currents.

Record every capacity multiplier and the resulting algebraic voltages.

## Stage D: WT verification and phenotype firewall

Before reading or using any AE4-loss or AE2-loss result:

- verify the exact inherited WT conserved state byte-for-byte;
- evaluate the full production resting RHS under the new parameterization;
- verify charge, current, carbon, alkalinity, water and lumen accounting;
- verify the target chloride share and the preserved positive loading pool;
- record the actual signed flux of every Cl pathway;
- record the complete capacity changes and voltage changes;
- retain all ten roots and both routing families;
- save every feasible and infeasible decision.

Write a frozen WT candidate manifest and commit and push it before genotype evaluation.

The known approximately 35% AE4-loss secretion reduction and the known Task 14B wrong-direction phenotype must not enter the WT inverse objective, feasibility constraints, root selection, parameter selection or ranking.

## Stage E: modified WT dynamics

For each frozen feasible modified parameterization, run the WT CCh + IPR dynamic protocol at calcium 0.10, 0.25 and 0.50 uM.

Use production Radau and inherited conservation/numerical checks. Reuse hash-valid baseline trajectories for the inherited condition.

Do not reject a modified model because its absolute whole-gland scale differs from 9 to 10 uL/min. The gland mapping remains nonblocking.

Do not impose any AE4-loss minute-by-minute time-shape criterion on WT.

If a modified WT trajectory is numerically valid but differs dynamically from the inherited WT trajectory, report the change. Do not retroactively alter the fixed-WT parameterization based on genotype outcomes.

## Stage F: post-freeze genotype evaluation

Only after the frozen WT checkpoint is pushed:

### AE4 near-loss

For every feasible modified model, continue AE4 expression from 1.0 to 0.05 with all Task 18 capacities otherwise frozen.

Do not refit or compensate after reducing AE4.

Do not require exact AE4 zero.

Run matched 600 s production trajectories at calcium 0.10, 0.25 and 0.50 uM and compute total secretion ratios relative to the corresponding modified WT model.

### AE2 loss

Evaluate matched AE2 deletion with all Task 18 capacities frozen.

Do not introduce genotype-specific scaling or compensation.

## Stage G: interpretation

For baseline, 0.10 and 0.30, report:

- target and realized AE4 positive Cl-loading share at WT;
- complete signed chloride flux ledger;
- all adjustable capacity multipliers;
- largest and RMS log-fold parameter changes;
- algebraic membrane voltages;
- whether the exact WT conserved state is preserved;
- WT dynamic secretion;
- AE4 5% total secretion ratio and reduction;
- AE2 loss ratio and reduction;
- AE4 Na/K branch cancellation/productive fraction;
- routing family and calcium dependence;
- numerical status.

Only after all frozen genotype results exist may the experimental approximately 35% AE4-loss reduction be shown as contextual magnitude. It is not an optimization target or acceptance band.

Detailed time-course shape remains diagnostic only.

Do not search for an intermediate AE4 share if 0.10 and 0.30 behave differently.

## Required classifications

Choose exactly one primary classification:

1. `FIXED-WT REBALANCING FEASIBLE; AE4 NEAR-LOSS ROBUSTLY REDUCES SECRETION`
2. `FIXED-WT REBALANCING FEASIBLE; AE4 NEAR-LOSS STILL INCREASES SECRETION`
3. `FIXED-WT REBALANCING FEASIBLE BUT GENOTYPE EFFECT DEPENDS ON SHARE OR ROUTING`
4. `FIXED-WT REBALANCING REQUIRES EXTREME PARAMETER SHIFTS`
5. `FIXED-WT TARGET STRUCTURALLY INFEASIBLE EVEN WITHOUT INHERITED CAPACITY BOUNDS`
6. `TASK 18 NUMERICALLY INCONCLUSIVE`

If some roots/shares are feasible and others are not, do not force a universal classification unsupported by the full ensemble. Use classification 3 or 4 where appropriate and report the complete pattern.

## Required artifacts

Create at minimum:

- `results/18_fixed_wt_inverse_rebalance/contract.json`
- `results/18_fixed_wt_inverse_rebalance/parameter_provenance.csv`
- `results/18_fixed_wt_inverse_rebalance/baseline_flux_ledger.csv`
- `results/18_fixed_wt_inverse_rebalance/source_rank.json`
- `results/18_fixed_wt_inverse_rebalance/wt_inverse_solutions.csv`
- `results/18_fixed_wt_inverse_rebalance/wt_parameter_payloads/`
- `results/18_fixed_wt_inverse_rebalance/frozen_wt_manifest.json`
- `results/18_fixed_wt_inverse_rebalance/wt_dynamic_results.csv`
- `results/18_fixed_wt_inverse_rebalance/ae4_5pct_results.csv`
- `results/18_fixed_wt_inverse_rebalance/ae2_results.csv`
- `results/18_fixed_wt_inverse_rebalance/final_comparison.csv`
- `results/18_fixed_wt_inverse_rebalance/test_results.txt`
- `analysis/18_fixed_wt_inverse_rebalance/method.md`
- `analysis/18_fixed_wt_inverse_rebalance/parameter_changes.md`
- `analysis/18_fixed_wt_inverse_rebalance/genotype_interpretation.md`
- `analysis/18_fixed_wt_inverse_rebalance/final_answer.md`

The final answer must say plainly:

- whether 10% and 30% AE4 chloride-loading shares can support the exact same WT conserved states once uncertain capacities are genuinely allowed to rebalance;
- which parameters had to move and by how much;
- whether the old 14 nS K and 31.4 nS CaCC values were the previous blockers;
- whether the required rebalancing is modest or extreme;
- whether AE4 5% then decreases or increases secretion;
- whether AE2 remains comparatively neutral;
- whether the result is robust across roots, routing families and calcium values;
- whether the current AE4 mixed-cation architecture remains the likely next mechanism problem.

## Final actions

Run focused tests proving:

- exact fixed-state preservation;
- complete coupled source scaling;
- target chloride-share enforcement;
- no hidden inherited capacity bounds;
- phenotype-firewall integrity;
- conservation and electrical closure;
- complete ensemble accounting.

Commit and push completed Task 18 work to `codex/task-18-fixed-wt-inverse-rebalance`.

Do not edit `archive/`.

Do not merge to `main`.

Do not draft manuscript text.
