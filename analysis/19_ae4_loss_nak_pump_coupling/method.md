# Task 19 method

This is a baseline-model factorial mechanism diagnostic. Repository inheritance is merged Task 18 (`52a56cc8e8788cdb5cae60dcbaa75257de2878dd`), while the scientific model and ten accepted baseline roots are unchanged. Neither Task 18 rebalanced allocation is used.

The machine-readable contract is committed before any new resting continuation or trajectory. Its exact panel is ten roots (five in each AE4 routing family), AE4 expression 1.0 or 0.05, pump capacity 1.0, 0.50 or 0.10, and calcium 0.10, 0.25 or 0.50 uM. No phenotype data enter the contract, root selection, parameter construction, or calculation of the results.

Only `membranes.nak_capacity_fmol_s` is multiplied by the declared scale. The apical fraction is retained exactly, hence both membrane capacities receive the same multiplier. The production pump law and 3Na:2K sources, AE4 architecture/routing, other transporters, bath, geometry, regulation, and numerical tolerances remain unchanged. AE4 expression acts through the inherited genotype object on the complete AE4 contribution.

For each reduced pump condition independently, continuation begins from the corresponding inherited WT root, lowers the pump with AE4 intact, and then lowers AE4 to 0.05 at that fixed pump scale. Three local starts and the inherited ten-coordinate root solver are used. Root selection uses distance to the preceding accepted root, with residual as tie breaker. Failed steps are halved deterministically, subject to the frozen budget and minimum step. All attempted steps and admissible alternatives are retained. A failed pump parent blocks its dependent AE4-low arm. No other parameters are fitted, and no disconnected root is substituted.

The inherited broad continuation search bounds are retained and explicitly recorded. A boundary failure describes the numerical search domain; it is not evidence that no physical equilibrium exists outside that domain. Resting residual, rank, charge, current, positivity, and conservation gates are inherited. Historical WT target concentrations and whole-gland flow are not imposed as genotype acceptance gates.

Hash-valid scale-1 WT trajectories from Task 14 and AE4-0.05 trajectories from Task 14B are reused with exact parameter, state, genotype, solver, time-grid and integral checks. New valid resting states receive the same 600 s CCh + IPR protocol under production Radau, with three calcium levels and the original time grid. Full state trajectories are retained. Unavailable conditions stay in the 180-row accounting table with missing secretion values.

For each root/calcium, Q(e,s) is integrated lumen outflow in pL. The five prespecified comparisons are:

| Comparison | Definition |
|---|---|
| Matched AE4 effect | R_AE4(s) = Q(0.05,s) / Q(1,s) |
| Pump-only effect | R_pump_WT(s) = Q(1,s) / Q(1,1) |
| Pump effect in AE4-low cells | R_pump_AE4(s) = Q(0.05,s) / Q(0.05,1) |
| Multiplicative interaction | I(s) = R_AE4(s) / R_AE4(1) |
| Coupled phenotype | R_coupled(s) = Q(0.05,s) / Q(1,1) |

Each ratio requires valid numerator and denominator trajectories; unrelated missing arms do not suppress a valid pump-only comparison. Reduction is reported as 100(1 − ratio) percent. Exact ratios are retained, with a prespecified 0.0001 numerical neutral band and a descriptive 1% interaction-materiality threshold. The latter is not a biological target. A reduced coupled ratio alone does not demonstrate an AE4-specific interaction.

Reproduction uses `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=src python -m modern_full_model.task19_pump_coupling STAGE --workers 9`, with stages `prepare`, `rest`, and `dynamic`; the contract must be committed between prepare and rest. Focused tests do not launch trajectories. The complete numerical results are frozen before experimental context is discussed.

After dynamics, `PYTHONPATH=src python -m modern_full_model.task19_reporting freeze` writes the diagnostic tables and result hashes. An explicit reporting correction is recorded in the final `results_frozen.json`, with the first receipt preserved as `results_frozen_v1.json`: the production AE4 evaluator exposes branch diagnostics before genotype scaling. The final resting and state/flux tables report actual branches as minus the conserved Na/K sources, retain the raw evaluator values in `before_expression` columns, and verify that the scaled branches sum to the net Cl source. Original continuation logs preserve the original pre-expression branch fields; their conserved source fields and all residuals were already correct. No resting solve, trajectory, Q value, factorial comparison or selection was changed for this correction.
