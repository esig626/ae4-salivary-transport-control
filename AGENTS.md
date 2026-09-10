# AGENTS.md

These instructions apply to all work in this repository.

## Current scientific objective

The active task is Task 15. The executable prompt is

`prompts/15_ae4_sign_and_mixed_cation_slip_audit.md`.

Work on branch

`codex/task-15-ae4-sign-slip-audit`.

The controlling completed Task 14B commit is

`4d4403f279fc36aec940c9e4759486d02f07e3f6`.

Task 15 asks a narrow diagnostic question:

> Are the modern AE4 and whole-cell transport signs correct, and does the current reversible shared Na/K AE4 architecture generate a large mixed-cation slip mode under the frozen physiological states?

This is not a calibration task, not a repair task, and not a new blind holdout.

## Scientific state entering Task 15

Task 14 retained all ten WT roots, calcium values 0.10, 0.25 and 0.50 uM, and regulatory member `R1_G125_P21_PROSE_TREFERENCE` with the model frozen. Exact AE4 continuation reached 5% expression in all ten roots but failed at zero; AE2 exact deletion succeeded.

Task 14B then used the already validated connected 5% AE4 resting states and the frozen Task 14 WT trajectories. Across all 30 root/calcium pairs, 5% AE4 increased total 0 to 600 s secretion by 8.23% to 21.66%, median 17.82%. The effect was opposite to the experimentally observed AE4-loss direction. All 30 trajectories passed numerical checks.

The AE4 phenotype is now known and cannot be treated as unseen again.

## Why Task 15 exists

Before adding compensation, changing regulation, changing transport laws, or fitting anything, audit the possibility that the wrong-direction phenotype arises from an incorrect sign convention or from an unsupported property of the current AE4 architecture.

The modern AE4 module declares positive branch flux as chloride loading:

`Cl_o + C_i + 2 HCO3_i <-> Cl_i + C_o + 2 HCO3_o`, with `C = Na or K`.

Preliminary inspection shows that representative frozen WT states have positive net AE4 chloride loading but large opposing Na and K branch currents. The current shared-carrier graph therefore permits a model-generated Na/K cation-slip cycle. Task 15 must quantify this rigorously across every retained case and determine whether any actual sign or source-mapping error exists.

Do not treat the preliminary inspection as the final result. Reproduce it independently in the Task 15 artifacts.

## Hard freeze

Base the Task 15 audit exactly on completed Task 14B commit `4d4403f279fc36aec940c9e4759486d02f07e3f6`.

Do not modify any scientific model equation or parameter.

Keep exactly as inherited:

- all ten retained WT roots;
- all whole-cell equations;
- all transporter capacities and source multipliers;
- AE4 carrier graph, stoichiometry, rates and Na/K routing;
- AE2, NKCC1 and NHE1 laws;
- Na/K pump and K-channel topology;
- CaCC and calcium gate;
- paracellular pathways;
- acid-base and carbon chemistry;
- hydraulic and lumen outflow laws;
- bath and geometry;
- beta/cAMP/PKA regulatory member and parameters;
- calcium inputs 0.10, 0.25 and 0.50 uM;
- the saved Task 14 WT trajectories;
- the saved Task 14B 5% AE4 trajectories;
- solver and conservation tolerances.

No root or calcium value may be dropped.

Do not attempt exact AE4 zero.

## Diagnostic-only rule

Task 15 may add audit code, tests, tables and reports only.

Do not:

- flip any sign in production code;
- replace the AE4 cycle;
- suppress Na/K slip;
- impose a no-slip constraint;
- retune the mixed-bath Na/K routing;
- change stoichiometry;
- fit to the known 35% secretion reduction;
- introduce genotype-specific compensation;
- rerun calibration;
- choose a root or calcium value because it agrees better with experiment.

If a clerical, algebraic, source-mapping or expression-scaling error is found, document it precisely and stop before repair. The repair belongs in a later task.

If signs are correct but large opposing Na/K branch currents are present, quantify them and classify them as a property of the current model architecture. Do not claim that such a slip cycle is established native physiology unless direct evidence supports it.

## Sign audit requirements

Construct a complete sign truth table for AE4, NKCC1, NHE1, AE2, CaCC, K channels, pumps, paracellular transport and water/outflow.

For each component identify the positive mathematical direction, biological direction, cell source sign, lumen source sign where applicable, and the exact source/RHS code path.

Conservation alone is not evidence of correct orientation. A reversed electroneutral cycle can conserve charge perfectly.

Independently derive the AE4 source vector from the declared reaction and verify the implementation and whole-cell insertion.

Explicitly audit genotype expression scaling. The raw AE4 diagnostic branch currents may be reported before expression scaling, whereas the whole-cell source is scaled by genotype expression. Distinguish these quantities and prove that the 5% state applies the 0.05 factor exactly once.

## Reversal and direction tests

Use constructed concentration states, independent of phenotype agreement, to verify:

- global AE4 equilibrium gives zero flux;
- positive Na-cycle affinity produces the declared positive Na branch and source signs;
- reversed Na gradients reverse those signs;
- the corresponding K tests;
- affinity, branch flux and entropy production are thermodynamically consistent;
- 5% expression preserves sign and scales the actual whole-cell AE4 source by 0.05;
- no second sign inversion occurs in the whole-cell RHS.

Add equivalent minimal direction tests for other secretion-relevant transport and water components if existing tests do not already prove them.

Do not change production equations to make a test pass.

## Full-ensemble slip audit

Evaluate the frozen model on every saved WT and 5% AE4 trajectory time point without regenerating the ODE trajectories.

For each root, calcium, genotype state and time record AE4 Na/K affinities, actual expression-scaled Na and K branch fluxes, net AE4 Cl source, source vector, entropy production, regulation gain, expression scale and relevant intracellular/bath state.

Use

`branch_turnover = abs(J_na) + abs(J_k)`

`productive_fraction = abs(J_na + J_k) / branch_turnover`

`cancellation_fraction = 1 - productive_fraction`

when branch turnover is nonzero.

When `J_na * J_k < 0`, record

`slip_pair_magnitude = min(abs(J_na), abs(J_k))`.

Positive `J_na + J_k` means net AE4 chloride loading under the declared sign convention. Negative means net chloride unloading.

Summarize the fraction of time with opposing branch directions, fraction of time with net loading/unloading, cancellation fraction, productive fraction, integrated branch turnover, and integrated net Cl source for every case and ensemble-wide.

## Evidence discipline

Use the retained repository evidence for the 2016 and 2025 AE4 experiments.

The data establish Na-supported and K-supported AE4 transport, reversibility in manipulated gradients, mutation-dependent cation specificity, and macroscopic electroneutrality. They do not automatically establish a large simultaneous opposing Na/K cycle in a physiological mixed bath.

Compare with the 2018 AE4 model only at the level of sign and cation-coupling structure. Determine whether it permitted independent opposite-sign Na and K cycles or instead partitioned one net AE4 exchange direction.

## Required output

Create the artifacts specified in the Task 15 prompt under:

`results/15_ae4_sign_slip_audit/`

and

`analysis/15_ae4_sign_slip_audit/`.

The final answer must state plainly:

1. whether any sign is actually wrong;
2. whether AE4 is net loading or unloading chloride in the frozen WT states;
3. whether Na and K branches oppose one another;
4. how large the cancellation/slip is;
5. whether that mixed-bath behavior is experimentally established or merely permitted by the model;
6. whether Task 16 should repair a sign error, test a no-slip/alternative AE4 architecture, or look elsewhere.

## Classification

Choose exactly one classification from the Task 15 prompt.

Do not blur a sign error and an unsupported architecture into the same conclusion.

## General discipline

Do not edit `archive/`.

Do not merge to `main`.

Do not draft manuscript text.

Do not perform new calibration, genotype fitting, model reduction, GSPT or a new identifiability study.

Commit and push the completed Task 15 audit to `codex/task-15-ae4-sign-slip-audit`.
