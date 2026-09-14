# Task 42: Catalán 2025 AE4 mechanism class discrimination

Work only on branch:

`codex/task-42-catalan-2025-ae4-mechanism-classes`

Read first:

1. `AGENTS.md`
2. `analysis/42_catalan_2025_ae4_mechanism_classes/source_contract.md`
3. `analysis/40_ae4_equal_cation_routing/final_answer.md`
4. `results/40_ae4_equal_cation_routing/`

Do not edit `main`.

Do not alter, merge, cherry pick, force push, or overwrite Task 41 draft PR #33 or commit `86f135e075aa91a70a752b650d7bb894353fee04`.

## Scientific question

The 2018 model assumed a 1 Cl inward : 1 monovalent cation outward : 2 HCO3 outward AE4 cycle. Task 40 retained the inherited scalar AE4 cycle law but imposed a 50:50 Na/K cation split. Under that model, AE4 loss caused only about a 3 to 4% cumulative secretion deficit because the network compensated strongly.

Catalán et al. 2025 now permits several more specific AE4 transport classes, including a proposed NaCl inward and K plus base outward cycle.

Test whether any of those source supported stoichiometric classes can reproduce a substantial AE4 loss secretion phenotype without the target selected AE4 dependent CaCC recruitment introduced in Task 41.

This is a mechanism discrimination task, not phenotype fitting.

## Immutable parent

Task 42 was branched directly from current `main` at:

`fe51473d7cc201e2e04007aa62fbc1b45098732e`

Use the completed Task 40 model as the scientific parent.

Keep fixed unless the source contract explicitly changes an AE4 source coefficient:

* the Task 40 scalar AE4 cycle law `J4(state)`;
* Palk/Benjamin NKCC1 law and inherited stimulation;
* Cha NHE1 implementation and stimulation;
* electrogenic NBC implementation and recruitment;
* Na/K pump;
* AE2;
* K channels;
* CaCC model and conductance;
* CO2 and acid base closure;
* paracellular transport;
* water transport;
* bath, geometry and stimulus;
* solver tolerances and production protocol.

The Task 41 AE4 dependent CaCC recruitment law is forbidden.

## Predeclared mechanism panel

Implement exactly the classes and source vectors in the source contract:

* C0 Task 40 legacy 50:50 control;
* C1 NaCl inward, K plus HCO3 outward, 1:1:1:1;
* C2 NaCl inward, 2K plus CO3 outward, 1:1:1:2;
* C3a K879 class with K plus 2 HCO3 outward;
* C3b K879 class with Na plus 2 HCO3 outward;
* C4a K879 class with K plus CO3 outward;
* C4b K879 class with Na plus CO3 outward.

Do not add intermediate Na/K mixtures. Do not sweep a cation fraction.

## Stage 1: source and thermodynamic audit

Before numerical solves:

1. Implement the source vectors in one small, testable mechanism layer without duplicating the full model.
2. Verify exact charge balance for every class.
3. Verify exact reversal when `J4` changes sign.
4. Derive and record the Na, K, Cl, TIC and TA balance identities for every class.
5. Calculate the forward chemical affinity at the accepted Task 40 WT state for every stoichiometry using the model concentrations and membrane potential. Treat carbonate explicitly through the model acid base state.
6. Record whether the source proposed forward direction is supported, opposed, or indeterminate under the frozen WT state.

No class may be discarded because its eventual secretion prediction looks inconvenient.

## Stage 2: WT resting state gate

For every mechanism class that is algebraically valid:

1. Start from the accepted Task 40 WT resting state.
2. Perform one intended WT stationary solve for that class.
3. Permit at most one purely numerical retry using the same equations and tolerances.
4. Do not optimise parameters, search roots, widen bounds, alter transport capacities or change the stimulus.
5. Apply the same resting physiology, charge, positivity and conservation gates used by Task 40.

If a class has no admissible WT rest under the frozen Task 40 model, record that as a result and do not run its genotype trajectories.

## Stage 3: frozen genotype predictions

For every class with an admissible WT rest:

1. Freeze that class specific WT rest.
2. From the same frozen WT rest run exactly three 600 second production trajectories:
   * WT, AE4 expression 1.0;
   * AE4 expression 0.05;
   * AE4 expression 0.00.
3. Use the inherited Task 40 production stimulus and numerical gates.
4. No genotype specific resting solve.
5. No genotype specific compensation, refitting or parameter adjustment.
6. Save compact time series and integrated flux summaries for all successful trajectories.

## Phenotype firewall

Do not inspect or use the approximately 35% experimental AE4 loss secretion deficit to make any modelling decision.

Do not compare candidate classes with Task 41 during execution.

After every predeclared class has either failed the WT gate or produced frozen genotype outputs, commit and push a checkpoint containing all predictions.

Only after that checkpoint may the held out phenotype and Task 41 be used for interpretation.

## Required outputs

For each class report:

* WT rest state and whether every gate passes;
* forward thermodynamic affinity at the Task 40 reference state;
* 600 second WT, 5% and null cumulative secretion where available;
* instantaneous secretion at 600 seconds;
* intracellular Na, K, Cl, pH and volume trajectories;
* AE4, NKCC1, AE2, NBC, NHE1, pump, CaCC and paracellular flux integrals;
* NKCC1 compensation relative to that class's WT;
* early versus late secretion deficit;
* exact source and current balance residuals;
* reason for any failed class.

Include a single comparison table across all mechanism classes.

After the prediction checkpoint, compare with:

* the original approximately 35% AE4 loss secretion phenotype;
* the experimental temporal pattern in which early secretion is relatively preserved and the sustained deficit develops later;
* Task 40;
* Task 41 as a target selected constructive comparison only.

## Interpretation rules

A class that lands near the experimental secretion deficit without phenotype based tuning is scientifically interesting, but not validated merely because the scalar deficit matches.

Its Na, Cl, pH, volume, NKCC1 response and temporal secretion pattern remain independent discriminators.

If none succeeds, do not tune the panel. Conclude only that Catalán 2025 stoichiometric source classes, with the inherited Task 40 scalar cycle law and all other mechanisms fixed, are insufficient to recover the phenotype.

Do not claim that Catalán 2025 is falsified, because that paper does not specify the complete salivary whole cell kinetic law or regulation.

## Compute discipline

This is a seven class finite panel, not an optimisation problem.

Maximum intended stationary solves: seven, plus at most one numerical retry per class.

Maximum intended production integrations: three per class that passes WT, plus at most one purely numerical retry per class if needed.

No parameter sweep.

No optimiser.

No phenotype driven root selection.

One worker and one BLAS thread.

## Publication

Keep all Task 42 source, tests, compact results and reports on the Task 42 branch.

Preserve every failed class and every successful class.

Do not delete or rewrite earlier task material.

Do not merge to `main`.

Open a draft PR only after the final prediction and interpretation report is complete.
