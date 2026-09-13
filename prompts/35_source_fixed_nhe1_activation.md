# Task 35: source-fixed stimulated NHE1 activation

Repository: `esig626/ae4-salivary-transport-control`
Branch: `codex/task-35-source-fixed-nhe1-activation`

## Mission

Task 34 stopped correctly because Peña-Münzenmayer et al. 2015 Fig. 2D reports the slope of normalized BCECF fluorescence rather than an absolute pH slope. Do not repeat that failed calibration attempt.

This task tests the same single missing mechanism using an independent primary source that reports calibrated intracellular pH and a quantitative stimulation fold-change for NHE1.

The hypothesis is:

> Muscarinic/Ca-dependent stimulation increases the effective activity of the already-implemented mechanistic NHE1 enough to replenish intracellular alkalinity/HCO3 during secretion, allowing the existing AE4 mechanism to contribute substantially during sustained stimulation even if AE4 flux is small at rest.

This is a source-fixed mechanism test, not an optimisation task.

## Primary source constraint

Before changing code, verify the following directly from:

Park et al., 1999, Journal of Biological Chemistry
"Targeted Disruption of the Nhe1 Gene Prevents Muscarinic Agonist-induced Up-regulation of Na+/H+ Exchange in Mouse Parotid Acinar Cells"
DOI: `10.1074/jbc.274.41.29025`
PMID: `10506152`

The source reports, in mouse parotid acinar cells:

- intracellular pH measured with BCECF and calibrated in situ with the high-K+/nigericin method;
- `10^-5 M` carbachol stimulation;
- a `2.3-fold` increase in the initial rate of NHE-mediated pH recovery in WT acini;
- an approximately `0.13 pH unit` stimulation-associated alkalinisation beyond the original resting pH in the relevant experiment;
- NHE1, rather than NHE2 or NHE3, as the exchanger responsible for muscarinic up-regulation;
- in HCO3-replete conditions, WT pH recovery after carbachol stimulation of about `0.15 pH units/min`, versus about `0.04 pH units/min` in Nhe1-null acini;
- the source describes both a rapid component associated with secretion/cell shrinkage and a slower Ca-dependent component of NHE1 up-regulation.

Also record that Peña-Münzenmayer et al. 2015 explicitly states that stimulation causes a dramatic increase in NHE activity, raises intracellular pH/HCO3, has a slow time course consistent with the delayed AE4 secretion phenotype, and cites prior NHE work in support of that interpretation.

The Park 1999 experiment is parotid and uses a higher carbachol concentration than the 2015 SMG protocol. Therefore the `2.3-fold` value is a source-supported mechanistic transfer for this bounded test, not a claim of an exact SMG dose-response measurement.

If the primary source does not support the statements above, stop before model modification and report:

`SOURCE_FIXED_NHE1_GAIN_NOT_REPRODUCED`

Do not replace it with another paper or another numerical gain inside this task.

## Fixed model state and architecture

Start from the completed Task 34 branch state. Task 34 made no model-mechanism changes.

Use the accepted Task 31 R09 WT resting state unchanged as the common initial condition.

Freeze all of the following:

- Cha et al. 2009 eight-state + Mod2 NHE1 kinetic mechanism;
- Task 31 effective resting NHE1 carrier amount: `2.339370005697548e-05 fmol`;
- all Cha kinetic constants, affinities, modifier exponents and stoichiometry;
- AE4 routed topology and 1:1:2 working stoichiometry;
- all AE4 resting parameters;
- existing beta/cAMP/PKA AE4 regulation, including its family, gain and time constants;
- NKCC1, AE2, pump, K channels, CaCC, CO2, paracellular and water functional forms and parameters;
- bath, geometry and membrane allocation fractions;
- the Task 31 R09 resting state itself.

Do not use Task 33 optimised parameter vectors.
Do not force a 70/30 chloride partition at rest.
Do not re-solve or re-optimise rest.

## The only new mechanism

Add exactly one task-specific stimulated NHE1 activity factor.

At rest:

`G_NHE = 1.0`

During the standard central muscarinic/Ca secretory stimulation used in this project:

`G_NHE = 2.3`

Implement this as a multiplicative factor on the effective active NHE1 carrier amount / cycle capacity consumed by the already-built Cha NHE1 evaluator.

Do not alter the Cha kinetic law itself.
Do not alter intracellular proton dependence.
Do not alter Na/H stoichiometry.
Do not add a new NHE state or a new adjustable time constant.
Do not invent an NHE dose-response curve.

This task tests only the source-supported stimulated magnitude. The existing whole-cell state dynamics determine how pH/HCO3 and AE4 respond over time.

The NHE1 stimulation factor must be driven by the muscarinic/Ca secretory condition, not by beta/PKA. Park 1999 reports muscarinic/Ca-linked NHE1 up-regulation and notes that the acinar response is not a PKA mechanism. The existing beta/PKA pathway remains reserved for the existing AE4 regulation.

Because this task permits only the standard resting and central stimulated conditions, no interpolation between `1.0` and `2.3` is required.

If implementing the source-fixed factor would require any additional free parameter, stop and report:

`SOURCE_FIXED_NHE1_TEST_NOT_PARAMETER_FREE`

## No calibration and no optimisation

There is **no parameter fitting in Task 35**.

The stimulated NHE1 factor is exactly `2.3`.

Do not evaluate 1.5, 2.0, 2.5, 3.0 or any other gain.
Do not bracket anything.
Do not optimise anything.
Do not perform a sensitivity sweep.

The `0.13 pH unit` and `0.15 pH units/min` source observations are validation context only because the source gland, agonist concentration and acid-load protocol differ from the present SMG whole-cell stimulation. Report the model's corresponding pH behaviour but do not fit to those values.

## Verify the existing AE4 stimulation path

Before dynamics, verify by code inspection/unit test that the existing beta/cAMP/PKA regulatory output reaches the routed AE4 evaluator during the standard IPR/beta stimulus.

Do not modify that regulation unless there is a literal implementation bug preventing the already-intended signal from reaching AE4. If such a bug is found, document it and make only the minimal wiring correction. Do not change regulatory parameters or choose another regulatory family.

## Dynamic protocol

Run exactly one standard WT stimulation first:

- initial condition: accepted Task 31 R09 WT resting state;
- AE4 expression: `1.00`;
- Ca = `0.25 uM`;
- existing standard beta/IPR input;
- time = `0--600 s`;
- NHE1 stimulated factor = exactly `2.3` under the stimulated condition.

Do not solve a new WT resting state.

### WT physical admissibility

Require throughout the WT trajectory:

- finite positive concentrations and volumes;
- `Na_i < 40 mM`;
- `30 <= Cl_i <= 80 mM`;
- `6.6 <= pH_i <= 7.3`;
- cell volume `< 3 pL`;
- no charge/current/carbon/water accounting violation beyond existing numerical tolerances.

### WT secretion activation gate

The WT stimulated trajectory must genuinely secrete more than the resting baseline.

Using the accepted Task 31 R09 resting outflow

`q_out_rest = 0.0010757073853493032 pL/s`,

require:

- mean `q_out` over `60--600 s` > `1.10 * q_out_rest`;
- cumulative outflow over `0--600 s` > `1.10 * (600 * q_out_rest)`;
- `q_out(600 s) >= q_out_rest`.

If WT fails physical admissibility or this secretion gate, stop immediately and report:

`WT_SOURCE_FIXED_NHE1_STIMULATION_FAILED`

Do not change any parameter and do not run the AE4=5% or AE4=0% trajectories.

## WT mechanistic readout

If WT passes, report at minimum:

- pH_i at 0, 60, 120, 180, 300 and 600 s;
- HCO3_i at the same times;
- Na_i and Cl_i at the same times;
- q_out at 60 s intervals;
- cumulative outflow over 600 s;
- integrated NHE1 Na-in/H-out flux over 0--600 and 60--600 s;
- integrated positive NKCC1 chloride loading over 60--600 s;
- integrated positive AE4 chloride loading over 60--600 s;
- signed integrated AE2 chloride flux over 60--600 s;
- integrated CaCC chloride export;
- integrated paracellular chloride return;
- lumen osmolarity trajectory;
- existing AE4 regulatory capacity multiplier at 0, 60, 180 and 600 s.

Define the stimulated positive basolateral chloride-loading pool over 60--600 s as the sum of positive integrated NKCC1, AE4 and AE2 loading contributions.

Report their shares.

Compare the NKCC1 share with the approximate `65--75%` experimental context.

This is validation only.
Do not tune or reject the trajectory solely because the share lies outside that range.

## Exactly two AE4 expression tests

Only after the WT trajectory passes, run exactly:

- AE4 expression = `0.05`;
- AE4 expression = `0.00`.

For both perturbations:

- use the identical accepted Task 31 R09 WT resting state as the initial condition;
- use identical Ca and beta/IPR stimulation;
- use the identical fixed `G_NHE = 2.3` stimulated NHE1 factor;
- change only AE4 expression.

Do not solve genotype-specific resting states in Task 35.

This remains an acute matched-initial-state stimulated-expression test. It is not a claim about constitutive knockout resting adaptation.

For each perturbation report:

- physical admissibility;
- cumulative 600 s secretion;
- secretion ratio to the frozen WT trajectory;
- q_out at 60 s intervals;
- pH_i, HCO3_i, Na_i and Cl_i trajectories;
- integrated NHE1 flux;
- integrated NKCC1 chloride loading;
- integrated AE4 chloride loading;
- signed integrated AE2 chloride flux;
- CaCC chloride export;
- paracellular chloride return;
- lumen osmolarity and water outflow;
- WT-minus-perturbation cumulative secretion as a function of time;
- whether the secretion deficit is small during the first `2--3 min` and becomes more pronounced thereafter.

The experimental ~35% AE4-null secretion deficit is held-out validation only.
Do not fit to it.

## Absolute anti-combinatorial restrictions

There is no optimisation, search or parameter sweep in Task 35.

Forbidden:

- any second NHE gain;
- any NHE gain sweep;
- Cartesian grids;
- one-at-a-time parameter sweeps;
- pairwise or factorial searches;
- random search;
- Latin hypercube;
- multistart;
- global optimisation;
- evolutionary optimisation;
- Bayesian optimisation;
- Shapley analysis;
- alternate parameter subsets;
- alternative NHE models;
- alternative AE4 regulatory families;
- NBC addition;
- carbonic-anhydrase activation addition;
- another bicarbonate source;
- changing NKCC1 capacity;
- changing AE4 capacity;
- changing AE2 capacity;
- changing water parameters;
- changing the calcium input;
- changing beta/IPR input;
- trying intermediate AE4 expression levels;
- R10;
- AE2 knockout;
- phenotype fitting;
- any automatic escalation after failure.

If this single fixed mechanism does not work, report the failure and stop.

## Hard numerical budget

- zero stationary solves;
- zero optimisation calls;
- exactly one intended WT integration;
- at most one WT numerical retry only if the first integration fails for a purely numerical solver reason;
- after WT passes, exactly one AE4=0.05 integration and exactly one AE4=0.00 integration;
- maximum 4 dynamic integration attempts total including the possible WT retry;
- one worker;
- one BLAS thread;
- maximum 10 minutes numerical execution.

Do not use the spare integration slot for another biological condition.

## Outputs

Write compact UTF-8 outputs under:

- `analysis/35_source_fixed_nhe1_activation/`
- `results/35_source_fixed_nhe1_activation/`

At minimum include:

- `source_record.md`;
- `implementation_note.md`;
- `budget.json`;
- `wt_summary.json`;
- `chloride_loading.csv` if WT passes;
- `genotype_comparison.csv` if WT passes;
- `timepoints.csv` if dynamics run;
- `verification.json`;
- `final_answer.md`.

Do not publish giant trajectory arrays or binary archives.

## GitHub publication

Use the connected GitHub integration for all remote writes.
Do not depend on shell Git authentication, PATs, SSH or `gh auth`.
Shell Git is only for local status/diff inspection.
Publication failure must not trigger scientific recomputation.
Do not merge to `main`.

The task ends after either the source gate, the WT gate, or the exactly two AE4 expression tests. Do not start another mechanism from inside Task 35.
