# Task 14: scale free AE4 and AE2 genotype holdout validation

## Mission

Task 13B produced a modern conservation explicit whole cell model family with ten retained native WT roots. Task 13C then tested whether the remaining absolute WT gland flow discrepancy could be removed by changing only stimulated calcium. It could not.

The controlling Task 13C commit is

`2f54e7c4f87b6da746d3a8427c7bdca40a77c3de`

Task 13C preserved every frozen model parameter, reproduced the saved 0.10 uM WT reference, newly simulated only 0.25 and 0.50 uM, and found that no retained root reaches the nominal 9 to 10 uL/min whole gland range under the inherited one SMG observation mapping. The best 0.50 uM case reaches 5.405993 uL/min under that mapping. All numerical and sustainment gates passed. Five roots also failed the inherited minute observation scale shape criterion.

Those findings remain recorded. They are not deleted or rewritten.

Task 14 changes the interpretation of one thing only:

> The absolute one SMG flow mapping is now a nonblocking observation scale diagnostic, not a prerequisite for testing genotype biology.

The cellular model and the whole gland measurement are connected by an uncertain multiplicative observation scale. A common multiplicative scale cancels from genotype to WT secretion ratios. Task 14 therefore asks the central held out question directly.

> With the model completely frozen, does exact AE4 deletion predict the observed WT relative secretion phenotype and its time structure, while exact AE2 deletion remains comparatively small?

This is a prediction and validation task. It is not a calibration task, not a reconstruction task and not a mechanism search.

Read `AGENTS.md` first and obey it throughout.

Do not edit anything under `archive/`.

Do not merge to `main`.

Do not draft a manuscript.


## 1. Scientific decision fixed before the run

The Task 13B and Task 13C absolute gland flow discrepancy remains an unresolved observation map issue.

For Task 14:

* Keep reporting absolute per cell and mapped whole gland flows as diagnostics.
* Keep reporting the old one SMG scale interval and the old minute observation shape result for provenance.
* Do not use either quantity to accept or reject a genotype prediction.
* Do not rescale AE4 null and WT separately.
* Do not estimate a genotype specific effective cell number.
* Do not change the gland scale to improve agreement with the knockout experiment.

The primary endpoints are scale free paired quantities formed from an AE4 null or AE2 null trajectory and its matching frozen WT trajectory.

For a common observation factor `S`,

`Y_WT(t) = S Q_WT(t)`

and

`Y_null(t) = S Q_null(t)`.

Therefore

`Y_null(t) / Y_WT(t) = Q_null(t) / Q_WT(t)`.

This cancellation is the reason Task 14 may proceed despite the unresolved absolute scale.

The observation scale is a nuisance quantity. It is not permitted to become a hidden fitting parameter in this task.


## 2. Hard model freeze

The branch starts exactly from Task 13C commit

`2f54e7c4f87b6da746d3a8427c7bdca40a77c3de`.

Before running a genotype trajectory, verify and record the inherited tree, Task 13B final manifest, Task 13C frozen manifest, all ten retained root payloads and the relevant equation source hashes.

Keep fixed exactly as inherited:

* all ten retained WT roots and their whole cell parameter payloads;
* the twelve core state equations and all regulatory state equations;
* NKCC1 family, source scale and kinetics;
* AE4 carrier model, total amount, rates, Na and K routing and stoichiometry;
* AE2 and NHE1;
* Na K pump capacities and membrane fractions;
* K channel conductances and membrane fractions;
* CaCC conductance and calcium gate;
* paracellular conductances;
* acid base and carbon chemistry;
* hydraulic coefficients, lumen parameters and outflow law;
* geometry and bath composition;
* beta input and cAMP PKA gains and kinetics;
* all solver tolerances and conservation tolerances;
* all resting OTHER osmole amounts;
* the exact genotype deletion rule.

No parameter optimization is permitted.

No parameter may depend on genotype unless the existing frozen genotype object already represents exact transporter expression deletion.

Do not add compensation, transporter recruitment, altered calcium, altered pump activity, altered conductance, altered regulation or altered geometry after deletion.


## 3. Calcium uncertainty is a frozen sensitivity panel, not a fitting knob

The exact matched SMG intracellular calcium amplitude under the dual agonist protocol remains uncertain.

Carry forward exactly these three already evaluated WT calcium amplitudes:

`0.10, 0.25, 0.50 uM`.

Do not introduce new calcium values.

Do not select one of the three after seeing the held out phenotype.

Use all three as a fixed sensitivity panel and report the genotype prediction separately at each value.

The preferred dynamic regulatory member is the Task 13C member

`R1_G125_P21_PROSE_TREFERENCE`.

The primary stimulus remains `CCH_IPR`, 0.3 uM CCh plus 5 uM IPR for 600 s, with the inherited beta input.

Resting calcium remains 0.058 uM.

The existing `genotype_evaluation.py` contains useful target free continuation and exact deletion logic, but its old frozen dynamic helper encodes an earlier R3 and 0.10 uM contract. Do not silently use that helper as the Task 14 dynamic protocol. Either add a thin Task 14 builder or minimally parameterise the builder so that the Task 13C R1 member and the three predeclared calcium values are explicit. Do not change the underlying model equations.


## 4. Resting genotype continuation

A null genotype must not be created by simply evaluating the null vector field at the WT resting state.

For each of the ten retained WT roots, independently continue transporter expression from 1 to exact 0 for

* AE4;
* AE2.

Use the existing target free 21 point expression grid unless a saved Task 13B contract specifies the same grid more exactly.

At every continuation step:

* solve only the resting state coordinates already declared by the genotype continuation code;
* fit zero model parameters;
* keep the intracellular OTHER osmole amount fixed;
* retain alternate admissible roots rather than hiding them;
* select the connected branch using only the predeclared distance rule;
* require exact deleted transporter flux at expression 0;
* require positivity, current closure, charge consistency, conservation and full numerical rank under the existing numerical gates.

If a branch fails before exact deletion, record that as a model prediction failure. Do not repair it using phenotype information.

Run the existing independent genotype continuation implementation as an audit on every branch that is numerically ambiguous and on at least one representative root from each retained AE4 routing family for both AE4 and AE2 deletion.

Resting genotype continuation is independent of stimulated calcium and must be performed only once per root and transporter.


## 5. Blind prediction firewall

This section is nonnegotiable.

Before the blind prediction checkpoint is committed and pushed:

* do not open the held out phenotype table;
* do not parse it through Python, pandas, shell tools or tests;
* do not inspect digitised AE4 null secretion values;
* do not inspect any stored KO to WT ratios or minute specific target values;
* do not use target values in a code comment, test fixture, threshold or branch selection rule;
* do not search repository history for the concealed values;
* do not compare a preliminary trajectory with the published knockout curve by eye.

The public literature context that AE4 deletion causes a substantial sustained secretion deficit and AE2 deletion has a much smaller secretion phenotype is already part of the project history. That qualitative knowledge cannot be undone. The purpose of the firewall is to prevent the exact hidden magnitude and time course from influencing the frozen numerical prediction.

Create before simulation:

`results/14_scale_free_genotype_holdout/prediction_contract.json`

It must state, before any genotype result exists:

* the exact ten root IDs;
* the two null genotypes;
* the three calcium amplitudes;
* the R1 regulatory member;
* the stimulus protocol;
* the solver specifications;
* the continuation rule;
* the reporting times;
* the scale free metrics below;
* the numerical failure rules;
* the no pruning rule;
* the rule that all roots and all calcium values remain in the final report regardless of agreement.

Hash the contract and include its hash in every blind output.


## 6. Blind dynamic prediction

For every retained root and every calcium value, generate matched trajectories for

* WT;
* AE4 null from its connected exact zero resting root;
* AE2 null from its connected exact zero resting root.

Use production Radau for the primary prediction.

This is 10 roots x 3 calcium values x 3 genotypes = 90 production trajectories if every continuation succeeds.

Parallel execution across roots is encouraged. Set numerical library thread counts so process level parallelism does not oversubscribe the machine.

Do not prune roots after a convenient prediction appears.

Use the exact same dynamic model, calcium input, regulatory member, solver and time grid within each WT null pair.

For each trajectory record at minimum:

* solver success;
* positivity;
* exact deletion check;
* conservation residuals;
* current and charge closure diagnostics available from the frozen model;
* flow at 60, 120, 180, 240, 300, 360, 420, 480, 540 and 600 s;
* cumulative flow at the same landmarks;
* total 0 to 600 s flow;
* mean flow;
* endpoint intracellular Na, K, Cl, pH and volume;
* endpoint luminal quantities already exposed by the model;
* apical and basolateral potentials where available.

For each matched null and WT pair compute, without any experimental target:

### Primary scale free secretion quantities

`R_total = cumulative_null(600) / cumulative_WT(600)`

`D_total = 1 - R_total`

minute specific instantaneous ratios

`R_flow(t) = Q_null(t) / Q_WT(t)`

minute specific cumulative ratios

`R_cum(t) = cumulative_null(t) / cumulative_WT(t)`.

### Prespecified early and sustained summaries

Define

* early interval: 0 to 180 s;
* sustained interval: 180 to 600 s.

Compute the integrated secretion ratio separately in those two intervals.

Also report

`delta_emergence = R_late - R_early`.

A negative value means the genotype deficit becomes stronger later. This is a descriptive prediction before reveal, not a target based pass rule.

### Genotype specificity

For each root and calcium value compute

`AE4_effect = 1 - R_total_AE4`

`AE2_effect = 1 - R_total_AE2`

and

`effect_contrast = AE4_effect - AE2_effect`.

Do not impose a hidden minimum effect contrast before reveal. Record the prediction as it is.

### Ionic phenotype

For AE4 null and AE2 null report WT relative differences in resting and endpoint

* intracellular Cl;
* intracellular pH;
* intracellular Na;
* intracellular K;
* cell volume.

Do not tune a resting root or dynamic parameter to improve those values.


## 7. Numerical confirmation before reveal

Before opening the holdout, perform BDF confirmation using a completely target independent selection rule.

At minimum crosscheck:

* one fixed representative root from each of the two retained AE4 routing families at each of the three calcium values;
* both AE4 null and AE2 null for those representatives;
* the matching WT denominators.

Choose representatives lexicographically within each routing family before running the primary predictions and record them in `prediction_contract.json`.

If Radau and BDF disagree beyond the inherited genotype solver tolerance, expand BDF confirmation to every root in that affected family and calcium condition. Do not use phenotype agreement to decide which case receives confirmation.

Any continuation anomaly or numerical anomaly must also receive the declared independent audit before reveal.


## 8. Blind checkpoint

Before any holdout value is opened, create all of the following:

### Analysis

`analysis/14_scale_free_genotype_holdout/freeze_and_design.md`

`analysis/14_scale_free_genotype_holdout/blind_prediction.md`

### Results

`results/14_scale_free_genotype_holdout/prediction_contract.json`

`results/14_scale_free_genotype_holdout/genotype_continuation.csv`

`results/14_scale_free_genotype_holdout/ae4_null_blind_predictions.csv`

`results/14_scale_free_genotype_holdout/ae2_null_blind_predictions.csv`

`results/14_scale_free_genotype_holdout/blind_pair_summary.csv`

`results/14_scale_free_genotype_holdout/blind_numerical_audit.json`

`results/14_scale_free_genotype_holdout/blind_artifact_hashes.csv`

The blind prediction report must state the complete model prediction without mentioning whether it agrees with experiment.

Then:

1. run focused Task 14 tests and relevant inherited genotype and Task 13C regression tests;
2. commit all blind implementation and prediction artifacts;
3. push that checkpoint to `codex/task-14-scale-free-genotype-holdout`;
4. record the exact blind checkpoint SHA in a machine readable file.

Only after that pushed commit exists may the holdout be opened.


## 9. One time holdout reveal

After the blind checkpoint is committed and pushed, locate the existing held out AE4 phenotype source already retained by the project.

Open it once for evaluation.

Do not modify the blind prediction files after reveal except to correct a demonstrated serialization or arithmetic bug that does not change any trajectory. Any such correction must preserve the original blind files and document both versions and hashes.

Record:

* the holdout source path;
* its file hash;
* the blind checkpoint SHA;
* the time of reveal;
* the exact fields used for comparison;
* any uncertainty or replicate information stored with the target;
* which target quantities are direct measurements and which are derived or digitised.

Do not convert experimental uncertainty into a biological law. Compare prediction and experiment at the resolution actually supported by the data.


## 10. Post reveal comparison

The comparison must use the already frozen blind predictions only.

### AE4 total secretion

Compare the predicted `R_total` and `D_total` with the held out 600 s WT relative secretion result.

Report for every root and every calcium value:

* signed error in reduction fraction;
* absolute error in reduction fraction;
* whether the prediction lies inside any directly reported experimental uncertainty interval, if such an interval is actually available;
* distance from that interval if outside it.

Do not select the best root or best calcium as the reported answer.

Report the complete range, median and routing family summaries.

### AE4 time structure

Compare the predicted minute specific ratios with the held out time course.

Use paired WT relative quantities so the unknown absolute scale cancels.

Where the holdout contains uncertainty intervals, compare against those intervals. Where it does not, report raw residuals without inventing confidence bands.

Report separately:

* early 0 to 180 s agreement;
* sustained 180 to 600 s agreement;
* whether the direction of `delta_emergence` agrees with the measured emergence of the deficit;
* a simple prespecified root mean square error over the available relative time points.

No time shift, gain, offset or smoothing parameter may be fitted after reveal.

### AE4 ionic phenotype

Compare predicted resting and stimulated intracellular Cl and pH with independent AE4 null measurements where protocol matching permits it.

Keep resting and stimulated measurements separate.

Do not force measurements from different gland, assay or stimulus contexts into one quantitative score.

### AE2 specificity

Compare the blind AE2 null predictions with the available AE2 phenotype evidence.

The central question is whether the same frozen model that gives an AE4 phenotype also predicts a substantially smaller AE2 secretion effect.

Do not use AE2 to rescale or repair AE4.


## 11. Interpretation of the absolute flow discrepancy

Task 14 must preserve the following distinction.

A failure of the inherited 9 to 10 uL/min one SMG mapping is not evidence, by itself, that the genotype ratio is invalid.

Conversely, a successful genotype ratio does not prove that the absolute observation map is correct.

The final report must therefore keep two statements separate:

1. Cellular genotype prediction under a common scale.
2. Unresolved mapping from cellular secretion to absolute whole gland output.

If genotype ratios agree well while the absolute scale remains low, record that as evidence that the transporter mechanism may be correct while the cell to gland observation map remains uncertain.

If genotype ratios fail, do not blame the observation scale, because a common multiplicative scale cancels from the primary ratios.


## 12. No post reveal repair

After holdout reveal, this task is analysis only.

Do not:

* change a parameter;
* choose a different root;
* choose a different calcium value;
* alter regulation;
* add compensation;
* alter transporter stoichiometry;
* alter the initial state;
* alter the solver to improve biological agreement;
* fit an observation scale;
* fit a time shift;
* fit a dynamic gain;
* add a new mechanism.

If the prediction fails, diagnose the failure using saved trajectories and end the task.

A repair, if scientifically justified, belongs in a later task with a new blinded target or an independent data source.


## 13. Required outputs after reveal

Create

### Analysis

`analysis/14_scale_free_genotype_holdout/holdout_reveal.md`

`analysis/14_scale_free_genotype_holdout/final_answer.md`

### Results

`results/14_scale_free_genotype_holdout/reveal_manifest.json`

`results/14_scale_free_genotype_holdout/ae4_holdout_comparison.csv`

`results/14_scale_free_genotype_holdout/ae2_validation_comparison.csv`

`results/14_scale_free_genotype_holdout/final_classification.json`

`results/14_scale_free_genotype_holdout/final_artifact_hashes.csv`

Machine readable outputs must preserve every root and every calcium value.

A figure is optional. The CSV and JSON outputs are mandatory.


## 14. Required tests

Add focused tests proving at minimum:

* the Task 14 branch freeze points to `2f54e7c4f87b6da746d3a8427c7bdca40a77c3de`;
* all ten root payloads are unchanged;
* no whole cell parameter differs among WT, AE4 null and AE2 null members except the intended expression field;
* the calcium panel contains exactly 0.10, 0.25 and 0.50 uM;
* the R1 regulatory member is exactly `R1_G125_P21_PROSE_TREFERENCE`;
* exact AE4 deletion makes every AE4 conserved source exactly zero;
* exact AE2 deletion makes the AE2 flux exactly zero;
* genotype resting continuation fits zero parameters;
* WT and null paired trajectories share the exact time grid and solver;
* paired secretion ratios are invariant under a common artificial multiplicative observation scale;
* blind prediction code contains no holdout reader or target interval;
* the reveal evaluator cannot alter blind trajectory files;
* all final summaries retain all roots and all calcium values;
* repeated evaluation of saved blind artifacts is deterministic.

Run the relevant inherited genotype continuation, conservation and Task 13C freeze tests as regression tests.


## 15. Final classification

End Task 14 with exactly one of the following classifications.

### Outcome 1

`AE4 PHENOTYPE AND AE2 SPECIFICITY RECOVERED WITHOUT RETUNING`

Use when the frozen blind predictions are quantitatively compatible with the held out AE4 total secretion and time structure across a scientifically coherent portion of the predeclared root and calcium panel, and AE2 remains comparatively small, with no target dependent selection or retuning.

### Outcome 2

`AE4 PHENOTYPE RECOVERED BUT AE2 SPECIFICITY FAILS`

Use when AE4 agreement is convincing but the same frozen model predicts an AE2 effect inconsistent with the available evidence.

### Outcome 3

`AE4 DIRECTION RECOVERED BUT MAGNITUDE OR TIMING DISAGREES`

Use when AE4 deletion lowers secretion in the correct qualitative direction but the held out magnitude, temporal structure or both are materially inconsistent with the blind prediction.

### Outcome 4

`GENOTYPE PREDICTION IS MATERIALLY ROOT OR CALCIUM DEPENDENT`

Use when the scientific conclusion changes across the predeclared retained roots or calcium panel strongly enough that no single mechanism level conclusion is justified without an additional independent measurement.

### Outcome 5

`AE4 PHENOTYPE NOT RECOVERED BY THE FROZEN MODEL`

Use when the frozen model gives little effect, the wrong direction or otherwise clearly fails the held out AE4 phenotype across the valid predeclared panel.

### Outcome 6

`GENOTYPE CONTINUATION OR NUMERICS PREVENT A VALID HOLDOUT TEST`

Use when exact null resting states or stable paired dynamics cannot be obtained under the predeclared numerical rules, so the biological comparison is not licensed.

Do not force Outcome 1 merely because some root or calcium value happens to match. The complete predeclared panel matters.


## 16. Final repository actions

After the one time reveal and evaluation:

* run all focused Task 14 tests and relevant inherited regressions;
* preserve the blind checkpoint SHA and hashes in the final artifacts;
* commit the reveal comparison and final report separately from the blind prediction commit;
* push the final branch;
* do not merge to `main`;
* do not edit `archive/`.

The final answer must state plainly what the model predicted before reveal, what the experiment showed after reveal, where they agree, where they disagree, whether AE2 specificity survives, and what scientific uncertainty remains.
