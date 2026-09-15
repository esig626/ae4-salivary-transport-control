# Task 46B addendum: independent model discovery and synthesis

This addendum is mandatory when executing Task 46 and Task 46A.

## Purpose

Task 46 is not an exercise in selecting among mechanisms already proposed by the user, the existing repository, Marty.tex, Tasks 37--44, or the current manuscript. Those materials define the failure to explain and the evidence that must be respected. They do **not** define the allowable model space.

Codex must act as a model developer, not a summariser.

The scientific objective is to discover or construct a better mechanistic model that explains the AE4 phenotype while satisfying conservation, thermodynamics, WT physiology and independent experimental holdouts. Existing equations and user-suggested mechanisms are hypotheses and baselines only.

## Independent discovery requirement

Before selecting any replacement architecture:

1. Search the primary literature independently for mechanistic information relevant to salivary acinar chloride loading, NKCC1 kinetics/regulation, AE4, NHE1, NBC-like bicarbonate transport, Na/K ATPase, CaCC/TMEM16A, cell-volume regulation and stimulus-dependent transporter trafficking/regulation.

2. Do not restrict the literature search to papers or mechanisms already cited in the repository. Follow useful references backward and forward where possible and record the search trail and source hierarchy.

3. Separate direct salivary-acinar evidence from evidence in other epithelia, neurons, erythrocytes, kidney, secretory cells or heterologous expression systems. Cross-tissue mechanisms may motivate a candidate, but they must not silently become salivary facts.

4. Produce a `MODEL_DISCOVERY_LEDGER` before phenotype fitting. For each newly identified mechanism or modelling principle record:
   - source and biological preparation;
   - what is directly measured;
   - what is inferred;
   - mathematical form suggested by the evidence;
   - whether it constrains kinetics, capacity, regulation, stoichiometry, localisation or coupling;
   - whether it is suitable for a candidate model and why.

5. Explicitly identify mechanisms and model structures not previously considered in Tasks 37--44. A candidate panel consisting only of variants of the existing Palk/Benjamin law, arbitrary capacity ceilings, or user-suggested couplings is insufficient.

Commit and push the discovery ledger before using the AE4-null secretion magnitude to compare candidate fits.

## First-principles model synthesis requirement

If no published kinetic model is demonstrably appropriate for salivary acinar NKCC1 or another critical transporter, construct a minimal mechanistic model from first principles rather than defaulting to the inherited law.

For NKCC1 specifically, derive at least one candidate from the following physical requirements:

- exact 1 Na : 1 K : 2 Cl stoichiometry;
- electroneutral transport;
- chemical-potential/thermodynamic driving force with the correct reversal condition;
- finite transporter abundance or finite turnover capacity;
- bounded saturation as substrate concentrations increase;
- physically meaningful behaviour as one substrate becomes limiting;
- no divergence or artificial increase in flux solely because intracellular chloride falls;
- detailed balance or a clearly stated thermodynamic consistency condition;
- a transparent distinction between transporter abundance, intrinsic turnover and stimulus-dependent regulation.

A carrier-state/Markov-cycle model, reversible binding-cycle model, or reduced thermodynamically consistent occupancy model may be used if justified. The number of states must be the minimum needed to express the relevant biology. Do not create kinetic complexity merely to gain fitting freedom.

Derive the equations in the repository. Do not simply cite an equation and paste it into the model.

## Genuine candidate diversity

Before model selection, evaluate at least three **structurally distinct** NKCC1 kinetic architectures unless the literature establishes one uniquely superior model. Parameter variants of the same equation do not count as distinct architectures.

The candidate set should contain, where scientifically defensible:

1. the current Task 44 law as the documented failure baseline;
2. a newly derived first-principles thermodynamic finite-capacity model;
3. at least one independently sourced mechanistic or regulatory architecture from the literature;
4. additional candidates discovered during the search if they materially differ in mechanism.

The same principle applies to any other subsystem found to be structurally responsible for the phenotype. Codex is allowed to replace an inherited law, introduce a justified regulatory state, alter coupling topology, or remove an unsupported phenomenological multiplier if the evidence supports doing so.

Codex is **not** required to preserve the current NKCC stimulation multiplier, current capacity calibration strategy, present state dependence, or current algebraic functional form.

## Model invention beyond the existing hypotheses

Codex may propose a new mechanistic architecture not explicitly present in the literature when all of the following hold:

- the architecture is assembled from experimentally supported biological components or physical principles;
- conservation, stoichiometry, charge and thermodynamics are exact;
- every new assumption is declared;
- the model is simpler than alternative constructions offering the same explanatory power;
- calibration and holdout evidence remain separated;
- the proposed mechanism yields testable predictions that distinguish it from alternatives.

Label such a construction as a **new model hypothesis**, not as established biology.

If a novel model is necessary, derive the mechanism before examining its detailed fit to holdout data. Preserve unsuccessful novel models in the rejected-candidate ledger.

## Do not optimise toward the answer

The target AE4 secretion deficit is not an objective function from which to reverse engineer arbitrary equations.

Candidate equations, state variables, regulatory topology and thermodynamic structure must be defined and remotely checkpointed before quantitative phenotype fitting within that candidate family.

Within a predeclared candidate, parameters may be estimated only from the declared calibration data and supported bounds. Do not mutate the mathematical form after observing a holdout failure without versioning the candidate and recording which holdout evidence has thereby been consumed.

A model that reaches the desired secretion deficit for the wrong biological reason must be rejected.

## What "better model" means

The accepted reconstruction must improve explanatory adequacy, not merely numerical fit. In order of priority, it should:

1. obey exact physical and mathematical constraints;
2. reproduce WT rest and stimulation;
3. avoid the unsupported NKCC1 recovery of missing AE4 chloride;
4. produce a substantial AE4-loss phenotype through a mechanistically intelligible pathway;
5. predict independent holdout observations better than Task 44;
6. use fewer or better-justified free assumptions than competing models;
7. expose identifiable experimental predictions that could falsify it.

A numerically poorer phenotype fit may be preferred to a closer fit if it is substantially better supported mechanistically and by holdouts.

## Required explanatory output

For every surviving architecture, Codex must answer in causal terms:

- Why does WT secretion work?
- What role does AE4 play in WT?
- What changes immediately when AE4 is lost?
- Why can NKCC1 not simply recover the lost chloride in this model?
- Which state or flux becomes limiting?
- Why does secretion fall?
- Why does the phenotype develop on the predicted timescale?
- Which independent experiments agree or disagree?
- What experiment most cleanly distinguishes this architecture from the next-best candidate?

If these questions cannot be answered from the equations and computed flux/state responses, the model is not sufficiently explanatory to accept.

## Mathematical reanalysis

After freezing the accepted model, repeat the full Task 44 programme on **that model**, not on a reduced surrogate. If the accepted architecture introduces new dynamic states, include them in the nondimensionalisation, charge-manifold analysis, equilibrium problem, Jacobian, eigenvalue spectrum, implicit sensitivities and compensation decomposition.

Use the mathematical analysis to explain why the new architecture succeeds where Task 44 failed. Do not merely report that its simulated curve is closer to experiment.

## Additional remote checkpoints

In addition to the Task 46/46A checkpoints, push immediately after:

1. independent literature/model discovery ledger;
2. first-principles NKCC1 derivation before simulation;
3. structurally distinct candidate equations frozen before phenotype comparison;
4. candidate discrimination results;
5. any newly invented model architecture before parameter fitting;
6. accepted mechanistic explanation before full Task 44-style reanalysis.

At every checkpoint fetch remote state, pull with `--ff-only`, verify, commit and push immediately. Never force-push.

## Failure is allowed

If independent model discovery does not yield an architecture that satisfies the physical, physiological and holdout constraints, publish that result. Do not retreat to the inherited model merely because it is available, and do not manufacture a coupling to obtain the desired number.

The final report must clearly distinguish:

- mechanisms inherited from earlier work;
- mechanisms discovered in the literature during Task 46;
- mechanisms derived from first principles;
- genuinely new hypotheses introduced by Task 46;
- evidence used for calibration;
- evidence retained as independent validation.
