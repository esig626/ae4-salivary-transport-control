# Task 45: agentic completion of parameter provenance and full-system mathematics

Repository: `esig626/ae4-salivary-transport-control`

This is a continuation and completion task. Do not restart Tasks 43 or 44 from scratch.

## Existing remote checkpoints

- Parameter provenance branch: `analysis/task-43-parameter-provenance`
  - known checkpoint: `693f44c0b05e92faec737aa7793c3a99ec76445c`
  - existing report: `analysis/43_parameter_provenance/report.tex`
  - existing compiled report: `analysis/43_parameter_provenance/AE4_parameter_provenance.pdf`
- Full-system mathematics branch: `analysis/task-44-full-system-mathematics`
  - known checkpoint: `b7b775a0277d93263c4de339745806e5658ec9c2`
  - existing protocol: `analysis/44_full_system_mathematics/PROTOCOL.md`
  - existing first numerical checkpoint: `analysis/44_full_system_mathematics/CHECKPOINT_01.json`
  - existing code includes `core.py`, `run_analysis.py`, `verify.py`, and `prepare_report.py`
  - this branch is not complete until it contains a final standalone TeX report and compiled PDF plus verified machine-readable results.

Read `AGENTS.md` first, then read the Task 43 and Task 44 protocols in full. Inspect the current remote refs before doing any work. If either branch has moved beyond the checkpoint above, use the newer remote head and do not overwrite or force-push it.

## Operating mode

Work agentically to completion. Do not stop for routine clarification. Make the most scientifically conservative defensible choice when a detail is underspecified, record that choice, and continue. Ask only if blocked by missing credentials, inaccessible source material, or a decision that would alter the scientific model rather than analyse it.

At every safe checkpoint:

1. `git fetch --all --prune` and inspect the current remote branch head.
2. `git pull --ff-only` on the branch before new work.
3. Run the relevant verification/tests.
4. Commit a coherent checkpoint with a descriptive message.
5. Push immediately.
6. Never force-push.

Keep terminal output compact. Do not dump large trajectory files or giant JSON objects to the terminal.

Do not modify `main`, model source, Tasks 1-42, frozen trajectories, frozen figure source values, or the current manuscript. Keep Tasks 43 and 44 separate and unmerged pending scientific review.

---

# Part A: finish and harden Task 43 on `analysis/task-43-parameter-provenance`

The purpose is to justify every active production parameter without pretending that an accepted physiological point is a measured uncertainty interval.

## Required audit standard

Audit the actual Task 40 production objects and selected laws, not dataclass defaults. For every active model setting record:

- parameter name and active value;
- units;
- where it enters the equations;
- whether the selected production law actually uses it;
- exact code/source location;
- provenance category;
- literature/model source where applicable;
- what mathematically or experimentally constrains it;
- what does **not** justify it;
- whether uncertainty is known, bounded only physiologically, or genuinely unknown;
- whether the conclusion of the current manuscript materially depends on it.

Maintain clear categories at least equivalent to:

1. physical constant;
2. published kinetic coefficient;
3. published model conversion/inherited law parameter;
4. resting-state or conservation-closure derived quantity;
5. architecture-derived quantity;
6. effective modelling assumption;
7. inherited effective/model setting;
8. protocol setting;
9. numerical bracket/regularisation/constructor seed;
10. inactive legacy field;
11. target-selected parameter.

Do not convert model assumptions into literature parameters by rhetoric. In particular:

- the proposed 70/30 loading partition used in deriving NBC capacity is an assumption, not a measurement;
- using a published kinetic law does not establish salivary carrier abundance;
- Task 41 parameters remain explicitly target-selected;
- an accepted physiological point is not an uncertainty interval.

## Parameter robustness work

Use the mathematical branch only for genuinely new sensitivity calculations. On Task 43 itself, report provenance and defensible admissible ranges only where they can be sourced or logically derived. Do not invent confidence intervals.

Cross-check the existing Task 43 inventory against the current production configuration and Task 44 sensitivity results. If Task 44 identifies a conclusion-sensitive effective parameter, ensure the Task 43 report flags that parameter prominently and states what future measurement would constrain it.

## Task 43 deliverables

Update as needed:

- machine-readable parameter inventory;
- audit summary;
- `report.tex`;
- compiled `AE4_parameter_provenance.pdf`;
- a concise `FINAL_STATUS.md` stating what is genuinely justified, what remains assumed, and which parameters are most important for the article.

Compile the report cleanly with no unresolved references/citations. Push the final verified Task 43 checkpoint.

---

# Part B: complete Task 44 on `analysis/task-44-full-system-mathematics`

This is the main mathematical analysis. It must remain biologically honest and must not repeat the reduction in `Marty.tex` that froze the very acid-base variables now central to the paper.

## Non-negotiable state-space requirement

Retain:

- intracellular sodium, potassium and chloride amounts;
- intracellular total inorganic carbon and total alkalinity;
- intracellular volume;
- luminal sodium, potassium and chloride amounts;
- luminal total inorganic carbon and total alkalinity;
- lumen volume;
- AE4 regulatory state(s);
- algebraic pH closure;
- algebraic electrical/current closure.

An exact charge-manifold coordinate elimination is allowed. Freezing carbon, alkalinity, either volume, lumen composition, or membrane voltage is not allowed unless a scaling argument is first derived and quantitatively demonstrates an asymptotically small parameter. Do not assume constant cell volume merely for convenience.

## 1. Verify the production system

Independently verify that the Task 44 implementation corresponds exactly to the current Task 40 production chassis and frozen parameters. Record hashes/source commits. Distinguish metadata changes from scientific source changes.

Reproduce key Task 40 reference observables before further analysis.

## 2. Nondimensionalisation

Use the existing protocol unless a demonstrable error is found:

- concentration scale: bath sodium concentration;
- cell volume scale: frozen WT resting cell volume;
- lumen volume scale: frozen WT resting lumen volume;
- flux scale: initial stimulated WT NKCC1 **cycle** flux;
- amount scales induced by concentration x volume;
- time scale induced by amount/flux.

Write the full dimensionless equations, not just a table of scales.

Preserve biologically meaningful groups. At minimum report and interpret groups equivalent to:

- lumen/cell volume scale ratio;
- AE4 regulatory timescale ratio `theta_4`;
- AE4/NKCC loading ratio(s);
- NBC/characteristic-flux capacity ratio;
- AE2/characteristic-flux capacity ratio;
- water/ionic timescale groups;
- buffer/carbonate groups;
- any electrical relaxation group that emerges;
- any genuinely small/large parameter revealed by the scaling.

Numerically verify dimensional/dimensionless trajectory equivalence to tight tolerance and show dimensional recovery.

Do not manufacture a singular perturbation argument if the computed groups are O(1).

## 3. Carbonate and electrical closure

Prove only what is supported.

For the carbonate/pH closure, establish existence/uniqueness and smooth dependence on the retained conserved variables over the stated physiological domain if valid. Provide derivative formulae needed later, including pH sensitivity to TIC, alkalinity, and volume.

For the electrical closure, establish the corresponding local uniqueness/smoothness on the analysed domain if valid.

If a global result is unavailable, state the precise local domain instead of exaggerating.

## 4. Formal alkalinity result

For the parent AE4 stoichiometry derive the exact intracellular alkalinity balance and the stimulated steady-state identity

`2 J_4 = J_H + 2 J_B - J_2`

with the repository's actual sign conventions.

Define a dimensionless alkalinity-support quantity such as

`A_support = (J_H + 2 J_B - J_2)/(2 J_4)`

only where the denominator is nonzero and define the degenerate case separately.

Derive the corresponding no-NBC quantity. Establish a **conditional capacity obstruction**: if the physiologically admissible NHE1/AE2 contribution cannot meet the required AE4 alkalinity demand, sustained AE4 loading at that level is impossible without an additional alkalinity source.

Do not claim a universal impossibility unless the admissible bounds are actually proved.

Also derive the finite-time alkalinity budget, so the argument remains valid away from equilibrium and explicitly accounts for storage.

## 5. Equal-routing structural proposition

Prove from the full conserved balances, under the exact equal-routing source convention, the transient identity already present in the manuscript:

`J_Cl,a = 6 P - J_H + 2 d n_Na,i/dt - d n_A,i/dt - d n_Cl,i/dt`

and hence at stationary intracellular state:

`J_Cl,a = 6 P - J_H`.

State this carefully:

- the explicit AE4 cycle disappears from this stationary balance;
- AE4 can still affect the state, pump, NHE1, NKCC1, pH, volume and therefore total secretion indirectly;
- this does **not** prove zero AE4 sensitivity, a small knockout phenotype for every parameterisation, or observational equivalence of mechanisms.

Turn this into a formal proposition and proof in the report.

## 6. Implicit-function sensitivity and network compensation

For constant stimulation define the exact equilibrium problem

`F(x, mu) = 0`

on the charge-constrained independent coordinates.

Where `D_x F` is nonsingular, compute

`d x*/d mu = -(D_x F)^(-1) F_mu`

and for stationary secretion `Q` compute

`d Q*/d mu = Q_mu - Q_x (D_x F)^(-1) F_mu`.

Interpret the decomposition conservatively as direct parameter dependence plus state-mediated network response. Do not call a term biologically direct if `Q_mu` is zero simply because AE4 acts through state variables.

Validate implicit sensitivities against independently solved nearby equilibria at at least two perturbation sizes.

For AE4 expression define NKCC chloride compensation correctly as

`C_N = -2 (d J_N*/d e_4)/(d J_4*/d e_4)`

because NKCC1 carries two chlorides per cycle.

Keep separate:

- local equilibrium derivative-based compensation;
- relative integrated NKCC increase;
- finite-time integrated chloride compensation fraction.

Do not conflate them.

## 7. Equilibria, stability and timescales

Under constant full stimulation:

- solve physiological equilibria for WT and AE4 expression continuation down to null;
- work on the exact neutral charge manifold;
- test local uniqueness by Jacobian nonsingularity only, not global uniqueness;
- compute eigenvalues on independent dynamic coordinates;
- verify Jacobian/eigenvalues at at least two numerical difference steps or by an independent differentiation method;
- report dominant decay times and modal structure;
- inspect nonnormality/input-output projection before interpreting transient observables.

Important existing checkpoint: preliminary calculations found slow local decay times of roughly 346 s in WT and 460 s in the null. Therefore **do not claim the model lacks a slow mode** unless later corrected calculations actually establish that. Investigate how those slow modes project onto secretion and genotype divergence instead.

Where continuation leaves the valid numerical domain, retry from documented nearby seeds/continuation steps. Retain failed numerical/physiological cases and do not reinterpret solver failure as nonexistence.

## 8. General Catalan source-vector analysis

Represent each AE4 source class as

`S^4 = J_4 nu`

with source vector

`nu = (nu_Na, nu_K, nu_Cl, nu_T, nu_A)`

using the repository's exact sign convention.

Derive the linear projections of `nu` into:

- chloride balance;
- alkalinity balance;
- carbon balance;
- current/charge closure;
- the equal-routing stationary identity or its appropriate generalisation.

Determine which source-vector directions are annihilated by particular conservation combinations and which remain observable in those balances.

Be precise: equal projected signatures in one balance do **not** imply global observational equivalence of two molecular mechanisms.

Evaluate all seven frozen Catalan classes under the same constant-stimulus equilibrium framework. Retain thermodynamically opposed, nonphysiological and unsuccessful cases as such. Do not rescue or refit them.

## 9. Task 41 inverse structural result

Keep the Task 41 family explicitly target-selected.

Use `rho = 1 - b` (or an equivalent clearly defined dimensionless recruitment fraction) for the AE4-dependent stimulated CaCC recruitment component.

Within the **one prescribed Task 41 family only**, search numerically for the first crossing at which the 600 s cumulative AE4-null deficit reaches the lower experimental comparison bound of 30.3%, while preserving:

- the inherited WT initial state;
- all inherited physiological gates;
- the original Task 41 mechanism family;
- the distinction between the adaptive integral and the original sampled trapezoidal observable.

Do not call this a global lower bound unless monotonicity plus global exclusion have actually been proved. Otherwise call it the first numerical crossing in the prescribed family.

Quantify the corresponding percentage of stimulated CaCC recruitment that must be AE4-dependent in that family and state exactly what experiment would test it.

Extend selected cases beyond 600 s as a separate diagnostic and report any loss of physiological validity without changing the declared 600 s Task 41 result.

## 10. Parameter sensitivity for the manuscript claims

Use the Task 43 classification to select uncertain/effective parameters rather than blindly perturbing everything.

For the main scientific claims, calculate local logarithmic elasticities and, where computationally feasible and scientifically justified, bounded robustness over defensible admissible ranges for:

- WT secretion;
- WT pH;
- AE4-null stationary secretion effect;
- alkalinity support;
- compensation gain;
- slowest relevant decay mode;
- Task 41 threshold/crossing.

Do not invent ranges. If a parameter has no defensible interval, report local sensitivity only and explicitly say that global robustness cannot yet be claimed.

Identify which uncertain effective parameters could overturn each principal conclusion and feed that result back into Task 43's final status.

## 11. Standalone mathematical report

Only after all calculations and verification are complete, write a self-contained TeX report under `analysis/44_full_system_mathematics/` and compile it to PDF.

The report should be suitable as the technical basis for a later BMB manuscript revision. It should contain:

1. relation to `Marty.tex` and why the old Na/K/Cl-only reduction is not reused;
2. full retained state and algebraic closures;
3. parameter/scaling choices;
4. nondimensional system;
5. dimensionless groups and their numerical values;
6. carbonate/electrical closure results;
7. formal alkalinity proposition/budget;
8. equal-routing proposition and proof;
9. equilibrium framework and local well-posedness statements;
10. IFT sensitivity formula and validation;
11. NKCC compensation analysis;
12. equilibrium continuation and local stability/eigenmodes;
13. general source-vector/Catalan mechanism analysis;
14. Task 41 inverse-family crossing and its limitations;
15. parameter sensitivity/robustness tied to Task 43;
16. direct comparison with literature experiments: what the model structurally supports, what is in tension, and what measurements would discriminate the unresolved mechanisms;
17. limitations and statements that must **not** be made.

Prefer propositions/lemmas with explicit assumptions over informal claims. Separate exact identities, local numerical results, finite-time simulations and biological interpretations.

The final report must not present Task 41 as independent validation, must not identify the structural NBC pathway with a specific isoform, and must preserve the experimental tension between strong modelled NKCC compensation and the isolated 2015 NKCC assay.

## 12. Verification and publication checkpoints

Before finalising Task 44:

- rerun all analysis from a clean checkout of the branch;
- verify saved machine-readable outputs against recomputation;
- verify exact identities numerically at random valid states where appropriate;
- verify dimensional/nondimensional equivalence;
- verify equilibrium residuals and Jacobian conditioning;
- verify IFT sensitivities by independent finite perturbations;
- compile the TeX report with all citations/references resolved;
- inspect all generated figures for clipping/misleading axes;
- write `FINAL_STATUS.md` with exact commit SHA, completed results, failures, caveats and files produced.

Commit and push each major scientific checkpoint separately rather than waiting until the end.

Suggested checkpoints:

A. recovered/reverified baseline + dimensional/nondimensional equivalence;
B. closure proofs + exact structural identities;
C. equilibria/IFT/compensation/stability results;
D. Catalan source projections + all mechanism equilibria;
E. Task 41 inverse-family crossing + long-time diagnostic;
F. parameter sensitivity crosswalk with Task 43;
G. final TeX/PDF report + verification + final status.

Do not merge either branch to `main`. Do not edit the main manuscript. Stop after both branches are complete, pushed, and independently reviewable.

## Final response

Report only after the remote branches contain the complete deliverables. Provide:

- final head SHA for Task 43;
- final head SHA for Task 44;
- exact paths to both TeX reports and PDFs;
- a compact table of the main mathematical conclusions and their status: exact/proved, local numerical, finite-time numerical, or unresolved;
- any corrected earlier claims;
- parameters/measurements most urgently needing experimental constraint;
- confirmation that no model source, Tasks 1-42, frozen outputs, `main`, or the manuscript were modified;
- confirmation that neither branch was merged.