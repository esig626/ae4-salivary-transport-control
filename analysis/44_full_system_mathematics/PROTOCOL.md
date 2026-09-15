# Task 44: full system mathematical analysis

Base main: c4d4f207f702d8eee09ab94d2d90ae90552dd641.
Companion branch: analysis/task-43-parameter-provenance.

Preserve every intracellular and luminal amount, both volumes, finite buffering, algebraic pH and electrical closure, and AE4 regulation. An exact charge-manifold coordinate change is allowed; freezing carbon, alkalinity, volume, lumen composition or voltage is not. No edits to frozen model code, existing numerical results, or manuscript text.

## Analysis

1. Independently verify Task 40 numerical sources, frozen parameters and resting state. Document subsequent metadata changes separately from numerical source changes.
2. Use bath sodium as concentration scale, the frozen cell and lumen resting volumes as their respective volume scales, and initial stimulated WT NKCC1 cycle flux as flux scale. Verify exact dimensional/dimensionless trajectory equivalence and recover all dimensional quantities.
3. Prove uniqueness and smoothness of the carbonate closure and the coupled electrical closure on their stated domains. Establish only local well posedness, not unproved global positivity or stability.
4. Prove the full alkalinity balance, its conditional capacity obstruction, the finite-time budget identity, the pH derivative including carbon and volume, and the general AE4 source-vector projection.
5. The cancellation of an explicit AE4 term does NOT prove zero total AE4 sensitivity, a small knockout, or impossibility of another parameterisation. Equal projected source signatures do NOT prove observational equivalence.
6. Extend the stimulus as a separate constant-input experiment. Solve equilibria on the exact neutral charge manifold and analyse its 11 independent dynamic coordinates. Compare Jacobians at two difference steps. Evaluate AE4-expression continuation and all seven predeclared source classes; retain nonphysical or unsuccessful cases. No global uniqueness assertion follows from local roots.
7. Compute implicit equilibrium sensitivities and validate them against independently solved nearby roots. Define chloride compensation as -2(dJ_N/de)/(dJ_4/de), not the cycle-flux ratio lacking the factor two. Distinguish total derivatives, relative increases, and integrated compensation fractions.
8. Interpret eigenvalues only locally. Existing slow modes, transient storage, nonnormality, input and output projection must be considered before claiming a missing regulatory timescale.
9. Keep Task 41 as an inverse construction. Search the one prescribed residual-recruitment family for the first numerical threshold reaching 30.3% cumulative loss at 600 s, preserving all inherited physiological gates and shared WT initial state. Do not claim a global lower bound without proving monotonicity/global exclusion. Keep adaptive integral and the original sampled trapezoidal observable distinct.
10. Write and compile a standalone TeX report, with scripts, compact machine-readable results, verification, assumptions and unresolved limitations. Preserve original results and leave both branches unmerged for later manuscript integration.

## Checkpoint discipline

The local workspace cannot resolve github.com for direct Git access. A GitHub Actions source bundle supplies exact history; connector reads and writes provide remote checkpoints. Check remote refs before publication and never force update a moved branch. Source and generated-result checkpoints must be remote before reporting completion.
