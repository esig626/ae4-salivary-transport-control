# Task 46D addendum: anti-combinatorial model-search discipline

This addendum is mandatory for Task 46 and supersedes any interpretation of Tasks 46/46A/46B/46C that would lead to a Cartesian-product search across transporters, regulatory mechanisms, parameterisations or model families.

## Core rule

Task 46 is a hypothesis-driven model reconstruction, not an exhaustive combinatorial search.

Do not enumerate all combinations of NKCC1, pump, NHE1, NBC, AE2, AE4, CaCC, water transport and regulatory variants. Do not create a grid of architecture permutations. Do not run large brute-force sweeps merely because the space is available.

Every new candidate must answer one specific mechanistic question raised by the previous checkpoint.

## Search tree, not Cartesian product

Use a sequential decision tree:

1. diagnose the dominant failure mechanism;
2. modify only the smallest mechanistic block needed to address that failure;
3. reject or accept that block using predeclared gates;
4. freeze the accepted block;
5. move to the next unresolved mechanism only if necessary.

Once a component is frozen, do not reopen it merely to explore combinations unless a later result gives specific evidence that the frozen choice is incompatible with the whole-cell model.

## NKCC1 stage

The NKCC1 reconstruction stage must be compact.

Use at most four primary NKCC1 model families unless a fifth is justified by a qualitatively distinct mechanism that cannot be represented by the first four.

The four preferred categories are:

1. current Task 44 law as failure baseline only;
2. one best-supported literature mechanistic model;
3. one independently derived thermodynamically reversible finite-capacity model;
4. one parsimonious salivary-specific effective model if needed.

Do not create multiple trivial variants of the same family and count them as distinct models.

Within each family, calibrate only the minimum number of free parameters needed to reproduce declared WT construction constraints. Do not sweep every parameter combination.

If a family fails a hard gate, stop exploring that family.

## Hard early rejection gates

Reject a candidate immediately if any of the following occurs:

- conservation or charge failure;
- wrong NKCC1 stoichiometry;
- wrong thermodynamic reversal;
- non-finite or unbounded compensatory flux;
- failure of WT resting physiology;
- failure of WT stimulated physiology;
- strong AE4-null NKCC compensation incompatible with the declared isolated-assay constraint;
- pathological Na, K, Cl, pH, volume or membrane potential;
- phenotype recovery only through an unsupported extreme coupling.

Do not continue tuning a rejected family.

## Parameter-search discipline

Do not use high-dimensional brute-force optimisation or exhaustive sweeps.

For each candidate model:

1. identify the smallest set of genuinely free parameters;
2. derive or fix all others from literature, conservation, equilibrium, closure or WT physiology where possible;
3. perform one low-dimensional calibration problem;
4. use local sensitivity/continuation to understand identifiability;
5. only expand the parameter set if a specific structural deficiency remains.

Prefer analytical constraints, continuation, monotonicity arguments and one-dimensional or two-dimensional searches over generic global optimisation.

Do not tune more than three new effective parameters simultaneously without first publishing a justification explaining why fewer cannot work.

If a three-parameter fit is non-identifiable, do not add more parameters. Reduce or redesign the model.

## Whole-cell reconstruction stage

After NKCC1 is frozen, do not branch simultaneously over pump, NHE1, NBC, AE4, CaCC and water transport.

Use the causal diagnosis to rank remaining deficiencies. Address only the highest-ranked unresolved mechanism first.

At each stage:

- one structural hypothesis;
- one candidate architecture;
- one checkpoint;
- one accept/reject decision.

Only move to a second mechanism if the first frozen repair still cannot satisfy the declared constraints.

Do not combine two unvalidated repairs in the same candidate unless the equations demonstrate that they are mathematically inseparable.

## Candidate budget

Before detailed simulation, maintain a candidate ledger with a hard budget:

- no more than 4 primary NKCC1 families;
- no more than 3 whole-cell architectural hypotheses active at one time;
- no more than 2 versions of any single architecture before it must be either accepted or abandoned.

Exceed these limits only if a checkpoint explicitly demonstrates that the surviving evidence cannot discriminate the existing candidates and a new qualitatively distinct mechanism is required.

## Escalation rule

If no compact candidate survives, do not respond by multiplying the search space.

Instead:

1. stop;
2. publish the negative result;
3. identify which observation cannot be reconciled;
4. identify the smallest missing biological mechanism or experiment needed;
5. only then define one new hypothesis.

## Checkpoint discipline

Every candidate family or architecture decision must be remotely published before another branch of the search tree is started.

The checkpoint must contain:

- the mechanistic hypothesis;
- frozen equations;
- free parameters;
- construction data used;
- rejection gates;
- result;
- accept/reject decision;
- reason for the next branch of the search tree.

Do not maintain unpublished parallel candidate workspaces.

## Final model-selection principle

The objective is not to find the numerically best member of a large search space.

The objective is to find the smallest mechanistically defensible model that satisfies the declared physical, physiological and experimental constraints and whose behaviour can be mathematically explained.

A compact model with a clear causal mechanism is preferred over a larger model with a slightly better numerical fit.
