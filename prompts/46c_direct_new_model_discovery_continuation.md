# Task 46C: direct new-model discovery continuation

This instruction supersedes any Task 46 step that would spend additional time reproducing, re-verifying, or forensically re-auditing the already preserved Task 43/44 baseline unless a new candidate model specifically requires a baseline comparator.

## Start point

Resume from the published Task 46 checkpoint 01 at:

`83e69048a329135b95ef00f007610610ff8afc6a`

Treat Task 43, Task 44 and Task 46 checkpoint 01 as authoritative failure evidence. Their purpose is to explain why the current architecture is inadequate, not to become the main object of further work.

Do not rerun the old model merely to confirm known values again.

## Immediate objective

The next scientific action is to **construct a genuinely new model**.

The existing Task 40/44 whole-cell model and its Palk/Benjamin-derived NKCC1 law are baselines only. They do not constrain the search space.

Use the preserved results to identify the failure that the new model must correct:

- AE4 loss is almost completely buffered by NKCC1 in the current architecture;
- stationary null secretion loss is about 0.946%;
- 600 s null cumulative loss is about 3.857%;
- local NKCC chloride compensation is about 92.9%;
- finite-time replacement of missing AE4 chloride is about 88.5%;
- isolated salivary NKCC evidence does not support a large genotype-dependent increase in NKCC1 after AE4 loss;
- the current model therefore has too much compensatory freedom and is biologically inadequate for the phenotype of interest.

These results are the **reason for replacing the model**, not a reason to spend further milestones reproducing it.

## Required work now

1. Freeze the calibration/holdout evidence split immediately.
2. Independently search the literature for salivary acinar chloride-loading, NKCC1 kinetics/regulation, Na/K pump coupling, NHE1, NBC, AE4, CaCC/TMEM16A, trafficking/phosphorylation, volume regulation and gland-level effects relevant to the phenotype.
3. Generate genuinely new mechanistic candidate models. Do not limit candidates to equations already present in this repository or suggested in earlier prompts.
4. For NKCC1, construct at least three materially distinct candidates, including:
   - one literature-derived mechanistic or carrier-state model if available;
   - one thermodynamically reversible finite-capacity saturating model derived from first principles;
   - one parsimonious salivary-specific effective model if the literature does not provide an adequate mechanistic law.
5. For each candidate, write the equations and provenance before simulation.
6. Evaluate one candidate family per remote checkpoint. Push after every candidate family.
7. Reject models whose AE4-null phenotype is achieved mainly by unsupported parameter tuning or by shifting the problem into another unphysical compensatory pathway.
8. Select the smallest biologically defensible whole-cell architecture that reproduces WT physiology and a substantial AE4-null secretion phenotype while improving independent holdouts.
9. Only after the new model is frozen, repeat the Task 44 mathematical analysis on that new model.

## Scientific selection rule

Use this order:

1. mathematical and conservation validity;
2. thermodynamic validity;
3. WT physiological validity;
4. qualitative agreement with holdout evidence;
5. absence of extreme unsupported compensation;
6. parsimony of assumptions and parameters;
7. quantitative phenotype agreement.

Do not optimise a single weighted objective across all data.

## Checkpoint discipline

At every scientific milestone:

1. fetch remote state;
2. fast-forward only;
3. verify the completed stage;
4. update `CURRENT_STATUS.md`;
5. create an immutable `CHECKPOINT_XX.md`;
6. commit;
7. push immediately;
8. verify the remote ref contains the commit;
9. only then continue.

Do not wait until several model families have been tested before publication.

If the run freezes, the next run must resume from the latest remote checkpoint and continue with the next candidate family.

## Final success condition

Success is not 'the old model was reproduced'. Success is one of:

- a new biologically defensible model is discovered, frozen, validated and mathematically analysed; or
- no defensible model survives, in which case publish the negative result and the minimum discriminating experiments needed.

Do not declare Task 46 complete until one of those outcomes is remotely published.