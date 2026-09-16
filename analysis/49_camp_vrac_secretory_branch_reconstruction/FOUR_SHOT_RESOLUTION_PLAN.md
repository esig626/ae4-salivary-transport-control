# Four-shot scientific resolution budget

There are four remaining Codex scientific execution runs. This file is binding project memory so later tasks do not reset the investigation or repeat excluded mechanisms.

## Shot 1 — Task 49

Use one session for two independent architectural failures:

1. prospectively reconstruct the missing TMEM16A-independent beta/cAMP volume-regulated apical chloride pathway using independent cAMP-pathway data, with **no simultaneous beta-NKCC modification**;
2. after publishing that result, solve the measured WT/AE4-KO/AE2-KO resting-state residual/source-signature inverse problem in conserved coordinates.

The resting-state inverse problem is **not a subset search**. Form the residual vectors at the measured states and the Jacobian/source-signature matrix for the already represented pathways. Use rank/SVD/projection and, where needed, one constrained continuous least-squares solve to ask whether the measured residual lies in the span of the existing source directions and which direction is identifiable. Do **not** enumerate transporter subsets, mechanism combinations, on/off patterns, parameter subsets, source-vector combinations, or candidate families. Do not run all pairs, all triples, best-subset regression, sparse combinatorial selection, Cartesian grids, or evolutionary/random mechanism searches. If the residual geometry does not identify a unique correction direction, publish non-identifiability and stop that correction stage.

Outcome required: either a validated stimulated branch plus a localized resting correction, or a mathematically explicit source/equation deficiency. Do not stop after the first negative subtest.

## Shot 2 — Task 50

Implement exactly the one resting/equation correction identified by Task 49, if one is uniquely supported. Combine it with any independently fixed Task 49 cAMP branch **without retuning either to saliva output**. Fit only admissible shared parameters to non-secretion cellular data. Freeze the entire model before the AE4 whole-gland phenotype reveal.

No new mechanism family is allowed unless Task 49 proves the existing source-signature basis cannot span the measured states and identifies the missing source direction.

## Shot 3 — Task 51

If the frozen Task 50 model still fails, perform one global residual/discrepancy analysis using the complete experimental ledger and all repository exclusions. Determine the minimum equation-level discrepancy required to satisfy the remaining constraints. Permit at most one new source-supported correction, selected from the residual geometry, not from a menu or parameter sweep. Freeze before whole-gland phenotype comparison.

This is the last scientific reconstruction shot.

## Shot 4 — Task 52

Final verification and article lock. No open-ended exploration. Reproduce all accepted predictions, run independent conservation/state checks, quantify uncertainty/sensitivity around the accepted reconstruction, merge the canonical results, and write the article-facing model narrative.

If Task 51 still leaves a contradiction, Task 52 must state the rigorous negative conclusion and exactly which experimental constraints prove the model class inadequate. Do not spend the final run inventing another mechanism.

## Global prohibitions

Across all four runs:

- no rerunning excluded Tasks 12-48 ideas;
- no broad mechanism searches or Cartesian sweeps;
- no combinatorial subset enumeration of transporters, mechanisms, source directions, parameters, or corrections;
- no all-pairs/all-triples/all-subsets search and no best-subset or evolutionary/random model selection;
- residual localisation must use direct algebra, rank/SVD/projection, sensitivities, or one continuous constrained solve, not discrete search;
- no fitting the ~35% AE4-null secretion deficit in isolation;
- no genotype-specific fudge multipliers;
- no post-reveal retuning;
- no unverified claims that a negative result identifies a mechanism;
- publish every dependency boundary immediately to GitHub with restart metadata.

The repository, not conversational memory, is the authoritative cumulative exclusion and decision record.