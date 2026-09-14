Task 41: algebraic design of an AE4-loss secretion deficit

The new user instruction, after completion of Task 40, is:

> your new task is to figure algebraically what do we need to get the result, dont be rigid a loss of 20-35% vs WT is enough. once you find it, then amend the model and simulate it.

This authorizes target-directed algebra, a model amendment, and simulations
as new work. It supersedes the completed Task 40 prohibition on further
mechanism changes and target-directed selection for this new task. Task 40
source, frozen states, results, instructions and completion ledger remain
unchanged. The inherited branch and GitHub publication restrictions remain:
codex/task-40-ae4-equal-cation-routing, connected GitHub remote writes, no main
merge. The scientific parent for this work is Task 40 commit
76f3144a4e662eb12d6da26c83ebf6f2ff47a59d.

Aim: a 20-35% reduction of 0-600 s cumulative secretion for AE4 null versus
the amended WT; also evaluate AE4=0.05 and report whether it falls in that
range. The target is now used for model design and is NOT held-out validation.

Plan before new simulations:

1. Derive exact Na/alkalinity and cell+lumen chloride balance identities.
2. Derive a shared finite NKCC1 capacity from a coupled algebraic stationary
   WT/null design problem, using the range midpoint (27.5%) as a guide.
   No trajectory is run during parameter selection. The algebraic guide
   concerns sustained flow; the requested finite-horizon cumulative outcome
   must still be checked dynamically.
3. If the algebraic solution is admissible, implement the shared capacity
   restriction on complete 1 Na:1 K:2 Cl NKCC1 cycles. It has no AE4-expression
   term. Retain equal AE4 routing and every other Task 40 mechanism.
4. Solve one amended WT rest from Task 40 WT rest if the cap changes rest.
   Freeze it and use the same state for WT, 5%, and null 600 s trajectories.
5. Keep the Task 40 physical/conservation and WT activation gates. Report the
   absolute WT change as well as relative deficits. Do not label target
   attainment as independent evidence for the added capacity restriction.

One scientific agent/worker and one BLAS thread. Algebraic stationary design
states are sustained-stimulus states; they are not genotype-specific resting
initial conditions for the trajectories. This task does not identify a unique
biological mechanism from one secretion ratio.
