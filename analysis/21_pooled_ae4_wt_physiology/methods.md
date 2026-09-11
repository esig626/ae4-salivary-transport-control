# Task 21 methods and numerical limitations

The Stage 0 physiology contract was pushed and remotely verified at
`9820a6b62b29cd9895b9020b85ae89ad7d44b3fb` before new WT optimization.
The contract JSON, its Git blob, provenance hashes, and the archive tree are
checked by the WT runner. No Task 21 genotype calculation occurs in the WT code.

The unchanged Task 20 pooled-cation evaluator supplies a single net cycle.
All ten inherited WT roots are separate numerical starts and capacity
references. The full neutral cell/lumen state, algebraic voltages, and the
15 declared capacities can move. Geometry **parameters**, buffer amounts,
bath, water coefficients/laws, chemistry and regulation remain fixed. There
is no state-distance penalty, chloride-share condition or legacy-flux target.

## Coordinates and production verification

The initial search uses chloride, pH, TIC, volume and Na/K partition in each
compartment. Electroneutrality determines the cation pool from chloride,
alkalinity and fixed charge. Fifteen independent log capacity folds and two
open voltages complete the variables. The fast residual transcribes the
production material and current balances; it is checked against production.
Every saved candidate is reconstructed and independently evaluated by the
unchanged production model. An optimizer success flag is never an acceptance
criterion. Both omitted alkalinity rows, all independent RHS rows, currents,
state charges and conservation identities must pass the inherited tolerances.
The final inventory uses JSON-reloaded Python float capacities. At these large
fluxes, one initial NumPy-scalar resting pass failed its conservation replay
after serialization. Both diagnostics are retained; the replay failure is
authoritative. Its already completed WT trajectories remain in the audit and
its state/parameter identity is verified when cached trajectories are reused.

Each of the ten roots receives a legacy start, an expanded-capacity start,
and two high-carbon lumen starts. Initial capacity exploration limits are
absolute log folds 8, 16 and 32; equality-constrained follow-ups expand to 40.
These are numerical search limits, not biological feasibility walls. Failed
starts and messages are preserved. The first equality-search batch encountered
an invalid trial state; the resumed search records/rejects invalid states.

## Independent flux search and lexicographic follow-up

At fixed cell concentrations and pH, electroneutrality determines cell volume.
The unchanged two water balances then determine lumen osmolarity, volume and
outflow. Lumen carbon/pH and Na/K partition determine its other concentrations.
At fixed state and voltages the remaining source equations are linear in
capacities. A nonnegative inner solve, normalized by outflow, avoids mistaking
uniformly vanishing rates for a steady solution.

The first exploratory `flux_search` outputs used subtraction of open residuals
to construct linear columns. Large outflow caused cancellation in that
construction. Those rejected exploratory outputs are retained as diagnostics,
but are superseded by analytic stoichiometric columns in `linear_system`.
`mean_flux_search` and all subsequent profiles use analytic columns. Production
verification, rather than the inner solver residual, decides acceptance.

Differential evolution searches the positive-carbon lumen branch independently
for all ten parameter backgrounds, first across the screening bands and then
at Cl = 50.10 mM and pH = 6.91. The latter reaches zero measured-error objective.
The subsequent joint state/voltage/minimax search profiles positive capacities
between exp(-t) and exp(t) times their references with a linear feasibility
solve. Its numerical range is 0 <= t <= 32. A capacity-only log-square
tie-break at the resulting state is followed by a joint tie-break with all
state and voltage coordinates free. The joint minimax search uses 100 population
members and 250 generations per seed; the mean-feasibility search uses up to
400 generations. Seeds and bounds are in the machine-readable outputs.

The optimization hierarchy uses WT measured error, then maximum absolute
capacity log fold, then summed squared log folds. Failed follow-ups do not
erase valid earlier WT candidates. We preserve and dynamically test every
distinct saved candidate that passes the same production resting gates,
including high-fold feasibility witnesses. Such witnesses are explicitly
not certificates of a global capacity optimum. Numerical search limits,
branch-focused robustness searches, optimizer nonconvergence and incomplete
tie-break convergence prevent a global optimality or infeasibility claim.

## Dynamic gate, retention and firewall

Deduplication compares complete states and complete parameter vectors using
the predeclared componentwise relative tolerance 1e-6. Distinct immutable
parameter backgrounds remain distinct. Every rest-valid candidate receives
600 s production Radau WT trajectories at Ca 0.10, 0.25 and 0.50 uM, with the
inherited right-limit time grid and 1 s sampling. All original solver and
conservation tolerances remain unchanged. Each full trajectory must have
positive finite core states, nonnegative flow, valid regulatory fractions,
charge closure and all conservation residuals within their native-unit limits.
Trajectory failure does not prove that the production equations have no exact
solution. In particular, large balanced currents can expose floating-point
and integration limitations under the inherited absolute tolerances.

The complete WT inventory, failed seed decisions, states, parameters, resting
ledgers, raw trajectories and gates are hashed in `wt_manifest.json`, committed,
pushed and remotely verified. Only a nonempty fully admissible frozen ensemble
can authorize genotype evaluation. If it is empty, AE4-loss and AE2-loss
tables contain explicit unavailable values, never inferred or fabricated
ratios. No loss magnitude is compared with experiment in that case.

**Ordering deviation:** a broad repository regression invocation executed
historical genotype tests before the Task 21 WT checkpoint was pushed. This
did not follow the user's requested ordering for *any* phenotype evaluation.
All Task 21 optimization had already completed; no Task 21 genotype was
calculated and no WT state, capacity, bound or selection was changed from
those historical tests. The deviation is recorded in
`validation_scope_deviation.json` and the broad test transcript is retained.
The broad suite also failed on unavailable historical files/Git objects in
this partial-object checkout. Focused Task 21 checks are reported separately;
the broad suite is not represented as passing.

## Reproduction

Run from this branch with `PYTHONPATH=src`. Stage 0 is immutable; do not rerun
its preparation after it exists. The WT commands, in order, are:

```text
python -m modern_full_model.task21_calibration --workers 4
python -m modern_full_model.task21_flux_search
python -m modern_full_model.task21_flux_search --means
python -m modern_full_model.task21_refine
python -m modern_full_model.task21_exact_calibration
python -m modern_full_model.task21_profile_calibration
python -m modern_full_model.task21_joint_tie_break
python -m modern_full_model.task21_wt_checkpoint dynamic --workers 4
python -m modern_full_model.task21_wt_checkpoint freeze
```

The present analytic flux-search implementation supersedes the documented
initial subtraction-column experiment, so replay need not duplicate those
superseded exploratory numbers. Production gates and conservation tolerances
are authoritative in either case. Frozen outputs must not be overwritten;
use an isolated, explicitly unfrozen reproduction checkout. The post-checkpoint
unavailable-result writer requires the remotely verified receipt and refuses
to replace genotype testing if any WT solution is admissible.
