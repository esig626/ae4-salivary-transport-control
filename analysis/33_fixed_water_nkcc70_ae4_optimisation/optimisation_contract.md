# Immutable Task 33 optimisation contract

Status: immutable once checkpointed. This contract was written before any Task 33 stationary solve, residual evaluation, WT parameter vector evaluation, or dynamic integration.

## Fixed starting point

* Repository: `esig626/ae4-salivary-transport-control`.
* Branch: `codex/task-33-fixed-water-nkcc70-ae4-optimisation`.
* Prepared head: `3a8560d39eb6771fce76895e285c682d5d3333ac`.
* Background: R09 only during WT optimisation.
* The Cha et al. 2009 eight state plus Mod2 NHE1 equations and the carrier amount `2.339370005697548e-05 fmol` are fixed.
* The routed AE4 topology, donor side Na/K partition, 1:1:2 working stoichiometry, membrane allocation fractions, water law, CO2 law, paracellular laws, bath, geometry, and stimulus and regulation equations are fixed.
* The Task 31 R09 WT state and flux ledger are the only numerical reference. Existing knockout records are excluded from the optimisation, objective, feasibility decisions, and initialisation.

## Free parameters and fixed bounds

Exactly the following six existing uncertain magnitudes are free. Their reference values are the Task 31 R09 values. No alternative subset will be tried.

| Parameter | Task 31 value | Multiplier bounds | Absolute bounds |
|---|---:|---:|---:|
| NKCC1 capacity | `0.32 fmol/s` | `[0.25, 2]` | `[0.08, 0.64] fmol/s` |
| Routed AE4 carrier amount | `0.1501264265070657 fmol` | `[0.1, 50]` | `[0.01501264265070657, 7.506321325353285] fmol` |
| AE2 capacity | `0.005 fmol/s` | `[0.25, 4]` | `[0.00125, 0.02] fmol/s` |
| Total Na/K pump capacity | `0.08 fmol/s` | `[0.5, 2]` | `[0.04, 0.16] fmol/s` |
| Total K conductance | `1.4e-08 S` | `[0.5, 2]` | `[7e-09, 2.8e-08] S` |
| Apical CaCC conductance | `3.14e-08 S` | `[0.5, 2]` | `[1.57e-08, 6.28e-08] S` |

For multiplier `m_i`, the optimiser co-ordinate is `z_i = log(m_i)`. Box bounds are the logs of the declared multiplier bounds. Values supplied to the model are `reference_i * exp(z_i)`.

## Algebraic elimination and deterministic initial point

Cell and lumen chloride are eliminated exactly with the existing electroneutral charge manifold, leaving the ten inherited physical stationary co-ordinates. Membrane voltages retain their existing exact linear current solution. No state, voltage, or transporter magnitude is introduced as another outer optimisation variable.

At the accepted Task 31 R09 WT state, the positive NKCC1 and AE4 chloride fluxes are `0.2562444047732706` and `0.004937067441546469 fmol/s`. The initial NKCC1 and AE4 multipliers preserve their sum at that fixed reference state while allocating 70 percent to NKCC1:

* NKCC1 multiplier: `0.7134869177422251`;
* AE4 multiplier: `15.870644383966862`.

The AE2, pump, K conductance and CaCC initial multipliers are exactly `1`. This algebraic initialisation is not a model evaluation and will not be revised after numerical results are seen.

## Stationary solver

Every distinct WT parameter vector receives at most one bounded stationary solve. The solver is the inherited ten co-ordinate `scipy.optimize.least_squares` solve with the existing co-ordinate bounds, residual rows, residual scaling, tolerances, charge manifold, and rank audit. The maximum actual stationary residual evaluations per solve, including final audit evaluations, is `1600`.

The first solve starts at the accepted Task 31 R09 WT co-ordinates. A later solve starts at the completed numerical stationary state belonging to the previously evaluated parameter vector nearest in Euclidean log multiplier distance, with ties resolved by earliest evaluation order. If no earlier numerical state exists, it uses the Task 31 R09 WT co-ordinates. This rule is deterministic and immutable.

Repeated requests for the identical six multiplier vector use an in-memory cached result and do not start another stationary solve or count as another distinct vector.

## Single local optimisation algorithm

The sole WT optimiser family is SciPy COBYLA applied to the six log multipliers with the declared box bounds and hard constraint margins below. The first run uses `rhobeg=0.35`, `tol=1e-4`, `catol=1e-8`, and a combined ceiling of `19` distinct WT parameter vectors. There is no grid, sweep, pairwise search, factorial search, random design, multistart cloud, global method, evolutionary method, Bayesian method, or change of optimiser family.

One deterministic restart is permitted only if the first COBYLA call terminates before the combined 19 vector ceiling and at least seven unused vector evaluations remain. It starts from that call's own least infeasible evaluated point, ranked by `(sum of squared negative normalised constraint margins, soft objective, evaluation order)`, uses the same bounds, `rhobeg=0.05`, `tol=1e-4`, and `catol=1e-8`, and may use only the remaining part of the same combined 19 vector allowance. It is not permitted after a hard limit or from any independently chosen point.

## Hard stationary constraint margins

A stationary candidate is feasible only when every item below passes. Optimiser margins are dimensionless and nonnegative on the admissible side. Final acceptance uses the unrounded values and the stated inequalities.

1. All inherited full stationary gates pass: maximum scaled independent RHS `<= 1e-7`, maximum omitted charge row `<= 1e-9 fmol/s`, maximum regulatory RHS `<= 1e-10 1/s`, maximum bulk charge `<= 1e-9 fmol`, maximum membrane current residual `<= 1e-18 A`, all conservation ratios `<= 1`, stationary Jacobian rank `10`, no inherited co-ordinate boundary hit, positive core amounts and volumes, time invariant REST RHS, and all finite capacity utilisation ratios `<= 1 + 1e-12`.
2. Resting outflow target `q_out = 0.0010757073853493032 pL/s` with absolute error `<= 1.0757073853493032e-09 pL/s`, equivalently relative error `<= 1e-6`.
3. NKCC1 share of the positive basolateral chloride loading pool lies in `[0.69, 0.71]`. `J_NKCC1,Cl = 2 J_NKCC1,cycle`. The positive pool is the sum of `max(J_NKCC1,Cl,0)`, `max(J_AE4,Cl,0)`, and `max(J_AE2,Cl,0)`; there is no other implemented basolateral chloride loader. Signed counterfluxes are retained separately.
4. `45 <= Cl_i <= 65 mM`.
5. `6.84 <= pH_i <= 6.98`.
6. AE4 chloride loading is strictly positive and its positive pool share is at least `0.05`.
7. `Na_i <= 30 mM`, all reported cell and lumen concentrations are finite and strictly positive, and both compartment volumes are finite and strictly positive.
8. Cell volume is strictly below `3 pL`.
9. The inherited Task 31 potassium and carbon sanity gates remain active: `60 <= K_i <= 210 mM`, `TIC_i < 100 mM`, and `HCO3_i < 100 mM`.

If a stationary solve cannot return finite diagnostics, its soft objective is `1e12` and every hard margin is `-1e6`. Otherwise the normalised margins are the signed distance to each stated threshold divided by that threshold or interval scale. A discrete gate has margin `+1` when true and `-1` when false.

## Soft objective and weights

Hard feasibility is never replaced by a penalty. Among feasible candidates, COBYLA minimises

`F = sum_i d_i^2 + 0.1 sum_j ((x_j - x31_j)/(upper_j-lower_j))^2 + 0.1 ((Cl_i-Cl31_i)/20)^2`,

where `x_j` are the ten inherited stationary physical co-ordinates, `x31_j` are the accepted Task 31 R09 WT co-ordinates, and `Cl31_i = 60.30496692587399 mM`. For a multiplier above one, `d_i = log(m_i)/log(upper_multiplier_i)`; below one, `d_i = log(m_i)/abs(log(lower_multiplier_i))`. Thus every parameter term has weight `1`, every state co-ordinate term has weight `0.1`, and the extra intracellular chloride term has weight `0.1`. No knockout quantity appears.

## Mandatory pre-freeze WT calcium gate

After the local optimisation, stationary feasible evaluated candidates are considered in ascending `(F, evaluation order)`. At most three are integrated, stopping at the first pass. Each uses only the existing central stimulus with `Ca = 0.25 uM` from `t = 0` to `600 s`, the production Radau method and tolerances, and the candidate's own stationary WT state. No separate resting integration is run.

A candidate passes only if all of the following hold:

1. The integration succeeds, reaches `600 s`, matches the REST right hand side at stimulus onset, and satisfies all existing dynamic conservation and finite capacity gates.
2. Every core amount and both compartment volumes remain finite and strictly positive; every evaluated concentration is finite and strictly positive; `q_out` remains finite and nonnegative.
3. The trapezoidal time mean of `q_out` on the fixed one second samples from `60` through `600 s` is at least `1.10 * 0.0010757073853493032 pL/s`.
4. Trapezoidal cumulative secretion on the fixed one second samples from `0` through `600 s` is at least `1.10 * (600 * 0.0010757073853493032) pL`.
5. `q_out(600 s) >= 0.0010757073853493032 pL/s`.

The calcium value, interval, sampling, activation factor and physical gates will not change after results are seen. A failed integration is consumed and will not be repeated with another solver.

## Freeze and held out genotypes

The first dynamic pass in the declared objective ordering is frozen. Before any genotype is evaluated, the complete WT parameter payload, state and flux ledger, dynamic activation summary, and SHA256 payload hash are written and checkpointed through the connected GitHub integration. The accepted pre-freeze WT trajectory is the WT denominator and is not rerun.

Only after this freeze are two genotypes released, in this order: AE4 expression `0.05`, then `0.00`. Every non-AE4 parameter remains frozen. Each receives one stationary solve from the frozen WT co-ordinates. It receives one `Ca = 0.25 uM`, `0` to `600 s` Radau integration only if its own REST passes the same numerical and physical admissibility gates. No intermediate expression, R10, AE2 loss, other calcium value, or refit is permitted.

## Counters and immediate stops

The entire task uses one scientific worker and one BLAS thread. Cumulative ceilings are:

* `20` distinct WT parameter vectors, with the declared optimiser voluntarily capped at `19`;
* `30` stationary solver calls including held out genotype solves;
* `40000` actual stationary residual evaluations;
* `3` pre-freeze WT dynamic integrations;
* `2` post-freeze genotype dynamic integrations, at most one per genotype;
* `1200 s` numerical time for WT optimisation and pre-freeze activation checks.

Counters are checked before an operation that would exceed a ceiling and after every operation. Reaching any Task 33 hard ceiling without an already accepted WT dynamic pass causes immediate termination with `WT_OPTIMISATION_FAILED_WITHIN_PREDECLARED_BOUNDS`. The contract, optimiser, free set, bounds, targets and gates cannot be changed after failure.
