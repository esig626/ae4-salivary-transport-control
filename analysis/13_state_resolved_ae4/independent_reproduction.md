# Independent numerical reproduction

## Outcome

The independent route reproduces the **unresolved mechanistic bottleneck**, not
a successful secretion phenotype. It confirms all of the numerical facts on
which that conclusion depends:

1. the frozen state-resolved CTMC occupancies and currents;
2. the bounded, candidate-independent fixed-chassis AE4-null root;
3. the wrong chloride and pH of that root after the Stage-A reveal;
4. the exact module-signature ranks and the conditional sparse projections;
5. the five-dimensional released target manifold and rank-five variation of
   its required source; and
6. the exact nesting, nonzero knockout root, and failed WT boundary attempts
   of the built Round-5 pump/K membrane split.

No admissible WT Gate-2 row exists, so there is no valid frozen WT denominator
or validated extension from which to reproduce a knockout/WT secretion
trajectory. Computing such a ratio would promote a failed chassis row and
violate the holdout discipline.

## Independence firewall

The reproducer is `src/state_resolved_ae4/independent.py`. It imports neither
Agent E's cycle/solver modules nor Agent F's localization module.

- CTMC stationarity is solved from the smallest right singular vector of
  `Q.T` and, separately, from principal generator cofactors using the
  matrix-tree theorem. Agent E used a row-replacement linear solve.
- The seven-state AE4-null equations are independently transcribed from the
  documented chassis equations. The two-voltage closure uses an explicit
  two-by-two inverse rather than the production linear solve.
- The nonlinear equilibrium uses `scipy.root` in positive coordinates that
  replace intracellular bicarbonate by the eliminated positive proton
  coordinate.
- Integer ranks use exact `Fraction` Gaussian elimination. Sparse projections
  use NumPy SVD least squares and the Lawson-Hanson NNLS algorithm.
- No AE4-null secretion magnitude or time-course value is read by a numerical
  calculation. All 17 frozen CTMC condition rows declare
  `heldout_ko_data_read=False`.

## State-resolved CTMC crosscheck

Agent E froze 17 conditions covering SR2, SR5, and SR6, with 288 state rows and
365 reversible-edge rows. The three input hashes match the manifest exactly.

| Comparison across all 17 conditions | Maximum absolute discrepancy |
| --- | ---: |
| SVD occupancy versus Agent E | `5.71e-12` |
| Matrix-tree occupancy versus Agent E | `2.90e-13` |
| SVD edge current versus Agent E | `1.32e-12` |
| SVD versus matrix-tree occupancy | `5.71e-12` |
| SVD stationarity residual, infinity norm | `1.26e-12` |
| Matrix-tree stationarity residual, infinity norm | `1.94e-11` |

Every condition passes the declared `1e-10` occupancy and `1e-9` current
tolerances. The slightly larger matrix-tree residual is absolute cancellation
against edge rates spanning many orders of magnitude; its probabilities agree
more closely with Agent E than the SVD route does.

The direct log-activity calculation also reproduces the frozen salivary-tuple
affinities: SR1-Na `-1.89404638`, SR1-K `+4.10064824`, and SR3 `+6.64341795`.
This confirms that SR3 is thermodynamically favorable under that tuple; it does
not remove the independent direct Na-direction conflict.

## Fixed-chassis AE4-null root

Starting from a deliberately perturbed positive state, the independent
equations and nonlinear solver return:

| Quantity | Agent E frozen row | Independent route | Absolute difference |
| --- | ---: | ---: | ---: |
| `Cl_i` (mM) | `48.7822150270` | `48.7822150255` | `1.55e-9` |
| pH | `7.37847391079` | `7.37847391079` | `<1e-14` |
| volume (pL) | `2.97417820581` | `2.97417820506` | `7.56e-10` |
| maximum state component | saved seven-vector | independent seven-vector | `1.67e-8` |

The independent raw residual is `2.38e-13`; apical and basolateral current
closure residuals are `-2.30e-17` and `+5.07e-17`. The small component
difference is solver tolerance on height, not a different equilibrium.

The bounded root remains `+7.676` reported SEM from the released knockout
chloride target and `+24.424` SEM from the released knockout pH target. The
independent residual function contains no AE4-family argument. Candidate
independence after deletion is therefore structural: every adapter's
genotype-zero AE4 source vanishes before the common non-AE4 field is evaluated.

### Root-domain caveat

A deterministic 30-start log-perturbation check (`seed=1313`) produced 19
convergences within the declared `2-200 um` height domain; all belong to one
cluster with maximum component spread `1.37e-7`. Ten starts failed or entered
nonphysical coordinates, and one converged to a remote huge-height branch.

Agent G found an exact positive root near height `1.25e14 um`. Starting from
that rounded row, the independent solver converged to the same ionic asymptote
(`Cl_i=55.31009 mM`, pH `7.568636`) but height `7.69e14 um`, with residual
`7.56e-16`. The large height drift is evidence of an ill-conditioned or
asymptotic remote branch. It is far outside the declared physiological domain
and cannot rescue the phenotype, but it forbids any claim of global root
uniqueness.

## Stage-B source geometry

The released means fix only `Cl_i=36.50 mM` and pH `6.89`. With cellular
electroneutrality, the seven-state target remains a five-dimensional positive
manifold parameterized by `(Na_l,K_l,height,Na_i,K_i)`. The independent
central-difference calculation reproduces the saved `12 x 5` conditional
source Jacobian to maximum absolute difference `7.41e-13`.

Its singular values are

`(0.4022354333, 0.0188763562, 0.00934813971, 0.00163359948, 0.000449923228)`.

All five are nonzero, so the required source changes in rank five across
states that are equally consistent with the released chloride and pH. This is
not an uncertainty in one scalar amplitude along a fixed vector.

Exact rational ranks are independently reproduced:

| Signature projection | Rank |
| --- | ---: |
| All 14 expanded module signatures | `12` |
| Six intracellular chemical rows | `6` |
| Twelve determined reduced rows | `11` |
| Pump/K membrane-split differences, expanded space | `2` |
| Pump/K membrane-split differences, intracellular chemistry | `0` |

The sparse projection coefficients and residuals in `module_projection.csv`
are reproduced to `2e-9` or better. In particular:

| Conditional slice | Rule | Best displayed pair | Residual fraction |
| --- | --- | --- | ---: |
| Model-anchored | signed | NKCC1 + carbon hydration | `0.120909` |
| Model-anchored | nonnegative | apical Cl + carbon hydration | `0.261994` |
| WT-reference cations | signed | NKCC1 + NHE1 | `0.297902` |
| WT-reference cations | nonnegative | apical Cl + carbon hydration | `0.656699` |

No displayed one- or two-module cell signature is exact, and the best pair
changes with the unmeasured state coordinates. Negative signed coefficients
also represent reversal or reduced activity, not addition of the named
physiological module.

## Minimal extension identifiability

### Independent Round-5 chassis reconstruction check

The split-cation evaluator was reimplemented locally from the published
equations, including an explicit two-by-two voltage inverse. The frozen
`SR2_SHARED_123/slow` WT source was also reconstructed as an eight-state CTMC
and solved by SVD; neither the production cycle adapter nor its optimizer was
called.

- At `f_P=f_K=0`, the split and fixed RHS agree to `2.11e-19`, and the M0
  root agrees with the Stage-A null state to `4.26e-12`.
- Starting at M0, an independent positive-coordinate `scipy.root` solve at
  `f_P=0.4,f_K=0.2` closes to `4.30e-14`. Its state agrees with the frozen M1
  row to `2.58e-9`; pH differs by `2.88e-10`, and both current balances close
  below `1.1e-17`.
- All seven saved WT equation values are independently reproduced for each of
  the five attempts. The largest component discrepancy is `4.98e-17`.

| WT attempt | Independent max raw RHS | Height (um) | Closure / boundary / WT gate |
| --- | ---: | ---: | --- |
| M0 | `3.05048e-5` | `199.999995` | fail / hit / fail |
| pump 20% | `3.86550e-5` | `199.999989` | fail / hit / fail |
| K 20% | `3.14718e-5` | `199.999987` | fail / hit / fail |
| coupled 20% | `4.37303e-5` | `199.999995` | fail / hit / fail |
| coupled 40% | `5.40240e-5` | `199.999984` | fail / hit / fail |

Each exceeds the `1e-8` raw-closure threshold and hits the declared upper
height boundary. Thus all five independent `root_success`,
`calibration_valid`, and `wt_gate2_pass` decisions are false. A separate
square-root route starting at these rows found no admissible solution: it
escaped toward height at least `7.65e11 um` and capacity tending to zero. That
behavior corroborates the boundary failure; it is not a global nonexistence
proof.

Finally, direct reevaluation of all eight M0/M1 topologies gives exactly zero
spread in the bicarbonate row on either named target slice. The invariant
values are `-0.006545010622425635` (model-anchored) and
`-0.0018064818538654252` (WT-reference cations). Since both are nonzero, cation
redistribution alone cannot make either full conditional slice stationary.

### What the extension does and does not identify

The coupled apical/basolateral Na/K-pump and K-conductance split remains the
smallest independently supported omitted topology to test. Its two difference
directions have exact rank two in lumen/current space but rank zero in
intracellular chemistry. Thus the released chloride/pH pair cannot estimate
the two split fractions or prove that this topology is the repair.

This is the irreducible missing information: the stimulus-dependent,
membrane-partitioned cation-current support tuple

`(P_a, P_b, J_K,a, J_K,b)`,

together with native knockout Na, K, and volume coordinates. The first decisive
cation-topology branch test is to measure time-resolved intracellular Na, K,
Cl, pH, cell volume, and luminal or effluent `Na_l` and `K_l` in WT and
`Ae4-/-` mouse submandibular acini under `0.3 uM CCh + 5 uM IPR`, while
separately estimating apical and basolateral ouabain-sensitive pump current and
Ca-activated K current. This is not a complete mechanism identifier: without
co-measured total inorganic carbon plus NHE and buffer boundary fluxes, normal
cation support would redirect the reconstruction to an acid-base branch but
could not select or parameterize that branch.

## Disposition

The independent calculations contain no unexplained numerical discrepancy.
The nonphysiological remote/asymptotic root behavior narrows all root claims
from global uniqueness to bounded physiological branches; it does not rescue
either the fixed null phenotype or the failed WT boundary attempts.

The currently frozen evidence supports:

`PHENOTYPE UNRESOLVED — SPECIFIC MISSING EXPERIMENT OR FLUX IDENTIFIED`

The machine-readable audit is in
`results/13_state_resolved_ae4/independent_reproduction.csv` and
`independent_reproduction.json`. Ten independent tests cover CTMC routes,
affinities, detailed-balance current signs, the nonlinear null root, exact
ranks, sparse projection, the rank-five Stage-B source family, exact M0
nesting, a nonzero M1 root, and the independently reconstructed WT SR2 source.
