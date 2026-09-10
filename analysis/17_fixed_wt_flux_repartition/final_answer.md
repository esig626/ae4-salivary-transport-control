# Task 17 completed: fixed-WT flux repartition

Primary classification:

**FIXED-WT FLUX REPARTITION STRUCTURALLY INCOMPATIBLE WITH THE EXISTING NETWORK**

This classification applies to the inherited native source constraints,
including its 14 nS total K maximum. The source equations have consistent
signed algebraic solutions. The obstruction is **capacity/source-bound
incompatibility**, not matrix rank deficiency or local optimizer failure.
It is not a claim that higher AE4 chloride loading is impossible in biology
or in a different source-supported channel model.

| Scientific condition | Fixed-WT feasible | Disposition |
|---|---:|---|
| Inherited baseline | 10/10 roots | Exact inherited parameters and states retained |
| AE4 positive Cl share 0.10 | 0/10 roots | Five negative-required-K cases; five native total-K bound conflicts |
| AE4 positive Cl share 0.30 | 0/10 roots | Ten negative-required-K cases |

No WT state coordinate was solved for, continued, rounded, or changed. All
30 records carry the exact inherited state and its hashes. Rejected records
do not masquerade as new resting states: their allocation-only RHS and
signed diagnostic residuals are distinguished from an accepted complete
physical RHS. Every accepted baseline satisfies the full production RHS at
that exact state. Maximum absolute resting RHS is `8.28e-13` (raw component
units), and maximum scaled independent residual is `1.30e-11`, versus the
inherited `1e-7` scaled tolerance. Charge, omitted alkalinity, current,
carbon, speciation and water checks pass their inherited tolerances.

## Why the higher shares fail at the unchanged WT state

AE4 net loads Cl and NKCC1 is the only other positive basolateral Cl loader
at every inherited root. AE2 is a negative counterflux and is excluded from
the positive-loading denominator. The complete AE4 carrier and complete
NKCC1 cycle were scaled exactly by the declared rule. Positive Cl loading
is conserved to floating-point precision; no Cl-only source adjustment was
made. The imposed AE4 carrier multipliers are 5.35–7.38 at share 0.10 and
16.06–22.14 at share 0.30, with a common proportional reduction of the
other positive loaders.

The fixed-state source matrix has rank 9 and nullity 1 at every root,
including the 14 nS total-K constraint and both current equations. A
smallest full-rank basis is recorded; feasibility retains the complete
remaining affine direction. Independent HiGHS dual-simplex and
interior-point routes both certify infeasibility at every modified endpoint.
Signed equality solutions have maximum residual `5.69e-14`, demonstrating
that the obstruction is physical bounds rather than unsolved equations.

An independent balance certificate makes the result explicit. Alkalinity
fixes the required NHE1 change. Na then fixes the total pump change, which
leaves a unique required **total outward K-channel flux**, irrespective of
pump localization or paracellular capacity. Fixed chloride transport fixes
both membrane potentials, so the two inherited K driving forces remain
outward. With total K conductance fixed at 14 nS, the total outward K rate
must lie between its all-apical and all-basolateral values.

| AE4 routing | Target | Required total outward K, fmol/s | Result |
|---|---:|---:|---|
| AE4NA05 | 0.10 | 0.06225–0.06709 | Below every root's minimum at 14 nS; required apical K fractions 1.20–5.67 exceed 1 |
| AE4NA05 | 0.30 | −0.20492 to −0.20383 | Requires inward K at outward driving forces |
| AE4NA20 | 0.10 | −0.24433 to −0.23813 | Requires inward K at outward driving forces |
| AE4NA20 | 0.30 | −1.19937 to −1.16343 | Requires inward K at outward driving forces |

The smallest source-bound violation is `0.02659 fmol/s`, far above inherited
numerical tolerances. These certificates allow all pump/K fractions in
`[0,1]`, arbitrary nonnegative existing pump and non-Cl paracellular
capacities, and relaxed NHE1 upper bounds. Hence neither narrow fraction
brackets nor an arbitrary upper-capacity cut-off causes rejection.

The five AE4NA05/share-0.10 cases are specifically limited by the retained
native total-K source constraint. They should not be described as failures
of the unconstrained stoichiometric network. Changing that sourced total
would define a different model assumption. The other fifteen cases require
negative total outward K even if nonnegative K capacities were otherwise
unrestricted. Neither conclusion is a proof that mixed-cation slip causes
the secretion phenotype.

Uniform complete-carrier scaling leaves the unchanged-state AE4 branch
cancellation fraction unchanged: 96.80–99.01%. It increases the large
opposing Na/K branch fluxes along with productive Cl loading; the other WT
balances cannot be presumed unchanged just because Cl loading is preserved.

## Pushed freeze and post-freeze genotype accounting

The complete WT-only checkpoint was committed, pushed and read back from
the requested remote branch before genotype payload access:

`0fe0d9a09020f7311c63a665b3f8a8adbed68ce6`

It includes every baseline/0.10/0.30 decision, all exact state payloads,
accepted full parameters, rank/source ledgers, certificates, methods and
focused pre-freeze tests. The remote verification receipt is
`results/17_fixed_wt_flux_repartition/pushed_fixed_wt_checkpoint.json`.

No modified endpoint qualified for new AE4 5% continuation, AE2 loss, or
matched WT/genotype trajectories. All 60 modified root/share/calcium
combinations are explicitly marked `NOT_EVALUATED_FIXED_WT_INFEASIBLE` in
each result table, without fabricated ratios. No exact-zero AE4 attempt
was made and no baseline ODE trajectory was recomputed.

After the pushed freeze, all 30 inherited root/calcium baseline pairs were
reused with their parameter, trajectory and denominator hashes checked:

| Baseline genotype comparison | Total secretion change from matched WT |
|---|---:|
| AE4 expression 0.05 | +8.23% to +21.66%; median +17.82% |
| AE2 exact loss | +0.0318% to +0.1601%; median +0.0868% |

Baseline AE4 loss therefore retains the inherited increased-secretion
direction while AE2 loss has a much smaller effect. There is no eligible
higher-share genotype result from which to infer a repaired direction or
new specificity. These outcomes did not construct, select or rank any WT
parameterization. Genotype result files were subsequently hashed and frozen.

## Verification and retained artifacts

`PYTHONPATH=src python -m unittest discover -s tests -p 'test_task17*.py' -v`
passes **11 tests, zero failures, zero skips**. The checks include exact
state bytes, all 187 inherited frozen inputs, full coupled transport
scaling, the fixed positive-loading pool, every accepted resting equation,
the complete affine null direction against the unchanged production model,
independent infeasibility certificates, pushed-freeze ordering, and complete
post-freeze accounting. No ODE or WT state solve occurs in these tests.

All 187 frozen inputs, the production scientific equations/parameters, and
`archive/` remain unchanged. Task-specific algebra, tests and reports are
added outside the archive. Work is confined to
`codex/task-17-fixed-wt-flux-repartition`; no merge to `main` was performed.

The required artifacts are in `results/17_fixed_wt_flux_repartition/` and
this analysis directory. Numerical source units, the eligibility audit and
reproduction steps are detailed in `method.md`. The post-freeze reuse
command is `PYTHONPATH=src python -m modern_full_model.task17_postfreeze`;
it refuses to overwrite frozen results or evaluate an unpushed checkpoint.
