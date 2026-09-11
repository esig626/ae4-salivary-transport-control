# Task 17 fixed-WT method

All ten WT amount/volume vectors are copied verbatim from the hash-valid
Task 14B scientific inheritance. Basal regulatory coordinates are appended
exactly as in Task 14. No WT root search, coordinate continuation, rounding,
or reconstruction through electroneutral coordinates is performed. Every
accepted parameterization must satisfy the complete production resting RHS,
including omitted alkalinity rows, conserved sources, current, charge and
water accounting at that exact state.

The three conditions are inherited baseline, AE4 positive Cl share 0.10,
and share 0.30. They are sensitivity conditions, not measured fractions.
At all inherited roots NKCC1 is the only other positive basolateral Cl
loader; AE2 is a negative counterflux. Thus `Lplus=A0+2*J_NKCC`,
`lambda=f*Lplus/A0`, and `rho=(1-f)*Lplus/(2*J_NKCC)`. Carrier amount scales
the complete AE4 Na/K/Cl/TIC/alkalinity vector, with unchanged rates,
occupancies, affinities and routing. NKCC1 scales as its full Na:K:2Cl cycle.
No chloride source is scaled independently of its coupled ions.

## Source eligibility and existing bounds

Task 17 reopens existing unmeasured WT capacities. It does not turn the
native source measurements into arbitrary nuisance parameters. The current
source records are `src/modern_full_model/native_source_panel.py`,
`analysis/13B_modern_full_model/cation_topology.md`, and the native whole-cell
channel section of `analysis/13B_modern_full_model/evidence_ledger.md`.
They explicitly supersede the early common-conductance calibration and fix
the CaCC maximum at 31.4 nS and total K maximum at 14 nS. Background K
conductances are zero. There is no declared uncertainty interval licensing
a change of the 14 nS total. This is an inherited source constraint, not a
new bound invented in Task 17. Its transfer from parotid evidence to this
SMG model remains a model assumption, not a claim of a universal biological
constant.

The maximal eligible balancing set includes NHE1, effective pump capacities
on both membranes, K partition at fixed total, neutral CO2 exchange on both
membranes, and existing non-Cl paracellular capacities. Independent membrane
capacities map back to the existing total-capacity and fraction parameters;
no transporter, source, background current or topology is added. For the
impossibility proof the narrower inherited fraction brackets and capacity
upper bounds are relaxed to all nonnegative capacities and all fractions
in [0,1]. Failure even in that superset proves failure within the original
bounds. The old AE4 carrier fit box and old common-conductance box are not
used to reject a target; they belong to superseded calibrations. NHE1's old
WT calibration bound would additionally be checked for any physical solution.

Changing negative AE2 would change cell chloride without being part of the
prescribed positive-loader reduction; offsetting it with CaCC would violate
the chloride constraint. CaCC, paracellular Cl, bath, geometry, chemistry,
hydraulics, outflow, calcium and regulation therefore stay fixed. At fixed
concentrations and fixed CaCC conductance, fixed cell Cl balance fixes the
apical potential. The fixed lumen Cl balance and fixed paracellular Cl
conductance then fix the transepithelial potential, hence the basolateral
potential. This is a deduction from the specified balances, not a voltage
fit or extra state search.

## Exact source matrix and robust feasibility

The baseline ledger derives individual source vectors from declared reaction
stoichiometry and current orientation and checks their sum against the
complete production RHS. AE4 total is a subtotal and is not double counted.
Water and outflow appear separately. Electrical rows use fmol charge/s;
multiply by Faraday/1e15 to recover amperes.

Before any target solve, the fixed-state matrix is formed for ten effective
capacity multipliers: NHE1; apical and basolateral pumps; apical and
basolateral K; basolateral and apical CO2; paracellular Na, K and HCO3.
The first 12 rows are core source rates, the next two are electrical closure,
and the last is `fK*mKa+(1-fK)*mKb=1`, the native 14 nS constraint.
The matrix has rank 9/nullity 1 at all roots (rank 8 without the total-K
constraint). The lexicographically first nine-column full-rank basis is
recorded. The full one-dimensional affine family is retained in feasibility
tests; setting the omitted basis coordinate arbitrarily cannot hide a
feasible solution. If a feasible family existed, the predeclared choice
would minimize squared distance of these multipliers from one. Baseline's
minimum-change representative is exactly the inherited parameter vector.

Writing `d` for the coupled AE4/other-loader source change, solve
`M*(m-1)=-d` with `m>=0` and the source-total equality. This is affine at
the inherited voltages. HiGHS dual-simplex and interior-point routes test
feasibility. Signed minimum-change solutions are retained only as algebraic
diagnostics, with their full equality residuals; they are not accepted
negative-capacity models. No local nonlinear optimizer determines rejection.

Independently, cell alkalinity fixes `dNHE=-dTA`. Cell Na then fixes total
pump change `dP=(dNa+dNHE)/3`. Cell K consequently requires

`K_required=K_baseline+dK+2*dP`.

With fixed total K conductance and fixed voltages, total outward K is a
convex combination of the all-apical and all-basolateral values. Therefore
it must lie in their closed interval for every partition in [0,1]. A
required value below the lower endpoint is an explicit source-bound
incompatibility certificate independent of pump and paracellular partition.
A negative required value is stronger: both channels' inherited driving
forces are outward, so even arbitrary nonnegative K capacities cannot
supply it. This distinguishes bound incompatibility from rank deficiency
and numerical failure.

## Freeze and genotype discipline

Save all 30 decisions, allocation-only residuals, signed diagnostics,
certificates and exact inherited states. Only physically feasible solutions
receive accepted full parameter payloads. All source files, tests, methods
and WT artifacts are hashed into the WT manifest, committed and pushed
before any genotype result is accessed. A rejection has no fabricated
complete physical RHS: its allocation-only RHS and signed affine residual
are explicitly separate. Baselines have the full production RHS.

The current algebra certifies no modified endpoint as feasible. Thus the
post-freeze action is to reuse hash-valid inherited baseline AE4 5% and AE2
loss results, and mark every modified endpoint not evaluated due to its
pre-frozen WT infeasibility. No new WT dynamics, genotype root continuation,
or expensive baseline simulations are justified. No genotype output can
alter these decisions.

## Reproduction

From the repository root, with `PYTHONPATH=src`:

1. `python -m modern_full_model.task17_fixed_wt prepare`
2. `python -m modern_full_model.task17_fixed_wt solve`
3. `python -m unittest discover -s tests -p test_task17_fixed_wt.py -v`
4. `python -m modern_full_model.task17_fixed_wt freeze`
5. Commit and push the fixed-WT checkpoint; record and verify the remote SHA.
6. Run the post-freeze report command after that remote verification.

Preparation/freeze refuse to overwrite their immutable artifacts. Use a
separate clean checkout/output directory for a full rerun; the tests can
verify the saved artifacts directly without any state solve.
