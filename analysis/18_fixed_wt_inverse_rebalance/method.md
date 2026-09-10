# Task 18 method

## Scientific contract

Task 18 asks whether the ten accepted WT conserved states can remain exactly unchanged when AE4 carries either 10% or 30% of the positive basolateral chloride load. The inherited baseline is retained as the third and only other condition. No state coordinate is solved, continued, rounded or fitted.

The signed production ledger identifies AE4 and NKCC1 as the only positive basolateral chloride loaders at every root. AE2 is a small negative counterflux and is excluded from the denominator. If the inherited positive pool is \(L_+\), the complete AE4 carrier is scaled to carry \(fL_+\), while the complete NKCC1 law is scaled to carry \((1-f)L_+\). This preserves the positive pool exactly and retains every coupled Na, K, Cl, carbon and alkalinity source.

No AE4 loss or AE2 loss result is read by the WT runner. The input ledger contains only the frozen Task 13B WT model, its WT validation artifacts, and the ten inherited Task 17 baseline state payloads.

## Fixed state source audit

At each inherited root, the source audit retains twelve independent positive capacity coordinates:

1. NHE1
2. apical and basolateral Na/K pump capacities
3. apical and basolateral K conductances
4. CaCC conductance
5. basolateral and apical neutral CO2 exchange capacities
6. paracellular Na, K, Cl and HCO3 conductances

AE4 and NKCC1 are additional positive coordinates whose values are fixed by the chloride allocation. Pump and K quantities are represented as independent membrane capacities, rather than totals and arbitrarily penalised fractions.

The augmented source matrix also retains apical and basolateral voltages and both current closure equations. It has 14 variables, rank 10 and nullity 4 at every root. The four dependent equation combinations arise from the two volume rows and conserved charge relations. The corresponding production sensitivity, with voltages eliminated through the exact two by two electrical closure, has rank 8 and nullity 4.

## Exact inverse reduction

Once the complete AE4 and NKCC1 scales are imposed, the fixed material balances determine NHE1 flux, total pump flux, total K channel flux, CaCC flux, both CO2 fluxes, and the paracellular Cl and HCO3 fluxes. Four freedoms remain:

1. the apical share of pump flux
2. the apical K flux
3. the apical voltage
4. the transepithelial voltage

For any choice of these four quantities, every required existing capacity follows directly from its production flux law. Only solutions with finite positive capacities are admitted. The resulting parameters are then placed back into the unchanged production model, which independently resolves both membrane voltages and evaluates the complete resting equations.

This reduction preserves AE2, CaCC, paracellular Cl and apical CO2 fluxes at their inherited WT values. The remaining capacity changes restore Na, K, carbon, alkalinity, lumen and current balance.

## Canonical solution

The selection rule was declared before any genotype calculation. For the fourteen positive coordinates comprising AE4, NKCC1 and the twelve rebalancing capacities, the solver first minimises the largest absolute natural log fold change from the inherited value. Subject to that minimum, it minimises the sum of squared natural log fold changes.

Each modified case uses deterministic broad random feasibility sampling, differential evolution, and independent SLSQP epigraph refinements. A second set of SLSQP refinements minimises the squared objective on the minimax face. The voltage and K flux search boxes are numerical exploration ranges only. Every accepted optimum lies far from those artificial limits. Several solutions approach genuine sign boundaries set by the paracellular reversal potentials, and these boundaries are retained as physical positivity constraints.

## Verification

Every accepted payload is checked by the complete production evaluator. The checks include:

* the exact Task 17 complete state byte hash
* the complete resting RHS
* both membrane current closures
* charge, carbon, alkalinity, water and lumen accounting
* complete AE4 coupled source scaling
* the preserved positive chloride loading pool
* the realised AE4 share
* the signed chloride ledger
* all capacity multipliers and algebraic voltages

All ten roots, both AE4 routing families and all three conditions are retained in the frozen WT manifest.
