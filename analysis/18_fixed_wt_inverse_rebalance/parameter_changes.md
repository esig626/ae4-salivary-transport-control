# Task 18 WT parameter changes

## Feasibility

All ten inherited WT roots support both requested chloride allocations while preserving the exact same conserved state. The 20 modified WT cases satisfy the complete production steady state equations and electrical closure. The positive loading pool changes only at floating point roundoff, with a largest absolute discrepancy of \(7.3\times10^{-16}\) fmol/s.

The 10% share is realised between 0.0999999999999996 and 0.1000000000000005. The 30% share is realised between 0.2999999999999987 and 0.3000000000000018.

## Departures at 10% AE4 share

The largest log fold change corresponds to a factor between 6.53 and 10.69 across the ten roots. AE4 carrier amount rises by factors 5.35 to 7.38, while NKCC1 capacity falls modestly to 0.912 to 0.917 of baseline.

The coupled AE4 carbon and cation sources require NHE1 to rise by factors 6.53 to 8.12 and basolateral CO2 exchange by factors 6.02 to 7.33. Basolateral pump capacity rises by factors 2.83 to 5.18. The apical pump coordinate ranges from 0.094 to 0.857 of baseline.

Total K conductance ranges from 10.79 to 149.63 nS. CaCC conductance falls from the inherited 31.4 nS reference to 14.23 to 24.24 nS so that its inherited WT chloride flux is preserved under the new voltage. Apical voltage ranges from minus 46.47 to minus 36.61 mV, basolateral voltage from minus 89.68 to minus 73.05 mV, and transepithelial voltage from minus 44.63 to minus 36.24 mV.

## Departures at 30% AE4 share

The largest log fold change corresponds to a factor between 20.62 and 416.05. AE4 carrier amount rises by factors 16.06 to 22.14, while NKCC1 capacity falls to 0.710 to 0.713 of baseline.

NHE1 rises by factors 20.14 to 24.59, basolateral CO2 exchange by factors 18.35 to 21.96, and basolateral pump capacity by factors 7.33 to 14.74. The apical pump coordinate ranges from 0.0024 to 0.622 of baseline.

Total K conductance ranges from 32.07 to 5824.65 nS. CaCC conductance falls to 5.48 to 16.97 nS while retaining the inherited CaCC chloride flux. Apical voltage ranges from minus 74.23 to minus 40.58 mV, basolateral voltage from minus 119.34 to minus 82.13 mV, and transepithelial voltage from minus 45.10 to minus 40.59 mV.

One AE4NA05 root, `N_ABS_NKCC_S4_PLOW_KHIGH_AE4NA05_P1_H1:B00`, requires a 416 fold minimax change and 5824.65 nS total K conductance. Its inherited paracellular HCO3 and K reversal potentials differ by only about 0.085 mV. The modified K balance forces the voltage into this narrow sign compatible interval, so the very large K and paracellular conductances are structural consequences of that root rather than inherited upper bounds or optimizer failure. Four independent searches returned the same minimax value.

## Interpretation before genotype evaluation

The old 14 nS total K equality was the direct blocker in Task 17. Once it is treated as a reference rather than a law, every target is algebraically feasible. The inherited 31.4 nS CaCC value also cannot remain fixed when voltage reclosure is allowed and the same chloride flux must be preserved.

Feasibility is therefore established, but the required changes are not modest. Even the 10% allocation needs several capacities to move by factors around seven, and the 30% allocation needs factors around 20 to 25 in most roots, with one root requiring a factor above 400. This conclusion uses WT information only.
