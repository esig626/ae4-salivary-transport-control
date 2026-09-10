# Task 18 final answer

## Primary classification

`FIXED-WT REBALANCING REQUIRES EXTREME PARAMETER SHIFTS`

## Answer

Yes. Both the 0.10 and 0.30 AE4 positive chloride loading shares can support every one of the ten exact inherited WT conserved states when the uncertain existing capacities are genuinely allowed to rebalance. All 30 WT cases retain the complete inherited state byte for byte. The largest complete resting RHS is (8.28\times10^{-13}), the largest membrane current residual is (1.93\times10^{-24}) A, the largest target share error is (1.84\times10^{-15}), and the largest positive chloride pool error is (7.22\times10^{-16}) fmol/s.

AE4 and NKCC1 are the only positive basolateral chloride loaders in every inherited signed ledger. AE4 is scaled through its complete coupled Na, K, Cl, HCO3, TIC, and alkalinity source. NKCC1 carries the proportional remainder. AE2 stays a negative counterflux and is not included in the positive pool. The accepted models also preserve the inherited AE2, CaCC, and paracellular chloride fluxes to numerical precision. Complete signed ledgers are stored in every frozen WT payload.

| Quantity | AE4 share 0.10 | AE4 share 0.30 |
| --- | ---: | ---: |
| Realised share | 0.100000000000000 within (5.3\times10^{-16}) | 0.300000000000000 within (1.8\times10^{-15}) |
| Largest capacity departure | 6.53 to 10.69 fold | 20.62 to 416.05 fold |
| AE4 carrier multiplier | 5.35 to 7.38 | 16.06 to 22.14 |
| NKCC1 multiplier | 0.912 to 0.917 | 0.710 to 0.713 |
| NHE1 multiplier | 6.53 to 8.12 | 20.14 to 24.59 |
| Basolateral CO2 multiplier | 6.02 to 7.33 | 18.35 to 21.96 |
| Basolateral pump multiplier | 2.83 to 5.18 | 7.33 to 14.74 |
| Total K conductance | 10.79 to 149.63 nS | 32.07 to 5824.65 nS |
| CaCC conductance | 14.23 to 24.24 nS | 5.48 to 16.97 nS |

The old 14 nS total K conductance was the direct Task 17 blocker. The inherited 31.4 nS CaCC conductance also cannot remain fixed while voltage recloses and the inherited CaCC chloride flux is preserved. These are reference values, not physical feasibility constraints. Once released, all requested WT allocations are feasible, but the required rebalancing is clearly extreme rather than modest. The 0.30 outlier needs a 416.05 fold minimax change and 5824.65 nS total K conductance because positivity confines its voltage close to neighbouring paracellular reversal potentials.

All modified WT models completed production Radau trajectories at calcium 0.10, 0.25, and 0.50 uM. Their 600 second secretion is 0.938 to 0.971 of inherited WT at share 0.10 and 0.865 to 0.962 at share 0.30.

After the WT checkpoint was pushed, every modified model was continued to AE4 expression 0.05 with all capacities frozen. Across all 90 baseline and modified comparisons, AE4 near loss increases secretion in 89 cases. Ratios are 1.082 to 1.217 at baseline, 1.186 to 1.490 at share 0.10, and 0.996 to 1.572 at share 0.30. The only decrease is 0.384 percent in one AE4NA20 root at calcium 0.10 uM; that root increases at the other calcium values. The approximately 35 percent experimental reduction is therefore not recovered.

AE2 loss remains comparatively neutral, with ratios from 1.00026 to 1.00160 across all conditions. Every one of the 270 trajectories passes its numerical gate. Feasibility is robust across roots, routing families, and calcium values. The AE4 phenotype direction is nearly universal but has the one negligible exception just described.

The AE4 productive fraction remains only 0.0099 to 0.0320, with 0.9680 to 0.9901 branch cancellation. Reallocating more WT chloride load to the unchanged complete transporter does not repair that cancellation. The current AE4 mixed cation architecture therefore remains the likely next mechanism problem.

The complete WT set was frozen and pushed at `2cdcfb3f3ab606b2c6cb84f3149f840bee1f72e0` before any genotype evaluation. No genotype result entered WT feasibility, optimisation, or selection, and no capacity was refitted after the checkpoint.
