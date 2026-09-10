# Task 18 genotype interpretation

## Firewall and numerical completion

Genotype evaluation began only after the complete WT solution set was committed and pushed at `2cdcfb3f3ab606b2c6cb84f3149f840bee1f72e0`. The checkpoint contains all ten roots and all three WT chloride allocation conditions. Its manifest confirms that no AE4 loss or AE2 loss result was read or used during WT construction or canonical selection.

The post freeze evaluation made no parameter refit. All 20 modified parameterisations completed connected AE4 expression continuation from 1.0 to 0.05 and AE2 expression continuation from 1.0 to 0.0. AE4 exact zero was never attempted. The old WT calibration coordinate boxes were not used as genotype feasibility walls. One valid AE4 0.05 state has a cell volume of 3.878 pL and would have been incorrectly rejected by the inherited 3 pL search ceiling. All terminal roots are full rank, lie strictly inside the Task 18 numerical enclosures, and have maximum scaled independent residuals of (2.36\times10^{-11}) for AE4 and (1.55\times10^{-11}) for AE2.

There are 90 rows in each dynamic result table: ten roots, three chloride allocation conditions, and three calcium values. All WT, AE4 5 percent, and AE2 trajectories complete the 600 second production Radau protocol and pass the inherited numerical and conservation gates.

## Dynamic results

| WT chloride allocation | Realised WT AE4 share | WT secretion, pL over 600 s | AE4 5 percent ratio | AE4 direction | AE2 loss ratio |
| --- | ---: | ---: | ---: | --- | ---: |
| Inherited baseline | 0.01355 to 0.01868 | 0.7144 to 1.0790 | 1.0823 to 1.2166 | Increase in 30 of 30 | 1.00032 to 1.00160 |
| 0.10 | 0.100000000000000 within (5.3\times10^{-16}) | 0.6887 to 1.0237 | 1.1862 to 1.4899 | Increase in 30 of 30 | 1.00026 to 1.00100 |
| 0.30 | 0.300000000000000 within (1.8\times10^{-15}) | 0.6689 to 1.0075 | 0.9962 to 1.5717 | Increase in 29 of 30 | 1.00027 to 1.00100 |

Relative to its corresponding inherited WT trajectory, modified WT secretion is 0.938 to 0.971 at share 0.10 and 0.865 to 0.962 at share 0.30. This dynamic change is reported rather than used to alter the frozen parameterisations.

At share 0.10, AE4 near loss increases secretion for every root, routing family, and calcium value. At share 0.30, 29 of 30 comparisons increase. The sole exception is `N_ABS_NKCC_S4_PLOW_KHIGH_AE4NA20_P2_H1:B00` at calcium 0.10 uM, where the ratio is 0.99616, a reduction of only 0.384 percent. The same parameterisation increases secretion at calcium 0.25 and 0.50 uM. This exception is neither robust nor close to the contextual experimental reduction of about 35 percent.

Across the modified models, the AE4NA05 routing family has ratios from 1.193 to 1.462 at share 0.10 and 1.297 to 1.572 at share 0.30. AE4NA20 has ratios from 1.186 to 1.490 at share 0.10 and 0.996 to 1.417 at share 0.30. The general wrong direction therefore spans both routing families and all calcium values, although it is not mathematically universal because of the single slight reduction.

AE2 remains comparatively neutral. Every ratio is just above one, and the complete range across baseline and modified conditions is 1.00026 to 1.00160.

## Mechanistic interpretation

The complete AE4 transporter retains a cancellation fraction of 0.9680 to 0.9901 and a productive fraction of only 0.0099 to 0.0320. These fractions are unchanged by the WT carrier scaling because routing and thermodynamic ratios are fixed. Increasing the WT share therefore raises the magnitude of two strongly cancelling AE4 branches without repairing the mixed cation mechanism.

The parameter rebalancing establishes fixed WT feasibility, but it does not recover a meaningful secretion reduction after AE4 near loss. Since most modified comparisons move farther in the wrong direction and the only decrease is negligible and calcium specific, the current AE4 mixed cation architecture remains the likely next mechanism problem.

## Primary classification

`FIXED-WT REBALANCING REQUIRES EXTREME PARAMETER SHIFTS`

This classification is primary because even the 0.10 allocation needs minimax capacity changes of 6.53 to 10.69 fold, while the 0.30 allocation needs 20.62 to 416.05 fold. The genotype result is an additional mechanistic finding and was not used to select any WT parameterisation.
