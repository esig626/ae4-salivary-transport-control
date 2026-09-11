# Five percent AE4 secretion check

Reducing AE4 activity to 5% **increases total secretion in all 30 cases**, by **8.23% to 21.66%** (median **17.82%**). The 5% AE4 to WT total secretion ratio ranges from **1.082287 to 1.216630** (median **1.178191**).

All **30/30** trajectories are numerically valid. The direction is consistent across all ten roots and all three calcium values. The requested percentage reduction, defined as `100 × (1 − ratio)`, is **−21.66% to −8.23%** (median **−17.82%**); its negative sign denotes increased secretion.

| Calcium (uM) | Ratio range | Median ratio | Reduction range | Median reduction |
| ---: | ---: | ---: | ---: | ---: |
| 0.10 | 1.082287 to 1.163180 | 1.132727 | -16.32% to -8.23% | -13.27% |
| 0.25 | 1.145429 to 1.202714 | 1.185277 | -20.27% to -14.54% | -18.53% |
| 0.50 | 1.166968 to 1.216630 | 1.200315 | -21.66% to -16.70% | -20.03% |

The saved Task 14 AE2 percentage reductions range from −0.160080% to −0.031795% (median −0.086836%). These are small increases of 0.0318% to 0.1601%. The matched AE4 effect magnitudes are **112.0 to 449.3 times** the AE2 effects (median **196.0 times**). AE4 has a much larger secretion effect, but it has the opposite direction to the experimental reduction. The frozen 5% functional loss panel therefore does **not** reproduce the general phenotype of reduced secretion after AE4 loss.

Both routing families predict an increase at every calcium value. There is no disagreement in direction. AE4NA05 has a larger median increase than AE4NA20 by **4.70 percentage points at 0.10 uM**, **1.99 at 0.25 uM** and **0.95 at 0.50 uM**. The family ranges are separated at 0.10 uM and overlap at the higher calcium values. Thus routing affects magnitude, especially at the lowest calcium, while both families retain the same qualitative disagreement with the experimental reduction. Medians describe the five retained roots in each family and are not uncertainty intervals.

| Routing family | Calcium (uM) | Ratio range | Median ratio | Reduction range (%) | Median reduction (%) |
| --- | ---: | ---: | ---: | ---: | ---: |
| AE4NA05 | 0.10 | 1.140521 to 1.163180 | 1.151855 | -16.32 to -14.05 | -15.19 |
| AE4NA20 | 0.10 | 1.082287 to 1.124933 | 1.104898 | -12.49 to -8.23 | -10.49 |
| AE4NA05 | 0.25 | 1.179915 to 1.202714 | 1.191545 | -20.27 to -17.99 | -19.15 |
| AE4NA20 | 0.25 | 1.145429 to 1.195148 | 1.171664 | -19.51 to -14.54 | -17.17 |
| AE4NA05 | 0.50 | 1.191979 to 1.213305 | 1.202563 | -21.33 to -19.20 | -20.26 |
| AE4NA20 | 0.50 | 1.166968 to 1.216630 | 1.193040 | -21.66 to -16.70 | -19.30 |

These are model predictions at 5% residual AE4 activity, a near complete functional loss state. They are not literal genetic zero. The experimental approximately 35% AE4 null reduction is contextual evidence, not a fitted target or an exact acceptance band. Detailed experimental time signatures are not acceptance criteria. The inherited absolute gland scale discrepancy remains recorded in Task 14 and is nonblocking for these ratios.

Scientific source: `a2fde97e30d79601e7455d9f1d77d8d67609cb49`. The saved selected connected 5% resting state was used without another resting solve for every root. Exactly 30 new production Radau trajectories use R1, frozen NKCC1 normalisation, CCh + IPR, the three declared calcium inputs, the inherited tolerances and the original 602 point WT grid over 0 to 600 s. Each denominator is the matching saved Task 14 WT production total, checked against the raw trajectory and the commit anchored artifact hashes. Model and root hashes were verified before and after execution. No model parameter, root, regulation law or observation scale was fitted or changed. No root or calcium was selected by phenotype agreement.

Maximum dimensionless conservation residual ratio: 0.002277289 (inherited numerical limit 1). Full solver messages, positivity, charge, current, carbon, water, ionic states, pH and flow/cumulative traces are saved alongside the pair table.

| Root | Calcium (uM) | AE4 5% / WT | Reduction (%) | AE2 reduction (%) | Numerical status |
| --- | ---: | ---: | ---: | ---: | --- |
| N_ABS_NKCC_S4_PHIGH_KHIGH_AE4NA05_P1_H1:B00 | 0.10 | 1.150954 | -15.095401 | -0.132090 | COMPLETE |
| N_ABS_NKCC_S4_PHIGH_KHIGH_AE4NA20_P2_H1:B00 | 0.10 | 1.092074 | -9.207378 | -0.034164 | COMPLETE |
| N_ABS_NKCC_S4_PHIGH_KLOW_AE4NA05_P1_H1:B00 | 0.10 | 1.140521 | -14.052130 | -0.125461 | COMPLETE |
| N_ABS_NKCC_S4_PHIGH_KLOW_AE4NA20_P2_H1:B00 | 0.10 | 1.082287 | -8.228692 | -0.031795 | COMPLETE |
| N_ABS_NKCC_S4_PLOW_KHIGH_AE4NA05_P1_H1:B00 | 0.10 | 1.163180 | -16.317970 | -0.139193 | COMPLETE |
| N_ABS_NKCC_S4_PLOW_KHIGH_AE4NA20_P2_H1:B00 | 0.10 | 1.124933 | -12.493289 | -0.038676 | COMPLETE |
| N_ABS_NKCC_S4_PLOW_KLOW_AE4NA05_P1_H1:B00 | 0.10 | 1.151855 | -15.185457 | -0.132237 | COMPLETE |
| N_ABS_NKCC_S4_PLOW_KLOW_AE4NA20_P2_H1:B00 | 0.10 | 1.112949 | -11.294851 | -0.036090 | COMPLETE |
| N_ABS_NKCC_S4_PNOM_KMID_AE4NA05_P1_H1:B00 | 0.10 | 1.153074 | -15.307414 | -0.133072 | COMPLETE |
| N_ABS_NKCC_S4_PNOM_KMID_AE4NA20_P2_H1:B00 | 0.10 | 1.104898 | -10.489850 | -0.035475 | COMPLETE |
| N_ABS_NKCC_S4_PHIGH_KHIGH_AE4NA05_P1_H1:B00 | 0.25 | 1.191545 | -19.154460 | -0.149313 | COMPLETE |
| N_ABS_NKCC_S4_PHIGH_KHIGH_AE4NA20_P2_H1:B00 | 0.25 | 1.162443 | -16.244273 | -0.041737 | COMPLETE |
| N_ABS_NKCC_S4_PHIGH_KLOW_AE4NA05_P1_H1:B00 | 0.25 | 1.179915 | -17.991463 | -0.143647 | COMPLETE |
| N_ABS_NKCC_S4_PHIGH_KLOW_AE4NA20_P2_H1:B00 | 0.25 | 1.145429 | -14.542950 | -0.039264 | COMPLETE |
| N_ABS_NKCC_S4_PLOW_KHIGH_AE4NA05_P1_H1:B00 | 0.25 | 1.202714 | -20.271414 | -0.155982 | COMPLETE |
| N_ABS_NKCC_S4_PLOW_KHIGH_AE4NA20_P2_H1:B00 | 0.25 | 1.195148 | -19.514795 | -0.046157 | COMPLETE |
| N_ABS_NKCC_S4_PLOW_KLOW_AE4NA05_P1_H1:B00 | 0.25 | 1.190639 | -19.063912 | -0.150203 | COMPLETE |
| N_ABS_NKCC_S4_PLOW_KLOW_AE4NA20_P2_H1:B00 | 0.25 | 1.176467 | -17.646658 | -0.043554 | COMPLETE |
| N_ABS_NKCC_S4_PNOM_KMID_AE4NA05_P1_H1:B00 | 0.25 | 1.192377 | -19.237654 | -0.150418 | COMPLETE |
| N_ABS_NKCC_S4_PNOM_KMID_AE4NA20_P2_H1:B00 | 0.25 | 1.171664 | -17.166380 | -0.042955 | COMPLETE |
| N_ABS_NKCC_S4_PHIGH_KHIGH_AE4NA05_P1_H1:B00 | 0.50 | 1.202403 | -20.240286 | -0.153520 | COMPLETE |
| N_ABS_NKCC_S4_PHIGH_KHIGH_AE4NA20_P2_H1:B00 | 0.50 | 1.184001 | -18.400108 | -0.043811 | COMPLETE |
| N_ABS_NKCC_S4_PHIGH_KLOW_AE4NA05_P1_H1:B00 | 0.50 | 1.191979 | -19.197946 | -0.148882 | COMPLETE |
| N_ABS_NKCC_S4_PHIGH_KLOW_AE4NA20_P2_H1:B00 | 0.50 | 1.166968 | -16.696831 | -0.041646 | COMPLETE |
| N_ABS_NKCC_S4_PLOW_KHIGH_AE4NA05_P1_H1:B00 | 0.50 | 1.213305 | -21.330519 | -0.160080 | COMPLETE |
| N_ABS_NKCC_S4_PLOW_KHIGH_AE4NA20_P2_H1:B00 | 0.50 | 1.216630 | -21.663007 | -0.048210 | COMPLETE |
| N_ABS_NKCC_S4_PLOW_KLOW_AE4NA05_P1_H1:B00 | 0.50 | 1.202563 | -20.256288 | -0.155382 | COMPLETE |
| N_ABS_NKCC_S4_PLOW_KLOW_AE4NA20_P2_H1:B00 | 0.50 | 1.198226 | -19.822616 | -0.045947 | COMPLETE |
| N_ABS_NKCC_S4_PNOM_KMID_AE4NA05_P1_H1:B00 | 0.50 | 1.203481 | -20.348075 | -0.154958 | COMPLETE |
| N_ABS_NKCC_S4_PNOM_KMID_AE4NA20_P2_H1:B00 | 0.50 | 1.193040 | -19.304026 | -0.045142 | COMPLETE |

All 11 focused tests pass, including direct integration of the saved flow traces, exact WT denominator reuse, all 30 pair identities, saved initial states, unchanged model parameter hashes and ensemble arithmetic. The tests run no ODEs or resting solves. Post execution interpretation uses saved results only; the model and production trajectories remain unchanged.

Machine readable results: `results/14B_five_percent_ae4_check/five_percent_ae4_flow.csv`, `summary.json`, `interpretation.json`, `frozen_inputs.json`, `postrun_integrity.json` and `trajectories/`. Execution and test logs and the artifact hash ledger are saved in the same results directory.
