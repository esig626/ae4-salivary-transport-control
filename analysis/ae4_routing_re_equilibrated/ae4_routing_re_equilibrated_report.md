# AE4 routing experiment after resting equilibration

## Primary result

5% AE4 increases total secretion in every matched comparison. All 30 WT and 5% AE4 pairs passed the inherited production and conservation gates.

The resting equations were solved directly under an explicit REST protocol. No resting integration, fitting, recalibration, parameter change, or phenotype guided root selection was used. For each root, expression and calcium comparison, the inherited saved rest and the matching earlier routed 600 s endpoint were used as independent warm starts.

| Calcium (µM) | Re equilibrated 5% AE4 change across roots | Previous non equilibrated change across roots |
| ---: | ---: | ---: |
| 0.10 | +0.373% to +0.686% | +3.488% to +6.178% |
| 0.25 | +0.515% to +0.946% | +4.959% to +7.831% |
| 0.50 | +0.581% to +1.037% | +5.499% to +8.335% |

## Secretion comparison for every root and calcium

Totals are model secretion per cell over 0 to 600 s, in pL. Positive percentages mean that 5% AE4 increased secretion relative to matched WT.

| Root | Calcium (µM) | WT total (pL) | 5% AE4 total (pL) | Ratio | Change | Previous change | Direction changed? |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| R01 | 0.10 | 0.856485950 | 0.862193098 | 1.006663446 | +0.6663% | +5.7939% | no |
| R01 | 0.25 | 1.008588949 | 1.017900104 | 1.009231863 | +0.9232% | +7.5134% | no |
| R01 | 0.50 | 1.051178073 | 1.061830637 | 1.010133929 | +1.0134% | +8.0361% | no |
| R02 | 0.10 | 0.989759194 | 0.993931014 | 1.004214985 | +0.4215% | +4.0162% | no |
| R02 | 0.25 | 1.211747201 | 1.218814076 | 1.005831971 | +0.5832% | +5.7030% | no |
| R02 | 0.50 | 1.280989011 | 1.289332582 | 1.006513382 | +0.6513% | +6.2478% | no |
| R03 | 0.10 | 0.841334891 | 0.846422623 | 1.006047214 | +0.6047% | +5.3998% | no |
| R03 | 0.25 | 0.989614672 | 0.997965837 | 1.008438805 | +0.8439% | +7.0336% | no |
| R03 | 0.50 | 1.033443642 | 1.043116223 | 1.009359563 | +0.9360% | +7.5896% | no |
| R04 | 0.10 | 0.960732286 | 0.964313930 | 1.003728035 | +0.3728% | +3.4883% | no |
| R04 | 0.25 | 1.173608877 | 1.179651000 | 1.005148328 | +0.5148% | +4.9588% | no |
| R04 | 0.50 | 1.244892053 | 1.252127374 | 1.005812006 | +0.5812% | +5.4989% | no |
| R05 | 0.10 | 0.836584406 | 0.842322273 | 1.006858684 | +0.6859% | +6.1782% | no |
| R05 | 0.25 | 0.985729555 | 0.995054603 | 1.009460046 | +0.9460% | +7.8310% | no |
| R05 | 0.50 | 1.027438018 | 1.038091322 | 1.010368804 | +1.0369% | +8.3350% | no |
| R06 | 0.10 | 0.965840318 | 0.970010845 | 1.004318030 | +0.4318% | +4.6257% | no |
| R06 | 0.25 | 1.182353766 | 1.189372225 | 1.005936006 | +0.5936% | +6.2360% | no |
| R06 | 0.50 | 1.249578625 | 1.257850677 | 1.006619874 | +0.6620% | +6.7600% | no |
| R07 | 0.10 | 0.822811075 | 0.827936355 | 1.006228988 | +0.6229% | +5.7660% | no |
| R07 | 0.25 | 0.967707558 | 0.976080465 | 1.008652311 | +0.8652% | +7.3494% | no |
| R07 | 0.50 | 1.010362350 | 1.020040404 | 1.009578795 | +0.9579% | +7.8901% | no |
| R08 | 0.10 | 0.938825325 | 0.942422387 | 1.003831449 | +0.3831% | +4.0533% | no |
| R08 | 0.25 | 1.146114954 | 1.152136235 | 1.005253645 | +0.5254% | +5.4682% | no |
| R08 | 0.50 | 1.215077256 | 1.222269515 | 1.005919178 | +0.5919% | +5.9932% | no |
| R09 | 0.10 | 0.843062901 | 0.848541644 | 1.006498616 | +0.6499% | +5.8295% | no |
| R09 | 0.25 | 0.990271902 | 0.999178350 | 1.008993943 | +0.8994% | +7.4704% | no |
| R09 | 0.50 | 1.032105815 | 1.042318071 | 1.009894583 | +0.9895% | +7.9915% | no |
| R10 | 0.10 | 0.969497209 | 0.973435424 | 1.004062121 | +0.4062% | +4.0666% | no |
| R10 | 0.25 | 1.183500704 | 1.190118640 | 1.005591831 | +0.5592% | +5.6060% | no |
| R10 | 0.50 | 1.251734199 | 1.259572185 | 1.006261702 | +0.6262% | +6.1328% | no |

## Resting states

The resting protocol is independent of the later stimulation calcium. The three calcium specific solves agree to a maximum normalised coordinate distance of 1.362e-11. The table therefore lists each common resting state once.

| Root | AE4 | Cell Na | Cell K | Cell Cl | Cell pH | Cell volume | Lumen Na | Lumen K | Lumen Cl | Lumen pH | Lumen volume |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| R01 | 100% | 11.473515 | 117.084552 | 60.927671 | 6.889590 | 1.453807 | 122.924389 | 27.424767 | 145.721110 | 6.907399 | 0.105422 |
| R01 | 5% | 11.481708 | 119.714171 | 60.165180 | 7.145212 | 1.540698 | 122.758202 | 27.586004 | 145.738754 | 6.898588 | 0.105451 |
| R02 | 100% | 7.882512 | 127.841101 | 67.147723 | 7.149021 | 1.603649 | 119.968082 | 30.398947 | 146.176989 | 6.868799 | 0.106029 |
| R02 | 5% | 7.882791 | 129.896347 | 66.558417 | 7.350373 | 1.740618 | 119.838039 | 30.524457 | 146.187169 | 6.861182 | 0.106051 |
| R03 | 100% | 11.708547 | 116.041990 | 59.863376 | 6.886039 | 1.446044 | 127.933517 | 22.414728 | 146.030393 | 6.876965 | 0.105401 |
| R03 | 5% | 11.716834 | 118.691956 | 59.103517 | 7.140557 | 1.530397 | 127.812563 | 22.530630 | 146.048320 | 6.867900 | 0.105427 |
| R04 | 100% | 7.965992 | 127.309123 | 66.136059 | 7.145050 | 1.585506 | 125.852273 | 24.512495 | 146.491348 | 6.833985 | 0.105963 |
| R04 | 5% | 7.966039 | 129.385395 | 65.550981 | 7.345863 | 1.717592 | 125.758454 | 24.601717 | 146.501692 | 6.826206 | 0.105983 |
| R05 | 100% | 11.404110 | 117.591503 | 61.298532 | 6.890138 | 1.452651 | 121.999145 | 28.346037 | 145.725550 | 6.905290 | 0.105300 |
| R05 | 5% | 11.412071 | 120.210640 | 60.538228 | 7.146216 | 1.540335 | 121.829556 | 28.510672 | 145.743126 | 6.896496 | 0.105329 |
| R06 | 100% | 7.854235 | 128.449428 | 67.504647 | 7.149517 | 1.598303 | 118.907863 | 31.454935 | 146.181367 | 6.866500 | 0.105899 |
| R06 | 5% | 7.854374 | 130.474672 | 66.924581 | 7.351205 | 1.736471 | 118.776945 | 31.581302 | 146.191333 | 6.858900 | 0.105921 |
| R07 | 100% | 11.646068 | 116.496172 | 60.193439 | 6.886515 | 1.444970 | 126.873795 | 23.470749 | 146.038401 | 6.874557 | 0.105287 |
| R07 | 5% | 11.654183 | 119.137286 | 59.435170 | 7.141438 | 1.529995 | 126.749577 | 23.589912 | 146.056246 | 6.865508 | 0.105314 |
| R08 | 100% | 7.941683 | 127.857984 | 66.450751 | 7.145458 | 1.580534 | 124.656021 | 25.704836 | 146.499086 | 6.831375 | 0.105843 |
| R08 | 5% | 7.941621 | 129.908014 | 65.873531 | 7.346569 | 1.713599 | 124.561198 | 25.795046 | 146.509224 | 6.823614 | 0.105863 |
| R09 | 100% | 11.594292 | 116.646263 | 60.398947 | 6.887592 | 1.447917 | 124.650264 | 25.697373 | 145.897346 | 6.889884 | 0.105379 |
| R09 | 5% | 11.602835 | 119.283515 | 59.637226 | 7.142690 | 1.533431 | 124.502552 | 25.840089 | 145.915099 | 6.880945 | 0.105407 |
| R10 | 100% | 7.928091 | 127.741284 | 66.607118 | 7.146633 | 1.588815 | 122.000969 | 28.363988 | 146.356237 | 6.848933 | 0.105967 |
| R10 | 5% | 7.928317 | 129.797530 | 66.022510 | 7.347746 | 1.722833 | 121.886563 | 28.473833 | 146.366413 | 6.841242 | 0.105989 |

Concentrations are mM and volumes are pL.

With no permitted recalibration, the new routed WT equilibria do not preserve the inherited WT resting calibration. WT cell chloride is 59.863 to 67.505 mM, cell pH is 6.886 to 7.150, and cell volume is 1.445 to 1.604 pL. These are the genuine mathematical resting states of the routing intervention with frozen parameters, rather than a preservation of the original fitted WT resting physiology.

## Resting residuals and gates

Every value below is from the complete resting evaluation, not only the ten solved rows.

| Root | AE4 | Scaled independent RHS | Amount RHS (fmol/s) | Volume RHS (pL/s) | Omitted RHS (fmol/s) | Regulatory RHS (/s) | Current (A) | Conservation ratio | Gate |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| R01 | 100% | 2.754e-12 | 2.845e-16 | 2.754e-16 | 2.845e-16 | 0.000e+00 | 9.694e-27 | 9.694e-07 | pass |
| R01 | 5% | 5.161e-13 | 2.220e-16 | 5.161e-17 | 2.151e-16 | 0.000e+00 | 4.524e-26 | 4.840e-06 | pass |
| R02 | 100% | 5.651e-12 | 2.880e-16 | 5.651e-16 | 2.880e-16 | 0.000e+00 | 9.694e-27 | 1.214e-06 | pass |
| R02 | 5% | 8.116e-12 | 2.498e-16 | 8.116e-16 | 1.084e-16 | 0.000e+00 | 9.694e-27 | 1.350e-04 | pass |
| R03 | 100% | 2.548e-12 | 2.733e-14 | 2.548e-16 | 2.733e-14 | 0.000e+00 | 3.231e-27 | 8.882e-07 | pass |
| R03 | 5% | 1.193e-11 | 7.216e-16 | 1.193e-15 | 1.023e-16 | 0.000e+00 | 4.524e-26 | 4.797e-06 | pass |
| R04 | 100% | 1.599e-11 | 3.053e-16 | 1.599e-15 | 2.151e-16 | 0.000e+00 | 1.616e-26 | 1.776e-06 | pass |
| R04 | 5% | 6.332e-12 | 2.776e-16 | 6.332e-16 | 2.299e-16 | 0.000e+00 | 3.231e-27 | 2.203e-04 | pass |
| R05 | 100% | 7.730e-12 | 9.714e-14 | 7.730e-16 | 9.714e-14 | 0.000e+00 | 6.462e-27 | 1.776e-06 | pass |
| R05 | 5% | 1.360e-11 | 2.220e-16 | 1.360e-15 | 2.056e-16 | 0.000e+00 | 1.939e-26 | 7.105e-06 | pass |
| R06 | 100% | 1.546e-12 | 1.388e-16 | 1.546e-16 | 5.117e-17 | 0.000e+00 | 1.292e-26 | 2.359e-06 | pass |
| R06 | 5% | 4.218e-12 | 2.220e-16 | 4.218e-16 | 5.551e-17 | 0.000e+00 | 1.939e-26 | 1.421e-04 | pass |
| R07 | 100% | 1.816e-11 | 8.323e-15 | 1.816e-15 | 8.323e-15 | 0.000e+00 | 1.292e-26 | 1.292e-06 | pass |
| R07 | 5% | 9.571e-12 | 3.331e-16 | 9.571e-16 | 9.801e-17 | 0.000e+00 | 2.585e-26 | 3.010e-06 | pass |
| R08 | 100% | 1.086e-11 | 2.046e-14 | 1.086e-15 | 2.046e-14 | 0.000e+00 | 2.262e-26 | 3.261e-06 | pass |
| R08 | 5% | 5.573e-13 | 2.073e-16 | 5.573e-17 | 2.073e-16 | 0.000e+00 | 2.262e-26 | 2.203e-04 | pass |
| R09 | 100% | 6.037e-12 | 3.331e-16 | 6.037e-16 | 1.180e-16 | 0.000e+00 | 2.262e-26 | 4.424e-06 | pass |
| R09 | 5% | 2.277e-13 | 1.431e-16 | 2.277e-17 | 1.431e-16 | 0.000e+00 | 1.616e-26 | 1.616e-06 | pass |
| R10 | 100% | 6.731e-12 | 1.518e-16 | 6.731e-16 | 1.518e-16 | 0.000e+00 | 9.694e-27 | 1.648e-06 | pass |
| R10 | 5% | 3.344e-12 | 4.586e-15 | 3.344e-16 | 4.586e-15 | 0.000e+00 | 9.694e-27 | 1.812e-04 | pass |

The complete calcium specific values, both warm start attempts, all conservation residuals, and both compartment states are in `resting_states.csv` and `warm_start_attempts.csv`. The full 30 row comparison, including both WT and 5% resting diagnostics on every row, is in `comparison.csv`.

## Intervention and verification

The existing net AE4 chloride source, bicarbonate source, and total cation source were unchanged at every sampled state. Only the Na and K shares of that total cation source were reassigned. No opposing Na and K AE4 sources occurred. The net cycle did not reverse at any sampled state in this re equilibrated run.

All 60 resting cases and 60 stimulated trajectories passed. The maximum warm start root distance was 2.160e-11, below the inherited clustering tolerance 2.000e-05. The maximum stimulation conservation ratio was 1.364e-03.

| Root | Inherited root identifier |
| --- | --- |
| R01 | `N_ABS_NKCC_S4_PHIGH_KHIGH_AE4NA05_P1_H1:B00` |
| R02 | `N_ABS_NKCC_S4_PHIGH_KHIGH_AE4NA20_P2_H1:B00` |
| R03 | `N_ABS_NKCC_S4_PHIGH_KLOW_AE4NA05_P1_H1:B00` |
| R04 | `N_ABS_NKCC_S4_PHIGH_KLOW_AE4NA20_P2_H1:B00` |
| R05 | `N_ABS_NKCC_S4_PLOW_KHIGH_AE4NA05_P1_H1:B00` |
| R06 | `N_ABS_NKCC_S4_PLOW_KHIGH_AE4NA20_P2_H1:B00` |
| R07 | `N_ABS_NKCC_S4_PLOW_KLOW_AE4NA05_P1_H1:B00` |
| R08 | `N_ABS_NKCC_S4_PLOW_KLOW_AE4NA20_P2_H1:B00` |
| R09 | `N_ABS_NKCC_S4_PNOM_KMID_AE4NA05_P1_H1:B00` |
| R10 | `N_ABS_NKCC_S4_PNOM_KMID_AE4NA20_P2_H1:B00` |
