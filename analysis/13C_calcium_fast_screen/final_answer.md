Increasing calcium to 0.25 or 0.50 µM does not close the absolute WT secretion deficit in any of the ten frozen roots. All twenty loose trajectories and both endpoint confirmation trajectories pass the inherited numerical and sustainment gates.

The best 0.50 µM root, R02, reaches only 60.07% of the required minimum minute flow. Even at the unchanged generous one SMG size ceiling, its limiting minute gives 5.405993 µL/min, against the 9 to 10 µL/min requirement.

The comparison below uses the ten saved Task 13B production rows at 0.10 µM, the twenty new loose screen rows, and the production Radau result for the single confirmed 0.50 µM root. No 0.10 µM trajectory was recomputed. Each range is across roots. Means and cumulative flow use the inherited 122 sample grid and trapezoidal output calculation.

| Ca (µM) | Ca gate | Mean per cell flow (pL/s) | 600 s total (pL/cell) | Paired mean fold increase | Minute shape passes / 10 | Absolute passes / 10 |
| --- | --- | --- | --- | --- | --- | --- |
| 0.10 | 0.198603 | 0.001190427 to 0.001510129 | 0.714256 to 0.906077 | 1.000000 | 10 | 0 |
| 0.25 | 0.485688 | 0.001354316 to 0.001732109 | 0.812590 to 1.039265 | 1.135672 to 1.152270 | 5 | 0 |
| 0.50 | 0.722066 | 0.001399910 to 0.001796221 | 0.839946 to 1.077733 | 1.172592 to 1.198540 | 5 | 0 |


Thus 0.25 µM fails the absolute gate, and 0.50 µM also fails it. The negative absolute conclusion is consistent across all ten roots, but their flow shapes differ: the five AE4NA05_P1 roots have empty minute scale intervals at both new amplitudes; the five AE4NA20_P2 roots retain nonempty minute intervals. An empty interval is a WT flow shape failure even before applying the gland size limit. It is distinct from the broader sustainment test, which all roots pass.

There is no new numerical or inherited physical domain failure. The model exposes no separate source based stimulated Na, K, Cl or pH target bands, so this statement does not establish broader physiological validity. Resting calibration ranges have not been reused as new stimulated acceptance thresholds. Endpoint concentrations, pH, volume and membrane potentials are reported below.

At each calcium input the following table gives the per cell flow range across roots at the unchanged minute landmarks. It retains the decline within each root instead of hiding it in an average.

| Time (s) | 0.10 µM (pL/s) | 0.25 µM (pL/s) | 0.50 µM (pL/s) |
| --- | --- | --- | --- |
| 60 | 0.001222593 to 0.001511396 | 0.001431325 to 0.001766925 | 0.001493169 to 0.001845421 |
| 120 | 0.001213106 to 0.001516345 | 0.001404901 to 0.001753937 | 0.001460020 to 0.001824019 |
| 180 | 0.001204224 to 0.001514550 | 0.001379350 to 0.001737830 | 0.001428066 to 0.001802035 |
| 240 | 0.001195012 to 0.001511401 | 0.001356462 to 0.001724762 | 0.001400236 to 0.001785123 |
| 300 | 0.001186522 to 0.001508578 | 0.001336924 to 0.001715184 | 0.001376881 to 0.001773072 |
| 360 | 0.001179065 to 0.001506373 | 0.001320533 to 0.001708414 | 0.001357515 to 0.001764720 |
| 420 | 0.001172647 to 0.001504730 | 0.001306850 to 0.001703695 | 0.001341481 to 0.001758991 |
| 480 | 0.001167170 to 0.001503528 | 0.001295421 to 0.001700424 | 0.001328167 to 0.001755074 |
| 540 | 0.001162510 to 0.001502654 | 0.001285848 to 0.001698158 | 0.001317063 to 0.001752394 |
| 600 | 0.001158548 to 0.001502020 | 0.001277799 to 0.001696586 | 0.001307758 to 0.001750558 |


Root identifiers used in the remaining tables:

| Root | Frozen root ID |
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


The observation interval is the intersection of [9/q, 10/q] over all ten minute flows within each root. Alternative roots do not share a fitted scale. Every interval fails the fixed ceiling 3088.15386. The required cell count below is the lower scale divided by 0.00006, even when infeasible; it is a requirement for that model, not an estimate of the actual gland count. The ceiling is 51,469,231 cells.

| Ca (µM) | Root | Scale lower | Scale upper | Required cells (millions) | q600/q60 |
| --- | --- | --- | --- | --- | --- |
| 0.10 | R01 | 7459.737 | 7863.894 | 124.329 | 0.948761 |
| 0.10 | R02 | 5991.930 | 6594.806 | 99.866 | 0.993796 |
| 0.10 | R03 | 7512.164 | 7935.464 | 125.203 | 0.950714 |
| 0.10 | R04 | 6110.226 | 6735.877 | 101.837 | 0.993831 |
| 0.10 | R05 | 7729.968 | 8122.270 | 128.833 | 0.945676 |
| 0.10 | R06 | 6343.314 | 6969.860 | 105.722 | 0.989674 |
| 0.10 | R07 | 7768.343 | 8179.337 | 129.472 | 0.947616 |
| 0.10 | R08 | 6446.737 | 7088.068 | 107.446 | 0.989533 |
| 0.10 | R09 | 7593.934 | 7998.931 | 126.566 | 0.947999 |
| 0.10 | R10 | 6194.720 | 6817.239 | 103.245 | 0.991808 |
| 0.25 | R01 | 6774.654 | 6687.901 | 112.911 | 0.888475 |
| 0.25 | R02 | 5304.770 | 5659.550 | 88.413 | 0.960192 |
| 0.25 | R03 | 6813.895 | 6787.085 | 113.565 | 0.896459 |
| 0.25 | R04 | 5385.487 | 5769.640 | 89.758 | 0.964198 |
| 0.25 | R05 | 7010.070 | 6889.800 | 116.834 | 0.884559 |
| 0.25 | R06 | 5607.956 | 5952.104 | 93.466 | 0.955231 |
| 0.25 | R07 | 7043.361 | 6986.532 | 117.389 | 0.892738 |
| 0.25 | R08 | 5681.157 | 6054.839 | 94.686 | 0.959198 |
| 0.25 | R09 | 6900.451 | 6828.608 | 115.008 | 0.890630 |
| 0.25 | R10 | 5477.968 | 5842.737 | 91.299 | 0.959930 |
| 0.50 | R01 | 6625.534 | 6404.180 | 110.426 | 0.869932 |
| 0.50 | R02 | 5141.217 | 5418.819 | 85.687 | 0.948596 |
| 0.50 | R03 | 6657.665 | 6506.101 | 110.961 | 0.879511 |
| 0.50 | R04 | 5204.021 | 5516.012 | 86.734 | 0.953957 |
| 0.50 | R05 | 6853.659 | 6593.688 | 114.228 | 0.865861 |
| 0.50 | R06 | 5433.767 | 5696.550 | 90.563 | 0.943525 |
| 0.50 | R07 | 6882.009 | 6697.166 | 114.700 | 0.875827 |
| 0.50 | R08 | 5491.283 | 5790.964 | 91.521 | 0.949117 |
| 0.50 | R09 | 6748.610 | 6545.468 | 112.477 | 0.872909 |
| 0.50 | R10 | 5305.130 | 5593.537 | 88.419 | 0.948927 |


Endpoint intracellular states at 600 s. Potentials are apical and basolateral, in that order. The frozen reference profile did not save endpoint potentials; those entries remain unavailable rather than being recomputed.

| Ca (µM) | Root | Na (mM) | K (mM) | Cl (mM) | pH | Volume (pL) | Va (mV) | Vb (mV) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0.10 | R01 | 26.703772 | 98.849953 | 51.393320 | 6.837909 | 1.294310 | Unavailable | Unavailable |
| 0.10 | R02 | 19.213798 | 112.574113 | 53.520262 | 7.036759 | 1.312552 | Unavailable | Unavailable |
| 0.10 | R03 | 27.150653 | 97.767067 | 51.101160 | 6.836895 | 1.300479 | Unavailable | Unavailable |
| 0.10 | R04 | 19.217896 | 112.324814 | 53.562625 | 7.036227 | 1.317848 | Unavailable | Unavailable |
| 0.10 | R05 | 26.262947 | 99.781737 | 51.815830 | 6.838384 | 1.293299 | Unavailable | Unavailable |
| 0.10 | R06 | 18.965515 | 113.582070 | 54.177933 | 7.038163 | 1.311766 | Unavailable | Unavailable |
| 0.10 | R07 | 26.740421 | 98.616105 | 51.487360 | 6.837303 | 1.299724 | Unavailable | Unavailable |
| 0.10 | R08 | 18.994489 | 113.231297 | 54.164768 | 7.037480 | 1.317342 | Unavailable | Unavailable |
| 0.10 | R09 | 26.803581 | 98.563732 | 51.379488 | 6.837476 | 1.297468 | Unavailable | Unavailable |
| 0.10 | R10 | 19.128642 | 112.809174 | 53.791428 | 7.036946 | 1.315171 | Unavailable | Unavailable |
| 0.25 | R01 | 30.841965 | 93.684371 | 47.796226 | 6.830499 | 1.244676 | -31.985098 | -73.793919 |
| 0.25 | R02 | 20.804223 | 110.110890 | 49.925920 | 7.026524 | 1.257137 | -31.394522 | -78.603373 |
| 0.25 | R03 | 31.077137 | 92.870311 | 47.767029 | 6.830210 | 1.254308 | -32.020488 | -74.608697 |
| 0.25 | R04 | 20.750063 | 109.939460 | 50.085650 | 7.026479 | 1.263947 | -31.312510 | -79.470449 |
| 0.25 | R05 | 30.266662 | 94.762654 | 48.210983 | 6.830906 | 1.243348 | -31.724744 | -74.073057 |
| 0.25 | R06 | 20.462679 | 111.260178 | 50.666081 | 7.028290 | 1.257222 | -30.928024 | -78.861929 |
| 0.25 | R07 | 30.526825 | 93.874497 | 48.165933 | 6.830591 | 1.253533 | -31.769869 | -74.878921 |
| 0.25 | R08 | 20.432367 | 110.989432 | 50.790357 | 7.028125 | 1.264651 | -30.869788 | -79.715800 |
| 0.25 | R09 | 30.741021 | 93.643302 | 47.958625 | 6.830499 | 1.250080 | -31.891678 | -74.369703 |
| 0.25 | R10 | 20.636395 | 110.465979 | 50.329071 | 7.027202 | 1.261457 | -31.150050 | -79.202181 |
| 0.50 | R01 | 32.097864 | 92.170433 | 46.895514 | 6.828582 | 1.232859 | -31.897572 | -74.740182 |
| 0.50 | R02 | 21.302848 | 109.368909 | 48.938633 | 7.023476 | 1.242666 | -31.205472 | -79.631306 |
| 0.50 | R03 | 32.260025 | 91.443344 | 46.930151 | 6.828487 | 1.243249 | -31.888411 | -75.275845 |
| 0.50 | R04 | 21.242508 | 109.202524 | 49.103058 | 7.023506 | 1.249447 | -31.118134 | -80.247941 |
| 0.50 | R05 | 31.482840 | 93.290276 | 47.304627 | 6.828958 | 1.231407 | -31.646611 | -75.043874 |
| 0.50 | R06 | 20.928932 | 110.564235 | 49.702458 | 7.025358 | 1.242983 | -30.741180 | -79.907946 |
| 0.50 | R07 | 31.664493 | 92.496093 | 47.331106 | 6.828852 | 1.242452 | -31.642663 | -75.568973 |
| 0.50 | R08 | 20.888737 | 110.304137 | 49.842638 | 7.025294 | 1.250548 | -30.670221 | -80.508762 |
| 0.50 | R09 | 31.921654 | 92.216652 | 47.104836 | 6.828704 | 1.238806 | -31.775404 | -75.172768 |
| 0.50 | R10 | 21.108218 | 109.760454 | 49.374245 | 7.024299 | 1.247351 | -30.947291 | -80.097605 |


The machine readable screen also retains sampled post stimulus flow extrema, all endpoint and trajectory domain diagnostics, each accounting residual in its inherited units, and solver counters. The minimum over t > 0 includes the inherited right limit at 0.000001 s, whose state is still the resting state; it is not the minimum of the ten minute measurements.

The predeclared 90% proximity threshold was 0.002622926307 pL/s. No root reached it at either new amplitude. Accordingly only R02 at 0.50 µM was confirmed with production Radau and BDF. The maximum relative Radau versus BDF differences were 3.37e-07 for states, 1.1e-06 for flow and 9.13e-09 for cumulative flow, all below the inherited 0.0001 gate.

Measured numerical phase time was 1.436 s for the twenty case primary screen including its pilot, and 2.436 s including the two confirmations, on nine available worker processes. These timings exclude preparation, interpreter startup, analysis and testing. The original ninety trajectory primary screen was not executed; its runtime can only be estimated from the observed production timing. The exact primary trajectory reduction is 77.8%. The runtime report gives the assumptions and the estimated elapsed saving.

A matched SMG calcium measurement remains useful for choosing the physiological input, but it is no longer the deciding measurement for the proposed amplitude only rescue within the tested range: even 0.50 µM remains far below the required flow. These two tests do not prove failure at every untested amplitude or outside the tested range. The unresolved distinction between cellular secretion capacity and the cell to gland observation map still calls for matched absolute secretion per cell. No calcium threshold was refined.

This result concerns WT calcium only. It does not identify an AE4 phenotype, establish AE4 insufficiency, or assess the sealed target.

CALCIUM 0.50 UM REMAINS INSUFFICIENT TO CLOSE WT ABSOLUTE FLOW DEFICIT
