# Task 38: frozen AE4 perturbation predictions

`AE4_PERTURBATION_VALIDATION_COMPLETE`

Both authorized trajectories passed every inherited Task 37 physical and conservation gate. Cumulative secretion was **98.4315% of frozen WT at AE4=0.05** and **98.1605% at AE4=0.00**. The null prediction is a small secretion loss, with substantial NKCC1 compensation. It does not reproduce the held-out approximately 35% deficit; nothing was fitted or changed in response.

Work and publication are confined to `codex/task-38-ae4-perturbation-validation`, prepared head `7aaa53414714bae459227b4b972c63d49e13b79e`. Frozen WT reference: `2ba0e4e0c2a9e940aa96ed86c3c72456afd0e3aa`. `AGENTS.md` was read first. All 97 prepared files matched their GitHub blob hashes and remained unchanged; the prepared Git tree was `d3e17198e16ff0c3a8cd06b1472a23d621880370`. Only the Task 38 instructions differ from the frozen WT tree.

## Frozen setup and execution

The Task 37 factory loaded the accepted Task 31 R09 WT resting checkpoint directly. All 12 core coordinates and the resting regulatory coordinate (zero) were copied bitwise. Full initial-state SHA-256: `879a650761e9044691ab6b2a9b1bfc3be8f951ece176925117596cc2baa31abe`. No genotype-specific resting state was solved. The only changed scientific input was `genotype.ae4_expression`, first 0.05 and then 0.00; all other genotype fields were identical to WT.

The NHE1 carrier amount remained `2.339370005697548e-05 fmol`. The unchanged electrogenic NBC was 1 Na : 2 HCO3, with capacity `0.11570913197464398 fmol/s`, voltage-dependent affinity/current closure and stimulus recruitment. The stimulated NHE1 multiplier remained 2.3. Routed AE4 retained the 1 Cl : 1 monovalent cation : 2 HCO3 stoichiometry and inherited beta/PKA regulation. NKCC1, AE2, pump, channels, CO2, paracellular transport, water, bath and geometry were unchanged.

Both trajectories used CCh 0.3 uM, isoproterenol 5 uM, stimulated Ca 0.25 uM and beta occupancy 1.0 for 600 s. The inherited onset convention is REST at exactly zero and stimulation for 0<t<=600; integration starts at 1e-6 s with the unchanged initial state. Each used one continuous production Radau integration: rtol 1e-7; amount/regulatory atol 1e-10; volume atol 1e-12 pL; max step 2 s. NumPy 2.3.5 and SciPy 1.17.0 match Task 37. No numerical retry was used.

## Secretion

| Readout | Frozen WT | AE4=0.05 | AE4=0.00 |
| --- | --- | --- | --- |
| Cumulative 0-600 s (pL) | 0.986995055625 | 0.97151438779 | 0.968838809526 |
| Ratio to frozen WT | 1 | 0.984315354219 | 0.981604521729 |
| Cumulative deficit (%) | 0 | 1.568465 | 1.839548 |
| Mean q_out 60-600 s (pL/s) | 0.00162624074525 | 0.00159981579494 | 0.00159501996561 |
| q_out at 600 s (pL/s) | 0.00156842215184 | 0.00158503240022 | 0.00158379970045 |

Headline totals and flux integrals use the same integer-second trapezoidal quadrature as Task 37. Both perturbations also exceed the original WT resting secretion thresholds; Task 38's decision gates are the specified physiology/conservation checks.

| Time (s) | WT q_out (pL/s) | 5% q_out (pL/s) | 0% q_out (pL/s) | 5% flow deficit (%) | 0% flow deficit (%) |
| --- | --- | --- | --- | --- | --- |
| 0 | 0.00107570739 | 0.00107570739 | 0.00107570739 | 0.0000 | 0.0000 |
| 60 | 0.00171328884 | 0.0016705883 | 0.00166725406 | 2.4923 | 2.6869 |
| 120 | 0.00169082188 | 0.00163070743 | 0.0016247764 | 3.5553 | 3.9061 |
| 180 | 0.00166686496 | 0.00160951507 | 0.00160264402 | 3.4406 | 3.8528 |
| 240 | 0.00164485084 | 0.00159801174 | 0.00159121814 | 2.8476 | 3.2606 |
| 300 | 0.00162588646 | 0.00159171798 | 0.00158556467 | 2.1015 | 2.4800 |
| 360 | 0.00160992395 | 0.00158827114 | 0.00158303847 | 1.3450 | 1.6700 |
| 420 | 0.0015965924 | 0.00158641911 | 0.00158221731 | 0.6372 | 0.9004 |
| 480 | 0.00158547436 | 0.0015854854 | 0.0015823257 | -0.0007 | 0.1986 |
| 540 | 0.00157619107 | 0.00158509482 | 0.0015829339 | -0.5649 | -0.4278 |
| 600 | 0.00156842215 | 0.0015850324 | 0.0015837997 | -1.0590 | -0.9804 |

## Intracellular time courses

All concentrations are mM; cell volume is pL. The CSV files additionally include lumen osmolarity, voltages, NBC activation/current, NHE1 multiplier, and every requested instantaneous transporter flux.

### AE4=0.05

| Time (s) | Na_i | K_i | Cl_i | pH_i | HCO3_i | TIC_i | Volume |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 11.63613 | 116.7632 | 60.30497 | 6.91 | 4.909375 | 5.671746 | 1.448426 |
| 60 | 15.99241 | 113.3805 | 57.20142 | 6.998899 | 9.989549 | 11.25534 | 1.496538 |
| 120 | 17.37999 | 112.4655 | 55.57127 | 7.051323 | 12.17586 | 13.54473 | 1.51321 |
| 180 | 17.94568 | 112.1785 | 54.63449 | 7.088089 | 13.18353 | 14.54663 | 1.518349 |
| 240 | 18.20402 | 112.1163 | 54.09689 | 7.115894 | 13.6731 | 15.00023 | 1.520035 |
| 300 | 18.3283 | 112.1427 | 53.78389 | 7.137661 | 13.92109 | 15.20715 | 1.520705 |
| 360 | 18.38795 | 112.2039 | 53.59826 | 7.155028 | 14.05301 | 15.30113 | 1.521143 |
| 420 | 18.41411 | 112.277 | 53.48646 | 7.169048 | 14.12814 | 15.34373 | 1.521593 |
| 480 | 18.42197 | 112.3518 | 53.41872 | 7.18046 | 14.175 | 15.36354 | 1.522112 |
| 540 | 18.41946 | 112.4237 | 53.37819 | 7.189808 | 14.20741 | 15.37377 | 1.522698 |
| 600 | 18.41098 | 112.4909 | 53.35498 | 7.197501 | 14.23208 | 15.38034 | 1.523329 |

### AE4=0.00

| Time (s) | Na_i | K_i | Cl_i | pH_i | HCO3_i | TIC_i | Volume |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 11.63613 | 116.7632 | 60.30497 | 6.91 | 4.909375 | 5.671746 | 1.448426 |
| 60 | 15.98565 | 113.4085 | 57.07041 | 7.000904 | 10.17242 | 11.4555 | 1.497839 |
| 120 | 17.34937 | 112.5555 | 55.30911 | 7.056072 | 12.58983 | 13.98999 | 1.516798 |
| 180 | 17.8833 | 112.343 | 54.30675 | 7.095454 | 13.77824 | 15.17917 | 1.524428 |
| 240 | 18.10547 | 112.359 | 53.75045 | 7.125672 | 14.40116 | 15.76825 | 1.528522 |
| 300 | 18.19111 | 112.463 | 53.44626 | 7.149642 | 14.74675 | 16.07257 | 1.531398 |
| 360 | 18.21159 | 112.5984 | 53.28409 | 7.169001 | 14.9508 | 16.23731 | 1.533806 |
| 420 | 18.19971 | 112.7406 | 53.20256 | 7.184814 | 15.08035 | 16.33242 | 1.535999 |
| 480 | 18.17182 | 112.8786 | 53.16732 | 7.197835 | 15.1694 | 16.39236 | 1.538058 |
| 540 | 18.13648 | 113.0077 | 53.15893 | 7.208621 | 15.23549 | 16.43425 | 1.540006 |
| 600 | 18.09832 | 113.126 | 53.16615 | 7.217602 | 15.2878 | 16.46657 | 1.541844 |

## Integrated fluxes and chloride partition

Window: 60-600 s. Positive chloride loading integrates max(J_Cl,0) for each transporter separately. Signed AE2 is positive for chloride entry to the cell. NBC cycles and NHE1 are signed inward fluxes; NBC bicarbonate-equivalent influx is twice the cycle integral. CaCC export is cell-to-lumen; paracellular chloride return is lumen-to-bath.

| Quantity | Unit | Frozen WT | AE4=0.05 | AE4=0.00 |
| --- | --- | --- | --- | --- |
| positive_nkcc1_cl_loading | fmol | 155.723019 | 193.2992595 | 197.6037211 |
| positive_ae4_cl_loading | fmol | 39.11392936 | 4.183594922 | 0 |
| positive_ae2_cl_loading | fmol | 0 | 0.3773563236 | 0.4560458999 |
| signed_ae2_cl_flux | fmol | -0.4649947574 | 0.3760292927 | 0.4553929111 |
| nbc_cycles | fmol | 36.97651379 | 8.477231656 | 5.352988627 |
| nbc_bicarbonate_equivalent_influx | fmol | 73.95302758 | 16.95446331 | 10.70597725 |
| nhe1_flux | fmol | 3.904657014 | 1.91085779 | 1.773666246 |
| cacc_cl_export | fmol | 205.9230127 | 202.1857717 | 201.5674776 |
| paracellular_cl_return | fmol | 76.31423568 | 74.84928993 | 74.63690175 |
| water_outflow | pL | 0.8781700024 | 0.8639005293 | 0.8613107814 |

| Positive chloride pool | Frozen WT | AE4=0.05 | AE4=0.00 |
| --- | --- | --- | --- |
| NKCC1 share (%) | 79.924789 | 97.694862 | 99.769743 |
| AE4 share (%) | 20.075211 | 2.114420 | 0.000000 |
| AE2 share (%) | 0.000000 | 0.190719 | 0.230257 |

## Cumulative secretion deficit versus time

Task 37 retained its 0-600 and 60-600 s integrals and q_out at 60-second intervals, but did not save its integer-second trajectory or intermediate cumulative values. Its exact stored-quadrature 0-60 s total is therefore 0.108825053189 pL by subtraction. Values at 120-540 s below use trapezoidal integration of the frozen 60-second q_out samples anchored at that 60 s total. They are explicitly estimates. At 600 s the original frozen total is used directly. Perturbation cumulative values use the original one-second grid throughout. The uncorrected coarse WT quadrature would exceed the stored 600 s total by 7.76826767e-05 pL (0.007871%). This is a resolution diagnostic, not a bound on intermediate errors. No WT trajectory was rerun.

| Time (s) | WT cumulative (pL) | WT basis | 5% cumulative (pL) | 0% cumulative (pL) | WT minus 5% (pL) | WT minus 0% (pL) |
| --- | --- | --- | --- | --- | --- | --- |
| 0 | 0.000000000 | stored/derived | 0.000000000 | 0.000000000 | 0.000000000 | 0.000000000 |
| 60 | 0.108825053 | stored/derived | 0.107613859 | 0.107528028 | 0.001211195 | 0.001297025 |
| 120 | 0.210948375 | estimated | 0.206512411 | 0.206138775 | 0.004435964 | 0.004809599 |
| 180 | 0.311678980 | estimated | 0.303654495 | 0.302890232 | 0.008024485 | 0.008788748 |
| 240 | 0.411030454 | estimated | 0.399845629 | 0.398667665 | 0.011184825 | 0.012362789 |
| 300 | 0.509152573 | estimated | 0.495518751 | 0.493950400 | 0.013633822 | 0.015202173 |
| 360 | 0.606226885 | estimated | 0.590908030 | 0.588997204 | 0.015318855 | 0.017229681 |
| 420 | 0.702422376 | estimated | 0.686142824 | 0.683948711 | 0.016279551 | 0.018473664 |
| 480 | 0.797884379 | estimated | 0.781296501 | 0.778881655 | 0.016587877 | 0.019002724 |
| 540 | 0.892734342 | estimated | 0.876411835 | 0.873837673 | 0.016322507 | 0.018896668 |
| 600 | 0.986995056 | stored/derived | 0.971514388 | 0.968838810 | 0.015480668 | 0.018156246 |

## Interpretation

**5% versus null:** the 5% trajectory secretes 0.00267557826 pL more, a difference of 0.271083 percentage points of WT. Null/5% secretion is 0.997245972. These are distinct numerical predictions with a small secretion difference. No materiality or equivalence threshold was declared, so this is not a formal equivalence claim. Their largest chloride difference on the reported grid is 0.346443 mM, and terminal pH differs by 0.020100.

**Timing:** the null's cumulative shortfall is 0.00129703 pL at 60 s; it grows to approximately 0.00480960 pL at 120 s and 0.00878875 pL at 180 s. The largest tabulated cumulative shortfall is around 480 s, about 0.0190027 pL, then it narrows to the stored-reference deficit of 0.01815625 pL at 600 s. Thus an initially small cumulative deficit develops over the first minutes, but instantaneous flow does not progressively diverge after 2-3 minutes: its largest reported deficit is 3.9061% at 120 s, 3.8528% at 180 s, then shrinks. Null flow is 0.9804% above WT at 600 s. The 5% case has the same qualitative time pattern and is 1.0590% above WT at 600 s. Cumulative loss remains positive because earlier deficits outweigh the late recovery.

**Chloride and pH:** loss of AE4 lowers intracellular chloride at every reported post-onset time. At 600 s, Cl_i is 54.07952 mM in WT, 53.35498 mM at 5%, and 53.16615 mM in null. pH stays within 6.6-7.3 throughout, but rises to 7.19750 and 7.21760 versus WT 7.03483. It is controlled within the declared gate while showing substantial relative alkalinization; the finite trajectory does not establish a new resting state or long-time stability.

**NKCC1 and AE2:** positive NKCC1 chloride loading rises by 24.1302% at 5% AE4 and 26.8944% in null despite identical capacity, expression and stimulation parameters. It supplies 97.6949% and 99.7697% of the respective positive pools. AE2 reverses from net WT chloride export (-0.464995 fmol) to small net import (+0.376029 and +0.455393 fmol). It is positive by the 120 s readout in both perturbations but contributes only 0.1907% and 0.2303% of the positive pools. The dominant compensatory loading response is NKCC1.

**NBC and NHE1:** NBC recruitment remains one after onset, and the NHE1 multiplier remains 2.3; changes in realized flux arise within the frozen coupled model as concentrations and membrane voltage evolve. NBC cycle influx falls by 77.0740% and 85.5233% over 60-600 s, and NHE1 flux falls by 51.0621% and 54.5756%. Both remain inward on monitored stimulated samples. Bicarbonate accumulates to 14.23208 and 15.28780 mM versus WT 7.27856 mM at 600 s. This behavior is consistent with reduced AE4 bicarbonate export and altered acid-base/voltage driving forces; it is an interpretation of the coupled trajectory, not an additional intervention.

**Held-out validation:** the predicted null cumulative deficit is 1.83955%, versus the supplied experimental context of approximately 35%. The frozen model therefore misses the magnitude of that held-out phenotype. No parameter was tuned, mechanism added, or rescue attempted.

## Gate evidence

Task 37's unmodified `diagnose` function was applied to each perturbation evaluation. Each trajectory had 333 accepted solver steps, 2336 RHS evaluations and 934 checked states, including every accepted endpoint and every integer-second dense-output point. All tracked core amounts/volumes and derived physiological readouts passed the finite/positive and specified range gates. These are monitored numerical checks, not a formal enclosure between samples.

| Readout | AE4=0.05 min to max | AE4=0.00 min to max | Gate |
| --- | --- | --- | --- |
| na_i_mM | 11.6361257 to 18.4221174 | 11.6361257 to 18.2115987 | <40 |
| k_i_mM | 112.11617 to 116.763209 | 112.332603 to 116.763209 | 50-200 |
| cl_i_mM | 53.3549827 to 60.3049669 | 53.1589057 to 60.3049669 | 30-80 |
| ph_i | 6.91000022 to 7.19750138 | 6.91000022 to 7.2176015 | 6.6-7.3 |
| hco3_i_mM | 4.90937455 to 14.2320805 | 4.90937455 to 15.2878014 | 0<HCO3<100 |
| tic_i_mM | 5.67174623 to 15.3803377 | 5.67174623 to 16.4665679 | >0 |
| volume_i_pL | 1.44835391 to 1.52332888 | 1.44835416 to 1.54184372 | 0<V<3 |

| Conservation diagnostic | Tolerance | 5% max absolute | 0% max absolute |
| --- | --- | --- | --- |
| ae4_charge_fmol_s | 1e-10 | 3.46944695e-18 | 0 |
| apical_current_A | 1e-20 | 2.52031593e-25 | 2.52031593e-25 |
| basolateral_current_A | 1e-20 | 2.48735795e-23 | 2.46344726e-23 |
| buffer_site_accounting_fmol_s | 1e-10 | 0 | 0 |
| carbon_accounting_fmol_s | 1e-10 | 5.55111512e-17 | 5.55111512e-17 |
| cell_bulk_charge_rate_fmol_s | 1e-10 | 2.58598698e-13 | 2.55281907e-13 |
| cell_speciation_alkalinity_mM | 1e-09 | 5.92592642e-12 | 5.88684657e-12 |
| homeostasis_charge_fmol_s | 1e-10 | 6.9388939e-17 | 6.59194921e-17 |
| lumen_bulk_charge_rate_minus_outflow_fmol_s | 1e-10 | 2.52402266e-15 | 2.64371858e-15 |
| lumen_speciation_alkalinity_mM | 1e-09 | 5.70654635e-13 | 5.99520433e-13 |
| water_volume_accounting_pL_s | 1e-12 | 2.16840434e-19 | 4.33680869e-19 |

The largest tolerance-normalized residual was 0.00592593 at 5% and 0.00588685 in null (cell speciation alkalinity). The physical NBC source fields were checked as charge=-J_NBC and carbon=2*J_NBC, rather than treated as zero residuals. No scientific gate failed. The first case's PASS was verified before the null integration started.

## Budget, provenance and files

Exactly two intended and completed perturbation integrations; zero WT reruns, numerical retries, stationary solves, optimisation calls, sweeps or additional scientific parameter evaluations. One scientific agent/worker and one thread per BLAS pool. The numerical runner used 2.872439 s of the 600 s limit. The persisted budget prohibits further scientific execution.

[5% readouts](../../results/38_ae4_perturbation_validation/ae4_5pct_timeseries.csv), [null readouts](../../results/38_ae4_perturbation_validation/ae4_null_timeseries.csv), [flow comparison](../../results/38_ae4_perturbation_validation/flow_comparison.csv), [intracellular comparison](../../results/38_ae4_perturbation_validation/intracellular_comparison.csv), [cumulative comparison](../../results/38_ae4_perturbation_validation/cumulative_secretion_comparison.csv), [integrated flux comparison](../../results/38_ae4_perturbation_validation/integrated_flux_comparison.csv), [summary](../../results/38_ae4_perturbation_validation/comparison_summary.json), [5% verification](../../results/38_ae4_perturbation_validation/ae4_5pct_verification.json), [null verification](../../results/38_ae4_perturbation_validation/ae4_null_verification.json), [frozen input manifest](../../results/38_ae4_perturbation_validation/frozen_inputs.json), [budget](../../results/38_ae4_perturbation_validation/budget.json).

The minimal runner is [run_perturbation.py](run_perturbation.py); [summarize_results.py](summarize_results.py) produces the comparisons from saved data only. Publication uses the connected GitHub integration on the specified branch. No merge to main or further scientific task is performed.
