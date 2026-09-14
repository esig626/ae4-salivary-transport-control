Task 40 — fixed 50:50 AE4 cation routing

**TASK40_EQUAL_ROUTING_VALIDATION_COMPLETE**


The single WT resting-state solve and all three production trajectories passed. Equal routing reduced AE4-loss NKCC1 compensation, but did not preserve higher intracellular Na concentration at 600 s or produce a large secretion deficit. This is a completed structural test with a partial compensation reduction; the proposed sustained-Na outcome is unsupported in this frozen architecture.

| Measure | AE4 = 0.05 | AE4 = 0.00 |
| --- | --- | --- |
| 600 s Na_i minus WT (mM) | -0.137993 | -0.092517 |
| 60–600 s NKCC1 compensation | +20.81935% | +23.16345% |
| Task 39 compensation | +29.65434% | +33.10131% |
| Compensation change (percentage points) | -8.83499 | -9.93786 |
| 600 s Cl_i minus WT (mM) | -1.894288 | -2.195703 |
| Cumulative secretion / Task 40 WT | 0.9658242239 | 0.9614251370 |
| Cumulative secretion deficit | 3.417578% | 3.857486% |
| Task 39 cumulative deficit | 0.482162% | 0.627320% |

The binding inputs were read in order: [AGENTS.md](../../AGENTS.md), [source contract](source_contract.md), then [Task 40 prompt](../../prompts/40_ae4_equal_cation_routing.md). The prepared head was `a6a58ebabcd2aedb1702b52dfdd86359c60cec32`, on `codex/task-40-ae4-equal-cation-routing`. The scientific parent was merged Task 39, `afd101448763439f369f2682c467f069ea18a442`. The preparation changed only these three instruction files relative to that parent. The downloaded scientific inputs were verified against the prepared Git tree; the manifest contains 150 unchanged files. The old donor-weighted evaluator remains selectable.

The new named evaluator calls the same inherited QSS evaluator for J4 and then returns the fixed source vector. There is no routing-fraction parameter or donor-dependent choice. For each evaluated state, expression and regulation, the J4 law is identical to Task 39. J4 can differ between evolved trajectories because their states differ; it was never clamped to its reference value.

`(S_Na, S_K, S_Cl, S_TIC, S_TA) = (-J4/2, -J4/2, J4, -2J4, -2J4)`. Thus `S_Na + S_K = -J4`, and `S_Na + S_K - S_Cl - S_TA = -J4/2 - J4/2 - J4 + 2J4 = 0`. Negative J4 reverses every source automatically. The focused tests verify actual forward and reversed environments, zero cycle, exact inherited J4 and anion sources, expression scaling, and unchanged nonrouting diagnostics. This is the stipulated ensemble source allocation, not a separately derived microscopic equal-routing carrier.

Independent 50-digit decimal arithmetic and full-model evaluation reproduced the preflight:

| Reference quantity | Reproduced value |
| --- | --- |
| Na_i / (Na_i + K_i) | 0.09062450161183984 |
| Inherited J4 (fmol/s) | 0.004937067441546469 |
| Expected delta S_Na (fmol/s) | −0.0020211144444590447 |
| Expected delta S_K (fmol/s) | +0.0020211144444590447 |
| Actual full-RHS Na/K shifts (fmol/s) | −0.0020211144444590434 / +0.0020211144444590434 |
| Direct Cl/TIC/TA shifts | Exactly zero |
| Charge perturbation (fmol/s) | −4.34e−19 (round-off) |

All five focused tests passed after one documented literal correction to a floating-point scaling assertion. The original error was 2.0817e−17 fmol/s against a 2e−17 allowance; a scale-aware 32-machine-epsilon comparison replaced that assertion. Exact identity to inherited J4 remained required. No evaluator, parameter, production tolerance or scientific gate was changed. Both logs and [the correction record](source_correction.md) are retained. [Source verification](../../results/40_ae4_equal_cation_routing/source_verification.json) records the algebra and parameter hashes.

Exact Task 31 REST nesting was not expected and was not imposed. The WT rest used only the accepted Task 31 R09 state. A single square ten-equation `root/hybr` solve used the inherited charge manifold, residual scales, coordinate bounds and acceptance tolerances. A logit coordinate transform enforced those numerical bounds. No least-squares fit, calibration, multistart, continuation, other root family, or genotype rest was used. There were 86 root residual evaluations and 12 final audit RHS evaluations, with no retry. The sole initial residual used the literal accepted state; its coordinate reconstruction differs by less than 1e−12.

The accepted core state and its shift are below. Amounts are fmol and volumes pL. The sole regulatory coordinate is zero at rest. The frozen full-state SHA-256 is `a0b96ef16e5967de4dfe2bdabb882f2a28d478e211f920bc19cd735a5a143acb`.

| Core state | Task 31 | Task 40 | Change |
| --- | --- | --- | --- |
| na_i_fmol | 16.8540680421 | 16.7203587327 | -0.13370931 |
| k_i_fmol | 169.122877953 | 170.002008522 | +0.87913057 |
| cl_i_fmol | 87.3472870434 | 88.0318019636 | +0.68451492 |
| tic_i_fmol | 8.2151051813 | 8.2635118084 | +0.048406627 |
| alk_i_fmol | 20.5678967454 | 20.6288030851 | +0.06090634 |
| volume_i_pL | 1.44842608322 | 1.45333439259 | +0.0049083094 |
| na_l_fmol | 13.1391270128 | 13.1277686682 | -0.011358345 |
| k_l_fmol | 2.70611880568 | 2.71897260899 | +0.012853803 |
| cl_l_fmol | 15.3762475046 | 15.3786878545 | +0.0024403499 |
| tic_l_fmol | 0.541174028869 | 0.540262316406 | -0.00091171246 |
| alk_l_fmol | 0.468998313825 | 0.468053422678 | -0.00094489115 |
| volume_l_pL | 0.105378536927 | 0.105388407354 | +9.8704269e-06 |

| Rest observable | Task 31 | Task 40 | Change |
| --- | --- | --- | --- |
| na_i_mM | 11.63612575 | 11.50482561 | -0.1313001 |
| k_i_mM | 116.7632093 | 116.9737738 | +0.2105645 |
| cl_i_mM | 60.30496693 | 60.57229665 | +0.2673297 |
| ph_i | 6.910000223 | 6.910983068 | +0.0009828443 |
| hco3_i_mM | 4.909374552 | 4.923112541 | +0.01373799 |
| tic_i_mM | 5.671746233 | 5.685898476 | +0.01415224 |
| volume_i_pL | 1.448426083 | 1.453334393 | +0.004908309 |
| q_out_pL_s | 0.001075707385 | 0.001077681471 | +1.974085e-06 |
| v_apical_V | -0.03153129513 | -0.03142872408 | +0.0001025711 |
| v_basolateral_V | -0.06710109366 | -0.06705742727 | +4.366639e-05 |
| v_transepithelial_V | -0.03556979853 | -0.0356287032 | -5.890467e-05 |
| osm_lumen_mOsm | 301.4149587 | 301.4154236 | +0.0004649646 |

The complete rest flux ledger, water fluxes, voltages, acid/base speciation, and current/charge/carbon residuals are in [wt_rest.json](../../results/40_ae4_equal_cation_routing/wt_rest.json); every ledger shift is in [wt_rest_flux_comparison.csv](../../results/40_ae4_equal_cation_routing/wt_rest_flux_comparison.csv). Rest J4 is 0.004933930166208064 fmol/s; each cation source is −0.002466965083104032 fmol/s. Rest NBC cycle and current are zero. The small change from the old J4 is a state consequence, not a rate-law change.

Rest passed finite positive tracked amounts/volumes, Na_i < 40 mM, K_i 50–200 mM, Cl_i 30–80 mM, pH_i 6.6–7.3, positive TIC/HCO3 with HCO3 < 100 mM, volume < 3 pL, and all inherited production conservation tolerances. The maximum scaled stationary residual was 1.464e-11 against 1e−7; maximum amount RHS 4.996e-16 fmol/s; maximum volume RHS 1.464e-15 pL/s. The independent Jacobian had rank 10, no coordinate boundary was hit, the omitted charge rows and regulatory RHS passed, and the resting RHS was time invariant.

The rest was frozen before the three trajectories. WT, 5%, and null have bitwise identical initial vectors and identical scientific parameter hashes. Only `genotype.ae4_expression` differs. All used the inherited 600 s standard stimulus, Radau rtol 1e−7, amount atol 1e−10 fmol, volume atol 1e−12 pL, regulatory atol 1e−10, and maximum step 2 s. The inherited onset convention keeps the exact t=0 resting-stimulus readout and begins integration at 1e−6 s from the same vector. Every accepted endpoint and integer-second dense sample was checked: 934 states per case. All gates passed, including WT secretion activation; none was rescued.

All requested observables at 0, 60, …, 600 s are in [WT](../../results/40_ae4_equal_cation_routing/wt_timeseries.csv), [5%](../../results/40_ae4_equal_cation_routing/ae4_5pct_timeseries.csv), and [null](../../results/40_ae4_equal_cation_routing/ae4_null_timeseries.csv). They include q_out, cumulative secretion, intracellular composition/volume, all requested transporter fluxes, separate AE4 Na/K sources and export branches, NBC current, voltages, lumen osmolarity, and Palk X. Integrals use the inherited one-second trapezoidal grid, with onset samples; they were not recomputed from the coarser report rows.

| 600 s observable | WT | 5% | Null |
| --- | --- | --- | --- |
| na_i_mM | 17.5542362 | 17.4162429 | 17.4617192 |
| k_i_mM | 110.925054 | 113.301933 | 113.531289 |
| cl_i_mM | 53.9799192 | 52.0856307 | 51.7842163 |
| ph_i | 7.04783303 | 7.20427148 | 7.22291962 |
| hco3_i_mM | 7.50255519 | 14.522195 | 15.5201338 |
| tic_i_mM | 8.35276449 | 15.6760934 | 16.7025643 |
| volume_i_pL | 1.40154084 | 1.50616891 | 1.52097397 |
| nkcc1_x_i_mM4 | 5673826.46 | 5353374.76 | 5316151.94 |
| q_out_pL_s | 0.00160212486 | 0.00156355656 | 0.00155681142 |
| cumulative_outflow_0_t_pL | 0.992524544 | 0.958604248 | 0.954238046 |
| v_apical_V | -0.0292891667 | -0.0301353779 | -0.0302743831 |
| v_basolateral_V | -0.0792037827 | -0.0786748046 | -0.078602441 |
| v_transepithelial_V | -0.049914616 | -0.0485394267 | -0.0483280579 |
| osm_lumen_mOsm | 301.538171 | 301.530164 | 301.528805 |

Na_i is higher after AE4 loss at the reported 60–480 s samples for 5%, and 60–540 s samples for null. The largest tabulated differences occur at 180 s: +0.9186 and +1.0317 mM. The 5% difference changes sign between 480 and 540 s; null between 540 and 600 s. At 600 s, both are below WT. These intervals use the saved 60 s readouts; no finer crossing time is claimed.

Positive NKCC1 chloride loading integrated over 60–600 s rises from 156.8722 fmol in WT to 189.5319 / 193.2092 fmol. Compensation is smaller than Task 39 by 8.8350 / 9.9379 percentage points, a roughly 30% relative reduction in the compensation percentage. It neither disappears nor reverses. The NKCC1 increase replaces 88.86% / 88.51% of lost AE4 chloride loading. Endpoint Palk X is lower by 5.648% / 6.304%; this is consistent with increased loading under the unchanged decreasing Palk concentration-response law, despite the transient Na elevation. These comparisons use each task's own frozen WT; Task 39 trajectories were not rerun.

Cl_i remains lower at every positive-time reported checkpoint. At 600 s the deficits are 1.8943 / 2.1957 mM (3.509% / 4.068%). They are smaller than the 180 s differences of 3.1704 / 3.5016 mM, so compensation partly closes the concentration gap without removing it within 600 s. This is a modest persistent concentration deficit, not evidence of a large depletion of cellular chloride content. Larger cell volumes matter: final Cl amounts are 75.6551 / 78.4498 / 78.7624 fmol (WT/5%/null), while Na amounts are 24.6030 / 26.2318 / 26.5588 fmol. Concentration and amount therefore have different comparisons.

Integrated fluxes over 60–600 s are shown below. Percentages for positive AE2 loading are undefined relative to its zero WT value; signed AE2 changes are described by direction. NBC and NHE1 stimulation settings remain unchanged despite their flux responses.

| Integrated quantity | Unit | WT | 5% | Null | 5% change | Null change |
| --- | --- | --- | --- | --- | --- | --- |
| positive_nkcc1_cl_loading | fmol | 156.87218 | 189.53195 | 193.20918 | +20.819% | +23.163% |
| positive_ae4_cl_loading | fmol | 41.055502 | 4.3016803 | 0 | -89.522% | -100.000% |
| positive_ae2_cl_loading | fmol | 0 | 0.43509709 | 0.51147964 | undefined | undefined |
| signed_ae2_cl_flux | fmol | -0.41841964 | 0.43457519 | 0.51130754 | sign reversal | sign reversal |
| nkcc1_cycles | fmol | 78.436089 | 94.765973 | 96.604592 | +20.819% | +23.163% |
| ae4_signed_cl_flux | fmol | 41.055502 | 4.3016803 | 0 | -89.522% | -100.000% |
| ae4_na_cell_source | fmol | -20.527751 | -2.1508401 | 0 | -89.522% | -100.000% |
| ae4_k_cell_source | fmol | -20.527751 | -2.1508401 | 0 | -89.522% | -100.000% |
| ae4_tic_cell_source | fmol | -82.111004 | -8.6033605 | 0 | -89.522% | -100.000% |
| ae4_ta_cell_source | fmol | -82.111004 | -8.6033605 | 0 | -89.522% | -100.000% |
| nbc_cycles | fmol | 39.110273 | 8.7281381 | 5.3939272 | -77.683% | -86.208% |
| nbc_bicarbonate_equivalent_influx | fmol | 78.220546 | 17.456276 | 10.787854 | -77.683% | -86.208% |
| nbc_current_charge | C | 3.7735677e-09 | 8.4213731e-10 | 5.2043486e-10 | -77.683% | -86.208% |
| nhe1_flux | fmol | 4.0565646 | 1.9026821 | 1.7621799 | -53.096% | -56.560% |
| pump_cycles | fmol | 32.721489 | 33.241195 | 33.329779 | +1.588% | +1.859% |
| cacc_cl_export | fmol | 207.32726 | 199.6824 | 198.69931 | -3.687% | -4.162% |
| paracellular_cl_return | fmol | 76.823455 | 74.261359 | 73.942305 | -3.335% | -3.750% |
| water_outflow | pL | 0.88377136 | 0.85124324 | 0.84697528 | -3.681% | -4.164% |

AE2 switches from signed export in WT to a small positive loading route after AE4 loss. Integrated positive loading is 0 / 0.43510 / 0.51148 fmol, only 0 / 0.22397% / 0.26403% of the positive chloride-loading pool. It remains negligible in that quantitative sense, rather than exactly absent. WT's NKCC1 positive-loading share is 79.257%, outside the inherited 65–75% contextual band; that band is explicitly validation context and was not an acceptance gate.

Loss of AE4 base export is accompanied by higher TIC, bicarbonate and pH. Integrated NBC cycles, bicarbonate influx and current fall by 77.683% / 86.208%; NHE1 influx falls by 53.096% / 56.560%. Integrated pump cycles rise modestly, by 1.588% / 1.859%, while their endpoint rates are slightly below WT (0.062119 / 0.062196 versus 0.062352 fmol/s). Thus a transient rise in Na concentration does not imply a sustained rise: reduced Na export through AE4 coexists with markedly reduced NBC/NHE1 Na entry, pump feedback, and volume change. This interpretation follows the recorded coupled balances and is not an additional intervention.

CaCC chloride export falls by 3.687% / 4.162%, and paracellular chloride return falls by 3.335% / 3.750%. At 600 s, NBC cycle flux is 0.071918 / 0.009343 / 0.001679 fmol/s and current is 6.939 / 0.901 / 0.162 pA (WT/5%/null); NHE1 is 0.006245 / 0.001836 / 0.001567 fmol/s. CaCC export is 0.376856 / 0.367142 / 0.365543 fmol/s, and paracellular return is 0.140254 / 0.136569 / 0.136006 fmol/s. The complete endpoint fluxes are in the trajectory CSVs.

Cumulative secretion deficits are `100 × (1 − Q_case(0,t)/Q_WT(0,t))`; absolute loss is `Q_WT − Q_case`. At t=0, both cumulatives are zero and their ratio/percentage deficit is undefined.

| Time (s) | WT cumulative (pL) | 5% deficit | Null deficit | 5% lost (pL) | Null lost (pL) |
| --- | --- | --- | --- | --- | --- |
| 0 | 0.00000000 | undefined | undefined | 0.00000000 | 0.00000000 |
| 60 | 0.10875318 | 1.28012% | 1.37046% | 0.00139217 | 0.00149041 |
| 120 | 0.21019117 | 2.49034% | 2.69605% | 0.00523447 | 0.00566685 |
| 180 | 0.31037630 | 3.16154% | 3.45199% | 0.00981267 | 0.01071416 |
| 240 | 0.40955916 | 3.50123% | 3.85019% | 0.01433960 | 0.01576882 |
| 300 | 0.50794464 | 3.64954% | 4.03777% | 0.01853767 | 0.02050964 |
| 360 | 0.60569525 | 3.68712% | 4.10071% | 0.02233273 | 0.02483783 |
| 420 | 0.70293841 | 3.66040% | 4.08937% | 0.02573037 | 0.02874572 |
| 480 | 0.79977429 | 3.59670% | 4.03389% | 0.02876548 | 0.03226202 |
| 540 | 0.89628205 | 3.51242% | 3.95276% | 0.03148122 | 0.03542790 |
| 600 | 0.99252454 | 3.41758% | 3.85749% | 0.03392030 | 0.03828650 |

The relative cumulative deficits grow early, reach their largest reported values at 360 s (3.6871% / 4.1007%), and then decline modestly to 3.4176% / 3.8575%. They persist through 600 s rather than collapsing toward the small Task 39 endpoints. Absolute cumulative lost volume continues to increase across the reported checkpoints. Final flow rates remain 2.4073% / 2.8283% below WT, so the late percentage decline is not recovery of the lost absolute secretion. No infinite-time or steady-state conclusion is claimed.

Final cumulative secretion is 0.9925245441 pL for WT, 0.9586042475 pL at 5%, and 0.9542380458 pL in null. The ratios are 0.9658242239 and 0.9614251370, with losses of 0.0339202966 and 0.0382864983 pL. Task 39's deficits were only 0.48216% / 0.62732%, but its WT cumulative was 0.9586989000 pL. The larger Task 40 deficit partly reflects its changed WT reference; the reported relative effects do not imply the same change in absolute perturbation secretion between tasks.

The held-out approximately 35% experimental null deficit is much larger than the predicted 3.8575% cumulative deficit, a descriptive gap of 31.1425 percentage points. The endpoint flow deficit is only 2.8283%. The experimental number was never a fit target, gate, root selector, or trigger for another mechanism; exact matching of experimental protocol and observation model is not asserted.

Residual maxima and unchanged tolerances for rest and all trajectories are in [conservation_summary.csv](../../results/40_ae4_equal_cation_routing/conservation_summary.csv). Across the trajectories, the largest cell charge-rate residual is 2.572e−13 fmol/s against 1e−10; carbon residual 5.551e−17 fmol/s against 1e−10; water residual 4.337e−19 pL/s against 1e−12; basolateral current residual 2.478e−23 A against 1e−20; and speciation residual 6.128e−12 mM against 1e−9. AE4 charge residual is exactly zero. These are sampled production checks, not a continuous-time mathematical proof.

Execution used one scientific agent/worker, one BLAS thread, one WT stationary solve, exactly three production integrations, zero numerical retries, zero calibration/optimisation/sweeps, and no genotype-specific rest. Recorded numerical execution was 6.943 seconds, below 900 seconds. The one source-test assertion correction is recorded separately from numerical retries. No split other than 50:50 was tested; no scientific parameter or second mechanism was changed. Scientific execution is locked complete. Compact UTF-8 source, tests and outputs are published only on the designated branch through the connected GitHub integration; main is not modified or merged.

Reproduction order on a fresh checkout without Task 40 outputs: `python analysis/40_ae4_equal_cation_routing/run_preflight.py`, then `run_rest.py`, followed by `run_case.py wt`, `run_case.py ae4_5pct`, and `run_case.py ae4_null` in the same directory. A saved completion ledger prevents reruns. The report-only scripts `summarize_results.py` and `write_report.py` read saved outputs and do not import or execute the model.
