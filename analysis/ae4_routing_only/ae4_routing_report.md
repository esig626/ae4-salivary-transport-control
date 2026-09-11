**No. Removing the Na/K countercycle does not change AE4 loss from increasing secretion to decreasing secretion in this controlled rerun.** At 5% AE4, total secretion over 600 s remains above WT in all 30 root and calcium pairs. The increase is smaller: 3.49% to 8.34%, compared with 8.23% to 21.66% before the routing change.

The calculation uses all ten original modern model parameter sets and their saved WT and 5% AE4 starting states. The calcium inputs, R1 regulation, CCh plus IPR protocol, NKCC1 normalisation, 602 point time grid and production Radau settings are inherited from Tasks 14 and 14B. No parameters or capacities were fitted or recalibrated.

| Calcium (µM) | Original secretion change at 5% AE4 | Change after cation routing intervention |
| ---: | ---: | ---: |
| 0.10 | +8.23% to +16.32% | +3.49% to +6.18% |
| 0.25 | +14.54% to +20.27% | +4.96% to +7.83% |
| 0.50 | +16.70% to +21.66% | +5.50% to +8.34% |

Every percentage compares 5% AE4 with WT under the same routing law. Ranges describe the ten retained parameter sets; they are not confidence intervals.

Only the Na and K sources were changed. The inherited evaluator still supplies the same net AE4 cycle J, equal to its net chloride source. The total cation source is minus J. For chloride entry, intracellular Na and K fractions divide this source; for chloride exit, extracellular Na and K fractions divide it. Equal concentration weighting introduces no fitted selectivity. Both components therefore reverse together and vanish when J vanishes. The existing 1:1:2 chloride:cation:bicarbonate stoichiometry and bicarbonate, carbon and alkalinity sources remain unchanged.

Preserving the chloride law means preserving its value at the same state. The chloride trajectory can subsequently change through the resulting changes in cellular Na and K. Pumps, acid and base chemistry, chloride allocation rules, transport capacities and all other equations retain their original definitions and parameters.

**Initial state qualification.** The saved starting states were retained exactly and were not solved again. After changing cation routing, they are no longer exact resting states: the largest initial amount derivative ranges from 0.0670 to 0.1908 fmol/s for WT, and from 0.00412 to 0.01152 fmol/s for 5% AE4. The result therefore answers the requested isolated 600 s rerun from the saved states. It does not establish a comparison from newly equilibrated resting states. No physiological retuning was used to suppress this drift.

Most of the smaller genotype contrast comes from higher WT secretion after the routing change. WT totals increase by 4.66% to 14.58% compared with their original values, while 5% AE4 totals increase by 0.079% to 0.556%. The wrong direction of the AE4 loss effect persists.

All 120 production trajectories passed the existing solver, positivity, nonnegative flow and conservation checks: 60 original trajectories and 60 with modified routing. The 30 original secretion ratios reproduce the stored Task 14B ratios exactly. Across 36,120 sampled states from the modified trajectories, the pointwise chloride law, bicarbonate law and every amount derivative outside Na and K agree exactly with the original evaluator at those same states. No opposing AE4 cation fluxes occur. There are 3,252 sampled states with a reversed net AE4 cycle; both cation sources reverse together there. The maximum AE4 charge residual is 1.74 × 10⁻¹⁸ fmol/s, consistent with rounding.

The smallest and largest modified secretion effects were checked using BDF for both WT and 5% AE4. All four confirmation trajectories pass. The largest difference in a total secretion integral between Radau and BDF is 5.60 × 10⁻⁸ in relative terms; the largest relative difference over the state and flow traces is 2.62 × 10⁻⁶, below the inherited 10⁻⁴ comparison limit. Twelve focused tests pass: three tests of the routing intervention and nine existing AE4 tests.

The table below gives every root. Each entry is the original percentage change followed by the percentage change with modified cation routing. Positive values mean increased secretion after AE4 loss.

| Root | Calcium 0.10 µM | Calcium 0.25 µM | Calcium 0.50 µM |
| --- | ---: | ---: | ---: |
| R01 | +15.095% → +5.794% | +19.154% → +7.513% | +20.240% → +8.036% |
| R02 | +9.207% → +4.016% | +16.244% → +5.703% | +18.400% → +6.248% |
| R03 | +14.052% → +5.400% | +17.991% → +7.034% | +19.198% → +7.590% |
| R04 | +8.229% → +3.488% | +14.543% → +4.959% | +16.697% → +5.499% |
| R05 | +16.318% → +6.178% | +20.271% → +7.831% | +21.331% → +8.335% |
| R06 | +12.493% → +4.626% | +19.515% → +6.236% | +21.663% → +6.760% |
| R07 | +15.185% → +5.766% | +19.064% → +7.349% | +20.256% → +7.890% |
| R08 | +11.295% → +4.053% | +17.647% → +5.468% | +19.823% → +5.993% |
| R09 | +15.307% → +5.829% | +19.238% → +7.470% | +20.348% → +7.992% |
| R10 | +10.490% → +4.067% | +17.166% → +5.606% | +19.304% → +6.133% |

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

Source repository: `esig626/ae4-salivary-transport-control`, source commit `3f41579ef81f1e73cc37dc1ce4ab3380f98c4659`. The ten root manifest SHA256 is `6d283961421a6718a54fd7c86e6575bd5c032a93e649923777b02295cee97f29`. The source and inherited result files used here are identical on the existing Task 20 branch. All input hashes remained unchanged during the runs.

This experiment is separate from the earlier Task 20 pooled affinity and capacity recalibration. It neither adopts that new affinity nor uses the recalibrated parameter sets. The inherited branch thermodynamic diagnostics are not reused as a validation of the new cation assignment.

Reproduce from the repository root with Python, NumPy and SciPy available:

```bash
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 -m modern_full_model.run_ae4_routing_only --workers 4
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 -m modern_full_model.run_ae4_routing_only --confirm --workers 4
```

The runner writes the complete state and flux traces as compressed NumPy files. `comparison.csv` contains the 30 primary comparisons; `trajectory_summary.csv` contains the 120 production trajectory summaries; `confirmation_trajectories.csv` records the four BDF runs; `verification.json` records the invariant and solver checks; `contract.json` records the frozen inputs and settings.
