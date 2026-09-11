# Task 19 mechanism interpretation

**Primary classification: `TASK 19 NUMERICALLY INCONCLUSIVE`.**

The full prespecified mechanism test is incomplete because the inherited resting-state search domain blocks all ten 0.10-pump paths and three 0.50-pump paths. Within the seven valid 0.50-pump roots, the result consistently opposes the proposed mechanism: AE4 near-loss increases secretion relative to its pump-matched control, and pump impairment amplifies that relative increase. This is a material interaction in the direction opposite to the hypothesis.

## Coverage and numerical limitation

| Pump capacity | Valid resting states, out of 20 | Valid dynamics, out of 60 | Complete root/calcium comparisons, out of 30 |
|---|---:|---:|---:|
| 1.0 | 20 | 60, all reused | 30 |
| 0.50 | 14 | 42, all newly computed | 21 |
| 0.10 | 0 | 0 | 0 |

All 180 dynamic conditions remain in the result table: 102 have valid trajectories and 78 explicitly have no trajectory because of a failed resting parent. There are no failed ODE integrations among the 42 new trajectories. Every available trajectory passes production Radau and inherited conservation, charge, positivity, and finite-flow gates.

At pump 0.50, both expression arms are available in all five AE4NA20 roots and two of five AE4NA05 roots. The unavailable AE4NA05 roots are `PHIGH_KLOW`, `PLOW_KLOW`, and `PNOM_KMID` (with the full IDs retained in the tables). The valid AE4NA05 roots are `PHIGH_KHIGH` and `PLOW_KHIGH`. Each valid root has all three calcium values and both mandatory expression arms.

Every failed pump path reaches the inherited upper intracellular-Na search boundary of 60 mM; the final rejected candidate also fails the resting residual gate. None is repaired by retuning a transporter, widening a bound, or selecting a disconnected state. The AE4-low continuation cannot start when its matched pump-reduced WT parent is unavailable. A continuation boundary is a limitation of this numerical protocol, not proof of physiological collapse or nonexistence of an equilibrium. In particular, the model's accepted stimulated trajectories can exceed 60 mM Na: that resting-search bound is not a physical positivity or dynamic conservation gate. Consequently, the stronger classification that severe pump impairment *breaks physiology* is not justified.

## Five matched comparisons

Ranges below contain every available root/calcium pair at the indicated pump scale. Ratios are calculated per pair, never as a ratio of pooled means.

| Metric | Pump 1.0, 30 pairs | Pump 0.50, 21 pairs | Pump 0.10 |
|---|---:|---:|---|
| Matched AE4 effect, R_AE4 | 1.082287–1.216630 | 1.391742–1.627816 | unavailable |
| Pump-only effect, R_pump_WT | 1 | 0.568690–0.654919 | unavailable |
| Pump effect in AE4-low cells, R_pump_AE4 | 1 | 0.670101–0.877638 | unavailable |
| Interaction, I | 1 | 1.178323–1.358398 | unavailable |
| Coupled phenotype, R_coupled | 1.082287–1.216630 | 0.805731–1.001802 | unavailable |

At pump 0.50, pump impairment alone reduces secretion by **34.51–43.13%**. It reduces secretion in AE4-low cells by **12.24–32.99%** relative to the corresponding AE4-low, pump-normal cells. Thus AE4-low cells are relatively less impaired by the pump reduction. Matched AE4 near-loss increases secretion by **39.17–62.78%**, rather than reducing it.

All 21 interaction ratios exceed one and the prespecified 1% descriptive materiality threshold. The interaction is therefore neither absent nor supportive of the proposed AE4-specific harm: the relative secretion gain from AE4 near-loss is amplified by **17.83–35.84%**. This also follows independently from I = R_pump_AE4 / R_pump_WT.

The coupled AE4-low/pump-low state lies below normal WT in 20 of 21 valid comparisons. Its largest reduction is 19.43%; one `PLOW_KHIGH_AE4NA20` case at calcium 0.50 uM is 0.18% above normal WT. The coupled phenotype is therefore not uniformly in the experimental direction even within the available subset. Its reduction is smaller than the matched pump-only reduction; it cannot be interpreted as an AE4-specific success.

## Routing and calcium

| Routing family at pump 0.50 | Available pairs | R_AE4 range | I range | R_coupled range |
|---|---:|---:|---:|---:|
| AE4NA05 | 6 / 15 | 1.391742–1.430615 | 1.178323–1.209507 | 0.805731–0.827456 |
| AE4NA20 | 15 / 15 | 1.457643–1.627816 | 1.337971–1.358398 | 0.945545–1.001802 |

The wrong-direction matched effect and opposite-direction interaction occur in both routing families at every available root and calcium value. Effect magnitudes and the coupled phenotype depend on routing; availability also depends strongly on routing. These observations establish consistency within the evaluated subset, not robustness across the complete ten-root ensemble or the 0.10 pump condition. The per-calcium and per-family breakdown is in `factorial_summary.csv`; individual pairs are in `factorial_comparisons.csv`.

## State and transporter interpretation

These are resting-state ranges for the **same seven roots** that have valid 0.50 pump comparisons.

| Quantity | AE4 1, pump 1 | AE4 1, pump 0.50 | AE4 0.05, pump 0.50 |
|---|---:|---:|---:|
| Intracellular Na, mM | 14.33–18.60 | 40.26–59.86 | 12.14–29.27 |
| Intracellular K, mM | 107.07–118.12 | 63.57–88.59 | 97.46–123.82 |
| Intracellular Cl, mM | 51.73–53.02 | 41.19–43.58 | 46.34–57.89 |
| Intracellular pH | 6.841–7.048 | 6.823–6.929 | 7.043–7.277 |
| Cell volume, pL | 1.300 | 1.130–1.197 | 1.284–1.494 |
| Total pump cycles, fmol/s | 0.06356–0.11120 | 0.03668–0.07292 | 0.03537–0.04882 |
| Actual AE4 Na branch, fmol/s | −0.19034 to −0.06739 | −0.11647 to −0.03081 | −0.00988 to −0.00243 |
| Actual AE4 K branch, fmol/s | 0.07179–0.19415 | 0.03449–0.11994 | 0.00574–0.01300 |
| AE4 net Cl loading, fmol/s | 0.003801–0.004408 | 0.003443–0.003704 | 0.003104–0.003331 |
| NKCC1 inward cycles, fmol/s | 0.11576–0.13833 | 0.07280–0.09918 | 0.09741–0.13014 |

The branch convention is inherited: a positive AE4 branch loads Cl and exports its coupled cation. The negative Na branch therefore imports Na, while the positive K branch exports K. Actual branches in the final tables include the expression multiplier; their sum equals the AE4 net Cl source.

At fixed pump 0.50, reducing AE4 expression lowers Na by 28.12–30.59 mM and raises K by 33.84–36.29 mM in every available matched resting pair. It also raises Cl by 2.77–16.04 mM, pH by 0.217–0.348, and volume by 0.089–0.358 pL. The AE4 counter-directed branch magnitudes become much smaller. This is consistent with relief of the AE4-associated Na-loading/K-loss burden under limited pump capacity.

Despite a small decrease in resting net AE4 Cl loading (0.000339–0.000373 fmol/s), NKCC1 inward cycles rise by 0.02444–0.03163 fmol/s in the same pairs. NHE1 increases by 0.000257–0.000418 fmol/s. AE2 shifts by 0.000999–0.001104 fmol/s toward intracellular Cl loading; it is less outward and in some cases becomes inward. Total pump flux falls further as the state-dependent Na drive decreases: scaling capacity does not require steady flux to retain the same proportional change after the state adjusts. The apical/basolateral **capacity** partition remains exactly fixed; the flux partition may change because the two extracellular K environments differ.

The same state contrast persists during stimulation. Across the 21 valid endpoints at 600 s, pump-reduced AE4-intact cells have Na 55.38–75.40 mM, K 49.19–73.94 mM, and Cl 42.53–46.96 mM; AE4-low cells at the same pump scale have Na 15.72–47.97 mM, K 78.03–120.18 mM, and Cl 44.05–58.39 mM. Complete apical/basolateral pump, voltage, AE4 branch, NKCC1, NHE1, AE2, and conservation diagnostics are retained at rest and 60, 180, and 600 s.

## Hypothesis decision

The valid 50% pump results argue against this particular coupling mechanism as an explanation of an AE4-specific secretion deficit. They do not motivate adding a biologically explicit AE4-to-pump regulatory law on the strength of this diagnostic. The untested 10% condition and three unavailable 50% roots prevent a complete exclusion across the prescribed panel; no conclusion about experimental AE4 regulation of Na/K-ATPase follows.

Only after all secretion outputs and comparisons were frozen was the user-supplied experimental context considered: an approximately 35% reduction corresponds to a ratio near 0.65. No valid coupled arm reaches that magnitude, and matching the direction in 20 cases is explained by the substantial pump-only effect alongside a protective relative AE4 interaction. No pump scale was selected, added, or fitted against that value. No continuous coupling law was fitted.

## Verification and provenance

The experiment contract was committed and pushed before new numerical work at `74159c95bedc739d63c0399c5f24166ad65de11d`. All 187 frozen baseline inputs, the production equation files, `model/`, and `archive/` remain unchanged. A reporting-only AE4 branch scaling correction is explicitly documented in `results_frozen.json`; its preceding receipt remains available. Every secretion, trajectory, and factorial-comparison hash is unchanged by that correction.

The final focused and inherited suite passes **72 tests**. Missing historical test fixtures in the local checkout were restored byte-for-byte from their repository blobs and checked against the inherited hashes; they were not recomputed or altered. The suite covers the exact intervention, both membrane pump sources, baseline reuse, complete failed/valid accounting, resting and dynamic numerical gates, actual expression-scaled AE4 branch accounting, and all five comparison formulas.
