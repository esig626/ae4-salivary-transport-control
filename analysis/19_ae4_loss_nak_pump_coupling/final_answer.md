# Task 19 final answer

**Primary classification: `TASK 19 NUMERICALLY INCONCLUSIVE`.** The available 50% pump results oppose the proposed mechanism, while resting-state continuation prevents the full prescribed panel from being evaluated.

The baseline chloride-allocation model was used throughout, with all ten inherited roots, both routing families, both AE4 expression arms, and exactly the three pump and calcium values specified. Pump scaling preserves the total-capacity partition and 3Na:2K law. No other transporter or parameter was changed or fitted.

| Pump capacity | Valid resting/dynamic solutions | Matched AE4 ratio | Pump-only ratio | Pump ratio in AE4-low cells | Interaction I | Coupled ratio versus normal WT |
|---|---|---:|---:|---:|---:|---:|
| 1.0 | All ten roots; 60 trajectories reused | 1.0823–1.2166 | 1 | 1 | 1 | 1.0823–1.2166 |
| 0.50 | Seven roots in both expression arms; 42 new trajectories | 1.3917–1.6278 | 0.5687–0.6549 | 0.6701–0.8776 | 1.1783–1.3584 | 0.8057–1.0018 |
| 0.10 | No valid target resting state under the inherited continuation bounds | unavailable | unavailable | unavailable | unavailable | unavailable |

At 50% pump capacity, pump inhibition alone reduces secretion by **34.5–43.1%**. Matched AE4 near-loss instead **increases** secretion by **39.2–62.8%**. The interaction is substantial and opposite to the hypothesis: AE4-low cells are relatively less affected by the pump reduction. The coupled phenotype falls below normal WT in **20 of 21** valid root/calcium comparisons, by at most **19.4%**; this is not an AE4-specific success.

The effect direction is consistent across both routing families and all calcium values among available cases. All five AE4NA20 roots and two AE4NA05 roots reach pump 0.50; three AE4NA05 roots do not. All ten pump 0.10 paths fail. These failures reach the inherited **60 mM intracellular-Na search boundary**, so they do not establish physical nonexistence or physiological collapse. Failed parents and their dependent arms remain explicit in all 180-condition accounting.

At pump 0.50, AE4 near-loss lowers resting Na by 28.1–30.6 mM and raises K by 33.8–36.3 mM and Cl by 2.8–16.0 mM in the valid matched pairs. Actual AE4 inward-Na and outward-K branch magnitudes shrink; NKCC1 chloride loading rises, NHE1 rises, and AE2 shifts toward chloride loading. These state changes are consistent with relief of the inherited AE4 Na/K transport burden under pump limitation. Full resting and time-resolved flux/conservation records are supplied.

This diagnostic does **not** justify adding an explicit biological AE4–pump coupling law. The valid 50% results argue against the proposed harmful interaction; the 10% condition remains unresolved. The experimental approximately 35% reduction was considered only after the numerical results were frozen, and was never used to select or fit a pump scale.

**Verification:** 72 focused/inherited tests passed; all 187 frozen inputs, production equations, `model/`, and `archive/` are unchanged. Sixty baseline trajectories were reused; 42 new trajectories passed; 78 conditions remain unavailable due to resting continuation failure. The branch-flux reporting correction is explicit and changes no simulated state, secretion value, or factorial comparison.
