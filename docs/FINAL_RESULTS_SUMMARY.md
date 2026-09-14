# Final scientific results used by the manuscript

This document is the concise public summary of Tasks 37-42. The numbered analysis and result directories remain the authoritative audit trail.

## Wild-type acid-base reconstruction

The conservation-explicit model evolves sodium, potassium, chloride, total inorganic carbon, total alkalinity and volume in the cell and lumen. A stimulus-recruited electrogenic NBC-like pathway imports one sodium and two bicarbonate equivalents per inward cycle.

The validated 600 s wild-type trajectory passes the declared physical, physiological and conservation gates. Positive chloride loading over 60-600 s is:

- NKCC1: 79.9248%
- AE4: 20.0752%
- AE2: 0%

This establishes the central acid-base result: an inward alkalinity source is required if AE4 is to carry a substantial stimulated chloride load while pH remains physiological.

## Task 40: equal AE4 cation routing

The total AE4 scalar cycle was retained, while its monovalent-cation source was divided equally between sodium and potassium.

| Quantity | AE4 5% | AE4 null |
| --- | ---: | ---: |
| Cumulative secretion deficit | 3.4176% | 3.8575% |
| NKCC1 compensation | +20.82% | +23.16% |

Equal routing alone is insufficient. NKCC1 replaces most of the missing AE4 chloride loading.

## Task 41: target-selected constructive coupling

A model was constructed in which stimulated CaCC recruitment depends on AE4 expression. The coupling strength was selected against a requested 20-35% secretion-loss range and is therefore not independent validation.

| Quantity | AE4 5% | AE4 null |
| --- | ---: | ---: |
| Cumulative secretion deficit | 23.16% | 30.26% |
| NKCC1 compensation | +5.68% | +2.79% |
| Full-stimulus CaCC conductance retained | 14.99% | 10.51% |

The model passes the declared 600 s gates and retains more intracellular chloride, but intracellular sodium remains below wild type and the early-time secretion deficit is too large. Long-time validity is not established.

## Task 42: Catalán 2025 source classes

Seven predeclared source classes were executed with no fitting, optimisation or parameter sweep. Their predictions were frozen at commit:

`4cb75308bd81d06f064bb80b11831985fa4bb219`

All seven wild-type resting states passed. Eighteen of 21 production trajectories passed. The wild-type trajectories for C2, C4a and C4b failed the pH gate during stimulation.

Among classes with a valid wild-type trajectory:

- C1, the proposed inward NaCl / outward K-HCO3 cycle, predicts 10.76% more secretion after AE4 deletion;
- C3a predicts essentially no phenotype;
- C3b gives the largest valid null deficit, 7.67%, but its forward affinity is opposed at the reference state;
- no class recovers the approximately 35% experimental phenotype.

The narrow conclusion is that source stoichiometry alone, under the inherited Task 40 scalar kinetic law and fixed whole-cell network, is insufficient. This does not falsify the molecular experiments or every kinetic realisation compatible with them.

## Manuscript interpretation

The biological importance of AE4 survives the decade. The original network explanation does not. The remaining unresolved mechanism lies in the interaction between bicarbonate supply, NKCC1 compensation and sustained regulation of apical chloride exit.
