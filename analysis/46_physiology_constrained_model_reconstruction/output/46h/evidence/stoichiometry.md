# Task 46H00 stoichiometry evidence review

Exact source checkouts: Task 43 `547f113d9123ab1976774639a69483faf40cef34`; Task 44 `e5fa9840bf96be3147c6118daee4942468c78f8b`. The relevant production source blobs match the frozen Task 46 base. No simulations or optimisation were run.

Cell vectors below use `(Na, K, Cl, TIC, TA)`. Positive directions are conventions and need separate thermodynamic justification before imposing one way bounds.

| Process | Cell vector | Lumen vector | Source |
| --- | --- | --- | --- |
| NKCC1 | (1, 1, 2, 0, 0) | (0, 0, 0, 0, 0) | membranes:215-249,286-304 |
| AE4 | (-0.5, -0.5, 1, -2, -2) | (0, 0, 0, 0, 0) | ae4:24-63 |
| AE2 | (0, 0, 1, -1, -1) | (0, 0, 0, 0, 0) | membranes:215-249,286-304 |
| NHE1 | (1, 0, 0, 0, 1) | (0, 0, 0, 0, 0) | membranes:124-126,215-249,286-304 |
| NBC | (1, 0, 0, 2, 2) | (0, 0, 0, 0, 0) | nbc:1-28,220-239 |
| pump_b | (-3, 2, 0, 0, 0) | (0, 0, 0, 0, 0) | membranes:380-396,462-470 |
| pump_a | (-3, 2, 0, 0, 0) | (3, -2, 0, 0, 0) | membranes:380-396,462-476 |
| CaCC | (0, 0, -1, 0, 0) | (0, 0, 1, 0, 0) | membranes:400-405,452-476 |
| K_b | (0, -1, 0, 0, 0) | (0, 0, 0, 0, 0) | membranes:398-404,452-470 |
| K_a | (0, -1, 0, 0, 0) | (0, 1, 0, 0, 0) | membranes:398-404,452-476 |
| CO2_b | (0, 0, 0, 1, 0) | (0, 0, 0, 0, 0) | model:624-655 |
| CO2_a | (0, 0, 0, 1, 0) | (0, 0, 0, -1, 0) | model:628-631,649-665 |
| para_Na | (0, 0, 0, 0, 0) | (-1, 0, 0, 0, 0) | membranes:406-411,455-476 |
| para_K | (0, 0, 0, 0, 0) | (0, -1, 0, 0, 0) | membranes:406-411,455-476 |
| para_Cl | (0, 0, 0, 0, 0) | (0, 0, -1, 0, 0) | membranes:406-411,455-476 |
| para_HCO3 | (0, 0, 0, 0, 0) | (0, 0, 0, -1, -1) | membranes:406-411,455-476 |
| out_Na | (0, 0, 0, 0, 0) | (-1, 0, 0, 0, 0) | model:636-642,662-666 |
| out_K | (0, 0, 0, 0, 0) | (0, -1, 0, 0, 0) | model:636-642,662-666 |
| out_Cl | (0, 0, 0, 0, 0) | (0, 0, -1, 0, 0) | model:636-642,662-666 |
| out_TIC | (0, 0, 0, 0, 0) | (0, 0, 0, -1, 0) | model:636-642,662-666 |
| out_TA | (0, 0, 0, 0, 0) | (0, 0, 0, 0, -1) | model:636-642,662-666 |

AE4 equal Na/K routing gives `(-1/2,-1/2,+1,-2,-2)`. It exports both cations and bicarbonate while importing chloride. This is the frozen ensemble source allocation, not an experimentally established microscopic stoichiometry. NBC imports one Na and two bicarbonates and therefore brings one net negative charge into the cell. Neither NBC nor the pump should be forced individually electroneutral.

The inherited network has both apical and basal pump and K routes. Their numerical partition parameters describe capacities or conductances, not measured localisation fractions or necessarily realised flux ratios. Luminal outflow in the source is one water flux times composition; independent CBM component sinks are a relaxation and must retain charge and any declared composition restrictions.

The exact source identities are `2 J4 = H + 2 B - E` and `CaCC = 6(P_a+P_b) - H` at intracellular steady state for the frozen equal routing model. The second does not mean AE4 has no causal effect: the feasible pump/NHE combinations and carbon support may still change.

The finite NBC capacity in the old model was constructed from an assumed WT chloride partition. NKCC amount came from an inherited model resting flux. Neither is directly measured. No numerical bounds are proposed by this review.

Source identities and Git blobs are in the companion JSON.
